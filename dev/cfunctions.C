#include "Farray.h"
#include "Sarray.h"
#include "cfunctions.h"
void interpolation_restriction(Farray &P, Farray &Rop, Farray &RPop) {
  P(-1) = -1.0 / 16.0;
  P(0) = 9.0 / 16.0;
  P(1) = 9.0 / 16.0;
  P(2) = -1.0 / 16.0;

  Rop(-4) = -1.0 / 32.0;
  Rop(-3) = 0.0;
  Rop(-2) = 9.0 / 32.0;
  Rop(-1) = 1.0 / 2.0;
  Rop(0) = 9.0 / 32.0;
  Rop(1) = 0.0;
  Rop(2) = -1.0 / 32.0;

  RPop(-3) = 1.0 / 512.0;
  RPop(-2) = -9.0 / 256.0;
  RPop(-1) = 63.0 / 512.0;
  RPop(0) = 105.0 / 128.0;
  RPop(1) = 63.0 / 512.0;
  RPop(2) = -9.0 / 256.0;
  RPop(3) = 1.0 / 512.0;
}

void central_difference(Farray &ux_cof, Farray &uxx_cof) {
  ux_cof(-2) = 1.0 / 12.0;
  ux_cof(-1) = -2.0 / 3.0;
  ux_cof(0) = 0.0;
  ux_cof(1) = 2.0 / 3.0;
  ux_cof(2) = -1.0 / 12.0;

  uxx_cof(-2) = -1.0 / 12.0;
  uxx_cof(-1) = 4.0 / 3.0;
  uxx_cof(0) = -5.0 / 2.0;
  uxx_cof(1) = 4.0 / 3.0;
  uxx_cof(2) = -1.0 / 12.0;
}

void varcoeffs4(Farray &acof, Farray &ghcof, Farray &Sb) {
  // acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j);
  // ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.;
  ghcof(1) = 12.0 / 17.0;
  ghcof(2) = 0;
  ghcof(3) = 0;
  ghcof(4) = 0;
  ghcof(5) = 0;
  ghcof(6) = 0;
  acof(1, 1, 1) = 104.0 / 289.0;
  acof(1, 1, 2) = -2476335.0 / 2435692.0;
  acof(1, 1, 3) = -16189.0 / 84966.0;
  acof(1, 1, 4) = -9.0 / 3332.0;
  acof(1, 1, 5) = 0;
  acof(1, 1, 6) = 0;
  acof(1, 1, 7) = 0;
  acof(1, 1, 8) = 0;
  acof(1, 2, 1) = -516.0 / 289.0;
  acof(1, 2, 2) = 544521.0 / 1217846.0;
  acof(1, 2, 3) = 2509879.0 / 3653538.0;
  acof(1, 2, 4) = 0;
  acof(1, 2, 5) = 0;
  acof(1, 2, 6) = 0;
  acof(1, 2, 7) = 0;
  acof(1, 2, 8) = 0;
  acof(1, 3, 1) = 312.0 / 289.0;
  acof(1, 3, 2) = 1024279.0 / 2435692.0;
  acof(1, 3, 3) = -687797.0 / 1217846.0;
  acof(1, 3, 4) = 177.0 / 3332.0;
  acof(1, 3, 5) = 0;
  acof(1, 3, 6) = 0;
  acof(1, 3, 7) = 0;
  acof(1, 3, 8) = 0;
  acof(1, 4, 1) = -104.0 / 289.0;
  acof(1, 4, 2) = 181507.0 / 1217846.0;
  acof(1, 4, 3) = 241309.0 / 3653538.0;
  acof(1, 4, 4) = 0;
  acof(1, 4, 5) = 0;
  acof(1, 4, 6) = 0;
  acof(1, 4, 7) = 0;
  acof(1, 4, 8) = 0;
  acof(1, 5, 1) = 0;
  acof(1, 5, 2) = 0;
  acof(1, 5, 3) = 5.0 / 2193.0;
  acof(1, 5, 4) = -48.0 / 833.0;
  acof(1, 5, 5) = 0;
  acof(1, 5, 6) = 0;
  acof(1, 5, 7) = 0;
  acof(1, 5, 8) = 0;
  acof(1, 6, 1) = 0;
  acof(1, 6, 2) = 0;
  acof(1, 6, 3) = 0;
  acof(1, 6, 4) = 6.0 / 833.0;
  acof(1, 6, 5) = 0;
  acof(1, 6, 6) = 0;
  acof(1, 6, 7) = 0;
  acof(1, 6, 8) = 0;
  acof(1, 7, 1) = 0;
  acof(1, 7, 2) = 0;
  acof(1, 7, 3) = 0;
  acof(1, 7, 4) = 0;
  acof(1, 7, 5) = 0;
  acof(1, 7, 6) = 0;
  acof(1, 7, 7) = 0;
  acof(1, 7, 8) = 0;
  acof(1, 8, 1) = 0;
  acof(1, 8, 2) = 0;
  acof(1, 8, 3) = 0;
  acof(1, 8, 4) = 0;
  acof(1, 8, 5) = 0;
  acof(1, 8, 6) = 0;
  acof(1, 8, 7) = 0;
  acof(1, 8, 8) = 0;
  acof(2, 1, 1) = 12.0 / 17.0;
  acof(2, 1, 2) = 544521.0 / 4226642.0;
  acof(2, 1, 3) = 2509879.0 / 12679926.0;
  acof(2, 1, 4) = 0;
  acof(2, 1, 5) = 0;
  acof(2, 1, 6) = 0;
  acof(2, 1, 7) = 0;
  acof(2, 1, 8) = 0;
  acof(2, 2, 1) = -59.0 / 68.0;
  acof(2, 2, 2) = -1633563.0 / 4226642.0;
  acof(2, 2, 3) = -21510077.0 / 25359852.0;
  acof(2, 2, 4) = -12655.0 / 372939.0;
  acof(2, 2, 5) = 0;
  acof(2, 2, 6) = 0;
  acof(2, 2, 7) = 0;
  acof(2, 2, 8) = 0;
  acof(2, 3, 1) = 2.0 / 17.0;
  acof(2, 3, 2) = 1633563.0 / 4226642.0;
  acof(2, 3, 3) = 2565299.0 / 4226642.0;
  acof(2, 3, 4) = 40072.0 / 372939.0;
  acof(2, 3, 5) = 0;
  acof(2, 3, 6) = 0;
  acof(2, 3, 7) = 0;
  acof(2, 3, 8) = 0;
  acof(2, 4, 1) = 3.0 / 68.0;
  acof(2, 4, 2) = -544521.0 / 4226642.0;
  acof(2, 4, 3) = 987685.0 / 25359852.0;
  acof(2, 4, 4) = -14762.0 / 124313.0;
  acof(2, 4, 5) = 0;
  acof(2, 4, 6) = 0;
  acof(2, 4, 7) = 0;
  acof(2, 4, 8) = 0;
  acof(2, 5, 1) = 0;
  acof(2, 5, 2) = 0;
  acof(2, 5, 3) = 1630.0 / 372939.0;
  acof(2, 5, 4) = 18976.0 / 372939.0;
  acof(2, 5, 5) = 0;
  acof(2, 5, 6) = 0;
  acof(2, 5, 7) = 0;
  acof(2, 5, 8) = 0;
  acof(2, 6, 1) = 0;
  acof(2, 6, 2) = 0;
  acof(2, 6, 3) = 0;
  acof(2, 6, 4) = -1.0 / 177.0;
  acof(2, 6, 5) = 0;
  acof(2, 6, 6) = 0;
  acof(2, 6, 7) = 0;
  acof(2, 6, 8) = 0;
  acof(2, 7, 1) = 0;
  acof(2, 7, 2) = 0;
  acof(2, 7, 3) = 0;
  acof(2, 7, 4) = 0;
  acof(2, 7, 5) = 0;
  acof(2, 7, 6) = 0;
  acof(2, 7, 7) = 0;
  acof(2, 7, 8) = 0;
  acof(2, 8, 1) = 0;
  acof(2, 8, 2) = 0;
  acof(2, 8, 3) = 0;
  acof(2, 8, 4) = 0;
  acof(2, 8, 5) = 0;
  acof(2, 8, 6) = 0;
  acof(2, 8, 7) = 0;
  acof(2, 8, 8) = 0;
  acof(3, 1, 1) = -96.0 / 731.0;
  acof(3, 1, 2) = 1024279.0 / 6160868.0;
  acof(3, 1, 3) = -687797.0 / 3080434.0;
  acof(3, 1, 4) = 177.0 / 8428.0;
  acof(3, 1, 5) = 0;
  acof(3, 1, 6) = 0;
  acof(3, 1, 7) = 0;
  acof(3, 1, 8) = 0;
  acof(3, 2, 1) = 118.0 / 731.0;
  acof(3, 2, 2) = 1633563.0 / 3080434.0;
  acof(3, 2, 3) = 2565299.0 / 3080434.0;
  acof(3, 2, 4) = 40072.0 / 271803.0;
  acof(3, 2, 5) = 0;
  acof(3, 2, 6) = 0;
  acof(3, 2, 7) = 0;
  acof(3, 2, 8) = 0;
  acof(3, 3, 1) = -16.0 / 731.0;
  acof(3, 3, 2) = -5380447.0 / 6160868.0;
  acof(3, 3, 3) = -3569115.0 / 3080434.0;
  acof(3, 3, 4) = -331815.0 / 362404.0;
  acof(3, 3, 5) = -283.0 / 6321.0;
  acof(3, 3, 6) = 0;
  acof(3, 3, 7) = 0;
  acof(3, 3, 8) = 0;
  acof(3, 4, 1) = -6.0 / 731.0;
  acof(3, 4, 2) = 544521.0 / 3080434.0;
  acof(3, 4, 3) = 2193521.0 / 3080434.0;
  acof(3, 4, 4) = 8065.0 / 12943.0;
  acof(3, 4, 5) = 381.0 / 2107.0;
  acof(3, 4, 6) = 0;
  acof(3, 4, 7) = 0;
  acof(3, 4, 8) = 0;
  acof(3, 5, 1) = 0;
  acof(3, 5, 2) = 0;
  acof(3, 5, 3) = -14762.0 / 90601.0;
  acof(3, 5, 4) = 32555.0 / 271803.0;
  acof(3, 5, 5) = -283.0 / 2107.0;
  acof(3, 5, 6) = 0;
  acof(3, 5, 7) = 0;
  acof(3, 5, 8) = 0;
  acof(3, 6, 1) = 0;
  acof(3, 6, 2) = 0;
  acof(3, 6, 3) = 0;
  acof(3, 6, 4) = 9.0 / 2107.0;
  acof(3, 6, 5) = -11.0 / 6321.0;
  acof(3, 6, 6) = 0;
  acof(3, 6, 7) = 0;
  acof(3, 6, 8) = 0;
  acof(3, 7, 1) = 0;
  acof(3, 7, 2) = 0;
  acof(3, 7, 3) = 0;
  acof(3, 7, 4) = 0;
  acof(3, 7, 5) = 0;
  acof(3, 7, 6) = 0;
  acof(3, 7, 7) = 0;
  acof(3, 7, 8) = 0;
  acof(3, 8, 1) = 0;
  acof(3, 8, 2) = 0;
  acof(3, 8, 3) = 0;
  acof(3, 8, 4) = 0;
  acof(3, 8, 5) = 0;
  acof(3, 8, 6) = 0;
  acof(3, 8, 7) = 0;
  acof(3, 8, 8) = 0;
  acof(4, 1, 1) = -36.0 / 833.0;
  acof(4, 1, 2) = 181507.0 / 3510262.0;
  acof(4, 1, 3) = 241309.0 / 10530786.0;
  acof(4, 1, 4) = 0;
  acof(4, 1, 5) = 0;
  acof(4, 1, 6) = 0;
  acof(4, 1, 7) = 0;
  acof(4, 1, 8) = 0;
  acof(4, 2, 1) = 177.0 / 3332.0;
  acof(4, 2, 2) = -544521.0 / 3510262.0;
  acof(4, 2, 3) = 987685.0 / 21061572.0;
  acof(4, 2, 4) = -14762.0 / 103243.0;
  acof(4, 2, 5) = 0;
  acof(4, 2, 6) = 0;
  acof(4, 2, 7) = 0;
  acof(4, 2, 8) = 0;
  acof(4, 3, 1) = -6.0 / 833.0;
  acof(4, 3, 2) = 544521.0 / 3510262.0;
  acof(4, 3, 3) = 2193521.0 / 3510262.0;
  acof(4, 3, 4) = 8065.0 / 14749.0;
  acof(4, 3, 5) = 381.0 / 2401.0;
  acof(4, 3, 6) = 0;
  acof(4, 3, 7) = 0;
  acof(4, 3, 8) = 0;
  acof(4, 4, 1) = -9.0 / 3332.0;
  acof(4, 4, 2) = -181507.0 / 3510262.0;
  acof(4, 4, 3) = -2647979.0 / 3008796.0;
  acof(4, 4, 4) = -80793.0 / 103243.0;
  acof(4, 4, 5) = -1927.0 / 2401.0;
  acof(4, 4, 6) = -2.0 / 49.0;
  acof(4, 4, 7) = 0;
  acof(4, 4, 8) = 0;
  acof(4, 5, 1) = 0;
  acof(4, 5, 2) = 0;
  acof(4, 5, 3) = 57418.0 / 309729.0;
  acof(4, 5, 4) = 51269.0 / 103243.0;
  acof(4, 5, 5) = 1143.0 / 2401.0;
  acof(4, 5, 6) = 8.0 / 49.0;
  acof(4, 5, 7) = 0;
  acof(4, 5, 8) = 0;
  acof(4, 6, 1) = 0;
  acof(4, 6, 2) = 0;
  acof(4, 6, 3) = 0;
  acof(4, 6, 4) = -283.0 / 2401.0;
  acof(4, 6, 5) = 403.0 / 2401.0;
  acof(4, 6, 6) = -6.0 / 49.0;
  acof(4, 6, 7) = 0;
  acof(4, 6, 8) = 0;
  acof(4, 7, 1) = 0;
  acof(4, 7, 2) = 0;
  acof(4, 7, 3) = 0;
  acof(4, 7, 4) = 0;
  acof(4, 7, 5) = 0;
  acof(4, 7, 6) = 0;
  acof(4, 7, 7) = 0;
  acof(4, 7, 8) = 0;
  acof(4, 8, 1) = 0;
  acof(4, 8, 2) = 0;
  acof(4, 8, 3) = 0;
  acof(4, 8, 4) = 0;
  acof(4, 8, 5) = 0;
  acof(4, 8, 6) = 0;
  acof(4, 8, 7) = 0;
  acof(4, 8, 8) = 0;
  acof(5, 1, 1) = 0;
  acof(5, 1, 2) = 0;
  acof(5, 1, 3) = 5.0 / 6192.0;
  acof(5, 1, 4) = -1.0 / 49.0;
  acof(5, 1, 5) = 0;
  acof(5, 1, 6) = 0;
  acof(5, 1, 7) = 0;
  acof(5, 1, 8) = 0;
  acof(5, 2, 1) = 0;
  acof(5, 2, 2) = 0;
  acof(5, 2, 3) = 815.0 / 151704.0;
  acof(5, 2, 4) = 1186.0 / 18963.0;
  acof(5, 2, 5) = 0;
  acof(5, 2, 6) = 0;
  acof(5, 2, 7) = 0;
  acof(5, 2, 8) = 0;
  acof(5, 3, 1) = 0;
  acof(5, 3, 2) = 0;
  acof(5, 3, 3) = -7381.0 / 50568.0;
  acof(5, 3, 4) = 32555.0 / 303408.0;
  acof(5, 3, 5) = -283.0 / 2352.0;
  acof(5, 3, 6) = 0;
  acof(5, 3, 7) = 0;
  acof(5, 3, 8) = 0;
  acof(5, 4, 1) = 0;
  acof(5, 4, 2) = 0;
  acof(5, 4, 3) = 28709.0 / 151704.0;
  acof(5, 4, 4) = 51269.0 / 101136.0;
  acof(5, 4, 5) = 381.0 / 784.0;
  acof(5, 4, 6) = 1.0 / 6.0;
  acof(5, 4, 7) = 0;
  acof(5, 4, 8) = 0;
  acof(5, 5, 1) = 0;
  acof(5, 5, 2) = 0;
  acof(5, 5, 3) = -349.0 / 7056.0;
  acof(5, 5, 4) = -247951.0 / 303408.0;
  acof(5, 5, 5) = -577.0 / 784.0;
  acof(5, 5, 6) = -5.0 / 6.0;
  acof(5, 5, 7) = -1.0 / 24.0;
  acof(5, 5, 8) = 0;
  acof(5, 6, 1) = 0;
  acof(5, 6, 2) = 0;
  acof(5, 6, 3) = 0;
  acof(5, 6, 4) = 1135.0 / 7056.0;
  acof(5, 6, 5) = 1165.0 / 2352.0;
  acof(5, 6, 6) = 1.0 / 2.0;
  acof(5, 6, 7) = 1.0 / 6.0;
  acof(5, 6, 8) = 0;
  acof(5, 7, 1) = 0;
  acof(5, 7, 2) = 0;
  acof(5, 7, 3) = 0;
  acof(5, 7, 4) = 0;
  acof(5, 7, 5) = -1.0 / 8.0;
  acof(5, 7, 6) = 1.0 / 6.0;
  acof(5, 7, 7) = -1.0 / 8.0;
  acof(5, 7, 8) = 0;
  acof(5, 8, 1) = 0;
  acof(5, 8, 2) = 0;
  acof(5, 8, 3) = 0;
  acof(5, 8, 4) = 0;
  acof(5, 8, 5) = 0;
  acof(5, 8, 6) = 0;
  acof(5, 8, 7) = 0;
  acof(5, 8, 8) = 0;
  acof(6, 1, 1) = 0;
  acof(6, 1, 2) = 0;
  acof(6, 1, 3) = 0;
  acof(6, 1, 4) = 1.0 / 392.0;
  acof(6, 1, 5) = 0;
  acof(6, 1, 6) = 0;
  acof(6, 1, 7) = 0;
  acof(6, 1, 8) = 0;
  acof(6, 2, 1) = 0;
  acof(6, 2, 2) = 0;
  acof(6, 2, 3) = 0;
  acof(6, 2, 4) = -1.0 / 144.0;
  acof(6, 2, 5) = 0;
  acof(6, 2, 6) = 0;
  acof(6, 2, 7) = 0;
  acof(6, 2, 8) = 0;
  acof(6, 3, 1) = 0;
  acof(6, 3, 2) = 0;
  acof(6, 3, 3) = 0;
  acof(6, 3, 4) = 3.0 / 784.0;
  acof(6, 3, 5) = -11.0 / 7056.0;
  acof(6, 3, 6) = 0;
  acof(6, 3, 7) = 0;
  acof(6, 3, 8) = 0;
  acof(6, 4, 1) = 0;
  acof(6, 4, 2) = 0;
  acof(6, 4, 3) = 0;
  acof(6, 4, 4) = -283.0 / 2352.0;
  acof(6, 4, 5) = 403.0 / 2352.0;
  acof(6, 4, 6) = -1.0 / 8.0;
  acof(6, 4, 7) = 0;
  acof(6, 4, 8) = 0;
  acof(6, 5, 1) = 0;
  acof(6, 5, 2) = 0;
  acof(6, 5, 3) = 0;
  acof(6, 5, 4) = 1135.0 / 7056.0;
  acof(6, 5, 5) = 1165.0 / 2352.0;
  acof(6, 5, 6) = 1.0 / 2.0;
  acof(6, 5, 7) = 1.0 / 6.0;
  acof(6, 5, 8) = 0;
  acof(6, 6, 1) = 0;
  acof(6, 6, 2) = 0;
  acof(6, 6, 3) = 0;
  acof(6, 6, 4) = -47.0 / 1176.0;
  acof(6, 6, 5) = -5869.0 / 7056.0;
  acof(6, 6, 6) = -3.0 / 4.0;
  acof(6, 6, 7) = -5.0 / 6.0;
  acof(6, 6, 8) = -1.0 / 24.0;
  acof(6, 7, 1) = 0;
  acof(6, 7, 2) = 0;
  acof(6, 7, 3) = 0;
  acof(6, 7, 4) = 0;
  acof(6, 7, 5) = 1.0 / 6.0;
  acof(6, 7, 6) = 1.0 / 2.0;
  acof(6, 7, 7) = 1.0 / 2.0;
  acof(6, 7, 8) = 1.0 / 6.0;
  acof(6, 8, 1) = 0;
  acof(6, 8, 2) = 0;
  acof(6, 8, 3) = 0;
  acof(6, 8, 4) = 0;
  acof(6, 8, 5) = 0;
  acof(6, 8, 6) = -1.0 / 8.0;
  acof(6, 8, 7) = 1.0 / 6.0;
  acof(6, 8, 8) = -1.0 / 8.0;
  // 129 non-zero out of 384.;
  ;
  Sb(0) = -3.0 / 12.0;
  Sb(1) = -10.0 / 12.0;
  Sb(2) = 18.0 / 12.0;
  Sb(3) = -6.0 / 12.0;
  Sb(4) = 1.0 / 12.0;
}

void varcoeff_noghost(Farray &acof_no_gp, Farray &ghcof_no_gp,
                      Farray &sbop_no_gp) {
  Farray d5(0, 8);
  float_sw4 w0;
  Farray acof(1, 6, 1, 8, 1, 8), ghcof(1, 6);
  int i, j, k;
  varcoeffs4(acof, ghcof, sbop_no_gp);
  acof_no_gp = acof;
  d5 = 0.0;
  // Use 5th divided difference to cancel the ghost point contribution
  d5(0) = -1.00;
  d5(1) = 5.00;
  d5(2) = -10.00;
  d5(3) = 10.00;
  d5(4) = -5.00;
  d5(5) = 1.00;
  w0 = 17.0 / 48.0;
  i = 1;
  k = 1;  // only depends on the coefficient a(1)
  for (j = 1; j <= 8; j++)
    acof_no_gp(i, j, k) = acof(i, j, k) + d5(j) / (4.0 * w0);

  // the coeff for all ghost points are zero (don't divided by them!)
  ghcof_no_gp = 0.0;

  // boundary normal derivative, not using ghost points
  // sb = (-25*f(1)/12 + 4*f(2) - 3*f(3) + 4*f(4)/3 - f(5)/4)/h(q);
  sbop_no_gp(0) = 0.0;
  sbop_no_gp(1) = -25.0 / 12.0;
  sbop_no_gp(2) = 4.0;
  sbop_no_gp(3) = -3.0;
  sbop_no_gp(4) = 4.0 / 3.0;
  sbop_no_gp(5) = -1.0 / 4.0;
}

void dx_46(Farray &bop) {
  bop(1, 1) = -24.0 / 17.0;
  bop(1, 2) = 59.0 / 34.0;
  bop(1, 3) = -4.0 / 17.0;
  bop(1, 4) = -3.0 / 34.0;
  bop(1, 5) = 0;
  bop(1, 6) = 0;
  bop(2, 1) = -1.0 / 2.0;
  bop(2, 2) = 0;
  bop(2, 3) = 1.0 / 2.0;
  bop(2, 4) = 0;
  bop(2, 5) = 0;
  bop(2, 6) = 0;
  bop(3, 1) = 4.0 / 43.0;
  bop(3, 2) = -59.0 / 86.0;
  bop(3, 3) = 0;
  bop(3, 4) = 59.0 / 86.0;
  bop(3, 5) = -4.0 / 43.0;
  bop(3, 6) = 0;
  bop(4, 1) = 3.0 / 98.0;
  bop(4, 2) = 0;
  bop(4, 3) = -59.0 / 98.0;
  bop(4, 4) = 0;
  bop(4, 5) = 32.0 / 49.0;
  bop(4, 6) = -4.0 / 49.0;
}

// describe the interface of geometry
inline float_sw4 interface_cf(float_sw4 x, float_sw4 y, PackArgs &a) {
  return a.int_pos * a.l3 + a.amp * sin(4.0 * a.pi * x) +
         a.amp * cos(4.0 * a.pi * y);
}

inline float_sw4 interface_cfx(float_sw4 x, float_sw4 y, PackArgs &a) {
  return a.amp * 4.0 * a.pi * cos(4.0 * a.pi * x);
}

inline float_sw4 interface_cfy(float_sw4 x, float_sw4 y, PackArgs &a) {
  return -a.amp * 4.0 * a.pi * sin(4.0 * a.pi * y);
}

// describing the top geometry
inline float_sw4 top(float_sw4 x, float_sw4 y, PackArgs &a) {
  return a.l3 + a.amp * exp(-pow((x - 0.50), 2) / a.peak) +
         a.amp * exp(-pow((y - 0.50), 2) / a.peak);
}
inline float_sw4 topx(float_sw4 x, float_sw4 y, PackArgs &a) {
  return -a.amp * exp(-pow((x - 0.50), 2) / a.peak) * (2.0 * x - 1.0) / a.peak;
}
inline float_sw4 topy(float_sw4 x, float_sw4 y, PackArgs &a) {
  return -a.amp * exp(-pow((y - 0.50), 2) / a.peak) * (2.0 * y - 1.0) / a.peak;
}

// describe the bottom geometry
inline float_sw4 bottom(float_sw4 x, float_sw4 y, PackArgs &a) {
  return a.amp * exp(-pow((x - 0.60), 2) / a.peak) +
         a.amp * exp(-pow((y - 0.60), 2) / a.peak);
}
inline float_sw4 bottomx(float_sw4 x, float_sw4 y, PackArgs &a) {
  return -a.amp * exp(-pow((x - 0.60), 2) / a.peak) * (2.0 * x - 1.20) / a.peak;
}
inline float_sw4 bottomy(float_sw4 x, float_sw4 y, PackArgs &a) {
  return -a.amp * exp(-pow((y - 0.60), 2) / a.peak) * (2.0 * y - 1.20) / a.peak;
}

void generate_grid(Farray &X1_c, Farray &X2_c, Farray &X3_c, Farray &X1_f,
                   Farray &X2_f, Farray &X3_f, PackArgs &a) {
  int i, j, k;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  float_sw4 l3 = a.l3;

  int n1_c = a.n1_c;
  int n2_c = a.n2_c;
  int n3_c = a.n3_c;

  int n1_f = a.n1_f;
  int n2_f = a.n2_f;
  int n3_f = a.n3_f;

  auto nrg = a.nrg;

  // coarse domain
  for (i = X1_c.debut(); i <= X1_c.fin(); i++)
    X1_c(i) = float_sw4(i - 1) / (n1_c - 1) * l1;
  for (i = 1 - nrg; i <= n2_c + nrg; i++)
    X2_c(i) = float_sw4(i - 1) / (n2_c - 1) * l2;
  for (i = 1 - nrg; i <= n3_c + nrg; i++) {
    for (j = 1 - nrg; j <= n2_c + nrg; j++) {
      for (k = 1 - nrg; k <= n1_c + nrg; k++) {
        X3_c(k, j, i) = float_sw4(i - 1) / (n3_c - 1) *
                            interface_cf(float_sw4(k - 1) / (n1_c - 1),
                                         float_sw4(j - 1) / (n2_c - 1), a) +
                        (1.0 - float_sw4(i - 1) / (n3_c - 1)) *
                            bottom(float_sw4(k - 1) / (n1_c - 1),
                                   float_sw4(j - 1) / (n2_c - 1), a);
      }
    }
  }

  // fine domain
  for (i = 1 - nrg; i <= n1_f + nrg; i++)
    X1_f(i) = float_sw4(i - 1) / (n1_f - 1) * l1;

  for (i = 1 - nrg; i <= n2_f + nrg; i++)
    X2_f(i) = float_sw4(i - 1) / (n2_f - 1) * l2;
  for (i = 1 - nrg; i <= n3_f + nrg; i++) {
    for (j = 1 - nrg; j <= n2_f + nrg; j++) {
      for (k = 1 - nrg; k <= n1_f + nrg; k++) {
        X3_f(k, j, i) = float_sw4(i - 1) / (n3_f - 1) *
                            top(float_sw4(k - 1) / (n1_f - 1),
                                float_sw4(j - 1) / (n2_f - 1), a) +
                        (1.0 - float_sw4(i - 1) / (n3_f - 1)) *
                            interface_cf(float_sw4(k - 1) / (n1_f - 1),
                                         float_sw4(j - 1) / (n2_f - 1), a);
      }
    }
  }
}

void metric_derivative(float_sw4 r1, float_sw4 r2, float_sw4 r3,
                       float_sw4 &xi13, float_sw4 &xi23, float_sw4 &xi33,
                       float_sw4 &J, int flag, PackArgs &a) {
  float_sw4 x3r1, x3r2, x3r3;
  // forward derivatives
  switch (flag) {
    case 0:  // coarse
      x3r1 = r3 * interface_cfx(r1, r2, a) + (1.0 - r3) * bottomx(r1, r2, a);
      x3r2 = r3 * interface_cfy(r1, r2, a) + (1.0 - r3) * bottomy(r1, r2, a);
      x3r3 = interface_cf(r1, r2, a) - bottom(r1, r2, a);
      break;
    case 1:  // fine
      x3r1 = r3 * topx(r1, r2, a) + (1.0 - r3) * interface_cfx(r1, r2, a);
      x3r2 = r3 * topy(r1, r2, a) + (1.0 - r3) * interface_cfy(r1, r2, a);
      x3r3 = top(r1, r2, a) - interface_cf(r1, r2, a);
      break;
    default:
      std::cerr << "This should never happen\n";
      abort();
      break;
  }
  J = a.l1 * a.l2 * x3r3;
  // backward derivative
  xi13 = -a.l2 * x3r1 / J;
  xi23 = -a.l1 * x3r2 / J;
  xi33 = a.l1 * a.l2 / J;
}

void equation_cof(Farray &Xgrid_c_1, Farray &Xgrid_c_2, Farray &Xgrid_c_3,
                  Farray &Xgrid_f_1, Farray &Xgrid_f_2, Farray &Xgrid_f_3,
                  Sarray &XI13_c, Sarray &XI23_c, Sarray &XI33_c,
                  Sarray &Jacobian_c, Sarray &rho_c, Sarray &rho_f,
                  Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f,
                  Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                  Sarray &lambda_c, Sarray &lambda_f, PackArgs &a)

{
  int i, j, k;
  int nrg = a.nrg;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  float_sw4 l3 = a.l3;

  int n1_c = a.n1_c;
  int n2_c = a.n2_c;
  int n3_c = a.n3_c;

  int n1_f = a.n1_f;
  int n2_f = a.n2_f;
  int n3_f = a.n3_f;
  generate_grid(Xgrid_c_1, Xgrid_c_2, Xgrid_c_3, Xgrid_f_1, Xgrid_f_2,
                Xgrid_f_3, a);

  // compute metric derivatives, coarse
  for (i = 1 - nrg; i <= n3_c + nrg; i++) {
    for (j = 1 - nrg; j <= n2_c + nrg; j++) {
      for (k = 1 - nrg; k <= n1_c + nrg; k++) {
        metric_derivative(
            float_sw4(k - 1) / (n1_c - 1), float_sw4(j - 1) / (n2_c - 1),
            float_sw4(i - 1) / (n3_c - 1), XI13_c(k, j, i), XI23_c(k, j, i),
            XI33_c(k, j, i), Jacobian_c(k, j, i), 0, a);
      }
    }
  }
  // compute metric derivatives, fine
  for (i = 1 - nrg; i <= n3_f + nrg; i++) {
    for (j = 1 - nrg; j <= n2_f + nrg; j++) {
      for (k = 1 - nrg; k <= n1_f + nrg; k++) {
        metric_derivative(
            float_sw4(k - 1) / (n1_f - 1), float_sw4(j - 1) / (n2_f - 1),
            float_sw4(i - 1) / (n3_f - 1), XI13_f(k, j, i), XI23_f(k, j, i),
            XI33_f(k, j, i), Jacobian_f(k, j, i), 1, a);
      }
    }
  }
  // variable coefficients
  for (i = 1 - nrg; i <= n3_c + nrg; i++) {
    for (j = 1 - nrg; j <= n2_c + nrg; j++) {
      for (k = 1 - nrg; k <= n1_c + nrg; k++) {
        mu_c(k, j, i) = 3.0 + sin(3.0 * Xgrid_c_1(k) + 0.10) *
                                  sin(3.0 * Xgrid_c_2(j) + 0.10) *
                                  sin(Xgrid_c_3(k, j, i));
        lambda_c(k, j, i) = 21.0 + cos(Xgrid_c_1(k) + 0.10) *
                                       cos(Xgrid_c_2(j) + 0.10) *
                                       pow(sin(3.0 * Xgrid_c_3(k, j, i)), 2);
        rho_c(k, j, i) = 2.0 + sin(Xgrid_c_1(k) + 0.30) *
                                   sin(Xgrid_c_2(j) + 0.30) *
                                   sin(Xgrid_c_3(k, j, i) - 0.20);
      }
    }
  }
  for (i = 1 - nrg; i <= n3_f + nrg; i++) {
    for (j = 1 - nrg; j <= n2_f + nrg; j++) {
      for (k = 1 - nrg; k <= n1_f + nrg; k++) {
        mu_f(k, j, i) = 3.0 + sin(3.0 * Xgrid_f_1(k) + 0.10) *
                                  sin(3.0 * Xgrid_f_2(j) + 0.10) *
                                  sin(Xgrid_f_3(k, j, i));
        lambda_f(k, j, i) = 21.0 + cos(Xgrid_f_1(k) + 0.10) *
                                       cos(Xgrid_f_2(j) + 0.10) *
                                       pow(sin(3.0 * Xgrid_f_3(k, j, i)), 2);
        rho_f(k, j, i) = 2.0 + sin(Xgrid_f_1(k) + 0.30) *
                                   sin(Xgrid_f_2(j) + 0.30) *
                                   sin(Xgrid_f_3(k, j, i) - 0.20);
      }
    }
  }

  rho_c *= Jacobian_c;
  rho_f *= Jacobian_f;
}
extern "C" {
void dsyev_(char &, char &, int &, float_sw4 *, int &, float_sw4 *, float_sw4 *,
            int &, int &);
}

void kappa(Sarray &mu_f, Sarray &lambda_f, Sarray &XI13_f, Sarray &XI23_f,
           Sarray &XI33_f, Sarray &Jacobian_f, Sarray &rho_f, float_sw4 &kappa2,
           PackArgs &a) {
  Farray kappa1(1, 3);
  Farray mat_t(1, 3, 1, 3);
  Farray work(1, 8);
  int info = -1;
  int i, j, k;

  kappa2 = 0.0;
  for (k = 1; k <= a.n3_f; k++) {
    for (j = 1; j <= a.n2_f; j++) {
      for (i = 1; i <= a.n1_f; i++) {
        mat_t(1, 1) = (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) / pow(a.l1, 2);
        mat_t(1, 2) = 0.0;
        mat_t(1, 3) =
            (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) / a.l1 * XI13_f(i, j, k);
        mat_t(2, 1) = 0.0;
        mat_t(2, 2) = (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) / pow(a.l2, 2);
        mat_t(2, 3) =
            (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) / a.l2 * XI23_f(i, j, k);
        mat_t(3, 1) = mat_t(1, 3);
        mat_t(3, 2) = mat_t(2, 3);
        mat_t(3, 3) = (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) *
                          XI13_f(i, j, k) * XI13_f(i, j, k) +
                      (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) *
                          XI23_f(i, j, k) * XI23_f(i, j, k) +
                      (4.0 * mu_f(i, j, k) + lambda_f(i, j, k)) *
                          XI33_f(i, j, k) * XI33_f(i, j, k);
        // std::cout<<"AND HERE "<<mat_t(1,1)<<"\n";
        for (int l = 0; l < 9; l++)
          mat_t[l] = Jacobian_f(i, j, k) * mat_t[l] / rho_f(i, j, k);
        char N = 'N';
        char U = 'U';
        int three = 3;
        int eight = 8;
        // std::cout<<"mat_t before "; for (int l=0;l<9;l++)
        // std::cout<<mat_t[l]<<" ";std::cout<<"\n";
        // std::cout<<"INPUTS"<<Jacobian_f(i,j,k)<<" "<<rho_f(i,j,k)<<"
        // "<<a.l1<<"\n";
        dsyev_(N, U, three, mat_t.get(), three, kappa1.get(), work.get(), eight,
               info);
        if (info != 0) {
          std::cerr << "DSYEV FAILED in KAPPA " << info << mat_t(info, info)
                    << "\n";
          abort();
        }
        // std::cout<<"kappa "<<kappa1(1)<<" "<<kappa1(2)<<" "<<kappa1(3)<<"\n";
        // std::cout<<"mat_t after "; for (int l=0;l<9;l++)
        // std::cout<<mat_t[l]<<" ";std::cout<<"\n";
        if (kappa1(3) > kappa2) {
          kappa2 = kappa1(3);
        }
      }
    }
  }
}

// ! exact solution
//  ! We want the time-term to be the same in u1,u2 and u3, so that we don't
//  need to call this function ! in every time step. We just need to scale the
//  initial condition.
void exact_solution(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
                    float_sw4 &u1, float_sw4 &u2, float_sw4 &u3, int flag) {
  u1 = cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  u2 = sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  u3 = sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
}

void interface_block(Farray &Rop, Farray &ghcof, Farray &Sb, Sarray &rho_c,
                     Sarray &lambda_c, Sarray &rho_f, Sarray &Jacobian_c,
                     Sarray &mu_c, Farray &P, Sarray &XI13_c, Sarray &XI23_c,
                     Sarray &XI33_c, Farray &Mass_block, PackArgs &a) {
  float_sw4 int_cof;
  int i, j, k, l;
  auto h3_c = a.h3_c;
  auto h3_f = a.h3_f;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  //
  int_cof = 17.0 / 48.0 * h3_f * ghcof(1) / pow(h3_c, 2);
  //
  Mass_block = 0.0;
  for (l = 1; l <= n2_c; l++) {
    for (k = 1; k <= n1_c; k++) {
      //
      for (j = -4; j <= 2; j += 2) {
        for (i = -4; i <= 2; i += 2) {
          // first set equation w.r.t the first component
          Mass_block(1, 1, k, l) =
              Mass_block(1, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                       pow(XI13_c(k, l, n3_c), 2) +
                   mu_c(k, l, n3_c) * (pow(XI23_c(k, l, n3_c), 2) +
                                       pow(XI33_c(k, l, n3_c), 2))) /
                  rho_c(k, l, n3_c) * int_cof;
          // first set equation w.r.t the second component
          Mass_block(1, 2, k, l) =
              Mass_block(1, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI13_c(k, l, n3_c) * XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // first set equation w.r.t the third component
          Mass_block(1, 3, k, l) =
              Mass_block(1, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI13_c(k, l, n3_c) * XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // second set equation w.r.t the first component
          Mass_block(2, 1, k, l) =
              Mass_block(2, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI13_c(k, l, n3_c) * XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // second set equation w.r.t the second component
          Mass_block(2, 2, k, l) =
              Mass_block(2, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                       pow(XI23_c(k, l, n3_c), 2) +
                   mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                       pow(XI33_c(k, l, n3_c), 2))) /
                  rho_c(k, l, n3_c) * int_cof;
          // second set equation w.r.t the third component
          Mass_block(2, 3, k, l) =
              Mass_block(2, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI23_c(k, l, n3_c) * XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // third set equation w.r.t the first component
          Mass_block(3, 1, k, l) =
              Mass_block(3, 1, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI13_c(k, l, n3_c) * XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // third set equation w.r.t the second component
          Mass_block(3, 2, k, l) =
              Mass_block(3, 2, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                  XI23_c(k, l, n3_c) * XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) *
                  int_cof;
          // third set equation w.r.t the third component
          Mass_block(3, 3, k, l) =
              Mass_block(3, 3, k, l) +
              Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(-j / 2) *
                  P(-i / 2) * Jacobian_c(k, l, n3_c) *
                  ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                       pow(XI33_c(k, l, n3_c), 2) +
                   mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                       pow(XI23_c(k, l, n3_c), 2))) /
                  rho_c(k, l, n3_c) * int_cof;
        }
      }
      //
      for (j = -4; j <= 2; j += 2) {
        // first set equation w.r.t the first component
        Mass_block(1, 1, k, l) =
            Mass_block(1, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI13_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI23_c(k, l, n3_c), 2) +
                                     pow(XI33_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
        // first set equation w.r.t the second component
        Mass_block(1, 2, k, l) =
            Mass_block(1, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // first set equation w.r.t the third component
        Mass_block(1, 3, k, l) =
            Mass_block(1, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the first component
        Mass_block(2, 1, k, l) =
            Mass_block(2, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the second component
        Mass_block(2, 2, k, l) =
            Mass_block(2, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI23_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                     pow(XI33_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the third component
        Mass_block(2, 3, k, l) =
            Mass_block(2, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the first component
        Mass_block(3, 1, k, l) =
            Mass_block(3, 1, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the second component
        Mass_block(3, 2, k, l) =
            Mass_block(3, 2, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the third component
        Mass_block(3, 3, k, l) =
            Mass_block(3, 3, k, l) +
            Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(-j / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI33_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                     pow(XI23_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
      }
      //
      for (i = -4; i <= 2; i += 2) {
        // first set equation w.r.t the first component
        Mass_block(1, 1, k, l) =
            Mass_block(1, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI13_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI23_c(k, l, n3_c), 2) +
                                     pow(XI33_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
        // first set equation w.r.t the second component
        Mass_block(1, 2, k, l) =
            Mass_block(1, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // first set equation w.r.t the third component
        Mass_block(1, 3, k, l) =
            Mass_block(1, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the first component
        Mass_block(2, 1, k, l) =
            Mass_block(2, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the second component
        Mass_block(2, 2, k, l) =
            Mass_block(2, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI23_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                     pow(XI33_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
        // second set equation w.r.t the third component
        Mass_block(2, 3, k, l) =
            Mass_block(2, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the first component
        Mass_block(3, 1, k, l) =
            Mass_block(3, 1, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the second component
        Mass_block(3, 2, k, l) =
            Mass_block(3, 2, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
                XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
        // third set equation w.r.t the third component
        Mass_block(3, 3, k, l) =
            Mass_block(3, 3, k, l) +
            Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(-i / 2) *
                Jacobian_c(k, l, n3_c) *
                ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                     pow(XI33_c(k, l, n3_c), 2) +
                 mu_c(k, l, n3_c) * (pow(XI13_c(k, l, n3_c), 2) +
                                     pow(XI23_c(k, l, n3_c), 2))) /
                rho_c(k, l, n3_c) * int_cof;
      }
      //
      // first set equation w.r.t the first component
      Mass_block(1, 1, k, l) =
          Mass_block(1, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI13_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI23_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * int_cof;
      // first set equation w.r.t the second component
      Mass_block(1, 2, k, l) =
          Mass_block(1, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // first set equation w.r.t the third component
      Mass_block(1, 3, k, l) =
          Mass_block(1, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // second set equation w.r.t the first component
      Mass_block(2, 1, k, l) =
          Mass_block(2, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // second set equation w.r.t the second component
      Mass_block(2, 2, k, l) =
          Mass_block(2, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI23_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * int_cof;
      // second set equation w.r.t the third component
      Mass_block(2, 3, k, l) =
          Mass_block(2, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // third set equation w.r.t the first component
      Mass_block(3, 1, k, l) =
          Mass_block(3, 1, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // third set equation w.r.t the second component
      Mass_block(3, 2, k, l) =
          Mass_block(3, 2, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * int_cof;
      // third set equation w.r.t the third component
      Mass_block(3, 3, k, l) =
          Mass_block(3, 3, k, l) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI33_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI23_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * int_cof;
      // from the norm derivative
      // first set equation w.r.t the first component
      Mass_block(1, 1, k, l) =
          Mass_block(1, 1, k, l) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI13_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI23_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              h3_c;
      // first set equation w.r.t the second component
      Mass_block(1, 2, k, l) = Mass_block(1, 2, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI13_c(k, l, n3_c) * XI23_c(k, l, n3_c) /
                                   h3_c;
      // first set equation w.r.t the third component
      Mass_block(1, 3, k, l) = Mass_block(1, 3, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI13_c(k, l, n3_c) * XI33_c(k, l, n3_c) /
                                   h3_c;
      // second set equation w.r.t the first component
      Mass_block(2, 1, k, l) = Mass_block(2, 1, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI13_c(k, l, n3_c) * XI23_c(k, l, n3_c) /
                                   h3_c;
      // second set equation w.r.t the second component
      Mass_block(2, 2, k, l) =
          Mass_block(2, 2, k, l) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI23_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              h3_c;
      // second set equation w.r.t the third component
      Mass_block(2, 3, k, l) = Mass_block(2, 3, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI23_c(k, l, n3_c) * XI33_c(k, l, n3_c) /
                                   h3_c;
      // third set equation w.r.t the first component
      Mass_block(3, 1, k, l) = Mass_block(3, 1, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI13_c(k, l, n3_c) * XI33_c(k, l, n3_c) /
                                   h3_c;
      // third set equation w.r.t the second component
      Mass_block(3, 2, k, l) = Mass_block(3, 2, k, l) -
                               Sb(0) * Jacobian_c(k, l, n3_c) *
                                   (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) *
                                   XI23_c(k, l, n3_c) * XI33_c(k, l, n3_c) /
                                   h3_c;
      // third set equation w.r.t the third component
      Mass_block(3, 3, k, l) =
          Mass_block(3, 3, k, l) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI33_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI23_c(k, l, n3_c), 2))) /
              h3_c;
    }
  }
}

void interface_lhs(Farray &LHS, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                   Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                   Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                   Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                   Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                   Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                   Farray &acof, Farray &bof, PackArgs &a, int index) {
  float_sw4 int_cof;

  int i, j, k, i1, j1, k1, l;

  auto dim = a.dim;
  auto h3_c = a.h3_c;
  auto h3_f = a.h3_f;

  auto nrg = a.nrg;
  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;
  Farray u_temp(1 - nrg, n1_c + nrg, 1 - nrg, n2_c + nrg, 1, dim);
  u_temp = 0.0;
  for (j = 1; j <= n2_c; j++) {
    for (i = 1; i <= n1_c; i++) {
      u_temp(i, j, 1) = u_c(1, i, j, n3_c + 1);
      u_temp(i, j, 2) = u_c(2, i, j, n3_c + 1);
      u_temp(i, j, 3) = u_c(3, i, j, n3_c + 1);
    }
  }

  int_cof = 17.0 / 48.0 * h3_f * ghcof(1) / pow(h3_c, 2);

  LHS = 0.0;
  for (l = 1; l <= n2_c; l++) {
    for (k = 1; k <= n1_c; k++) {
      for (j = -4; j <= 2; j += 2) {
        for (i = -4; i <= 2; i += 2) {
          for (j1 = -1; j1 <= 2; j1++) {
            for (i1 = -1; i1 <= 2; i1++) {
              // first set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                           pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                           (pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2) +
                            pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
              // second set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                           pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                           (pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2) +
                            pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
              // third set equation
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
                  LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 1) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      (lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                      XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 2) * int_cof +
                  Rop(j) * Rop(i) * rho_f(2 * k + i, 2 * l + j, 1) * P(j1) *
                      P(i1) * Jacobian_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      ((2.0 * mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) +
                        lambda_c(k + i / 2 + i1, l + j / 2 + j1, n3_c)) *
                           pow(XI33_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                               2) +
                       mu_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                           (pow(XI13_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2) +
                            pow(XI23_c(k + i / 2 + i1, l + j / 2 + j1, n3_c),
                                2))) /
                      rho_c(k + i / 2 + i1, l + j / 2 + j1, n3_c) *
                      u_temp(k + i / 2 + i1, l + j / 2 + j1, 3) * int_cof;
            }
          }
        }
      }
      //
      for (j = -4; j <= 2; j += 2) {
        for (j1 = -1; j1 <= 2; j1++) {
          // first set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, n3_c) +
                    lambda_c(k, l + j / 2 + j1, n3_c)) *
                       pow(XI13_c(k, l + j / 2 + j1, n3_c), 2) +
                   mu_c(k, l + j / 2 + j1, n3_c) *
                       (pow(XI23_c(k, l + j / 2 + j1, n3_c), 2) +
                        pow(XI33_c(k, l + j / 2 + j1, n3_c), 2))) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI13_c(k, l + j / 2 + j1, n3_c) *
                  XI23_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI13_c(k, l + j / 2 + j1, n3_c) *
                  XI33_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
          // second set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI13_c(k, l + j / 2 + j1, n3_c) *
                  XI23_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, n3_c) +
                    lambda_c(k, l + j / 2 + j1, n3_c)) *
                       pow(XI23_c(k, l + j / 2 + j1, n3_c), 2) +
                   mu_c(k, l + j / 2 + j1, n3_c) *
                       (pow(XI13_c(k, l + j / 2 + j1, n3_c), 2) +
                        pow(XI33_c(k, l + j / 2 + j1, n3_c), 2))) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI23_c(k, l + j / 2 + j1, n3_c) *
                  XI33_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
          // third set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI13_c(k, l + j / 2 + j1, n3_c) *
                  XI33_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 1) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  (lambda_c(k, l + j / 2 + j1, n3_c) +
                   mu_c(k, l + j / 2 + j1, n3_c)) *
                  XI23_c(k, l + j / 2 + j1, n3_c) *
                  XI33_c(k, l + j / 2 + j1, n3_c) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 2) * int_cof +
              Rop(j) * Rop(-1) * rho_f(2 * k - 1, 2 * l + j, 1) * P(j1) *
                  Jacobian_c(k, l + j / 2 + j1, n3_c) *
                  ((2.0 * mu_c(k, l + j / 2 + j1, n3_c) +
                    lambda_c(k, l + j / 2 + j1, n3_c)) *
                       pow(XI33_c(k, l + j / 2 + j1, n3_c), 2) +
                   mu_c(k, l + j / 2 + j1, n3_c) *
                       (pow(XI13_c(k, l + j / 2 + j1, n3_c), 2) +
                        pow(XI23_c(k, l + j / 2 + j1, n3_c), 2))) /
                  rho_c(k, l + j / 2 + j1, n3_c) *
                  u_temp(k, l + j / 2 + j1, 3) * int_cof;
        }
      }
      //
      for (i = -4; i <= 2; i += 2) {
        for (i1 = -1; i1 <= 2; i1++) {
          // first set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, n3_c) +
                    lambda_c(k + i / 2 + i1, l, n3_c)) *
                       pow(XI13_c(k + i / 2 + i1, l, n3_c), 2) +
                   mu_c(k + i / 2 + i1, l, n3_c) *
                       (pow(XI23_c(k + i / 2 + i1, l, n3_c), 2) +
                        pow(XI33_c(k + i / 2 + i1, l, n3_c), 2))) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI13_c(k + i / 2 + i1, l, n3_c) *
                  XI23_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI13_c(k + i / 2 + i1, l, n3_c) *
                  XI33_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
          // second set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI13_c(k + i / 2 + i1, l, n3_c) *
                  XI23_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, n3_c) +
                    lambda_c(k + i / 2 + i1, l, n3_c)) *
                       pow(XI23_c(k + i / 2 + i1, l, n3_c), 2) +
                   mu_c(k + i / 2 + i1, l, n3_c) *
                       (pow(XI13_c(k + i / 2 + i1, l, n3_c), 2) +
                        pow(XI33_c(k + i / 2 + i1, l, n3_c), 2))) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI23_c(k + i / 2 + i1, l, n3_c) *
                  XI33_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
          // third set equation
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
              LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI13_c(k + i / 2 + i1, l, n3_c) *
                  XI33_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 1) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  (lambda_c(k + i / 2 + i1, l, n3_c) +
                   mu_c(k + i / 2 + i1, l, n3_c)) *
                  XI23_c(k + i / 2 + i1, l, n3_c) *
                  XI33_c(k + i / 2 + i1, l, n3_c) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 2) * int_cof +
              Rop(-1) * Rop(i) * rho_f(2 * k + i, 2 * l - 1, 1) * P(i1) *
                  Jacobian_c(k + i / 2 + i1, l, n3_c) *
                  ((2.0 * mu_c(k + i / 2 + i1, l, n3_c) +
                    lambda_c(k + i / 2 + i1, l, n3_c)) *
                       pow(XI33_c(k + i / 2 + i1, l, n3_c), 2) +
                   mu_c(k + i / 2 + i1, l, n3_c) *
                       (pow(XI13_c(k + i / 2 + i1, l, n3_c), 2) +
                        pow(XI23_c(k + i / 2 + i1, l, n3_c), 2))) /
                  rho_c(k + i / 2 + i1, l, n3_c) *
                  u_temp(k + i / 2 + i1, l, 3) * int_cof;
        }
      }
      //
      // first set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI13_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI23_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * u_temp(k, l, 1) * int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 2) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 3) *
              int_cof;
      // second set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 1) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI23_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * u_temp(k, l, 2) * int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 3) *
              int_cof;
      // third set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 1) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / rho_c(k, l, n3_c) * u_temp(k, l, 2) *
              int_cof +
          Rop(-1) * Rop(-1) * rho_f(2 * k - 1, 2 * l - 1, 1) *
              Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI33_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI23_c(k, l, n3_c), 2))) /
              rho_c(k, l, n3_c) * u_temp(k, l, 3) * int_cof;
      //
      // first set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 1) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI13_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI23_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / h3_c * u_temp(k, l, 3);
      // second set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 2) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI23_c(k, l, n3_c) / h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI23_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI33_c(k, l, n3_c), 2))) /
              h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / h3_c * u_temp(k, l, 3);
      // third set equation
      LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) =
          LHS((l - 1) * 3 * n1_c + (k - 1) * 3 + 3) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI13_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / h3_c * u_temp(k, l, 1) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              (lambda_c(k, l, n3_c) + mu_c(k, l, n3_c)) * XI23_c(k, l, n3_c) *
              XI33_c(k, l, n3_c) / h3_c * u_temp(k, l, 2) -
          Sb(0) * Jacobian_c(k, l, n3_c) *
              ((2.0 * mu_c(k, l, n3_c) + lambda_c(k, l, n3_c)) *
                   pow(XI33_c(k, l, n3_c), 2) +
               mu_c(k, l, n3_c) *
                   (pow(XI13_c(k, l, n3_c), 2) + pow(XI23_c(k, l, n3_c), 2))) /
              h3_c * u_temp(k, l, 3);
    }
  }
}  // END INTERFACE_LHS

void injection(Sarray &u_f, Sarray &u_c, Farray &P, PackArgs &a, int index) {
  int i, j, k;
  // Injection at the interface

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;
  auto nrg = a.nrg;
  for (int l = 1; l <= a.dim; l++) {
    for (j = 1; j <= n2_c; j++) {
      for (i = 1; i <= n1_c; i++) {
        u_f(l, 2 * i - 1, 2 * j - 1, 1) = u_c(l, i, j, n3_c);
        u_f(l, 2 * i, 2 * j - 1, 1) =
            P(-1) * u_c(l, i - 1, j, n3_c) + P(0) * u_c(l, i, j, n3_c) +
            P(1) * u_c(l, i + 1, j, n3_c) + P(2) * u_c(l, i + 2, j, n3_c);
        u_f(l, 2 * i - 1, 2 * j, 1) =
            P(-1) * u_c(l, i, j - 1, n3_c) + P(0) * u_c(l, i, j, n3_c) +
            P(1) * u_c(l, i, j + 1, n3_c) + P(2) * u_c(l, i, j + 2, n3_c);
        u_f(l, 2 * i, 2 * j, 1) =
            P(-1) * (P(-1) * u_c(l, i - 1, j - 1, n3_c) +
                     P(0) * u_c(l, i, j - 1, n3_c) +
                     P(1) * u_c(l, i + 1, j - 1, n3_c) +
                     P(2) * u_c(l, i + 2, j - 1, n3_c)) +
            P(0) * (P(-1) * u_c(l, i - 1, j, n3_c) + P(0) * u_c(l, i, j, n3_c) +
                    P(1) * u_c(l, i + 1, j, n3_c) +
                    P(2) * u_c(l, i + 2, j, n3_c)) +
            P(1) * (P(-1) * u_c(l, i - 1, j + 1, n3_c) +
                    P(0) * u_c(l, i, j + 1, n3_c) +
                    P(1) * u_c(l, i + 1, j + 1, n3_c) +
                    P(2) * u_c(l, i + 2, j + 1, n3_c)) +
            P(2) * (P(-1) * u_c(l, i - 1, j + 2, n3_c) +
                    P(0) * u_c(l, i, j + 2, n3_c) +
                    P(1) * u_c(l, i + 1, j + 2, n3_c) +
                    P(2) * u_c(l, i + 2, j + 2, n3_c));
      }
    }
  }
}

void interface_rhs(Farray &Vass, Farray &lh_c, Farray &lh_f, Sarray &Jacobian_c,
                   Sarray &Jacobian_f, Sarray &mu_c, Sarray &mu_f,
                   Sarray &lambda_c, Sarray &lambda_f, Sarray &rho_c,
                   Sarray &rho_f, Sarray &XI13_c, Sarray &XI23_c,
                   Sarray &XI33_c, Sarray &XI13_f, Sarray &XI23_f,
                   Sarray &XI33_f, Farray &P, Farray &Sb, Farray &Rop,
                   Farray &sbop_no_gp, Farray &acof_no_gp, Sarray &u_c,
                   Sarray &u_f, Farray &Mass_f1, Farray &ux_cof, Farray &ghcof,
                   Farray &acof, Farray &bof, PackArgs &a, int index) {
  int i, j, k, k1, l, m;
  //
  Vass = 0.0;
  lh_c = 0.0;
  lh_f = 0.0;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  float_sw4 l3 = a.l3;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n1_f = a.n1_f;
  auto n2_f = a.n2_f;
  auto n3_f = a.n3_f;

  auto h1_c = a.h1_c;
  auto h2_c = a.h2_c;
  auto h3_c = a.h3_c;

  auto h1_f = a.h1_f;
  auto h2_f = a.h2_f;
  auto h3_f = a.h3_f;

  auto nrg = a.nrg;

  //
  // term 1
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = 1; j <= 4; j++) {
        // 33
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) +
            Jacobian_c(i, k, n3_c) *
                ((2.0 * mu_c(i, k, n3_c) + lambda_c(i, k, n3_c)) *
                     pow(XI13_c(i, k, n3_c), 2) +
                 mu_c(i, k, n3_c) * (pow(XI23_c(i, k, n3_c), 2) +
                                     pow(XI33_c(i, k, n3_c), 2))) *
                Sb(j) * u_c(1, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI13_c(i, k, n3_c) * XI23_c(i, k, n3_c) * Sb(j) *
                u_c(2, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI13_c(i, k, n3_c) * XI33_c(i, k, n3_c) * Sb(j) *
                u_c(3, i, k, n3_c + 1 - j) / h3_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI13_c(i, k, n3_c) * XI23_c(i, k, n3_c) * Sb(j) *
                u_c(1, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) *
                ((2.0 * mu_c(i, k, n3_c) + lambda_c(i, k, n3_c)) *
                     pow(XI23_c(i, k, n3_c), 2) +
                 mu_c(i, k, n3_c) * (pow(XI13_c(i, k, n3_c), 2) +
                                     pow(XI33_c(i, k, n3_c), 2))) *
                Sb(j) * u_c(2, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI23_c(i, k, n3_c) * XI33_c(i, k, n3_c) * Sb(j) *
                u_c(3, i, k, n3_c + 1 - j) / h3_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI13_c(i, k, n3_c) * XI33_c(i, k, n3_c) * Sb(j) *
                u_c(1, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) * (lambda_c(i, k, n3_c) + mu_c(i, k, n3_c)) *
                XI23_c(i, k, n3_c) * XI33_c(i, k, n3_c) * Sb(j) *
                u_c(2, i, k, n3_c + 1 - j) / h3_c +
            Jacobian_c(i, k, n3_c) *
                ((2.0 * mu_c(i, k, n3_c) + lambda_c(i, k, n3_c)) *
                     pow(XI33_c(i, k, n3_c), 2) +
                 mu_c(i, k, n3_c) * (pow(XI13_c(i, k, n3_c), 2) +
                                     pow(XI23_c(i, k, n3_c), 2))) *
                Sb(j) * u_c(3, i, k, n3_c + 1 - j) / h3_c;
      }
    }
  }

  //
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -2; j <= 2; j++) {
        // 31
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
            Jacobian_c(i, k, n3_c) *
                (2.0 * mu_c(i, k, n3_c) + lambda_c(i, k, n3_c)) / l1 *
                XI13_c(i, k, n3_c) * ux_cof(j) * u_c(1, i + j, k, n3_c) / h1_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l1 *
                XI23_c(i, k, n3_c) * ux_cof(j) * u_c(2, i + j, k, n3_c) / h1_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l1 *
                XI33_c(i, k, n3_c) * ux_cof(j) * u_c(3, i + j, k, n3_c) / h1_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
            Jacobian_c(i, k, n3_c) * lambda_c(i, k, n3_c) / l1 *
                XI23_c(i, k, n3_c) * ux_cof(j) * u_c(1, i + j, k, n3_c) / h1_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l1 *
                XI13_c(i, k, n3_c) * ux_cof(j) * u_c(2, i + j, k, n3_c) / h1_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
            Jacobian_c(i, k, n3_c) * lambda_c(i, k, n3_c) / l1 *
                XI33_c(i, k, n3_c) * ux_cof(j) * u_c(1, i + j, k, n3_c) / h1_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l1 *
                XI13_c(i, k, n3_c) * ux_cof(j) * u_c(3, i + j, k, n3_c) / h1_c;
        // 32
        // first set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l2 *
                XI23_c(i, k, n3_c) * ux_cof(j) * u_c(1, i, k + j, n3_c) / h2_c -
            Jacobian_c(i, k, n3_c) * lambda_c(i, k, n3_c) / l2 *
                XI13_c(i, k, n3_c) * ux_cof(j) * u_c(2, i, k + j, n3_c) / h2_c;
        // second set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l2 *
                XI13_c(i, k, n3_c) * ux_cof(j) * u_c(1, i, k + j, n3_c) / h2_c -
            Jacobian_c(i, k, n3_c) *
                (2.0 * mu_c(i, k, n3_c) + lambda_c(i, k, n3_c)) / l2 *
                XI23_c(i, k, n3_c) * ux_cof(j) * u_c(2, i, k + j, n3_c) / h2_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l2 *
                XI33_c(i, k, n3_c) * ux_cof(j) * u_c(3, i, k + j, n3_c) / h2_c;
        // third set equation
        Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
            Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
            Jacobian_c(i, k, n3_c) * lambda_c(i, k, n3_c) / l2 *
                XI33_c(i, k, n3_c) * ux_cof(j) * u_c(2, i, k + j, n3_c) / h2_c -
            Jacobian_c(i, k, n3_c) * mu_c(i, k, n3_c) / l2 *
                XI23_c(i, k, n3_c) * ux_cof(j) * u_c(3, i, k + j, n3_c) / h2_c;
      }
    }
  }

  // term 2
  // interior
  for (j = -2; j <= n2_c + 3; j++) {
    for (i = -2; i <= n1_c + 3; i++) {
      // second derivative 11 & 22 & 12 & 21
      // first set
      lh_c(i, j, n3_c, 1) =
          lh_c(i, j, n3_c, 1) +
          ((-Jacobian_c(i - 2, j, n3_c) *
                (2.0 * mu_c(i - 2, j, n3_c) + lambda_c(i - 2, j, n3_c)) / 8.0 +
            Jacobian_c(i - 1, j, n3_c) *
                (2.0 * mu_c(i - 1, j, n3_c) + lambda_c(i - 1, j, n3_c)) / 6.0 -
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 8.0) *
               u_c(1, i - 2, j, n3_c) +
           (Jacobian_c(i - 2, j, n3_c) *
                (2.0 * mu_c(i - 2, j, n3_c) + lambda_c(i - 2, j, n3_c)) / 6.0 +
            Jacobian_c(i - 1, j, n3_c) *
                (2.0 * mu_c(i - 1, j, n3_c) + lambda_c(i - 1, j, n3_c)) / 2.0 +
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) *
                (2.0 * mu_c(i + 1, j, n3_c) + lambda_c(i + 1, j, n3_c)) / 6.0) *
               u_c(1, i - 1, j, n3_c) +
           (-Jacobian_c(i - 2, j, n3_c) *
                (2.0 * mu_c(i - 2, j, n3_c) + lambda_c(i - 2, j, n3_c)) / 24.0 -
            Jacobian_c(i - 1, j, n3_c) *
                (2.0 * mu_c(i - 1, j, n3_c) + lambda_c(i - 1, j, n3_c)) * 5.0 /
                6.0 -
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, n3_c) *
                (2.0 * mu_c(i + 1, j, n3_c) + lambda_c(i + 1, j, n3_c)) * 5.0 /
                6.0 -
            Jacobian_c(i + 2, j, n3_c) *
                (2.0 * mu_c(i + 2, j, n3_c) + lambda_c(i + 2, j, n3_c)) /
                24.0) *
               u_c(1, i - 0, j, n3_c) +
           (Jacobian_c(i - 1, j, n3_c) *
                (2.0 * mu_c(i - 1, j, n3_c) + lambda_c(i - 1, j, n3_c)) / 6.0 +
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) *
                (2.0 * mu_c(i + 1, j, n3_c) + lambda_c(i + 1, j, n3_c)) / 2.0 +
            Jacobian_c(i + 2, j, n3_c) *
                (2.0 * mu_c(i + 2, j, n3_c) + lambda_c(i + 2, j, n3_c)) / 6.0) *
               u_c(1, i + 1, j, n3_c) +
           (-Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 8.0 +
            Jacobian_c(i + 1, j, n3_c) *
                (2.0 * mu_c(i + 1, j, n3_c) + lambda_c(i + 1, j, n3_c)) / 6.0 -
            Jacobian_c(i + 2, j, n3_c) *
                (2.0 * mu_c(i + 2, j, n3_c) + lambda_c(i + 2, j, n3_c)) / 8.0) *
               u_c(1, i + 2, j, n3_c)) /
              (h1_c * h1_c) / (l1 * l1) +
          ((-Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 8.0 +
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0) *
               u_c(1, i, j - 2, n3_c) +
           (Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 6.0 +
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 2.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 6.0) *
               u_c(1, i, j - 1, n3_c) +
           (-Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 24.0 -
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 24.0) *
               u_c(1, i, j - 0, n3_c) +
           (Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 6.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 2.0 +
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 6.0) *
               u_c(1, i, j + 1, n3_c) +
           (-Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 6.0 -
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 8.0) *
               u_c(1, i, j + 2, n3_c)) /
              (h2_c * h2_c) / (l2 * l2) +
          (Jacobian_c(i - 2, j, n3_c) * lambda_c(i - 2, j, n3_c) *
               (u_c(2, i - 2, j - 2, n3_c) / 12.0 -
                u_c(2, i - 2, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i - 2, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i - 2, j + 2, n3_c) / 12.0) /
               12.0 -
           Jacobian_c(i - 1, j, n3_c) * lambda_c(i - 1, j, n3_c) *
               (u_c(2, i - 1, j - 2, n3_c) / 12.0 -
                u_c(2, i - 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i - 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i - 1, j + 2, n3_c) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i + 1, j, n3_c) * lambda_c(i + 1, j, n3_c) *
               (u_c(2, i + 1, j - 2, n3_c) / 12.0 -
                u_c(2, i + 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 1, j + 2, n3_c) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i + 2, j, n3_c) * lambda_c(i + 2, j, n3_c) *
               (u_c(2, i + 2, j - 2, n3_c) / 12.0 -
                u_c(2, i + 2, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 2, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 2, n3_c) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c +
          (Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) *
               (u_c(2, i - 2, j - 2, n3_c) / 12.0 -
                u_c(2, i - 1, j - 2, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 1, j - 2, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 2, j - 2, n3_c) / 12.0) /
               12.0 -
           Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) *
               (u_c(2, i - 2, j - 1, n3_c) / 12.0 -
                u_c(2, i - 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 1, j - 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 2, j - 1, n3_c) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) *
               (u_c(2, i - 2, j + 1, n3_c) / 12.0 -
                u_c(2, i - 1, j + 1, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 1, n3_c) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) *
               (u_c(2, i - 2, j + 2, n3_c) / 12.0 -
                u_c(2, i - 1, j + 2, n3_c) * 2.0 / 3.0 +
                u_c(2, i + 1, j + 2, n3_c) * 2.0 / 3.0 -
                u_c(2, i + 2, j + 2, n3_c) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c;
      //
      // second set
      lh_c(i, j, n3_c, 2) =
          lh_c(i, j, n3_c, 2) +
          ((-Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 8.0 +
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0) *
               u_c(2, i - 2, j, n3_c) +
           (Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 6.0 +
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 2.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 6.0) *
               u_c(2, i - 1, j, n3_c) +
           (-Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 24.0 -
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 24.0) *
               u_c(2, i - 0, j, n3_c) +
           (Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 6.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 2.0 +
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 6.0) *
               u_c(2, i + 1, j, n3_c) +
           (-Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 6.0 -
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 8.0) *
               u_c(2, i + 2, j, n3_c)) /
              pow(h1_c, 2) / pow(l1, 2) +
          ((-Jacobian_c(i, j - 2, n3_c) *
                (2.0 * mu_c(i, j - 2, n3_c) + lambda_c(i, j - 2, n3_c)) / 8.0 +
            Jacobian_c(i, j - 1, n3_c) *
                (2.0 * mu_c(i, j - 1, n3_c) + lambda_c(i, j - 1, n3_c)) / 6.0 -
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 8.0) *
               u_c(2, i, j - 2, n3_c) +
           (Jacobian_c(i, j - 2, n3_c) *
                (2.0 * mu_c(i, j - 2, n3_c) + lambda_c(i, j - 2, n3_c)) / 6.0 +
            Jacobian_c(i, j - 1, n3_c) *
                (2.0 * mu_c(i, j - 1, n3_c) + lambda_c(i, j - 1, n3_c)) / 2.0 +
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) *
                (2.0 * mu_c(i, j + 1, n3_c) + lambda_c(i, j + 1, n3_c)) / 6.0) *
               u_c(2, i, j - 1, n3_c) +
           (-Jacobian_c(i, j - 2, n3_c) *
                (2.0 * mu_c(i, j - 2, n3_c) + lambda_c(i, j - 2, n3_c)) / 24.0 -
            Jacobian_c(i, j - 1, n3_c) *
                (2.0 * mu_c(i, j - 1, n3_c) + lambda_c(i, j - 1, n3_c)) * 5.0 /
                6.0 -
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, n3_c) *
                (2.0 * mu_c(i, j + 1, n3_c) + lambda_c(i, j + 1, n3_c)) * 5.0 /
                6.0 -
            Jacobian_c(i, j + 2, n3_c) *
                (2.0 * mu_c(i, j + 2, n3_c) + lambda_c(i, j + 2, n3_c)) /
                24.0) *
               u_c(2, i, j - 0, n3_c) +
           (Jacobian_c(i, j - 1, n3_c) *
                (2.0 * mu_c(i, j - 1, n3_c) + lambda_c(i, j - 1, n3_c)) / 6.0 +
            Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) *
                (2.0 * mu_c(i, j + 1, n3_c) + lambda_c(i, j + 1, n3_c)) / 2.0 +
            Jacobian_c(i, j + 2, n3_c) *
                (2.0 * mu_c(i, j + 2, n3_c) + lambda_c(i, j + 2, n3_c)) / 6.0) *
               u_c(2, i, j + 1, n3_c) +
           (-Jacobian_c(i, j, n3_c) *
                (2.0 * mu_c(i, j, n3_c) + lambda_c(i, j, n3_c)) / 8.0 +
            Jacobian_c(i, j + 1, n3_c) *
                (2.0 * mu_c(i, j + 1, n3_c) + lambda_c(i, j + 1, n3_c)) / 6.0 -
            Jacobian_c(i, j + 2, n3_c) *
                (2.0 * mu_c(i, j + 2, n3_c) + lambda_c(i, j + 2, n3_c)) / 8.0) *
               u_c(2, i, j + 2, n3_c)) /
              pow(h2_c, 2) / pow(l2, 2) +
          (Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) *
               (u_c(1, i - 2, j - 2, n3_c) / 12.0 -
                u_c(1, i - 2, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i - 2, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i - 2, j + 2, n3_c) / 12.0) /
               12.0 -
           Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) *
               (u_c(1, i - 1, j - 2, n3_c) / 12.0 -
                u_c(1, i - 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i - 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i - 1, j + 2, n3_c) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) *
               (u_c(1, i + 1, j - 2, n3_c) / 12.0 -
                u_c(1, i + 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 1, j + 2, n3_c) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) *
               (u_c(1, i + 2, j - 2, n3_c) / 12.0 -
                u_c(1, i + 2, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 2, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 2, n3_c) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c +
          (Jacobian_c(i, j - 2, n3_c) * lambda_c(i, j - 2, n3_c) *
               (u_c(1, i - 2, j - 2, n3_c) / 12.0 -
                u_c(1, i - 1, j - 2, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 1, j - 2, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 2, j - 2, n3_c) / 12.0) /
               12.0 -
           Jacobian_c(i, j - 1, n3_c) * lambda_c(i, j - 1, n3_c) *
               (u_c(1, i - 2, j - 1, n3_c) / 12.0 -
                u_c(1, i - 1, j - 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 1, j - 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 2, j - 1, n3_c) / 12.0) *
               2.0 / 3.0 +
           Jacobian_c(i, j + 1, n3_c) * lambda_c(i, j + 1, n3_c) *
               (u_c(1, i - 2, j + 1, n3_c) / 12.0 -
                u_c(1, i - 1, j + 1, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 1, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 1, n3_c) / 12.0) *
               2.0 / 3.0 -
           Jacobian_c(i, j + 2, n3_c) * lambda_c(i, j + 2, n3_c) *
               (u_c(1, i - 2, j + 2, n3_c) / 12.0 -
                u_c(1, i - 1, j + 2, n3_c) * 2.0 / 3.0 +
                u_c(1, i + 1, j + 2, n3_c) * 2.0 / 3.0 -
                u_c(1, i + 2, j + 2, n3_c) / 12.0) /
               12.0) /
              l1 / l2 / h1_c / h2_c;
      // third set
      lh_c(i, j, n3_c, 3) =
          lh_c(i, j, n3_c, 3) +
          ((-Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 8.0 +
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0) *
               u_c(3, i - 2, j, n3_c) +
           (Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 6.0 +
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 2.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 6.0) *
               u_c(3, i - 1, j, n3_c) +
           (-Jacobian_c(i - 2, j, n3_c) * mu_c(i - 2, j, n3_c) / 24.0 -
            Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) * 3.0 / 4.0 -
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 24.0) *
               u_c(3, i - 0, j, n3_c) +
           (Jacobian_c(i - 1, j, n3_c) * mu_c(i - 1, j, n3_c) / 6.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 2.0 +
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 6.0) *
               u_c(3, i + 1, j, n3_c) +
           (-Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0 +
            Jacobian_c(i + 1, j, n3_c) * mu_c(i + 1, j, n3_c) / 6.0 -
            Jacobian_c(i + 2, j, n3_c) * mu_c(i + 2, j, n3_c) / 8.0) *
               u_c(3, i + 2, j, n3_c)) /
              pow(h1_c, 2) / pow(l1, 2) +
          ((-Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 8.0 +
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0) *
               u_c(3, i, j - 2, n3_c) +
           (Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 6.0 +
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 2.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 6.0) *
               u_c(3, i, j - 1, n3_c) +
           (-Jacobian_c(i, j - 2, n3_c) * mu_c(i, j - 2, n3_c) / 24.0 -
            Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) * 3.0 / 4.0 -
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) * 5.0 / 6.0 -
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 24.0) *
               u_c(3, i, j - 0, n3_c) +
           (Jacobian_c(i, j - 1, n3_c) * mu_c(i, j - 1, n3_c) / 6.0 +
            Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 2.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 2.0 +
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 6.0) *
               u_c(3, i, j + 1, n3_c) +
           (-Jacobian_c(i, j, n3_c) * mu_c(i, j, n3_c) / 8.0 +
            Jacobian_c(i, j + 1, n3_c) * mu_c(i, j + 1, n3_c) / 6.0 -
            Jacobian_c(i, j + 2, n3_c) * mu_c(i, j + 2, n3_c) / 8.0) *
               u_c(3, i, j + 2, n3_c)) /
              pow(h2_c, 2) / pow(l2, 2);
    }
  }
  //
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = -2; k <= n1_c + 3; k++) {
      for (k1 = 1; k1 <= 8; k1++) {
        for (m = 1; m <= 8; m++) {
          // second derivative 33
          // first set equation
          lh_c(k, j, n3_c, 1) =
              lh_c(k, j, n3_c, 1) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                 lambda_c(k, j, n3_c + 1 - m)) *
                    pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                mu_c(k, j, n3_c + 1 - m) *
                    (pow(XI23_c(k, j, n3_c + 1 - m), 2) +
                     pow(XI33_c(k, j, n3_c + 1 - m), 2))) *
               u_c(1, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI13_c(k, j, n3_c + 1 - m) * XI23_c(k, j, n3_c + 1 - m) *
               u_c(2, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI13_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
               u_c(3, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2);
          // second set equation PROBLEM FIXED HERE
          lh_c(k, j, n3_c, 2) =
              lh_c(k, j, n3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI13_c(k, j, n3_c + 1 - m) * XI23_c(k, j, n3_c + 1 - m) *
               u_c(1, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                 lambda_c(k, j, n3_c + 1 - m)) *
                    pow(XI23_c(k, j, n3_c + 1 - m), 2) +
                mu_c(k, j, n3_c + 1 - m) *
                    (pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                     pow(XI33_c(k, j, n3_c + 1 - m), 2))) *
               u_c(2, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI23_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
               u_c(3, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2);
          // third set equation
          lh_c(k, j, n3_c, 3) =
              lh_c(k, j, n3_c, 3) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI13_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
               u_c(1, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
               XI23_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
               u_c(2, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2) +
              (acof(1, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
               ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                 lambda_c(k, j, n3_c + 1 - m)) *
                    pow(XI33_c(k, j, n3_c + 1 - m), 2) +
                mu_c(k, j, n3_c + 1 - m) *
                    (pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                     pow(XI23_c(k, j, n3_c + 1 - m), 2))) *
               u_c(3, k, j, n3_c + 1 - k1)) /
                  pow(h3_c, 2);
        }
      }
    }
  }
  // ghost points
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = -2; k <= 0; k++) {
      // first set equation
      lh_c(k, j, n3_c, 1) =
          lh_c(k, j, n3_c, 1) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI13_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI23_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, n3_c, 2) =
          lh_c(k, j, n3_c, 2) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI23_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, n3_c, 3) =
          lh_c(k, j, n3_c, 3) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI33_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI23_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2);
    }
  }

  //
  for (j = -2; j <= n2_c + 3; j++) {
    for (k = n1_c + 1; k <= n1_c + 3; k++) {
      // first set equation
      lh_c(k, j, n3_c, 1) =
          lh_c(k, j, n3_c, 1) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI13_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI23_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, n3_c, 2) =
          lh_c(k, j, n3_c, 2) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI23_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, n3_c, 3) =
          lh_c(k, j, n3_c, 3) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI33_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI23_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2);
    }
  }
  //

  for (j = -2; j <= 0; j++) {
    for (k = 1; k <= n1_c; k++) {
      // first set equation
      lh_c(k, j, n3_c, 1) =
          lh_c(k, j, n3_c, 1) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI13_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI23_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, n3_c, 2) =
          lh_c(k, j, n3_c, 2) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI23_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, n3_c, 3) =
          lh_c(k, j, n3_c, 3) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI33_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI23_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2);
    }
  }
  //
  for (j = n2_c + 1; j <= n2_c + 3; j++) {
    for (k = 1; k <= n1_c; k++) {
      // first set equation
      lh_c(k, j, n3_c, 1) =
          lh_c(k, j, n3_c, 1) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI13_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI23_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, n3_c, 2) =
          lh_c(k, j, n3_c, 2) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI23_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, n3_c, 3) =
          lh_c(k, j, n3_c, 3) +
          (u_c(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI33_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI23_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2);
    }
  }
  //
  for (i = -2; i <= n2_c + 3; i++) {
    for (j = -2; j <= n1_c + 3; j++) {
      for (k1 = 1; k1 <= 6; k1++) {
        // mixed derivative 13  23  31  32
        // first set equation
        lh_c(j, i, n3_c, 1) =
            lh_c(j, i, n3_c, 1) +
            (-Jacobian_c(j - 2, i, n3_c) *
                 (2.0 * mu_c(j - 2, i, n3_c) + lambda_c(j - 2, i, n3_c)) *
                 XI13_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) *
                 (2.0 * mu_c(j - 1, i, n3_c) + lambda_c(j - 1, i, n3_c)) *
                 XI13_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) *
                 (2.0 * mu_c(j + 1, i, n3_c) + lambda_c(j + 1, i, n3_c)) *
                 XI13_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) *
                 (2.0 * mu_c(j + 2, i, n3_c) + lambda_c(j + 2, i, n3_c)) *
                 XI13_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, n3_c) * lambda_c(j - 2, i, n3_c) *
                 XI23_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(2, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * lambda_c(j - 1, i, n3_c) *
                 XI23_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * lambda_c(j + 1, i, n3_c) *
                 XI23_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * lambda_c(j + 2, i, n3_c) *
                 XI23_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(2, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, n3_c) * lambda_c(j - 2, i, n3_c) *
                 XI33_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(3, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * lambda_c(j - 1, i, n3_c) *
                 XI33_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * lambda_c(j + 1, i, n3_c) *
                 XI33_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * lambda_c(j + 2, i, n3_c) *
                 XI33_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(3, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * mu_c(j, i - 2, n3_c) *
                 XI23_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(1, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * mu_c(j, i - 1, n3_c) *
                 XI23_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * mu_c(j, i + 1, n3_c) *
                 XI23_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * mu_c(j, i + 2, n3_c) *
                 XI23_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(1, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * mu_c(j, i - 2, n3_c) *
                 XI13_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * mu_c(j, i - 1, n3_c) *
                 XI13_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * mu_c(j, i + 1, n3_c) *
                 XI13_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * mu_c(j, i + 2, n3_c) *
                 XI13_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             (2.0 * mu_c(j, i, n3_c + 1 - k1) + lambda_c(j, i, n3_c + 1 - k1)) *
             XI13_c(j, i, n3_c + 1 - k1) *
             (u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h3_c / h1_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI23_c(j, i, n3_c + 1 - k1) *
             (u_c(2, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(2, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h3_c / h1_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI33_c(j, i, n3_c + 1 - k1) *
             (u_c(3, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(3, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h1_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI23_c(j, i, n3_c + 1 - k1) *
             (u_c(1, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(1, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h3_c / h2_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             lambda_c(j, i, n3_c + 1 - k1) * XI13_c(j, i, n3_c + 1 - k1) *
             (u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h3_c / h2_c;
        // second set equation
        lh_c(j, i, n3_c, 2) =
            lh_c(j, i, n3_c, 2) +
            (-Jacobian_c(j - 2, i, n3_c) * mu_c(j - 2, i, n3_c) *
                 XI23_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * mu_c(j - 1, i, n3_c) *
                 XI23_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * mu_c(j + 1, i, n3_c) *
                 XI23_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * mu_c(j + 2, i, n3_c) *
                 XI23_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, n3_c) * mu_c(j - 2, i, n3_c) *
                 XI13_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(2, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * mu_c(j - 1, i, n3_c) *
                 XI13_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * mu_c(j + 1, i, n3_c) *
                 XI13_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * mu_c(j + 2, i, n3_c) *
                 XI13_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(2, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * lambda_c(j, i - 2, n3_c) *
                 XI13_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(1, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * lambda_c(j, i - 1, n3_c) *
                 XI13_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * lambda_c(j, i + 1, n3_c) *
                 XI13_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * lambda_c(j, i + 2, n3_c) *
                 XI13_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(1, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) *
                 (2.0 * mu_c(j, i - 2, n3_c) + lambda_c(j, i - 2, n3_c)) *
                 XI23_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) *
                 (2.0 * mu_c(j, i - 1, n3_c) + lambda_c(j, i - 1, n3_c)) *
                 XI23_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) *
                 (2.0 * mu_c(j, i + 1, n3_c) + lambda_c(j, i + 1, n3_c)) *
                 XI23_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) *
                 (2.0 * mu_c(j, i + 2, n3_c) + lambda_c(j, i + 2, n3_c)) *
                 XI23_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * lambda_c(j, i - 2, n3_c) *
                 XI33_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(3, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * lambda_c(j, i - 1, n3_c) *
                 XI33_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * lambda_c(j, i + 1, n3_c) *
                 XI33_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * lambda_c(j, i + 2, n3_c) *
                 XI33_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(3, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             lambda_c(j, i, n3_c + 1 - k1) * XI23_c(j, i, n3_c + 1 - k1) *
             (u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h1_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI13_c(j, i, n3_c + 1 - k1) *
             (u_c(2, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(2, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h1_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI13_c(j, i, n3_c + 1 - k1) *
             (u_c(1, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(1, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h3_c / h2_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             (2.0 * mu_c(j, i, n3_c + 1 - k1) + lambda_c(j, i, n3_c + 1 - k1)) *
             XI23_c(j, i, n3_c + 1 - k1) *
             (u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h3_c / h2_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI33_c(j, i, n3_c + 1 - k1) *
             (u_c(3, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(3, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h3_c / h2_c;
        // third set equation
        lh_c(j, i, n3_c, 3) =
            lh_c(j, i, n3_c, 3) +
            (-Jacobian_c(j - 2, i, n3_c) * mu_c(j - 2, i, n3_c) *
                 XI33_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * mu_c(j - 1, i, n3_c) *
                 XI33_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * mu_c(j + 1, i, n3_c) *
                 XI33_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * mu_c(j + 2, i, n3_c) *
                 XI33_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j - 2, i, n3_c) * mu_c(j - 2, i, n3_c) *
                 XI13_c(j - 2, i, n3_c) * bof(1, k1) *
                 u_c(3, j - 2, i, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j - 1, i, n3_c) * mu_c(j - 1, i, n3_c) *
                 XI13_c(j - 1, i, n3_c) * bof(1, k1) *
                 u_c(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j + 1, i, n3_c) * mu_c(j + 1, i, n3_c) *
                 XI13_c(j + 1, i, n3_c) * bof(1, k1) *
                 u_c(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j + 2, i, n3_c) * mu_c(j + 2, i, n3_c) *
                 XI13_c(j + 2, i, n3_c) * bof(1, k1) *
                 u_c(3, j + 2, i, n3_c + 1 - k1) / 12.0) /
                l1 / h1_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * mu_c(j, i - 2, n3_c) *
                 XI33_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * mu_c(j, i - 1, n3_c) *
                 XI33_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * mu_c(j, i + 1, n3_c) *
                 XI33_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * mu_c(j, i + 2, n3_c) *
                 XI33_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-Jacobian_c(j, i - 2, n3_c) * mu_c(j, i - 2, n3_c) *
                 XI23_c(j, i - 2, n3_c) * bof(1, k1) *
                 u_c(3, j, i - 2, n3_c + 1 - k1) / 12.0 +
             Jacobian_c(j, i - 1, n3_c) * mu_c(j, i - 1, n3_c) *
                 XI23_c(j, i - 1, n3_c) * bof(1, k1) *
                 u_c(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
             Jacobian_c(j, i + 1, n3_c) * mu_c(j, i + 1, n3_c) *
                 XI23_c(j, i + 1, n3_c) * bof(1, k1) *
                 u_c(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
             Jacobian_c(j, i + 2, n3_c) * mu_c(j, i + 2, n3_c) *
                 XI23_c(j, i + 2, n3_c) * bof(1, k1) *
                 u_c(3, j, i + 2, n3_c + 1 - k1) / 12.0) /
                l2 / h2_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             lambda_c(j, i, n3_c + 1 - k1) * XI33_c(j, i, n3_c + 1 - k1) *
             (u_c(1, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(1, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h3_c / h1_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI13_c(j, i, n3_c + 1 - k1) *
             (u_c(3, j - 2, i, n3_c + 1 - k1) / 12.0 -
              u_c(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(3, j + 2, i, n3_c + 1 - k1) / 12.0)) /
                l1 / h3_c / h1_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             lambda_c(j, i, n3_c + 1 - k1) * XI33_c(j, i, n3_c + 1 - k1) *
             (u_c(2, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(2, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h2_c / h3_c +
            (-bof(1, k1) * Jacobian_c(j, i, n3_c + 1 - k1) *
             mu_c(j, i, n3_c + 1 - k1) * XI23_c(j, i, n3_c + 1 - k1) *
             (u_c(3, j, i - 2, n3_c + 1 - k1) / 12.0 -
              u_c(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
              u_c(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
              u_c(3, j, i + 2, n3_c + 1 - k1) / 12.0)) /
                l2 / h2_c / h3_c;
      }
    }
  }

  // project
  // first set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, n3_c, 1) /
                                          rho_c(i + l, k + j, n3_c));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, n3_c, 1) / rho_c(i, j + k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, n3_c, 1) / rho_c(i + j, k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, n3_c, 1) / rho_c(i, k, n3_c);
    }
  }
  // restrict
  // first set
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) -
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, 1) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // second set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, n3_c, 2) /
                                          rho_c(i + l, k + j, n3_c));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, n3_c, 2) / rho_c(i, j + k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, n3_c, 2) / rho_c(i + j, k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, n3_c, 2) / rho_c(i, k, n3_c);
    }
  }
  // restriction
  // second set
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) -
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, 1) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // third set
  Mass_f1 = 0.0;
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        for (l = -1; l <= 2; l++) {
          Mass_f1(2 * i, 2 * k) = Mass_f1(2 * i, 2 * k) +
                                  P(j) * (P(l) * lh_c(i + l, k + j, n3_c, 3) /
                                          rho_c(i + l, k + j, n3_c));
        }
      }
    }
  }
  //
  for (k = -1; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i - 1, 2 * k) =
            Mass_f1(2 * i - 1, 2 * k) +
            P(j) * lh_c(i, k + j, n3_c, 3) / rho_c(i, j + k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = -1; i <= n1_c + 1; i++) {
      for (j = -1; j <= 2; j++) {
        Mass_f1(2 * i, 2 * k - 1) =
            Mass_f1(2 * i, 2 * k - 1) +
            P(j) * lh_c(i + j, k, n3_c, 3) / rho_c(i + j, k, n3_c);
      }
    }
  }
  //
  for (k = 0; k <= n2_c + 1; k++) {
    for (i = 0; i <= n1_c + 1; i++) {
      Mass_f1(2 * i - 1, 2 * k - 1) = lh_c(i, k, n3_c, 3) / rho_c(i, k, n3_c);
    }
  }
  // restrict
  // third set
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) -
              17.0 / 48.0 * h3_f * Rop(j) *
                  (Rop(l) * rho_f(2 * i + l, 2 * k + j, 1) *
                   Mass_f1(2 * i + l, 2 * k + j) * 1.0);
        }
      }
    }
  }
  // term 3
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      // second derivative 11  22  12  21
      // first set
      lh_f(i, j, 1, 1) =
          lh_f(i, j, 1, 1) +
          ((-Jacobian_f(i - 2, j, 1) *
                (2.0 * mu_f(i - 2, j, 1) + lambda_f(i - 2, j, 1)) / 8.0 +
            Jacobian_f(i - 1, j, 1) *
                (2.0 * mu_f(i - 1, j, 1) + lambda_f(i - 1, j, 1)) / 6.0 -
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                8.0) *
               u_f(1, i - 2, j, 1) +
           (Jacobian_f(i - 2, j, 1) *
                (2.0 * mu_f(i - 2, j, 1) + lambda_f(i - 2, j, 1)) / 6.0 +
            Jacobian_f(i - 1, j, 1) *
                (2.0 * mu_f(i - 1, j, 1) + lambda_f(i - 1, j, 1)) / 2.0 +
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                2.0 +
            Jacobian_f(i + 1, j, 1) *
                (2.0 * mu_f(i + 1, j, 1) + lambda_f(i + 1, j, 1)) / 6.0) *
               u_f(1, i - 1, j, 1) +
           (-Jacobian_f(i - 2, j, 1) *
                (2.0 * mu_f(i - 2, j, 1) + lambda_f(i - 2, j, 1)) / 24.0 -
            Jacobian_f(i - 1, j, 1) *
                (2.0 * mu_f(i - 1, j, 1) + lambda_f(i - 1, j, 1)) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) *
                3.0 / 4.0 -
            Jacobian_f(i + 1, j, 1) *
                (2.0 * mu_f(i + 1, j, 1) + lambda_f(i + 1, j, 1)) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, 1) *
                (2.0 * mu_f(i + 2, j, 1) + lambda_f(i + 2, j, 1)) / 24.0) *
               u_f(1, i - 0, j, 1) +
           (Jacobian_f(i - 1, j, 1) *
                (2.0 * mu_f(i - 1, j, 1) + lambda_f(i - 1, j, 1)) / 6.0 +
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                2.0 +
            Jacobian_f(i + 1, j, 1) *
                (2.0 * mu_f(i + 1, j, 1) + lambda_f(i + 1, j, 1)) / 2.0 +
            Jacobian_f(i + 2, j, 1) *
                (2.0 * mu_f(i + 2, j, 1) + lambda_f(i + 2, j, 1)) / 6.0) *
               u_f(1, i + 1, j, 1) +
           (-Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                8.0 +
            Jacobian_f(i + 1, j, 1) *
                (2.0 * mu_f(i + 1, j, 1) + lambda_f(i + 1, j, 1)) / 6.0 -
            Jacobian_f(i + 2, j, 1) *
                (2.0 * mu_f(i + 2, j, 1) + lambda_f(i + 2, j, 1)) / 8.0) *
               u_f(1, i + 2, j, 1)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 8.0 +
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0) *
               u_f(1, i, j - 2, 1) +
           (Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 6.0 +
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 2.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 6.0) *
               u_f(1, i, j - 1, 1) +
           (-Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 24.0 -
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) * 3.0 / 4.0 -
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 24.0) *
               u_f(1, i, j - 0, 1) +
           (Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 6.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 2.0 +
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 6.0) *
               u_f(1, i, j + 1, 1) +
           (-Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 6.0 -
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 8.0) *
               u_f(1, i, j + 2, 1)) /
              pow(h2_f, 2) / pow(l2, 2) +
          (Jacobian_f(i - 2, j, 1) * lambda_f(i - 2, j, 1) *
               (u_f(2, i - 2, j - 2, 1) / 12.0 -
                u_f(2, i - 2, j - 1, 1) * 2.0 / 3.0 +
                u_f(2, i - 2, j + 1, 1) * 2.0 / 3.0 -
                u_f(2, i - 2, j + 2, 1) / 12.0) /
               12.0 -
           Jacobian_f(i - 1, j, 1) * lambda_f(i - 1, j, 1) *
               (u_f(2, i - 1, j - 2, 1) / 12.0 -
                u_f(2, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(2, i - 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(2, i - 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i + 1, j, 1) * lambda_f(i + 1, j, 1) *
               (u_f(2, i + 1, j - 2, 1) / 12.0 -
                u_f(2, i + 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(2, i + 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i + 2, j, 1) * lambda_f(i + 2, j, 1) *
               (u_f(2, i + 2, j - 2, 1) / 12.0 -
                u_f(2, i + 2, j - 1, 1) * 2.0 / 3.0 +
                u_f(2, i + 2, j + 1, 1) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f +
          (Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) *
               (u_f(2, i - 2, j - 2, 1) / 12.0 -
                u_f(2, i - 1, j - 2, 1) * 2.0 / 3.0 +
                u_f(2, i + 1, j - 2, 1) * 2.0 / 3.0 -
                u_f(2, i + 2, j - 2, 1) / 12.0) /
               12.0 -
           Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) *
               (u_f(2, i - 2, j - 1, 1) / 12.0 -
                u_f(2, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(2, i + 1, j - 1, 1) * 2.0 / 3.0 -
                u_f(2, i + 2, j - 1, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) *
               (u_f(2, i - 2, j + 1, 1) / 12.0 -
                u_f(2, i - 1, j + 1, 1) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 1, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) *
               (u_f(2, i - 2, j + 2, 1) / 12.0 -
                u_f(2, i - 1, j + 2, 1) * 2.0 / 3.0 +
                u_f(2, i + 1, j + 2, 1) * 2.0 / 3.0 -
                u_f(2, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f;
      // second set
      lh_f(i, j, 1, 2) =
          lh_f(i, j, 1, 2) +
          ((-Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 8.0 +
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0) *
               u_f(2, i - 2, j, 1) +
           (Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 6.0 +
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 2.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 6.0) *
               u_f(2, i - 1, j, 1) +
           (-Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 24.0 -
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) * 3.0 / 4.0 -
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 24.0) *
               u_f(2, i - 0, j, 1) +
           (Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 6.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 2.0 +
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 6.0) *
               u_f(2, i + 1, j, 1) +
           (-Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 6.0 -
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 8.0) *
               u_f(2, i + 2, j, 1)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, 1) *
                (2.0 * mu_f(i, j - 2, 1) + lambda_f(i, j - 2, 1)) / 8.0 +
            Jacobian_f(i, j - 1, 1) *
                (2.0 * mu_f(i, j - 1, 1) + lambda_f(i, j - 1, 1)) / 6.0 -
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                8.0) *
               u_f(2, i, j - 2, 1) +
           (Jacobian_f(i, j - 2, 1) *
                (2.0 * mu_f(i, j - 2, 1) + lambda_f(i, j - 2, 1)) / 6.0 +
            Jacobian_f(i, j - 1, 1) *
                (2.0 * mu_f(i, j - 1, 1) + lambda_f(i, j - 1, 1)) / 2.0 +
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                2.0 +
            Jacobian_f(i, j + 1, 1) *
                (2.0 * mu_f(i, j + 1, 1) + lambda_f(i, j + 1, 1)) / 6.0) *
               u_f(2, i, j - 1, 1) +
           (-Jacobian_f(i, j - 2, 1) *
                (2.0 * mu_f(i, j - 2, 1) + lambda_f(i, j - 2, 1)) / 24.0 -
            Jacobian_f(i, j - 1, 1) *
                (2.0 * mu_f(i, j - 1, 1) + lambda_f(i, j - 1, 1)) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) *
                3.0 / 4.0 -
            Jacobian_f(i, j + 1, 1) *
                (2.0 * mu_f(i, j + 1, 1) + lambda_f(i, j + 1, 1)) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, 1) *
                (2.0 * mu_f(i, j + 2, 1) + lambda_f(i, j + 2, 1)) / 24.0) *
               u_f(2, i, j - 0, 1) +
           (Jacobian_f(i, j - 1, 1) *
                (2.0 * mu_f(i, j - 1, 1) + lambda_f(i, j - 1, 1)) / 6.0 +
            Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                2.0 +
            Jacobian_f(i, j + 1, 1) *
                (2.0 * mu_f(i, j + 1, 1) + lambda_f(i, j + 1, 1)) / 2.0 +
            Jacobian_f(i, j + 2, 1) *
                (2.0 * mu_f(i, j + 2, 1) + lambda_f(i, j + 2, 1)) / 6.0) *
               u_f(2, i, j + 1, 1) +
           (-Jacobian_f(i, j, 1) * (2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) /
                8.0 +
            Jacobian_f(i, j + 1, 1) *
                (2.0 * mu_f(i, j + 1, 1) + lambda_f(i, j + 1, 1)) / 6.0 -
            Jacobian_f(i, j + 2, 1) *
                (2.0 * mu_f(i, j + 2, 1) + lambda_f(i, j + 2, 1)) / 8.0) *
               u_f(2, i, j + 2, 1)) /
              pow(h2_f, 2) / pow(l2, 2) +
          (Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) *
               (u_f(1, i - 2, j - 2, 1) / 12.0 -
                u_f(1, i - 2, j - 1, 1) * 2.0 / 3.0 +
                u_f(1, i - 2, j + 1, 1) * 2.0 / 3.0 -
                u_f(1, i - 2, j + 2, 1) / 12.0) /
               12.0 -
           Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) *
               (u_f(1, i - 1, j - 2, 1) / 12.0 -
                u_f(1, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(1, i - 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(1, i - 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) *
               (u_f(1, i + 1, j - 2, 1) / 12.0 -
                u_f(1, i + 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(1, i + 1, j + 2, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) *
               (u_f(1, i + 2, j - 2, 1) / 12.0 -
                u_f(1, i + 2, j - 1, 1) * 2.0 / 3.0 +
                u_f(1, i + 2, j + 1, 1) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f +
          (Jacobian_f(i, j - 2, 1) * lambda_f(i, j - 2, 1) *
               (u_f(1, i - 2, j - 2, 1) / 12.0 -
                u_f(1, i - 1, j - 2, 1) * 2.0 / 3.0 +
                u_f(1, i + 1, j - 2, 1) * 2.0 / 3.0 -
                u_f(1, i + 2, j - 2, 1) / 12.0) /
               12.0 -
           Jacobian_f(i, j - 1, 1) * lambda_f(i, j - 1, 1) *
               (u_f(1, i - 2, j - 1, 1) / 12.0 -
                u_f(1, i - 1, j - 1, 1) * 2.0 / 3.0 +
                u_f(1, i + 1, j - 1, 1) * 2.0 / 3.0 -
                u_f(1, i + 2, j - 1, 1) / 12.0) *
               2.0 / 3.0 +
           Jacobian_f(i, j + 1, 1) * lambda_f(i, j + 1, 1) *
               (u_f(1, i - 2, j + 1, 1) / 12.0 -
                u_f(1, i - 1, j + 1, 1) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 1, 1) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 1, 1) / 12.0) *
               2.0 / 3.0 -
           Jacobian_f(i, j + 2, 1) * lambda_f(i, j + 2, 1) *
               (u_f(1, i - 2, j + 2, 1) / 12.0 -
                u_f(1, i - 1, j + 2, 1) * 2.0 / 3.0 +
                u_f(1, i + 1, j + 2, 1) * 2.0 / 3.0 -
                u_f(1, i + 2, j + 2, 1) / 12.0) /
               12.0) /
              l1 / l2 / h1_f / h2_f;
      // third set
      lh_f(i, j, 1, 3) =
          lh_f(i, j, 1, 3) +
          ((-Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 8.0 +
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0) *
               u_f(3, i - 2, j, 1) +
           (Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 6.0 +
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 2.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 6.0) *
               u_f(3, i - 1, j, 1) +
           (-Jacobian_f(i - 2, j, 1) * mu_f(i - 2, j, 1) / 24.0 -
            Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) * 3.0 / 4.0 -
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) * 5.0 / 6.0 -
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 24.0) *
               u_f(3, i - 0, j, 1) +
           (Jacobian_f(i - 1, j, 1) * mu_f(i - 1, j, 1) / 6.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 2.0 +
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 6.0) *
               u_f(3, i + 1, j, 1) +
           (-Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0 +
            Jacobian_f(i + 1, j, 1) * mu_f(i + 1, j, 1) / 6.0 -
            Jacobian_f(i + 2, j, 1) * mu_f(i + 2, j, 1) / 8.0) *
               u_f(3, i + 2, j, 1)) /
              pow(h1_f, 2) / pow(l1, 2) +
          ((-Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 8.0 +
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0) *
               u_f(3, i, j - 2, 1) +
           (Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 6.0 +
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 2.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 6.0) *
               u_f(3, i, j - 1, 1) +
           (-Jacobian_f(i, j - 2, 1) * mu_f(i, j - 2, 1) / 24.0 -
            Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) * 3.0 / 4.0 -
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) * 5.0 / 6.0 -
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 24.0) *
               u_f(3, i, j - 0, 1) +
           (Jacobian_f(i, j - 1, 1) * mu_f(i, j - 1, 1) / 6.0 +
            Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 2.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 2.0 +
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 6.0) *
               u_f(3, i, j + 1, 1) +
           (-Jacobian_f(i, j, 1) * mu_f(i, j, 1) / 8.0 +
            Jacobian_f(i, j + 1, 1) * mu_f(i, j + 1, 1) / 6.0 -
            Jacobian_f(i, j + 2, 1) * mu_f(i, j + 2, 1) / 8.0) *
               u_f(3, i, j + 2, 1)) /
              pow(h2_f, 2) / pow(l2, 2);
    }
  }
  //
  for (j = -2; j <= n2_f + 3; j++) {
    for (k = -2; k <= n1_f + 3; k++) {
      for (k1 = 1; k1 <= 8; k1++) {
        for (m = 1; m <= 8; m++) {
          // second derivative 33
          // first set equation
          lh_f(k, j, 1, 1) =
              lh_f(k, j, 1, 1) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI13_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI23_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
               u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI23_f(k, j, m) * u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI33_f(k, j, m) * u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
          // second set equation
          lh_f(k, j, 1, 2) =
              lh_f(k, j, 1, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI23_f(k, j, m) * u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI23_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI13_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
               u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
               XI33_f(k, j, m) * u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
          // third set equation
          lh_f(k, j, 1, 3) =
              lh_f(k, j, 1, 3) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
               XI33_f(k, j, m) * u_f(1, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
               XI33_f(k, j, m) * u_f(2, k, j, k1)) /
                  pow(h3_f, 2) +
              (acof_no_gp(1, k1, m) * Jacobian_f(k, j, m) *
               ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                    pow(XI33_f(k, j, m), 2) +
                mu_f(k, j, m) *
                    (pow(XI13_f(k, j, m), 2) + pow(XI23_f(k, j, m), 2))) *
               u_f(3, k, j, k1)) /
                  pow(h3_f, 2);
        }
      }
    }
  }
  //
  for (i = -2; i <= n2_f + 3; i++) {
    for (j = -2; j <= n1_f + 3; j++) {
      for (k1 = 1; k1 <= 6; k1++) {
        // mixed derivative 13  23  31  32
        // first set equation
        lh_f(j, i, 1, 1) =
            lh_f(j, i, 1, 1) +
            (Jacobian_f(j - 2, i, 1) *
                 (2.0 * mu_f(j - 2, i, 1) + lambda_f(j - 2, i, 1)) *
                 XI13_f(j - 2, i, 1) * bof(1, k1) * u_f(1, j - 2, i, k1) /
                 12.0 -
             Jacobian_f(j - 1, i, 1) *
                 (2.0 * mu_f(j - 1, i, 1) + lambda_f(j - 1, i, 1)) *
                 XI13_f(j - 1, i, 1) * bof(1, k1) * u_f(1, j - 1, i, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, 1) *
                 (2.0 * mu_f(j + 1, i, 1) + lambda_f(j + 1, i, 1)) *
                 XI13_f(j + 1, i, 1) * bof(1, k1) * u_f(1, j + 1, i, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, 1) *
                 (2.0 * mu_f(j + 2, i, 1) + lambda_f(j + 2, i, 1)) *
                 XI13_f(j + 2, i, 1) * bof(1, k1) * u_f(1, j + 2, i, k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, 1) * lambda_f(j - 2, i, 1) *
                 XI23_f(j - 2, i, 1) * bof(1, k1) * u_f(2, j - 2, i, k1) /
                 12.0 -
             Jacobian_f(j - 1, i, 1) * lambda_f(j - 1, i, 1) *
                 XI23_f(j - 1, i, 1) * bof(1, k1) * u_f(2, j - 1, i, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, 1) * lambda_f(j + 1, i, 1) *
                 XI23_f(j + 1, i, 1) * bof(1, k1) * u_f(2, j + 1, i, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, 1) * lambda_f(j + 2, i, 1) *
                 XI23_f(j + 2, i, 1) * bof(1, k1) * u_f(2, j + 2, i, k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, 1) * lambda_f(j - 2, i, 1) *
                 XI33_f(j - 2, i, 1) * bof(1, k1) * u_f(3, j - 2, i, k1) /
                 12.0 -
             Jacobian_f(j - 1, i, 1) * lambda_f(j - 1, i, 1) *
                 XI33_f(j - 1, i, 1) * bof(1, k1) * u_f(3, j - 1, i, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j + 1, i, 1) * lambda_f(j + 1, i, 1) *
                 XI33_f(j + 1, i, 1) * bof(1, k1) * u_f(3, j + 1, i, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j + 2, i, 1) * lambda_f(j + 2, i, 1) *
                 XI33_f(j + 2, i, 1) * bof(1, k1) * u_f(3, j + 2, i, k1) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * mu_f(j, i - 2, 1) * XI23_f(j, i - 2, 1) *
                 bof(1, k1) * u_f(1, j, i - 2, k1) / 12.0 -
             Jacobian_f(j, i - 1, 1) * mu_f(j, i - 1, 1) * XI23_f(j, i - 1, 1) *
                 bof(1, k1) * u_f(1, j, i - 1, k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, 1) * mu_f(j, i + 1, 1) * XI23_f(j, i + 1, 1) *
                 bof(1, k1) * u_f(1, j, i + 1, k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, 1) * mu_f(j, i + 2, 1) * XI23_f(j, i + 2, 1) *
                 bof(1, k1) * u_f(1, j, i + 2, k1) / 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * mu_f(j, i - 2, 1) * XI13_f(j, i - 2, 1) *
                 bof(1, k1) * u_f(2, j, i - 2, k1) / 12.0 -
             Jacobian_f(j, i - 1, 1) * mu_f(j, i - 1, 1) * XI13_f(j, i - 1, 1) *
                 bof(1, k1) * u_f(2, j, i - 1, k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, 1) * mu_f(j, i + 1, 1) * XI13_f(j, i + 1, 1) *
                 bof(1, k1) * u_f(2, j, i + 1, k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, 1) * mu_f(j, i + 2, 1) * XI13_f(j, i + 2, 1) *
                 bof(1, k1) * u_f(2, j, i + 2, k1) / 12.0) /
                l2 / h2_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) *
             (2.0 * mu_f(j, i, k1) + lambda_f(j, i, k1)) * XI13_f(j, i, k1) *
             (u_f(1, j - 2, i, k1) / 12.0 - u_f(1, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, k1) * 2.0 / 3.0 - u_f(1, j + 2, i, k1) / 12.0)) /
                l1 / h3_f / h1_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI23_f(j, i, k1) *
             (u_f(2, j - 2, i, k1) / 12.0 - u_f(2, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(2, j + 1, i, k1) * 2.0 / 3.0 - u_f(2, j + 2, i, k1) / 12.0)) /
                l1 / h3_f / h1_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI33_f(j, i, k1) *
             (u_f(3, j - 2, i, k1) / 12.0 - u_f(3, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(3, j + 1, i, k1) * 2.0 / 3.0 - u_f(3, j + 2, i, k1) / 12.0)) /
                l1 / h1_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI23_f(j, i, k1) *
             (u_f(1, j, i - 2, k1) / 12.0 - u_f(1, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(1, j, i + 1, k1) * 2.0 / 3.0 - u_f(1, j, i + 2, k1) / 12.0)) /
                l2 / h3_f / h2_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * lambda_f(j, i, k1) *
             XI13_f(j, i, k1) *
             (u_f(2, j, i - 2, k1) / 12.0 - u_f(2, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, k1) * 2.0 / 3.0 - u_f(2, j, i + 2, k1) / 12.0)) /
                l2 / h3_f / h2_f;
        // second set equation
        lh_f(j, i, 1, 2) =
            lh_f(j, i, 1, 2) +
            (Jacobian_f(j - 2, i, 1) * mu_f(j - 2, i, 1) * XI23_f(j - 2, i, 1) *
                 bof(1, k1) * u_f(1, j - 2, i, k1) / 12.0 -
             Jacobian_f(j - 1, i, 1) * mu_f(j - 1, i, 1) * XI23_f(j - 1, i, 1) *
                 bof(1, k1) * u_f(1, j - 1, i, k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, 1) * mu_f(j + 1, i, 1) * XI23_f(j + 1, i, 1) *
                 bof(1, k1) * u_f(1, j + 1, i, k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, 1) * mu_f(j + 2, i, 1) * XI23_f(j + 2, i, 1) *
                 bof(1, k1) * u_f(1, j + 2, i, k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, 1) * mu_f(j - 2, i, 1) * XI13_f(j - 2, i, 1) *
                 bof(1, k1) * u_f(2, j - 2, i, k1) / 12.0 -
             Jacobian_f(j - 1, i, 1) * mu_f(j - 1, i, 1) * XI13_f(j - 1, i, 1) *
                 bof(1, k1) * u_f(2, j - 1, i, k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, 1) * mu_f(j + 1, i, 1) * XI13_f(j + 1, i, 1) *
                 bof(1, k1) * u_f(2, j + 1, i, k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, 1) * mu_f(j + 2, i, 1) * XI13_f(j + 2, i, 1) *
                 bof(1, k1) * u_f(2, j + 2, i, k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * lambda_f(j, i - 2, 1) *
                 XI13_f(j, i - 2, 1) * bof(1, k1) * u_f(1, j, i - 2, k1) /
                 12.0 -
             Jacobian_f(j, i - 1, 1) * lambda_f(j, i - 1, 1) *
                 XI13_f(j, i - 1, 1) * bof(1, k1) * u_f(1, j, i - 1, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, 1) * lambda_f(j, i + 1, 1) *
                 XI13_f(j, i + 1, 1) * bof(1, k1) * u_f(1, j, i + 1, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, 1) * lambda_f(j, i + 2, 1) *
                 XI13_f(j, i + 2, 1) * bof(1, k1) * u_f(1, j, i + 2, k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, 1) *
                 (2.0 * mu_f(j, i - 2, 1) + lambda_f(j, i - 2, 1)) *
                 XI23_f(j, i - 2, 1) * bof(1, k1) * u_f(2, j, i - 2, k1) /
                 12.0 -
             Jacobian_f(j, i - 1, 1) *
                 (2.0 * mu_f(j, i - 1, 1) + lambda_f(j, i - 1, 1)) *
                 XI23_f(j, i - 1, 1) * bof(1, k1) * u_f(2, j, i - 1, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, 1) *
                 (2.0 * mu_f(j, i + 1, 1) + lambda_f(j, i + 1, 1)) *
                 XI23_f(j, i + 1, 1) * bof(1, k1) * u_f(2, j, i + 1, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, 1) *
                 (2.0 * mu_f(j, i + 2, 1) + lambda_f(j, i + 2, 1)) *
                 XI23_f(j, i + 2, 1) * bof(1, k1) * u_f(2, j, i + 2, k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * lambda_f(j, i - 2, 1) *
                 XI33_f(j, i - 2, 1) * bof(1, k1) * u_f(3, j, i - 2, k1) /
                 12.0 -
             Jacobian_f(j, i - 1, 1) * lambda_f(j, i - 1, 1) *
                 XI33_f(j, i - 1, 1) * bof(1, k1) * u_f(3, j, i - 1, k1) * 2.0 /
                 3.0 +
             Jacobian_f(j, i + 1, 1) * lambda_f(j, i + 1, 1) *
                 XI33_f(j, i + 1, 1) * bof(1, k1) * u_f(3, j, i + 1, k1) * 2.0 /
                 3.0 -
             Jacobian_f(j, i + 2, 1) * lambda_f(j, i + 2, 1) *
                 XI33_f(j, i + 2, 1) * bof(1, k1) * u_f(3, j, i + 2, k1) /
                 12.0) /
                l2 / h2_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * lambda_f(j, i, k1) *
             XI23_f(j, i, k1) *
             (u_f(1, j - 2, i, k1) / 12.0 - u_f(1, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, k1) * 2.0 / 3.0 - u_f(1, j + 2, i, k1) / 12.0)) /
                l1 / h1_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI13_f(j, i, k1) *
             (u_f(2, j - 2, i, k1) / 12.0 - u_f(2, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(2, j + 1, i, k1) * 2.0 / 3.0 - u_f(2, j + 2, i, k1) / 12.0)) /
                l1 / h1_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI13_f(j, i, k1) *
             (u_f(1, j, i - 2, k1) / 12.0 - u_f(1, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(1, j, i + 1, k1) * 2.0 / 3.0 - u_f(1, j, i + 2, k1) / 12.0)) /
                l2 / h3_f / h2_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) *
             (2.0 * mu_f(j, i, k1) + lambda_f(j, i, k1)) * XI23_f(j, i, k1) *
             (u_f(2, j, i - 2, k1) / 12.0 - u_f(2, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, k1) * 2.0 / 3.0 - u_f(2, j, i + 2, k1) / 12.0)) /
                l2 / h3_f / h2_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI33_f(j, i, k1) *
             (u_f(3, j, i - 2, k1) / 12.0 - u_f(3, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(3, j, i + 1, k1) * 2.0 / 3.0 - u_f(3, j, i + 2, k1) / 12.0)) /
                l2 / h3_f / h2_f;
        // third set equation
        lh_f(j, i, 1, 3) =
            lh_f(j, i, 1, 3) +
            (Jacobian_f(j - 2, i, 1) * mu_f(j - 2, i, 1) * XI33_f(j - 2, i, 1) *
                 bof(1, k1) * u_f(1, j - 2, i, k1) / 12.0 -
             Jacobian_f(j - 1, i, 1) * mu_f(j - 1, i, 1) * XI33_f(j - 1, i, 1) *
                 bof(1, k1) * u_f(1, j - 1, i, k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, 1) * mu_f(j + 1, i, 1) * XI33_f(j + 1, i, 1) *
                 bof(1, k1) * u_f(1, j + 1, i, k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, 1) * mu_f(j + 2, i, 1) * XI33_f(j + 2, i, 1) *
                 bof(1, k1) * u_f(1, j + 2, i, k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, 1) * mu_f(j - 2, i, 1) * XI13_f(j - 2, i, 1) *
                 bof(1, k1) * u_f(3, j - 2, i, k1) / 12.0 -
             Jacobian_f(j - 1, i, 1) * mu_f(j - 1, i, 1) * XI13_f(j - 1, i, 1) *
                 bof(1, k1) * u_f(3, j - 1, i, k1) * 2.0 / 3.0 +
             Jacobian_f(j + 1, i, 1) * mu_f(j + 1, i, 1) * XI13_f(j + 1, i, 1) *
                 bof(1, k1) * u_f(3, j + 1, i, k1) * 2.0 / 3.0 -
             Jacobian_f(j + 2, i, 1) * mu_f(j + 2, i, 1) * XI13_f(j + 2, i, 1) *
                 bof(1, k1) * u_f(3, j + 2, i, k1) / 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * mu_f(j, i - 2, 1) * XI33_f(j, i - 2, 1) *
                 bof(1, k1) * u_f(2, j, i - 2, k1) / 12.0 -
             Jacobian_f(j, i - 1, 1) * mu_f(j, i - 1, 1) * XI33_f(j, i - 1, 1) *
                 bof(1, k1) * u_f(2, j, i - 1, k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, 1) * mu_f(j, i + 1, 1) * XI33_f(j, i + 1, 1) *
                 bof(1, k1) * u_f(2, j, i + 1, k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, 1) * mu_f(j, i + 2, 1) * XI33_f(j, i + 2, 1) *
                 bof(1, k1) * u_f(2, j, i + 2, k1) / 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, 1) * mu_f(j, i - 2, 1) * XI23_f(j, i - 2, 1) *
                 bof(1, k1) * u_f(3, j, i - 2, k1) / 12.0 -
             Jacobian_f(j, i - 1, 1) * mu_f(j, i - 1, 1) * XI23_f(j, i - 1, 1) *
                 bof(1, k1) * u_f(3, j, i - 1, k1) * 2.0 / 3.0 +
             Jacobian_f(j, i + 1, 1) * mu_f(j, i + 1, 1) * XI23_f(j, i + 1, 1) *
                 bof(1, k1) * u_f(3, j, i + 1, k1) * 2.0 / 3.0 -
             Jacobian_f(j, i + 2, 1) * mu_f(j, i + 2, 1) * XI23_f(j, i + 2, 1) *
                 bof(1, k1) * u_f(3, j, i + 2, k1) / 12.0) /
                l2 / h2_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * lambda_f(j, i, k1) *
             XI33_f(j, i, k1) *
             (u_f(1, j - 2, i, k1) / 12.0 - u_f(1, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(1, j + 1, i, k1) * 2.0 / 3.0 - u_f(1, j + 2, i, k1) / 12.0)) /
                l1 / h1_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI13_f(j, i, k1) *
             (u_f(3, j - 2, i, k1) / 12.0 - u_f(3, j - 1, i, k1) * 2.0 / 3.0 +
              u_f(3, j + 1, i, k1) * 2.0 / 3.0 - u_f(3, j + 2, i, k1) / 12.0)) /
                l1 / h3_f / h1_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * lambda_f(j, i, k1) *
             XI33_f(j, i, k1) *
             (u_f(2, j, i - 2, k1) / 12.0 - u_f(2, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(2, j, i + 1, k1) * 2.0 / 3.0 - u_f(2, j, i + 2, k1) / 12.0)) /
                l2 / h2_f / h3_f +
            (bof(1, k1) * Jacobian_f(j, i, k1) * mu_f(j, i, k1) *
             XI23_f(j, i, k1) *
             (u_f(3, j, i - 2, k1) / 12.0 - u_f(3, j, i - 1, k1) * 2.0 / 3.0 +
              u_f(3, j, i + 1, k1) * 2.0 / 3.0 - u_f(3, j, i + 2, k1) / 12.0)) /
                l2 / h2_f / h3_f;
      }
    }
  }

  // scale L_f
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      lh_f(i, j, 1, 1) = lh_f(i, j, 1, 1) * 17.0 / 48.0 * h3_f;
      lh_f(i, j, 1, 2) = lh_f(i, j, 1, 2) * 17.0 / 48.0 * h3_f;
      lh_f(i, j, 1, 3) = lh_f(i, j, 1, 3) * 17.0 / 48.0 * h3_f;
    }
  }
  //
  for (j = -2; j <= n2_f + 3; j++) {
    for (i = -2; i <= n1_f + 3; i++) {
      for (k = 1; k <= 5; k++) {
        // first set equation
        lh_f(i, j, 1, 1) =
            lh_f(i, j, 1, 1) +
            Jacobian_f(i, j, 1) *
                ((2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) *
                     pow(XI13_f(i, j, 1), 2) +
                 mu_f(i, j, 1) *
                     (pow(XI23_f(i, j, 1), 2) + pow(XI33_f(i, j, 1), 2))) *
                sbop_no_gp(k) * u_f(1, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI13_f(i, j, 1) * XI23_f(i, j, 1) * sbop_no_gp(k) *
                u_f(2, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI13_f(i, j, 1) * XI33_f(i, j, 1) * sbop_no_gp(k) *
                u_f(3, i, j, k) / h3_f;
        // second set equation
        lh_f(i, j, 1, 2) =
            lh_f(i, j, 1, 2) +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI13_f(i, j, 1) * XI23_f(i, j, 1) * sbop_no_gp(k) *
                u_f(1, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) *
                ((2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) *
                     pow(XI23_f(i, j, 1), 2) +
                 mu_f(i, j, 1) *
                     (pow(XI13_f(i, j, 1), 2) + pow(XI33_f(i, j, 1), 2))) *
                sbop_no_gp(k) * u_f(2, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI23_f(i, j, 1) * XI33_f(i, j, 1) * sbop_no_gp(k) *
                u_f(3, i, j, k) / h3_f;
        // third set equation
        lh_f(i, j, 1, 3) =
            lh_f(i, j, 1, 3) +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI13_f(i, j, 1) * XI33_f(i, j, 1) * sbop_no_gp(k) *
                u_f(1, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) * (lambda_f(i, j, 1) + mu_f(i, j, 1)) *
                XI23_f(i, j, 1) * XI33_f(i, j, 1) * sbop_no_gp(k) *
                u_f(2, i, j, k) / h3_f +
            Jacobian_f(i, j, 1) *
                ((2.0 * mu_f(i, j, 1) + lambda_f(i, j, 1)) *
                     pow(XI33_f(i, j, 1), 2) +
                 mu_f(i, j, 1) *
                     (pow(XI13_f(i, j, 1), 2) + pow(XI23_f(i, j, 1), 2))) *
                sbop_no_gp(k) * u_f(3, i, j, k) / h3_f;
      }
    }
  }
  //
  for (k = -2; k <= n2_f + 3; k++) {
    for (i = -2; i <= n1_f + 3; i++) {
      for (j = -2; j <= 2; j++) {
        // 31  32
        // first set equation
        lh_f(i, k, 1, 1) =
            lh_f(i, k, 1, 1) +
            Jacobian_f(i, k, 1) * (2.0 * mu_f(i, k, 1) + lambda_f(i, k, 1)) /
                l1 * XI13_f(i, k, 1) * ux_cof(j) * u_f(1, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l1 * XI23_f(i, k, 1) *
                ux_cof(j) * u_f(2, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l1 * XI33_f(i, k, 1) *
                ux_cof(j) * u_f(3, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l2 * XI23_f(i, k, 1) *
                ux_cof(j) * u_f(1, i, k + j, 1) / h2_f +
            Jacobian_f(i, k, 1) * lambda_f(i, k, 1) / l2 * XI13_f(i, k, 1) *
                ux_cof(j) * u_f(2, i, k + j, 1) / h2_f;
        // second set equation
        lh_f(i, k, 1, 2) =
            lh_f(i, k, 1, 2) +
            Jacobian_f(i, k, 1) * lambda_f(i, k, 1) / l1 * XI23_f(i, k, 1) *
                ux_cof(j) * u_f(1, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l1 * XI13_f(i, k, 1) *
                ux_cof(j) * u_f(2, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l2 * XI13_f(i, k, 1) *
                ux_cof(j) * u_f(1, i, k + j, 1) / h2_f +
            Jacobian_f(i, k, 1) * (2.0 * mu_f(i, k, 1) + lambda_f(i, k, 1)) /
                l2 * XI23_f(i, k, 1) * ux_cof(j) * u_f(2, i, k + j, 1) / h2_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l2 * XI33_f(i, k, 1) *
                ux_cof(j) * u_f(3, i, k + j, 1) / h2_f;
        // third set equation
        lh_f(i, k, 1, 3) =
            lh_f(i, k, 1, 3) +
            Jacobian_f(i, k, 1) * lambda_f(i, k, 1) / l1 * XI33_f(i, k, 1) *
                ux_cof(j) * u_f(1, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l1 * XI13_f(i, k, 1) *
                ux_cof(j) * u_f(3, i + j, k, 1) / h1_f +
            Jacobian_f(i, k, 1) * lambda_f(i, k, 1) / l2 * XI33_f(i, k, 1) *
                ux_cof(j) * u_f(2, i, k + j, 1) / h2_f +
            Jacobian_f(i, k, 1) * mu_f(i, k, 1) / l2 * XI23_f(i, k, 1) *
                ux_cof(j) * u_f(3, i, k + j, 1) / h2_f;
      }
    }
  }
  // now restrict it to the coarse grid
  for (k = 1; k <= n2_c; k++) {
    for (i = 1; i <= n1_c; i++) {
      for (j = -4; j <= 2; j++) {
        for (l = -4; l <= 2; l++) {
          // first set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 1) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, 1, 1));
          // second set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 2) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, 1, 2));
          // third set
          Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) =
              Vass((k - 1) * 3 * n1_c + (i - 1) * 3 + 3) +
              Rop(j) * (Rop(l) * lh_f(2 * i + l, 2 * k + j, 1, 3));
        }
      }
    }
  }
}  // END INTERFACE_RHS

void update_interior(Sarray &u_c_t, Sarray &u_f_t, Farray &bof, Farray &ghcof,
                     Farray &acof, Farray &acof_no_gp, Farray &lh_f,
                     Sarray &Jacobian_f, Sarray &mu_f, Sarray &lambda_f,
                     Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f,
                     Farray &lh_c, Sarray &Jacobian_c, Sarray &mu_c,
                     Sarray &lambda_c, Sarray &XI13_c, Sarray &XI23_c,
                     Sarray &XI33_c, PackArgs &a) {
  // std::cout<<"START...";
  // real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) ::
  // u_c_t,lh_c real(dp), dimension
  // (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: u_f_t,lh_f
  int k1, k, j, i, m;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  float_sw4 l3 = a.l3;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n1_f = a.n1_f;
  auto n2_f = a.n2_f;
  auto n3_f = a.n3_f;

  auto h1_c = a.h1_c;
  auto h2_c = a.h2_c;
  auto h3_c = a.h3_c;

  auto h1_f = a.h1_f;
  auto h2_f = a.h2_f;
  auto h3_f = a.h3_f;

  auto nrg = a.nrg;

  // Difference operators in the interior of the domains
  // fine mesh
  lh_f = 0.0;
  for (k = 1; k <= n3_f; k++) {
    for (j = 1; j <= n2_f; j++) {
      for (i = 1; i <= n1_f; i++) {
        // second derivative 11  22  12  21
        // first set
        lh_f(i, j, k, 1) =
            lh_f(i, j, k, 1) +
            ((-Jacobian_f(i - 2, j, k) *
                  (2.0 * mu_f(i - 2, j, k) + lambda_f(i - 2, j, k)) / 8.0 +
              Jacobian_f(i - 1, j, k) *
                  (2.0 * mu_f(i - 1, j, k) + lambda_f(i - 1, j, k)) / 6.0 -
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  8.0) *
                 u_f_t(1, i - 2, j, k) +
             (Jacobian_f(i - 2, j, k) *
                  (2.0 * mu_f(i - 2, j, k) + lambda_f(i - 2, j, k)) / 6.0 +
              Jacobian_f(i - 1, j, k) *
                  (2.0 * mu_f(i - 1, j, k) + lambda_f(i - 1, j, k)) / 2.0 +
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  2.0 +
              Jacobian_f(i + 1, j, k) *
                  (2.0 * mu_f(i + 1, j, k) + lambda_f(i + 1, j, k)) / 6.0) *
                 u_f_t(1, i - 1, j, k) +
             (-Jacobian_f(i - 2, j, k) *
                  (2.0 * mu_f(i - 2, j, k) + lambda_f(i - 2, j, k)) / 24.0 -
              Jacobian_f(i - 1, j, k) *
                  (2.0 * mu_f(i - 1, j, k) + lambda_f(i - 1, j, k)) * 5.0 /
                  6.0 -
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) *
                  3.0 / 4.0 -
              Jacobian_f(i + 1, j, k) *
                  (2.0 * mu_f(i + 1, j, k) + lambda_f(i + 1, j, k)) * 5.0 /
                  6.0 -
              Jacobian_f(i + 2, j, k) *
                  (2.0 * mu_f(i + 2, j, k) + lambda_f(i + 2, j, k)) / 24.0) *
                 u_f_t(1, i - 0, j, k) +
             (Jacobian_f(i - 1, j, k) *
                  (2.0 * mu_f(i - 1, j, k) + lambda_f(i - 1, j, k)) / 6.0 +
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  2.0 +
              Jacobian_f(i + 1, j, k) *
                  (2.0 * mu_f(i + 1, j, k) + lambda_f(i + 1, j, k)) / 2.0 +
              Jacobian_f(i + 2, j, k) *
                  (2.0 * mu_f(i + 2, j, k) + lambda_f(i + 2, j, k)) / 6.0) *
                 u_f_t(1, i + 1, j, k) +
             (-Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  8.0 +
              Jacobian_f(i + 1, j, k) *
                  (2.0 * mu_f(i + 1, j, k) + lambda_f(i + 1, j, k)) / 6.0 -
              Jacobian_f(i + 2, j, k) *
                  (2.0 * mu_f(i + 2, j, k) + lambda_f(i + 2, j, k)) / 8.0) *
                 u_f_t(1, i + 2, j, k)) /
                pow(h1_f, 2) / pow(l1, 2) +
            ((-Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 8.0 +
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0) *
                 u_f_t(1, i, j - 2, k) +
             (Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 6.0 +
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 2.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 6.0) *
                 u_f_t(1, i, j - 1, k) +
             (-Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 24.0 -
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) * 5.0 / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) * 3.0 / 4.0 -
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) * 5.0 / 6.0 -
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 24.0) *
                 u_f_t(1, i, j - 0, k) +
             (Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 6.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 2.0 +
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 6.0) *
                 u_f_t(1, i, j + 1, k) +
             (-Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 6.0 -
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 8.0) *
                 u_f_t(1, i, j + 2, k)) /
                pow(h2_f, 2) / pow(l2, 2) +
            (Jacobian_f(i - 2, j, k) * lambda_f(i - 2, j, k) *
                 (u_f_t(2, i - 2, j - 2, k) / 12.0 -
                  u_f_t(2, i - 2, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i - 2, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i - 2, j + 2, k) / 12.0) /
                 12.0 -
             Jacobian_f(i - 1, j, k) * lambda_f(i - 1, j, k) *
                 (u_f_t(2, i - 1, j - 2, k) / 12.0 -
                  u_f_t(2, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i - 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i - 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(i + 1, j, k) * lambda_f(i + 1, j, k) *
                 (u_f_t(2, i + 1, j - 2, k) / 12.0 -
                  u_f_t(2, i + 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(i + 2, j, k) * lambda_f(i + 2, j, k) *
                 (u_f_t(2, i + 2, j - 2, k) / 12.0 -
                  u_f_t(2, i + 2, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 2, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_f / h2_f +
            (Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) *
                 (u_f_t(2, i - 2, j - 2, k) / 12.0 -
                  u_f_t(2, i - 1, j - 2, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 1, j - 2, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 2, j - 2, k) / 12.0) /
                 12.0 -
             Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) *
                 (u_f_t(2, i - 2, j - 1, k) / 12.0 -
                  u_f_t(2, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 1, j - 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 2, j - 1, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) *
                 (u_f_t(2, i - 2, j + 1, k) / 12.0 -
                  u_f_t(2, i - 1, j + 1, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 2, j + 1, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) *
                 (u_f_t(2, i - 2, j + 2, k) / 12.0 -
                  u_f_t(2, i - 1, j + 2, k) * 2.0 / 3.0 +
                  u_f_t(2, i + 1, j + 2, k) * 2.0 / 3.0 -
                  u_f_t(2, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_f / h2_f;
        // second set
        lh_f(i, j, k, 2) =
            lh_f(i, j, k, 2) +
            ((-Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 8.0 +
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0) *
                 u_f_t(2, i - 2, j, k) +
             (Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 6.0 +
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 2.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 6.0) *
                 u_f_t(2, i - 1, j, k) +
             (-Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 24.0 -
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) * 5.0 / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) * 3.0 / 4.0 -
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) * 5.0 / 6.0 -
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 24.0) *
                 u_f_t(2, i - 0, j, k) +
             (Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 6.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 2.0 +
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 6.0) *
                 u_f_t(2, i + 1, j, k) +
             (-Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 6.0 -
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 8.0) *
                 u_f_t(2, i + 2, j, k)) /
                pow(h1_f, 2) / pow(l1, 2) +
            ((-Jacobian_f(i, j - 2, k) *
                  (2.0 * mu_f(i, j - 2, k) + lambda_f(i, j - 2, k)) / 8.0 +
              Jacobian_f(i, j - 1, k) *
                  (2.0 * mu_f(i, j - 1, k) + lambda_f(i, j - 1, k)) / 6.0 -
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  8.0) *
                 u_f_t(2, i, j - 2, k) +
             (Jacobian_f(i, j - 2, k) *
                  (2.0 * mu_f(i, j - 2, k) + lambda_f(i, j - 2, k)) / 6.0 +
              Jacobian_f(i, j - 1, k) *
                  (2.0 * mu_f(i, j - 1, k) + lambda_f(i, j - 1, k)) / 2.0 +
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  2.0 +
              Jacobian_f(i, j + 1, k) *
                  (2.0 * mu_f(i, j + 1, k) + lambda_f(i, j + 1, k)) / 6.0) *
                 u_f_t(2, i, j - 1, k) +
             (-Jacobian_f(i, j - 2, k) *
                  (2.0 * mu_f(i, j - 2, k) + lambda_f(i, j - 2, k)) / 24.0 -
              Jacobian_f(i, j - 1, k) *
                  (2.0 * mu_f(i, j - 1, k) + lambda_f(i, j - 1, k)) * 5.0 /
                  6.0 -
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) *
                  3.0 / 4.0 -
              Jacobian_f(i, j + 1, k) *
                  (2.0 * mu_f(i, j + 1, k) + lambda_f(i, j + 1, k)) * 5.0 /
                  6.0 -
              Jacobian_f(i, j + 2, k) *
                  (2.0 * mu_f(i, j + 2, k) + lambda_f(i, j + 2, k)) / 24.0) *
                 u_f_t(2, i, j - 0, k) +
             (Jacobian_f(i, j - 1, k) *
                  (2.0 * mu_f(i, j - 1, k) + lambda_f(i, j - 1, k)) / 6.0 +
              Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  2.0 +
              Jacobian_f(i, j + 1, k) *
                  (2.0 * mu_f(i, j + 1, k) + lambda_f(i, j + 1, k)) / 2.0 +
              Jacobian_f(i, j + 2, k) *
                  (2.0 * mu_f(i, j + 2, k) + lambda_f(i, j + 2, k)) / 6.0) *
                 u_f_t(2, i, j + 1, k) +
             (-Jacobian_f(i, j, k) * (2.0 * mu_f(i, j, k) + lambda_f(i, j, k)) /
                  8.0 +
              Jacobian_f(i, j + 1, k) *
                  (2.0 * mu_f(i, j + 1, k) + lambda_f(i, j + 1, k)) / 6.0 -
              Jacobian_f(i, j + 2, k) *
                  (2.0 * mu_f(i, j + 2, k) + lambda_f(i, j + 2, k)) / 8.0) *
                 u_f_t(2, i, j + 2, k)) /
                pow(h2_f, 2) / pow(l2, 2) +
            (Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) *
                 (u_f_t(1, i - 2, j - 2, k) / 12.0 -
                  u_f_t(1, i - 2, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i - 2, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i - 2, j + 2, k) / 12.0) /
                 12.0 -
             Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) *
                 (u_f_t(1, i - 1, j - 2, k) / 12.0 -
                  u_f_t(1, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i - 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i - 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) *
                 (u_f_t(1, i + 1, j - 2, k) / 12.0 -
                  u_f_t(1, i + 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) *
                 (u_f_t(1, i + 2, j - 2, k) / 12.0 -
                  u_f_t(1, i + 2, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 2, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_f / h2_f +
            (Jacobian_f(i, j - 2, k) * lambda_f(i, j - 2, k) *
                 (u_f_t(1, i - 2, j - 2, k) / 12.0 -
                  u_f_t(1, i - 1, j - 2, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 1, j - 2, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 2, j - 2, k) / 12.0) /
                 12.0 -
             Jacobian_f(i, j - 1, k) * lambda_f(i, j - 1, k) *
                 (u_f_t(1, i - 2, j - 1, k) / 12.0 -
                  u_f_t(1, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 1, j - 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 2, j - 1, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(i, j + 1, k) * lambda_f(i, j + 1, k) *
                 (u_f_t(1, i - 2, j + 1, k) / 12.0 -
                  u_f_t(1, i - 1, j + 1, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 2, j + 1, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(i, j + 2, k) * lambda_f(i, j + 2, k) *
                 (u_f_t(1, i - 2, j + 2, k) / 12.0 -
                  u_f_t(1, i - 1, j + 2, k) * 2.0 / 3.0 +
                  u_f_t(1, i + 1, j + 2, k) * 2.0 / 3.0 -
                  u_f_t(1, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_f / h2_f;
        // third set
        lh_f(i, j, k, 3) =
            lh_f(i, j, k, 3) +
            ((-Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 8.0 +
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0) *
                 u_f_t(3, i - 2, j, k) +
             (Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 6.0 +
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 2.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 6.0) *
                 u_f_t(3, i - 1, j, k) +
             (-Jacobian_f(i - 2, j, k) * mu_f(i - 2, j, k) / 24.0 -
              Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) * 5.0 / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) * 3.0 / 4.0 -
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) * 5.0 / 6.0 -
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 24.0) *
                 u_f_t(3, i - 0, j, k) +
             (Jacobian_f(i - 1, j, k) * mu_f(i - 1, j, k) / 6.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 2.0 +
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 6.0) *
                 u_f_t(3, i + 1, j, k) +
             (-Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0 +
              Jacobian_f(i + 1, j, k) * mu_f(i + 1, j, k) / 6.0 -
              Jacobian_f(i + 2, j, k) * mu_f(i + 2, j, k) / 8.0) *
                 u_f_t(3, i + 2, j, k)) /
                pow(h1_f, 2) / pow(l1, 2) +
            ((-Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 8.0 +
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0) *
                 u_f_t(3, i, j - 2, k) +
             (Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 6.0 +
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 2.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 6.0) *
                 u_f_t(3, i, j - 1, k) +
             (-Jacobian_f(i, j - 2, k) * mu_f(i, j - 2, k) / 24.0 -
              Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) * 5.0 / 6.0 -
              Jacobian_f(i, j, k) * mu_f(i, j, k) * 3.0 / 4.0 -
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) * 5.0 / 6.0 -
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 24.0) *
                 u_f_t(3, i, j - 0, k) +
             (Jacobian_f(i, j - 1, k) * mu_f(i, j - 1, k) / 6.0 +
              Jacobian_f(i, j, k) * mu_f(i, j, k) / 2.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 2.0 +
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 6.0) *
                 u_f_t(3, i, j + 1, k) +
             (-Jacobian_f(i, j, k) * mu_f(i, j, k) / 8.0 +
              Jacobian_f(i, j + 1, k) * mu_f(i, j + 1, k) / 6.0 -
              Jacobian_f(i, j + 2, k) * mu_f(i, j + 2, k) / 8.0) *
                 u_f_t(3, i, j + 2, k)) /
                pow(h2_f, 2) / pow(l2, 2);
      }
    }
  }
  //
  for (i = 7; i <= n3_f - 6; i++) {
    for (j = 1; j <= n2_f; j++) {
      for (k = 1; k <= n1_f; k++) {
        // second derivative 33
        // first set equation
        lh_f(k, j, i, 1) =
            lh_f(k, j, i, 1) +
            ((-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI13_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI23_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI13_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI23_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI13_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI23_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  8.0) *
                 u_f_t(1, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI13_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI23_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI13_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI23_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI13_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI23_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI13_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI23_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  6.0) *
                 u_f_t(1, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI13_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI23_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI13_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI23_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI13_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI23_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI13_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI23_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI13_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI23_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  24.0) *
                 u_f_t(1, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI13_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI23_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI13_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI23_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI13_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI23_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI13_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI23_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  6.0) *
                 u_f_t(1, k, j, i + 1) +
             (-Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI13_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI23_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  8.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI13_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI23_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI13_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI23_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  8.0) *
                 u_f_t(1, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 8.0) *
                 u_f_t(2, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 6.0) *
                 u_f_t(2, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 24.0) *
                 u_f_t(2, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 6.0) *
                 u_f_t(2, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 8.0) *
                 u_f_t(2, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 8.0) *
                 u_f_t(3, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0) *
                 u_f_t(3, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 24.0) *
                 u_f_t(3, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 6.0) *
                 u_f_t(3, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 8.0) *
                 u_f_t(3, k, j, i + 2)) /
                pow(h3_f, 2);
        // second set of equation
        lh_f(k, j, i, 2) =
            lh_f(k, j, i, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 8.0) *
                 u_f_t(1, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 6.0) *
                 u_f_t(1, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI23_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 24.0) *
                 u_f_t(1, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI23_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 6.0) *
                 u_f_t(1, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI23_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI23_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI23_f(k, j, i + 2) / 8.0) *
                 u_f_t(1, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI23_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI23_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI23_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  8.0) *
                 u_f_t(2, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI23_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI23_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI23_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI23_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  6.0) *
                 u_f_t(2, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI23_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI33_f(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI23_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI23_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI23_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI23_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  24.0) *
                 u_f_t(2, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI23_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI33_f(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI23_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI23_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI23_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  6.0) *
                 u_f_t(2, k, j, i + 1) +
             (-Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI23_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI33_f(k, j, i), 2))) /
                  8.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI23_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI33_f(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI23_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI33_f(k, j, i + 2), 2))) /
                  8.0) *
                 u_f_t(2, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 8.0) *
                 u_f_t(3, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0) *
                 u_f_t(3, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 24.0) *
                 u_f_t(3, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 6.0) *
                 u_f_t(3, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 8.0) *
                 u_f_t(3, k, j, i + 2)) /
                pow(h3_f, 2);
        // third set equation
        lh_f(k, j, i, 3) =
            lh_f(k, j, i, 3) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 8.0) *
                 u_f_t(1, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0) *
                 u_f_t(1, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI13_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 24.0) *
                 u_f_t(1, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI13_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 6.0) *
                 u_f_t(1, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI13_f(k, j, i) * XI33_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI13_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI13_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 8.0) *
                 u_f_t(1, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 8.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 8.0) *
                 u_f_t(2, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 6.0 +
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 2.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0) *
                 u_f_t(2, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  (mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                  XI23_f(k, j, i - 2) * XI33_f(k, j, i - 2) / 24.0 -
              Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) * 3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 24.0) *
                 u_f_t(2, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  (mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                  XI23_f(k, j, i - 1) * XI33_f(k, j, i - 1) / 6.0 +
              Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 2.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 2.0 +
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 6.0) *
                 u_f_t(2, k, j, i + 1) +
             (-Jacobian_f(k, j, i) * (mu_f(k, j, i) + lambda_f(k, j, i)) *
                  XI23_f(k, j, i) * XI33_f(k, j, i) / 8.0 +
              Jacobian_f(k, j, i + 1) *
                  (mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                  XI23_f(k, j, i + 1) * XI33_f(k, j, i + 1) / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  (mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                  XI23_f(k, j, i + 2) * XI33_f(k, j, i + 2) / 8.0) *
                 u_f_t(2, k, j, i + 2)) /
                pow(h3_f, 2) +
            ((-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI33_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI23_f(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI33_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI23_f(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI33_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI23_f(k, j, i), 2))) /
                  8.0) *
                 u_f_t(3, k, j, i - 2) +
             (Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI33_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI23_f(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI33_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI23_f(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI33_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI23_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI33_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI23_f(k, j, i + 1), 2))) /
                  6.0) *
                 u_f_t(3, k, j, i - 1) +
             (-Jacobian_f(k, j, i - 2) *
                  ((2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                       pow(XI33_f(k, j, i - 2), 2) +
                   mu_f(k, j, i - 2) * (pow(XI13_f(k, j, i - 2), 2) +
                                        pow(XI23_f(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI33_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI23_f(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI33_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI23_f(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI33_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI23_f(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI33_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI23_f(k, j, i + 2), 2))) /
                  24.0) *
                 u_f_t(3, k, j, i) +
             (Jacobian_f(k, j, i - 1) *
                  ((2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                       pow(XI33_f(k, j, i - 1), 2) +
                   mu_f(k, j, i - 1) * (pow(XI13_f(k, j, i - 1), 2) +
                                        pow(XI23_f(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI33_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI23_f(k, j, i), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI33_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI23_f(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI33_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI23_f(k, j, i + 2), 2))) /
                  6.0) *
                 u_f_t(3, k, j, i + 1) +
             (-Jacobian_f(k, j, i) *
                  ((2.0 * mu_f(k, j, i) + lambda_f(k, j, i)) *
                       pow(XI33_f(k, j, i), 2) +
                   mu_f(k, j, i) *
                       (pow(XI13_f(k, j, i), 2) + pow(XI23_f(k, j, i), 2))) /
                  8.0 +
              Jacobian_f(k, j, i + 1) *
                  ((2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                       pow(XI33_f(k, j, i + 1), 2) +
                   mu_f(k, j, i + 1) * (pow(XI13_f(k, j, i + 1), 2) +
                                        pow(XI23_f(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_f(k, j, i + 2) *
                  ((2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                       pow(XI33_f(k, j, i + 2), 2) +
                   mu_f(k, j, i + 2) * (pow(XI13_f(k, j, i + 2), 2) +
                                        pow(XI23_f(k, j, i + 2), 2))) /
                  8.0) *
                 u_f_t(3, k, j, i + 2)) /
                pow(h3_f, 2);
      }
    }
  }
  for (j = 1; j <= n2_f; j++) {
    for (k = 1; k <= n1_f; k++) {
      for (i = 1; i <= 6; i++) {
        for (k1 = 1; k1 <= 8; k1++) {
          for (m = 1; m <= 8; m++) {
            // second derivative 33
            // first set equation
            lh_f(k, j, i, 1) =
                lh_f(k, j, i, 1) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                      pow(XI13_f(k, j, m), 2) +
                  mu_f(k, j, m) *
                      (pow(XI23_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
                 u_f_t(1, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
                 XI23_f(k, j, m) * u_f_t(2, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
                 XI33_f(k, j, m) * u_f_t(3, k, j, k1)) /
                    pow(h3_f, 2);

            lh_f(k, j, n3_f + 1 - i, 1) =
                lh_f(k, j, n3_f + 1 - i, 1) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 ((2.0 * mu_f(k, j, n3_f + 1 - m) +
                   lambda_f(k, j, n3_f + 1 - m)) *
                      pow(XI13_f(k, j, n3_f + 1 - m), 2) +
                  mu_f(k, j, n3_f + 1 - m) *
                      (pow(XI23_f(k, j, n3_f + 1 - m), 2) +
                       pow(XI33_f(k, j, n3_f + 1 - m), 2))) *
                 u_f_t(1, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI13_f(k, j, n3_f + 1 - m) * XI23_f(k, j, n3_f + 1 - m) *
                 u_f_t(2, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI13_f(k, j, n3_f + 1 - m) * XI33_f(k, j, n3_f + 1 - m) *
                 u_f_t(3, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2);
            // second set equation
            lh_f(k, j, i, 2) =
                lh_f(k, j, i, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
                 XI23_f(k, j, m) * u_f_t(1, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                      pow(XI23_f(k, j, m), 2) +
                  mu_f(k, j, m) *
                      (pow(XI13_f(k, j, m), 2) + pow(XI33_f(k, j, m), 2))) *
                 u_f_t(2, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
                 XI33_f(k, j, m) * u_f_t(3, k, j, k1)) /
                    pow(h3_f, 2);

            lh_f(k, j, n3_f + 1 - i, 2) =
                lh_f(k, j, n3_f + 1 - i, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI13_f(k, j, n3_f + 1 - m) * XI23_f(k, j, n3_f + 1 - m) *
                 u_f_t(1, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 ((2.0 * mu_f(k, j, n3_f + 1 - m) +
                   lambda_f(k, j, n3_f + 1 - m)) *
                      pow(XI23_f(k, j, n3_f + 1 - m), 2) +
                  mu_f(k, j, n3_f + 1 - m) *
                      (pow(XI13_f(k, j, n3_f + 1 - m), 2) +
                       pow(XI33_f(k, j, n3_f + 1 - m), 2))) *
                 u_f_t(2, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI23_f(k, j, n3_f + 1 - m) * XI33_f(k, j, n3_f + 1 - m) *
                 u_f_t(3, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2);
            // third set equation
            lh_f(k, j, i, 3) =
                lh_f(k, j, i, 3) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI13_f(k, j, m) *
                 XI33_f(k, j, m) * u_f_t(1, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 (mu_f(k, j, m) + lambda_f(k, j, m)) * XI23_f(k, j, m) *
                 XI33_f(k, j, m) * u_f_t(2, k, j, k1)) /
                    pow(h3_f, 2) +
                (acof_no_gp(i, k1, m) * Jacobian_f(k, j, m) *
                 ((2.0 * mu_f(k, j, m) + lambda_f(k, j, m)) *
                      pow(XI33_f(k, j, m), 2) +
                  mu_f(k, j, m) *
                      (pow(XI13_f(k, j, m), 2) + pow(XI23_f(k, j, m), 2))) *
                 u_f_t(3, k, j, k1)) /
                    pow(h3_f, 2);

            lh_f(k, j, n3_f + 1 - i, 3) =
                lh_f(k, j, n3_f + 1 - i, 3) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI13_f(k, j, n3_f + 1 - m) * XI33_f(k, j, n3_f + 1 - m) *
                 u_f_t(1, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 (mu_f(k, j, n3_f + 1 - m) + lambda_f(k, j, n3_f + 1 - m)) *
                 XI23_f(k, j, n3_f + 1 - m) * XI33_f(k, j, n3_f + 1 - m) *
                 u_f_t(2, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2) +
                (acof(i, k1, m) * Jacobian_f(k, j, n3_f + 1 - m) *
                 ((2.0 * mu_f(k, j, n3_f + 1 - m) +
                   lambda_f(k, j, n3_f + 1 - m)) *
                      pow(XI33_f(k, j, n3_f + 1 - m), 2) +
                  mu_f(k, j, n3_f + 1 - m) *
                      (pow(XI13_f(k, j, n3_f + 1 - m), 2) +
                       pow(XI23_f(k, j, n3_f + 1 - m), 2))) *
                 u_f_t(3, k, j, n3_f + 1 - k1)) /
                    pow(h3_f, 2);
          }
        }
      }
      // first set equation
      lh_f(k, j, n3_f, 1) =
          lh_f(k, j, n3_f, 1) +
          (u_f_t(1, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           ((2.0 * mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) *
                pow(XI13_f(k, j, n3_f), 2) +
            mu_f(k, j, n3_f) *
                (pow(XI23_f(k, j, n3_f), 2) + pow(XI33_f(k, j, n3_f), 2)))) /
              pow(h3_f, 2) +
          (u_f_t(2, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI13_f(k, j, n3_f) *
           XI23_f(k, j, n3_f)) /
              pow(h3_f, 2) +
          (u_f_t(3, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI13_f(k, j, n3_f) *
           XI33_f(k, j, n3_f)) /
              pow(h3_f, 2);
      // second set equation
      lh_f(k, j, n3_f, 2) =
          lh_f(k, j, n3_f, 2) +
          (u_f_t(1, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI13_f(k, j, n3_f) *
           XI23_f(k, j, n3_f)) /
              pow(h3_f, 2) +
          (u_f_t(2, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           ((2.0 * mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) *
                pow(XI23_f(k, j, n3_f), 2) +
            mu_f(k, j, n3_f) *
                (pow(XI13_f(k, j, n3_f), 2) + pow(XI33_f(k, j, n3_f), 2)))) /
              pow(h3_f, 2) +
          (u_f_t(3, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI23_f(k, j, n3_f) *
           XI33_f(k, j, n3_f)) /
              pow(h3_f, 2);
      // third set equation
      lh_f(k, j, n3_f, 3) =
          lh_f(k, j, n3_f, 3) +
          (u_f_t(1, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI13_f(k, j, n3_f) *
           XI33_f(k, j, n3_f)) /
              pow(h3_f, 2) +
          (u_f_t(2, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           (mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) * XI23_f(k, j, n3_f) *
           XI33_f(k, j, n3_f)) /
              pow(h3_f, 2) +
          (u_f_t(3, k, j, n3_f + 1) * ghcof(1) * Jacobian_f(k, j, n3_f) *
           ((2.0 * mu_f(k, j, n3_f) + lambda_f(k, j, n3_f)) *
                pow(XI33_f(k, j, n3_f), 2) +
            mu_f(k, j, n3_f) *
                (pow(XI13_f(k, j, n3_f), 2) + pow(XI23_f(k, j, n3_f), 2)))) /
              pow(h3_f, 2);
    }
  }
  // mixed derivative 13
  for (k = 1; k <= 4; k++) {
    for (i = 1; i <= n2_f; i++) {
      for (j = 1; j <= n1_f; j++) {
        for (k1 = 1; k1 <= 6; k1++) {
          // mixed derivative 13  23
          // first set equation
          lh_f(j, i, k, 1) =
              lh_f(j, i, k, 1) +
              (Jacobian_f(j - 2, i, k) *
                   (2.0 * mu_f(j - 2, i, k) + lambda_f(j - 2, i, k)) *
                   XI13_f(j - 2, i, k) * bof(k, k1) * u_f_t(1, j - 2, i, k1) /
                   12.0 -
               Jacobian_f(j - 1, i, k) *
                   (2.0 * mu_f(j - 1, i, k) + lambda_f(j - 1, i, k)) *
                   XI13_f(j - 1, i, k) * bof(k, k1) * u_f_t(1, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j + 1, i, k) *
                   (2.0 * mu_f(j + 1, i, k) + lambda_f(j + 1, i, k)) *
                   XI13_f(j + 1, i, k) * bof(k, k1) * u_f_t(1, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j + 2, i, k) *
                   (2.0 * mu_f(j + 2, i, k) + lambda_f(j + 2, i, k)) *
                   XI13_f(j + 2, i, k) * bof(k, k1) * u_f_t(1, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_f / h3_f +
              (Jacobian_f(j - 2, i, k) * lambda_f(j - 2, i, k) *
                   XI23_f(j - 2, i, k) * bof(k, k1) * u_f_t(2, j - 2, i, k1) /
                   12.0 -
               Jacobian_f(j - 1, i, k) * lambda_f(j - 1, i, k) *
                   XI23_f(j - 1, i, k) * bof(k, k1) * u_f_t(2, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j + 1, i, k) * lambda_f(j + 1, i, k) *
                   XI23_f(j + 1, i, k) * bof(k, k1) * u_f_t(2, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j + 2, i, k) * lambda_f(j + 2, i, k) *
                   XI23_f(j + 2, i, k) * bof(k, k1) * u_f_t(2, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_f / h3_f +
              (Jacobian_f(j - 2, i, k) * lambda_f(j - 2, i, k) *
                   XI33_f(j - 2, i, k) * bof(k, k1) * u_f_t(3, j - 2, i, k1) /
                   12.0 -
               Jacobian_f(j - 1, i, k) * lambda_f(j - 1, i, k) *
                   XI33_f(j - 1, i, k) * bof(k, k1) * u_f_t(3, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j + 1, i, k) * lambda_f(j + 1, i, k) *
                   XI33_f(j + 1, i, k) * bof(k, k1) * u_f_t(3, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j + 2, i, k) * lambda_f(j + 2, i, k) *
                   XI33_f(j + 2, i, k) * bof(k, k1) * u_f_t(3, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_f / h3_f +
              (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) *
                   XI23_f(j, i - 2, k) * bof(k, k1) * u_f_t(1, j, i - 2, k1) /
                   12.0 -
               Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) *
                   XI23_f(j, i - 1, k) * bof(k, k1) * u_f_t(1, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) *
                   XI23_f(j, i + 1, k) * bof(k, k1) * u_f_t(1, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) *
                   XI23_f(j, i + 2, k) * bof(k, k1) * u_f_t(1, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_f / h3_f +
              (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) *
                   XI13_f(j, i - 2, k) * bof(k, k1) * u_f_t(2, j, i - 2, k1) /
                   12.0 -
               Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) *
                   XI13_f(j, i - 1, k) * bof(k, k1) * u_f_t(2, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) *
                   XI13_f(j, i + 1, k) * bof(k, k1) * u_f_t(2, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) *
                   XI13_f(j, i + 2, k) * bof(k, k1) * u_f_t(2, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_f / h3_f;

          lh_f(j, i, n3_f + 1 - k, 1) =
              lh_f(j, i, n3_f + 1 - k, 1) +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   (2.0 * mu_f(j - 2, i, n3_f + 1 - k) +
                    lambda_f(j - 2, i, n3_f + 1 - k)) *
                   XI13_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   (2.0 * mu_f(j - 1, i, n3_f + 1 - k) +
                    lambda_f(j - 1, i, n3_f + 1 - k)) *
                   XI13_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   (2.0 * mu_f(j + 1, i, n3_f + 1 - k) +
                    lambda_f(j + 1, i, n3_f + 1 - k)) *
                   XI13_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   (2.0 * mu_f(j + 2, i, n3_f + 1 - k) +
                    lambda_f(j + 2, i, n3_f + 1 - k)) *
                   XI13_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   lambda_f(j - 2, i, n3_f + 1 - k) *
                   XI23_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   lambda_f(j - 1, i, n3_f + 1 - k) *
                   XI23_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   lambda_f(j + 1, i, n3_f + 1 - k) *
                   XI23_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   lambda_f(j + 2, i, n3_f + 1 - k) *
                   XI23_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   lambda_f(j - 2, i, n3_f + 1 - k) *
                   XI33_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   lambda_f(j - 1, i, n3_f + 1 - k) *
                   XI33_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   lambda_f(j + 1, i, n3_f + 1 - k) *
                   XI33_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   lambda_f(j + 2, i, n3_f + 1 - k) *
                   XI33_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   mu_f(j, i - 2, n3_f + 1 - k) *
                   XI23_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   mu_f(j, i - 1, n3_f + 1 - k) *
                   XI23_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   mu_f(j, i + 1, n3_f + 1 - k) *
                   XI23_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   mu_f(j, i + 2, n3_f + 1 - k) *
                   XI23_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   mu_f(j, i - 2, n3_f + 1 - k) *
                   XI13_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   mu_f(j, i - 1, n3_f + 1 - k) *
                   XI13_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   mu_f(j, i + 1, n3_f + 1 - k) *
                   XI13_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   mu_f(j, i + 2, n3_f + 1 - k) *
                   XI13_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f;
          // second set equation
          lh_f(j, i, k, 2) =
              lh_f(j, i, k, 2) +
              (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) *
                   XI23_f(j - 2, i, k) * bof(k, k1) * u_f_t(1, j - 2, i, k1) /
                   12.0 -
               Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) *
                   XI23_f(j - 1, i, k) * bof(k, k1) * u_f_t(1, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) *
                   XI23_f(j + 1, i, k) * bof(k, k1) * u_f_t(1, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) *
                   XI23_f(j + 2, i, k) * bof(k, k1) * u_f_t(1, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_f / h3_f +
              (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) *
                   XI13_f(j - 2, i, k) * bof(k, k1) * u_f_t(2, j - 2, i, k1) /
                   12.0 -
               Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) *
                   XI13_f(j - 1, i, k) * bof(k, k1) * u_f_t(2, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) *
                   XI13_f(j + 1, i, k) * bof(k, k1) * u_f_t(2, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) *
                   XI13_f(j + 2, i, k) * bof(k, k1) * u_f_t(2, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_f / h3_f +
              (Jacobian_f(j, i - 2, k) * lambda_f(j, i - 2, k) *
                   XI13_f(j, i - 2, k) * bof(k, k1) * u_f_t(1, j, i - 2, k1) /
                   12.0 -
               Jacobian_f(j, i - 1, k) * lambda_f(j, i - 1, k) *
                   XI13_f(j, i - 1, k) * bof(k, k1) * u_f_t(1, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j, i + 1, k) * lambda_f(j, i + 1, k) *
                   XI13_f(j, i + 1, k) * bof(k, k1) * u_f_t(1, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j, i + 2, k) * lambda_f(j, i + 2, k) *
                   XI13_f(j, i + 2, k) * bof(k, k1) * u_f_t(1, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_f / h3_f +
              (Jacobian_f(j, i - 2, k) *
                   (2.0 * mu_f(j, i - 2, k) + lambda_f(j, i - 2, k)) *
                   XI23_f(j, i - 2, k) * bof(k, k1) * u_f_t(2, j, i - 2, k1) /
                   12.0 -
               Jacobian_f(j, i - 1, k) *
                   (2.0 * mu_f(j, i - 1, k) + lambda_f(j, i - 1, k)) *
                   XI23_f(j, i - 1, k) * bof(k, k1) * u_f_t(2, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j, i + 1, k) *
                   (2.0 * mu_f(j, i + 1, k) + lambda_f(j, i + 1, k)) *
                   XI23_f(j, i + 1, k) * bof(k, k1) * u_f_t(2, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j, i + 2, k) *
                   (2.0 * mu_f(j, i + 2, k) + lambda_f(j, i + 2, k)) *
                   XI23_f(j, i + 2, k) * bof(k, k1) * u_f_t(2, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_f / h3_f +
              (Jacobian_f(j, i - 2, k) * lambda_f(j, i - 2, k) *
                   XI33_f(j, i - 2, k) * bof(k, k1) * u_f_t(3, j, i - 2, k1) /
                   12.0 -
               Jacobian_f(j, i - 1, k) * lambda_f(j, i - 1, k) *
                   XI33_f(j, i - 1, k) * bof(k, k1) * u_f_t(3, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_f(j, i + 1, k) * lambda_f(j, i + 1, k) *
                   XI33_f(j, i + 1, k) * bof(k, k1) * u_f_t(3, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_f(j, i + 2, k) * lambda_f(j, i + 2, k) *
                   XI33_f(j, i + 2, k) * bof(k, k1) * u_f_t(3, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_f / h3_f;

          lh_f(j, i, n3_f + 1 - k, 2) =
              lh_f(j, i, n3_f + 1 - k, 2) +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   mu_f(j - 2, i, n3_f + 1 - k) *
                   XI23_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   mu_f(j - 1, i, n3_f + 1 - k) *
                   XI23_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   mu_f(j + 1, i, n3_f + 1 - k) *
                   XI23_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   mu_f(j + 2, i, n3_f + 1 - k) *
                   XI23_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   mu_f(j - 2, i, n3_f + 1 - k) *
                   XI13_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   mu_f(j - 1, i, n3_f + 1 - k) *
                   XI13_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   mu_f(j + 1, i, n3_f + 1 - k) *
                   XI13_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   mu_f(j + 2, i, n3_f + 1 - k) *
                   XI13_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   lambda_f(j, i - 2, n3_f + 1 - k) *
                   XI13_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   lambda_f(j, i - 1, n3_f + 1 - k) *
                   XI13_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   lambda_f(j, i + 1, n3_f + 1 - k) *
                   XI13_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   lambda_f(j, i + 2, n3_f + 1 - k) *
                   XI13_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   (2.0 * mu_f(j, i - 2, n3_f + 1 - k) +
                    lambda_f(j, i - 2, n3_f + 1 - k)) *
                   XI23_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   (2.0 * mu_f(j, i - 1, n3_f + 1 - k) +
                    lambda_f(j, i - 1, n3_f + 1 - k)) *
                   XI23_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   (2.0 * mu_f(j, i + 1, n3_f + 1 - k) +
                    lambda_f(j, i + 1, n3_f + 1 - k)) *
                   XI23_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   (2.0 * mu_f(j, i + 2, n3_f + 1 - k) +
                    lambda_f(j, i + 2, n3_f + 1 - k)) *
                   XI23_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   lambda_f(j, i - 2, n3_f + 1 - k) *
                   XI33_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   lambda_f(j, i - 1, n3_f + 1 - k) *
                   XI33_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   lambda_f(j, i + 1, n3_f + 1 - k) *
                   XI33_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   lambda_f(j, i + 2, n3_f + 1 - k) *
                   XI33_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f;
          // third set equation
          lh_f(j, i, k, 3) = lh_f(j, i, k, 3) +
                             (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) *
                                  XI33_f(j - 2, i, k) * bof(k, k1) *
                                  u_f_t(1, j - 2, i, k1) / 12.0 -
                              Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) *
                                  XI33_f(j - 1, i, k) * bof(k, k1) *
                                  u_f_t(1, j - 1, i, k1) * 2.0 / 3.0 +
                              Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) *
                                  XI33_f(j + 1, i, k) * bof(k, k1) *
                                  u_f_t(1, j + 1, i, k1) * 2.0 / 3.0 -
                              Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) *
                                  XI33_f(j + 2, i, k) * bof(k, k1) *
                                  u_f_t(1, j + 2, i, k1) / 12.0) /
                                 l1 / h1_f / h3_f +
                             (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) *
                                  XI13_f(j - 2, i, k) * bof(k, k1) *
                                  u_f_t(3, j - 2, i, k1) / 12.0 -
                              Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) *
                                  XI13_f(j - 1, i, k) * bof(k, k1) *
                                  u_f_t(3, j - 1, i, k1) * 2.0 / 3.0 +
                              Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) *
                                  XI13_f(j + 1, i, k) * bof(k, k1) *
                                  u_f_t(3, j + 1, i, k1) * 2.0 / 3.0 -
                              Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) *
                                  XI13_f(j + 2, i, k) * bof(k, k1) *
                                  u_f_t(3, j + 2, i, k1) / 12.0) /
                                 l1 / h1_f / h3_f +
                             (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) *
                                  XI33_f(j, i - 2, k) * bof(k, k1) *
                                  u_f_t(2, j, i - 2, k1) / 12.0 -
                              Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) *
                                  XI33_f(j, i - 1, k) * bof(k, k1) *
                                  u_f_t(2, j, i - 1, k1) * 2.0 / 3.0 +
                              Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) *
                                  XI33_f(j, i + 1, k) * bof(k, k1) *
                                  u_f_t(2, j, i + 1, k1) * 2.0 / 3.0 -
                              Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) *
                                  XI33_f(j, i + 2, k) * bof(k, k1) *
                                  u_f_t(2, j, i + 2, k1) / 12.0) /
                                 l2 / h2_f / h3_f +
                             (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) *
                                  XI23_f(j, i - 2, k) * bof(k, k1) *
                                  u_f_t(3, j, i - 2, k1) / 12.0 -
                              Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) *
                                  XI23_f(j, i - 1, k) * bof(k, k1) *
                                  u_f_t(3, j, i - 1, k1) * 2.0 / 3.0 +
                              Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) *
                                  XI23_f(j, i + 1, k) * bof(k, k1) *
                                  u_f_t(3, j, i + 1, k1) * 2.0 / 3.0 -
                              Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) *
                                  XI23_f(j, i + 2, k) * bof(k, k1) *
                                  u_f_t(3, j, i + 2, k1) / 12.0) /
                                 l2 / h2_f / h3_f;

          lh_f(j, i, n3_f + 1 - k, 3) =
              lh_f(j, i, n3_f + 1 - k, 3) +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   mu_f(j - 2, i, n3_f + 1 - k) *
                   XI33_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   mu_f(j - 1, i, n3_f + 1 - k) *
                   XI33_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   mu_f(j + 1, i, n3_f + 1 - k) *
                   XI33_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   mu_f(j + 2, i, n3_f + 1 - k) *
                   XI33_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(1, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j - 2, i, n3_f + 1 - k) *
                   mu_f(j - 2, i, n3_f + 1 - k) *
                   XI13_f(j - 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j - 2, i, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j - 1, i, n3_f + 1 - k) *
                   mu_f(j - 1, i, n3_f + 1 - k) *
                   XI13_f(j - 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j - 1, i, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j + 1, i, n3_f + 1 - k) *
                   mu_f(j + 1, i, n3_f + 1 - k) *
                   XI13_f(j + 1, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j + 1, i, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j + 2, i, n3_f + 1 - k) *
                   mu_f(j + 2, i, n3_f + 1 - k) *
                   XI13_f(j + 2, i, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j + 2, i, n3_f + 1 - k1) / 12.0) /
                  l1 / h1_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   mu_f(j, i - 2, n3_f + 1 - k) *
                   XI33_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   mu_f(j, i - 1, n3_f + 1 - k) *
                   XI33_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   mu_f(j, i + 1, n3_f + 1 - k) *
                   XI33_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   mu_f(j, i + 2, n3_f + 1 - k) *
                   XI33_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(2, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f +
              (-Jacobian_f(j, i - 2, n3_f + 1 - k) *
                   mu_f(j, i - 2, n3_f + 1 - k) *
                   XI23_f(j, i - 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i - 2, n3_f + 1 - k1) / 12.0 +
               Jacobian_f(j, i - 1, n3_f + 1 - k) *
                   mu_f(j, i - 1, n3_f + 1 - k) *
                   XI23_f(j, i - 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i - 1, n3_f + 1 - k1) * 2.0 / 3.0 -
               Jacobian_f(j, i + 1, n3_f + 1 - k) *
                   mu_f(j, i + 1, n3_f + 1 - k) *
                   XI23_f(j, i + 1, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i + 1, n3_f + 1 - k1) * 2.0 / 3.0 +
               Jacobian_f(j, i + 2, n3_f + 1 - k) *
                   mu_f(j, i + 2, n3_f + 1 - k) *
                   XI23_f(j, i + 2, n3_f + 1 - k) * bof(k, k1) *
                   u_f_t(3, j, i + 2, n3_f + 1 - k1) / 12.0) /
                  l2 / h2_f / h3_f;
        }
      }
    }
  }
  for (k = 5; k <= n3_f - 4; k++) {
    for (i = 1; i <= n2_f; i++) {
      for (j = 1; j <= n1_f; j++) {
        // mixed derivative 13  23
        // first set equation
        lh_f(j, i, k, 1) =
            lh_f(j, i, k, 1) +
            (Jacobian_f(j - 2, i, k) *
                 (2.0 * mu_f(j - 2, i, k) + lambda_f(j - 2, i, k)) *
                 XI13_f(j - 2, i, k) *
                 (u_f_t(1, j - 2, i, k - 2) / 12.0 -
                  u_f_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) *
                 (2.0 * mu_f(j - 1, i, k) + lambda_f(j - 1, i, k)) *
                 XI13_f(j - 1, i, k) *
                 (u_f_t(1, j - 1, i, k - 2) / 12.0 -
                  u_f_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) *
                 (2.0 * mu_f(j + 1, i, k) + lambda_f(j + 1, i, k)) *
                 XI13_f(j + 1, i, k) *
                 (u_f_t(1, j + 1, i, k - 2) / 12.0 -
                  u_f_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) *
                 (2.0 * mu_f(j + 2, i, k) + lambda_f(j + 2, i, k)) *
                 XI13_f(j + 2, i, k) *
                 (u_f_t(1, j + 2, i, k - 2) / 12.0 -
                  u_f_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, k) * lambda_f(j - 2, i, k) *
                 XI23_f(j - 2, i, k) *
                 (u_f_t(2, j - 2, i, k - 2) / 12.0 -
                  u_f_t(2, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * lambda_f(j - 1, i, k) *
                 XI23_f(j - 1, i, k) *
                 (u_f_t(2, j - 1, i, k - 2) / 12.0 -
                  u_f_t(2, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * lambda_f(j + 1, i, k) *
                 XI23_f(j + 1, i, k) *
                 (u_f_t(2, j + 1, i, k - 2) / 12.0 -
                  u_f_t(2, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * lambda_f(j + 2, i, k) *
                 XI23_f(j + 2, i, k) *
                 (u_f_t(2, j + 2, i, k - 2) / 12.0 -
                  u_f_t(2, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, k) * lambda_f(j - 2, i, k) *
                 XI33_f(j - 2, i, k) *
                 (u_f_t(3, j - 2, i, k - 2) / 12.0 -
                  u_f_t(3, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * lambda_f(j - 1, i, k) *
                 XI33_f(j - 1, i, k) *
                 (u_f_t(3, j - 1, i, k - 2) / 12.0 -
                  u_f_t(3, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * lambda_f(j + 1, i, k) *
                 XI33_f(j + 1, i, k) *
                 (u_f_t(3, j + 1, i, k - 2) / 12.0 -
                  u_f_t(3, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * lambda_f(j + 2, i, k) *
                 XI33_f(j + 2, i, k) *
                 (u_f_t(3, j + 2, i, k - 2) / 12.0 -
                  u_f_t(3, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) * XI23_f(j, i - 2, k) *
                 (u_f_t(1, j, i - 2, k - 2) / 12.0 -
                  u_f_t(1, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) * XI23_f(j, i - 1, k) *
                 (u_f_t(1, j, i - 1, k - 2) / 12.0 -
                  u_f_t(1, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) * XI23_f(j, i + 1, k) *
                 (u_f_t(1, j, i + 1, k - 2) / 12.0 -
                  u_f_t(1, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) * XI23_f(j, i + 2, k) *
                 (u_f_t(1, j, i + 2, k - 2) / 12.0 -
                  u_f_t(1, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) * XI13_f(j, i - 2, k) *
                 (u_f_t(2, j, i - 2, k - 2) / 12.0 -
                  u_f_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) * XI13_f(j, i - 1, k) *
                 (u_f_t(2, j, i - 1, k - 2) / 12.0 -
                  u_f_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) * XI13_f(j, i + 1, k) *
                 (u_f_t(2, j, i + 1, k - 2) / 12.0 -
                  u_f_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) * XI13_f(j, i + 2, k) *
                 (u_f_t(2, j, i + 2, k - 2) / 12.0 -
                  u_f_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f;
        // second set equation
        lh_f(j, i, k, 2) =
            lh_f(j, i, k, 2) +
            (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) * XI23_f(j - 2, i, k) *
                 (u_f_t(1, j - 2, i, k - 2) / 12.0 -
                  u_f_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) * XI23_f(j - 1, i, k) *
                 (u_f_t(1, j - 1, i, k - 2) / 12.0 -
                  u_f_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) * XI23_f(j + 1, i, k) *
                 (u_f_t(1, j + 1, i, k - 2) / 12.0 -
                  u_f_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) * XI23_f(j + 2, i, k) *
                 (u_f_t(1, j + 2, i, k - 2) / 12.0 -
                  u_f_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) * XI13_f(j - 2, i, k) *
                 (u_f_t(2, j - 2, i, k - 2) / 12.0 -
                  u_f_t(2, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) * XI13_f(j - 1, i, k) *
                 (u_f_t(2, j - 1, i, k - 2) / 12.0 -
                  u_f_t(2, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) * XI13_f(j + 1, i, k) *
                 (u_f_t(2, j + 1, i, k - 2) / 12.0 -
                  u_f_t(2, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) * XI13_f(j + 2, i, k) *
                 (u_f_t(2, j + 2, i, k - 2) / 12.0 -
                  u_f_t(2, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, k) * lambda_f(j, i - 2, k) *
                 XI13_f(j, i - 2, k) *
                 (u_f_t(1, j, i - 2, k - 2) / 12.0 -
                  u_f_t(1, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * lambda_f(j, i - 1, k) *
                 XI13_f(j, i - 1, k) *
                 (u_f_t(1, j, i - 1, k - 2) / 12.0 -
                  u_f_t(1, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * lambda_f(j, i + 1, k) *
                 XI13_f(j, i + 1, k) *
                 (u_f_t(1, j, i + 1, k - 2) / 12.0 -
                  u_f_t(1, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * lambda_f(j, i + 2, k) *
                 XI13_f(j, i + 2, k) *
                 (u_f_t(1, j, i + 2, k - 2) / 12.0 -
                  u_f_t(1, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, k) *
                 (2.0 * mu_f(j, i - 2, k) + lambda_f(j, i - 2, k)) *
                 XI23_f(j, i - 2, k) *
                 (u_f_t(2, j, i - 2, k - 2) / 12.0 -
                  u_f_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) *
                 (2.0 * mu_f(j, i - 1, k) + lambda_f(j, i - 1, k)) *
                 XI23_f(j, i - 1, k) *
                 (u_f_t(2, j, i - 1, k - 2) / 12.0 -
                  u_f_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) *
                 (2.0 * mu_f(j, i + 1, k) + lambda_f(j, i + 1, k)) *
                 XI23_f(j, i + 1, k) *
                 (u_f_t(2, j, i + 1, k - 2) / 12.0 -
                  u_f_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) *
                 (2.0 * mu_f(j, i + 2, k) + lambda_f(j, i + 2, k)) *
                 XI23_f(j, i + 2, k) *
                 (u_f_t(2, j, i + 2, k - 2) / 12.0 -
                  u_f_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, k) * lambda_f(j, i - 2, k) *
                 XI33_f(j, i - 2, k) *
                 (u_f_t(3, j, i - 2, k - 2) / 12.0 -
                  u_f_t(3, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * lambda_f(j, i - 1, k) *
                 XI33_f(j, i - 1, k) *
                 (u_f_t(3, j, i - 1, k - 2) / 12.0 -
                  u_f_t(3, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * lambda_f(j, i + 1, k) *
                 XI33_f(j, i + 1, k) *
                 (u_f_t(3, j, i + 1, k - 2) / 12.0 -
                  u_f_t(3, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * lambda_f(j, i + 2, k) *
                 XI33_f(j, i + 2, k) *
                 (u_f_t(3, j, i + 2, k - 2) / 12.0 -
                  u_f_t(3, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f;
        // third set equation
        lh_f(j, i, k, 3) =
            lh_f(j, i, k, 3) +
            (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) * XI33_f(j - 2, i, k) *
                 (u_f_t(1, j - 2, i, k - 2) / 12.0 -
                  u_f_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) * XI33_f(j - 1, i, k) *
                 (u_f_t(1, j - 1, i, k - 2) / 12.0 -
                  u_f_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) * XI33_f(j + 1, i, k) *
                 (u_f_t(1, j + 1, i, k - 2) / 12.0 -
                  u_f_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) * XI33_f(j + 2, i, k) *
                 (u_f_t(1, j + 2, i, k - 2) / 12.0 -
                  u_f_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j - 2, i, k) * mu_f(j - 2, i, k) * XI13_f(j - 2, i, k) *
                 (u_f_t(3, j - 2, i, k - 2) / 12.0 -
                  u_f_t(3, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j - 1, i, k) * mu_f(j - 1, i, k) * XI13_f(j - 1, i, k) *
                 (u_f_t(3, j - 1, i, k - 2) / 12.0 -
                  u_f_t(3, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j + 1, i, k) * mu_f(j + 1, i, k) * XI13_f(j + 1, i, k) *
                 (u_f_t(3, j + 1, i, k - 2) / 12.0 -
                  u_f_t(3, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j + 2, i, k) * mu_f(j + 2, i, k) * XI13_f(j + 2, i, k) *
                 (u_f_t(3, j + 2, i, k - 2) / 12.0 -
                  u_f_t(3, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) * XI33_f(j, i - 2, k) *
                 (u_f_t(2, j, i - 2, k - 2) / 12.0 -
                  u_f_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) * XI33_f(j, i - 1, k) *
                 (u_f_t(2, j, i - 1, k - 2) / 12.0 -
                  u_f_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) * XI33_f(j, i + 1, k) *
                 (u_f_t(2, j, i + 1, k - 2) / 12.0 -
                  u_f_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) * XI33_f(j, i + 2, k) *
                 (u_f_t(2, j, i + 2, k - 2) / 12.0 -
                  u_f_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(j, i - 2, k) * mu_f(j, i - 2, k) * XI23_f(j, i - 2, k) *
                 (u_f_t(3, j, i - 2, k - 2) / 12.0 -
                  u_f_t(3, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_f(j, i - 1, k) * mu_f(j, i - 1, k) * XI23_f(j, i - 1, k) *
                 (u_f_t(3, j, i - 1, k - 2) / 12.0 -
                  u_f_t(3, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(j, i + 1, k) * mu_f(j, i + 1, k) * XI23_f(j, i + 1, k) *
                 (u_f_t(3, j, i + 1, k - 2) / 12.0 -
                  u_f_t(3, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(j, i + 2, k) * mu_f(j, i + 2, k) * XI23_f(j, i + 2, k) *
                 (u_f_t(3, j, i + 2, k - 2) / 12.0 -
                  u_f_t(3, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_f_t(3, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_f_t(3, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f;
      }
    }
  }
  //
  for (i = 5; i <= n3_f - 4; i++) {
    for (j = 1; j <= n2_f; j++) {
      for (k = 1; k <= n1_f; k++) {
        // mixed derivative 31  32
        // first set equation
        lh_f(k, j, i, 1) =
            lh_f(k, j, i, 1) +
            (Jacobian_f(k, j, i - 2) *
                 (2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                 XI13_f(k, j, i - 2) *
                 (u_f_t(1, k - 2, j, i - 2) / 12.0 -
                  u_f_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) *
                 (2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                 XI13_f(k, j, i - 1) *
                 (u_f_t(1, k - 2, j, i - 1) / 12.0 -
                  u_f_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) *
                 (2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                 XI13_f(k, j, i + 1) *
                 (u_f_t(1, k - 2, j, i + 1) / 12.0 -
                  u_f_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) *
                 (2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                 XI13_f(k, j, i + 2) *
                 (u_f_t(1, k - 2, j, i + 2) / 12.0 -
                  u_f_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_f / h1_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI23_f(k, j, i - 2) *
                 (u_f_t(2, k - 2, j, i - 2) / 12.0 -
                  u_f_t(2, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI23_f(k, j, i - 1) *
                 (u_f_t(2, k - 2, j, i - 1) / 12.0 -
                  u_f_t(2, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI23_f(k, j, i + 1) *
                 (u_f_t(2, k - 2, j, i + 1) / 12.0 -
                  u_f_t(2, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI23_f(k, j, i + 2) *
                 (u_f_t(2, k - 2, j, i + 2) / 12.0 -
                  u_f_t(2, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_f / h1_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI33_f(k, j, i - 2) *
                 (u_f_t(3, k - 2, j, i - 2) / 12.0 -
                  u_f_t(3, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI33_f(k, j, i - 1) *
                 (u_f_t(3, k - 2, j, i - 1) / 12.0 -
                  u_f_t(3, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI33_f(k, j, i + 1) *
                 (u_f_t(3, k - 2, j, i + 1) / 12.0 -
                  u_f_t(3, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI33_f(k, j, i + 2) *
                 (u_f_t(3, k - 2, j, i + 2) / 12.0 -
                  u_f_t(3, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI23_f(k, j, i - 2) *
                 (u_f_t(1, k, j - 2, i - 2) / 12.0 -
                  u_f_t(1, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI23_f(k, j, i - 1) *
                 (u_f_t(1, k, j - 2, i - 1) / 12.0 -
                  u_f_t(1, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI23_f(k, j, i + 1) *
                 (u_f_t(1, k, j - 2, i + 1) / 12.0 -
                  u_f_t(1, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI23_f(k, j, i + 2) *
                 (u_f_t(1, k, j - 2, i + 2) / 12.0 -
                  u_f_t(1, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_f / h2_f +
            (Jacobian_f(k, j, i - 2) * lambda_f(k, j, i - 2) *
                 XI13_f(k, j, i - 2) *
                 (u_f_t(2, k, j - 2, i - 2) / 12.0 -
                  u_f_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * lambda_f(k, j, i - 1) *
                 XI13_f(k, j, i - 1) *
                 (u_f_t(2, k, j - 2, i - 1) / 12.0 -
                  u_f_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * lambda_f(k, j, i + 1) *
                 XI13_f(k, j, i + 1) *
                 (u_f_t(2, k, j - 2, i + 1) / 12.0 -
                  u_f_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * lambda_f(k, j, i + 2) *
                 XI13_f(k, j, i + 2) *
                 (u_f_t(2, k, j - 2, i + 2) / 12.0 -
                  u_f_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_f / h2_f;
        // second set equation
        lh_f(k, j, i, 2) =
            lh_f(k, j, i, 2) +
            (Jacobian_f(k, j, i - 2) * lambda_f(k, j, i - 2) *
                 XI23_f(k, j, i - 2) *
                 (u_f_t(1, k - 2, j, i - 2) / 12.0 -
                  u_f_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * lambda_f(k, j, i - 1) *
                 XI23_f(k, j, i - 1) *
                 (u_f_t(1, k - 2, j, i - 1) / 12.0 -
                  u_f_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * lambda_f(k, j, i + 1) *
                 XI23_f(k, j, i + 1) *
                 (u_f_t(1, k - 2, j, i + 1) / 12.0 -
                  u_f_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * lambda_f(k, j, i + 2) *
                 XI23_f(k, j, i + 2) *
                 (u_f_t(1, k - 2, j, i + 2) / 12.0 -
                  u_f_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI13_f(k, j, i - 2) *
                 (u_f_t(2, k - 2, j, i - 2) / 12.0 -
                  u_f_t(2, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI13_f(k, j, i - 1) *
                 (u_f_t(2, k - 2, j, i - 1) / 12.0 -
                  u_f_t(2, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI13_f(k, j, i + 1) *
                 (u_f_t(2, k - 2, j, i + 1) / 12.0 -
                  u_f_t(2, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI13_f(k, j, i + 2) *
                 (u_f_t(2, k - 2, j, i + 2) / 12.0 -
                  u_f_t(2, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(2, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(2, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_f / h3_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI13_f(k, j, i - 2) *
                 (u_f_t(1, k, j - 2, i - 2) / 12.0 -
                  u_f_t(1, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI13_f(k, j, i - 1) *
                 (u_f_t(1, k, j - 2, i - 1) / 12.0 -
                  u_f_t(1, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI13_f(k, j, i + 1) *
                 (u_f_t(1, k, j - 2, i + 1) / 12.0 -
                  u_f_t(1, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI13_f(k, j, i + 2) *
                 (u_f_t(1, k, j - 2, i + 2) / 12.0 -
                  u_f_t(1, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(1, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(1, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_f / h2_f +
            (Jacobian_f(k, j, i - 2) *
                 (2.0 * mu_f(k, j, i - 2) + lambda_f(k, j, i - 2)) *
                 XI23_f(k, j, i - 2) *
                 (u_f_t(2, k, j - 2, i - 2) / 12.0 -
                  u_f_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) *
                 (2.0 * mu_f(k, j, i - 1) + lambda_f(k, j, i - 1)) *
                 XI23_f(k, j, i - 1) *
                 (u_f_t(2, k, j - 2, i - 1) / 12.0 -
                  u_f_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) *
                 (2.0 * mu_f(k, j, i + 1) + lambda_f(k, j, i + 1)) *
                 XI23_f(k, j, i + 1) *
                 (u_f_t(2, k, j - 2, i + 1) / 12.0 -
                  u_f_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) *
                 (2.0 * mu_f(k, j, i + 2) + lambda_f(k, j, i + 2)) *
                 XI23_f(k, j, i + 2) *
                 (u_f_t(2, k, j - 2, i + 2) / 12.0 -
                  u_f_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_f / h2_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI33_f(k, j, i - 2) *
                 (u_f_t(3, k, j - 2, i - 2) / 12.0 -
                  u_f_t(3, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI33_f(k, j, i - 1) *
                 (u_f_t(3, k, j - 2, i - 1) / 12.0 -
                  u_f_t(3, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI33_f(k, j, i + 1) *
                 (u_f_t(3, k, j - 2, i + 1) / 12.0 -
                  u_f_t(3, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI33_f(k, j, i + 2) *
                 (u_f_t(3, k, j - 2, i + 2) / 12.0 -
                  u_f_t(3, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_f / h2_f;
        // third set equation
        lh_f(k, j, i, 3) =
            lh_f(k, j, i, 3) +
            (Jacobian_f(k, j, i - 2) * lambda_f(k, j, i - 2) *
                 XI33_f(k, j, i - 2) *
                 (u_f_t(1, k - 2, j, i - 2) / 12.0 -
                  u_f_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * lambda_f(k, j, i - 1) *
                 XI33_f(k, j, i - 1) *
                 (u_f_t(1, k - 2, j, i - 1) / 12.0 -
                  u_f_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * lambda_f(k, j, i + 1) *
                 XI33_f(k, j, i + 1) *
                 (u_f_t(1, k - 2, j, i + 1) / 12.0 -
                  u_f_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * lambda_f(k, j, i + 2) *
                 XI33_f(k, j, i + 2) *
                 (u_f_t(1, k - 2, j, i + 2) / 12.0 -
                  u_f_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_f / h1_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI13_f(k, j, i - 2) *
                 (u_f_t(3, k - 2, j, i - 2) / 12.0 -
                  u_f_t(3, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI13_f(k, j, i - 1) *
                 (u_f_t(3, k - 2, j, i - 1) / 12.0 -
                  u_f_t(3, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI13_f(k, j, i + 1) *
                 (u_f_t(3, k - 2, j, i + 1) / 12.0 -
                  u_f_t(3, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI13_f(k, j, i + 2) *
                 (u_f_t(3, k - 2, j, i + 2) / 12.0 -
                  u_f_t(3, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_f_t(3, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_f_t(3, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_f / h1_f +
            (Jacobian_f(k, j, i - 2) * lambda_f(k, j, i - 2) *
                 XI33_f(k, j, i - 2) *
                 (u_f_t(2, k, j - 2, i - 2) / 12.0 -
                  u_f_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * lambda_f(k, j, i - 1) *
                 XI33_f(k, j, i - 1) *
                 (u_f_t(2, k, j - 2, i - 1) / 12.0 -
                  u_f_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * lambda_f(k, j, i + 1) *
                 XI33_f(k, j, i + 1) *
                 (u_f_t(2, k, j - 2, i + 1) / 12.0 -
                  u_f_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * lambda_f(k, j, i + 2) *
                 XI33_f(k, j, i + 2) *
                 (u_f_t(2, k, j - 2, i + 2) / 12.0 -
                  u_f_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f +
            (Jacobian_f(k, j, i - 2) * mu_f(k, j, i - 2) * XI23_f(k, j, i - 2) *
                 (u_f_t(3, k, j - 2, i - 2) / 12.0 -
                  u_f_t(3, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_f(k, j, i - 1) * mu_f(k, j, i - 1) * XI23_f(k, j, i - 1) *
                 (u_f_t(3, k, j - 2, i - 1) / 12.0 -
                  u_f_t(3, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_f(k, j, i + 1) * mu_f(k, j, i + 1) * XI23_f(k, j, i + 1) *
                 (u_f_t(3, k, j - 2, i + 1) / 12.0 -
                  u_f_t(3, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_f(k, j, i + 2) * mu_f(k, j, i + 2) * XI23_f(k, j, i + 2) *
                 (u_f_t(3, k, j - 2, i + 2) / 12.0 -
                  u_f_t(3, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_f_t(3, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_f_t(3, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h2_f / h3_f;
      }
    }
  }
  for (i = 1; i <= 4; i++) {
    for (j = 1; j <= n2_f; j++) {
      for (k = 1; k <= n1_f; k++) {
        for (k1 = 1; k1 <= 6; k1++) {
          // mixed derivative 31  32
          // first set equation
          lh_f(k, j, i, 1) =
              lh_f(k, j, i, 1) +
              (bof(i, k1) * Jacobian_f(k, j, k1) *
               (2.0 * mu_f(k, j, k1) + lambda_f(k, j, k1)) * XI13_f(k, j, k1) *
               (u_f_t(1, k - 2, j, k1) / 12.0 -
                u_f_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                u_f_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                u_f_t(1, k + 2, j, k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI23_f(k, j, k1) *
               (u_f_t(2, k - 2, j, k1) / 12.0 -
                u_f_t(2, k - 1, j, k1) * 2.0 / 3.0 +
                u_f_t(2, k + 1, j, k1) * 2.0 / 3.0 -
                u_f_t(2, k + 2, j, k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI33_f(k, j, k1) *
               (u_f_t(3, k - 2, j, k1) / 12.0 -
                u_f_t(3, k - 1, j, k1) * 2.0 / 3.0 +
                u_f_t(3, k + 1, j, k1) * 2.0 / 3.0 -
                u_f_t(3, k + 2, j, k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI23_f(k, j, k1) *
               (u_f_t(1, k, j - 2, k1) / 12.0 -
                u_f_t(1, k, j - 1, k1) * 2.0 / 3.0 +
                u_f_t(1, k, j + 1, k1) * 2.0 / 3.0 -
                u_f_t(1, k, j + 2, k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * lambda_f(k, j, k1) *
               XI13_f(k, j, k1) *
               (u_f_t(2, k, j - 2, k1) / 12.0 -
                u_f_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                u_f_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                u_f_t(2, k, j + 2, k1) / 12.0)) /
                  l2 / h3_f / h2_f;

          lh_f(k, j, n3_f + 1 - i, 1) =
              lh_f(k, j, n3_f + 1 - i, 1) +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               (2.0 * mu_f(k, j, n3_f + 1 - k1) +
                lambda_f(k, j, n3_f + 1 - k1)) *
               XI13_f(k, j, n3_f + 1 - k1) *
               (u_f_t(1, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(1, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(1, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(1, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI23_f(k, j, n3_f + 1 - k1) *
               (u_f_t(2, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(2, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(2, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(2, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI33_f(k, j, n3_f + 1 - k1) *
               (u_f_t(3, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(3, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(3, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(3, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI23_f(k, j, n3_f + 1 - k1) *
               (u_f_t(1, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(1, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(1, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(1, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               lambda_f(k, j, n3_f + 1 - k1) * XI13_f(k, j, n3_f + 1 - k1) *
               (u_f_t(2, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(2, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(2, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(2, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h3_f / h2_f;
          // second set equation
          lh_f(k, j, i, 2) =
              lh_f(k, j, i, 2) +
              (bof(i, k1) * Jacobian_f(k, j, k1) * lambda_f(k, j, k1) *
               XI23_f(k, j, k1) *
               (u_f_t(1, k - 2, j, k1) / 12.0 -
                u_f_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                u_f_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                u_f_t(1, k + 2, j, k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI13_f(k, j, k1) *
               (u_f_t(2, k - 2, j, k1) / 12.0 -
                u_f_t(2, k - 1, j, k1) * 2.0 / 3.0 +
                u_f_t(2, k + 1, j, k1) * 2.0 / 3.0 -
                u_f_t(2, k + 2, j, k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI13_f(k, j, k1) *
               (u_f_t(1, k, j - 2, k1) / 12.0 -
                u_f_t(1, k, j - 1, k1) * 2.0 / 3.0 +
                u_f_t(1, k, j + 1, k1) * 2.0 / 3.0 -
                u_f_t(1, k, j + 2, k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) *
               (2.0 * mu_f(k, j, k1) + lambda_f(k, j, k1)) * XI23_f(k, j, k1) *
               (u_f_t(2, k, j - 2, k1) / 12.0 -
                u_f_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                u_f_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                u_f_t(2, k, j + 2, k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (bof(i, k1) * Jacobian_f(k, j, k1) * mu_f(k, j, k1) *
               XI33_f(k, j, k1) *
               (u_f_t(3, k, j - 2, k1) / 12.0 -
                u_f_t(3, k, j - 1, k1) * 2.0 / 3.0 +
                u_f_t(3, k, j + 1, k1) * 2.0 / 3.0 -
                u_f_t(3, k, j + 2, k1) / 12.0)) /
                  l2 / h3_f / h2_f;

          lh_f(k, j, n3_f + 1 - i, 2) =
              lh_f(k, j, n3_f + 1 - i, 2) +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               lambda_f(k, j, n3_f + 1 - k1) * XI23_f(k, j, n3_f + 1 - k1) *
               (u_f_t(1, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(1, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(1, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(1, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI13_f(k, j, n3_f + 1 - k1) *
               (u_f_t(2, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(2, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(2, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(2, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h1_f / h3_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI13_f(k, j, n3_f + 1 - k1) *
               (u_f_t(1, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(1, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(1, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(1, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               (2.0 * mu_f(k, j, n3_f + 1 - k1) +
                lambda_f(k, j, n3_f + 1 - k1)) *
               XI23_f(k, j, n3_f + 1 - k1) *
               (u_f_t(2, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(2, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(2, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(2, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h3_f / h2_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI33_f(k, j, n3_f + 1 - k1) *
               (u_f_t(3, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(3, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(3, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(3, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h3_f / h2_f;
          // third set equation
          lh_f(k, j, i, 3) = lh_f(k, j, i, 3) +
                             (bof(i, k1) * Jacobian_f(k, j, k1) *
                              lambda_f(k, j, k1) * XI33_f(k, j, k1) *
                              (u_f_t(1, k - 2, j, k1) / 12.0 -
                               u_f_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                               u_f_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                               u_f_t(1, k + 2, j, k1) / 12.0)) /
                                 l1 / h1_f / h3_f +
                             (bof(i, k1) * Jacobian_f(k, j, k1) *
                              mu_f(k, j, k1) * XI13_f(k, j, k1) *
                              (u_f_t(3, k - 2, j, k1) / 12.0 -
                               u_f_t(3, k - 1, j, k1) * 2.0 / 3.0 +
                               u_f_t(3, k + 1, j, k1) * 2.0 / 3.0 -
                               u_f_t(3, k + 2, j, k1) / 12.0)) /
                                 l1 / h3_f / h1_f +
                             (bof(i, k1) * Jacobian_f(k, j, k1) *
                              lambda_f(k, j, k1) * XI33_f(k, j, k1) *
                              (u_f_t(2, k, j - 2, k1) / 12.0 -
                               u_f_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                               u_f_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                               u_f_t(2, k, j + 2, k1) / 12.0)) /
                                 l2 / h2_f / h3_f +
                             (bof(i, k1) * Jacobian_f(k, j, k1) *
                              mu_f(k, j, k1) * XI23_f(k, j, k1) *
                              (u_f_t(3, k, j - 2, k1) / 12.0 -
                               u_f_t(3, k, j - 1, k1) * 2.0 / 3.0 +
                               u_f_t(3, k, j + 1, k1) * 2.0 / 3.0 -
                               u_f_t(3, k, j + 2, k1) / 12.0)) /
                                 l2 / h2_f / h3_f;

          lh_f(k, j, n3_f + 1 - i, 3) =
              lh_f(k, j, n3_f + 1 - i, 3) +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               lambda_f(k, j, n3_f + 1 - k1) * XI33_f(k, j, n3_f + 1 - k1) *
               (u_f_t(1, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(1, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(1, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(1, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI13_f(k, j, n3_f + 1 - k1) *
               (u_f_t(3, k - 2, j, n3_f + 1 - k1) / 12.0 -
                u_f_t(3, k - 1, j, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(3, k + 1, j, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(3, k + 2, j, n3_f + 1 - k1) / 12.0)) /
                  l1 / h3_f / h1_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               lambda_f(k, j, n3_f + 1 - k1) * XI33_f(k, j, n3_f + 1 - k1) *
               (u_f_t(2, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(2, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(2, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(2, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h2_f / h3_f +
              (-bof(i, k1) * Jacobian_f(k, j, n3_f + 1 - k1) *
               mu_f(k, j, n3_f + 1 - k1) * XI23_f(k, j, n3_f + 1 - k1) *
               (u_f_t(3, k, j - 2, n3_f + 1 - k1) / 12.0 -
                u_f_t(3, k, j - 1, n3_f + 1 - k1) * 2.0 / 3.0 +
                u_f_t(3, k, j + 1, n3_f + 1 - k1) * 2.0 / 3.0 -
                u_f_t(3, k, j + 2, n3_f + 1 - k1) / 12.0)) /
                  l2 / h2_f / h3_f;
        }
      }
    }
  }
  //
  // coarse mesh
  lh_c = 0.0;
  for (k = 1; k <= n3_c; k++) {
    for (j = 1; j <= n2_c; j++) {
      for (i = 1; i <= n1_c; i++) {
        // second derivative 11  22  12  21
        // first set
        lh_c(i, j, k, 1) =
            lh_c(i, j, k, 1) +
            ((-Jacobian_c(i - 2, j, k) *
                  (2.0 * mu_c(i - 2, j, k) + lambda_c(i - 2, j, k)) / 8.0 +
              Jacobian_c(i - 1, j, k) *
                  (2.0 * mu_c(i - 1, j, k) + lambda_c(i - 1, j, k)) / 6.0 -
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  8.0) *
                 u_c_t(1, i - 2, j, k) +
             (Jacobian_c(i - 2, j, k) *
                  (2.0 * mu_c(i - 2, j, k) + lambda_c(i - 2, j, k)) / 6.0 +
              Jacobian_c(i - 1, j, k) *
                  (2.0 * mu_c(i - 1, j, k) + lambda_c(i - 1, j, k)) / 2.0 +
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  2.0 +
              Jacobian_c(i + 1, j, k) *
                  (2.0 * mu_c(i + 1, j, k) + lambda_c(i + 1, j, k)) / 6.0) *
                 u_c_t(1, i - 1, j, k) +
             (-Jacobian_c(i - 2, j, k) *
                  (2.0 * mu_c(i - 2, j, k) + lambda_c(i - 2, j, k)) / 24.0 -
              Jacobian_c(i - 1, j, k) *
                  (2.0 * mu_c(i - 1, j, k) + lambda_c(i - 1, j, k)) * 5.0 /
                  6.0 -
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) *
                  3.0 / 4.0 -
              Jacobian_c(i + 1, j, k) *
                  (2.0 * mu_c(i + 1, j, k) + lambda_c(i + 1, j, k)) * 5.0 /
                  6.0 -
              Jacobian_c(i + 2, j, k) *
                  (2.0 * mu_c(i + 2, j, k) + lambda_c(i + 2, j, k)) / 24.0) *
                 u_c_t(1, i - 0, j, k) +
             (Jacobian_c(i - 1, j, k) *
                  (2.0 * mu_c(i - 1, j, k) + lambda_c(i - 1, j, k)) / 6.0 +
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  2.0 +
              Jacobian_c(i + 1, j, k) *
                  (2.0 * mu_c(i + 1, j, k) + lambda_c(i + 1, j, k)) / 2.0 +
              Jacobian_c(i + 2, j, k) *
                  (2.0 * mu_c(i + 2, j, k) + lambda_c(i + 2, j, k)) / 6.0) *
                 u_c_t(1, i + 1, j, k) +
             (-Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  8.0 +
              Jacobian_c(i + 1, j, k) *
                  (2.0 * mu_c(i + 1, j, k) + lambda_c(i + 1, j, k)) / 6.0 -
              Jacobian_c(i + 2, j, k) *
                  (2.0 * mu_c(i + 2, j, k) + lambda_c(i + 2, j, k)) / 8.0) *
                 u_c_t(1, i + 2, j, k)) /
                pow(h1_c, 2) / pow(l1, 2) +
            ((-Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 8.0 +
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0) *
                 u_c_t(1, i, j - 2, k) +
             (Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 6.0 +
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 2.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 6.0) *
                 u_c_t(1, i, j - 1, k) +
             (-Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 24.0 -
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) * 5.0 / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) * 3.0 / 4.0 -
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) * 5.0 / 6.0 -
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 24.0) *
                 u_c_t(1, i, j - 0, k) +
             (Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 6.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 2.0 +
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 6.0) *
                 u_c_t(1, i, j + 1, k) +
             (-Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 6.0 -
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 8.0) *
                 u_c_t(1, i, j + 2, k)) /
                pow(h2_c, 2) / pow(l2, 2) +
            (Jacobian_c(i - 2, j, k) * lambda_c(i - 2, j, k) *
                 (u_c_t(2, i - 2, j - 2, k) / 12.0 -
                  u_c_t(2, i - 2, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i - 2, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i - 2, j + 2, k) / 12.0) /
                 12.0 -
             Jacobian_c(i - 1, j, k) * lambda_c(i - 1, j, k) *
                 (u_c_t(2, i - 1, j - 2, k) / 12.0 -
                  u_c_t(2, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i - 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i - 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(i + 1, j, k) * lambda_c(i + 1, j, k) *
                 (u_c_t(2, i + 1, j - 2, k) / 12.0 -
                  u_c_t(2, i + 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(i + 2, j, k) * lambda_c(i + 2, j, k) *
                 (u_c_t(2, i + 2, j - 2, k) / 12.0 -
                  u_c_t(2, i + 2, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 2, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_c / h2_c +
            (Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) *
                 (u_c_t(2, i - 2, j - 2, k) / 12.0 -
                  u_c_t(2, i - 1, j - 2, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 1, j - 2, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 2, j - 2, k) / 12.0) /
                 12.0 -
             Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) *
                 (u_c_t(2, i - 2, j - 1, k) / 12.0 -
                  u_c_t(2, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 1, j - 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 2, j - 1, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) *
                 (u_c_t(2, i - 2, j + 1, k) / 12.0 -
                  u_c_t(2, i - 1, j + 1, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 2, j + 1, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) *
                 (u_c_t(2, i - 2, j + 2, k) / 12.0 -
                  u_c_t(2, i - 1, j + 2, k) * 2.0 / 3.0 +
                  u_c_t(2, i + 1, j + 2, k) * 2.0 / 3.0 -
                  u_c_t(2, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_c / h2_c;
        // second set
        lh_c(i, j, k, 2) =
            lh_c(i, j, k, 2) +
            ((-Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 8.0 +
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0) *
                 u_c_t(2, i - 2, j, k) +
             (Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 6.0 +
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 2.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 6.0) *
                 u_c_t(2, i - 1, j, k) +
             (-Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 24.0 -
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) * 5.0 / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) * 3.0 / 4.0 -
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) * 5.0 / 6.0 -
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 24.0) *
                 u_c_t(2, i - 0, j, k) +
             (Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 6.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 2.0 +
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 6.0) *
                 u_c_t(2, i + 1, j, k) +
             (-Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 6.0 -
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 8.0) *
                 u_c_t(2, i + 2, j, k)) /
                pow(h1_c, 2) / pow(l1, 2) +
            ((-Jacobian_c(i, j - 2, k) *
                  (2.0 * mu_c(i, j - 2, k) + lambda_c(i, j - 2, k)) / 8.0 +
              Jacobian_c(i, j - 1, k) *
                  (2.0 * mu_c(i, j - 1, k) + lambda_c(i, j - 1, k)) / 6.0 -
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  8.0) *
                 u_c_t(2, i, j - 2, k) +
             (Jacobian_c(i, j - 2, k) *
                  (2.0 * mu_c(i, j - 2, k) + lambda_c(i, j - 2, k)) / 6.0 +
              Jacobian_c(i, j - 1, k) *
                  (2.0 * mu_c(i, j - 1, k) + lambda_c(i, j - 1, k)) / 2.0 +
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  2.0 +
              Jacobian_c(i, j + 1, k) *
                  (2.0 * mu_c(i, j + 1, k) + lambda_c(i, j + 1, k)) / 6.0) *
                 u_c_t(2, i, j - 1, k) +
             (-Jacobian_c(i, j - 2, k) *
                  (2.0 * mu_c(i, j - 2, k) + lambda_c(i, j - 2, k)) / 24.0 -
              Jacobian_c(i, j - 1, k) *
                  (2.0 * mu_c(i, j - 1, k) + lambda_c(i, j - 1, k)) * 5.0 /
                  6.0 -
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) *
                  3.0 / 4.0 -
              Jacobian_c(i, j + 1, k) *
                  (2.0 * mu_c(i, j + 1, k) + lambda_c(i, j + 1, k)) * 5.0 /
                  6.0 -
              Jacobian_c(i, j + 2, k) *
                  (2.0 * mu_c(i, j + 2, k) + lambda_c(i, j + 2, k)) / 24.0) *
                 u_c_t(2, i, j - 0, k) +
             (Jacobian_c(i, j - 1, k) *
                  (2.0 * mu_c(i, j - 1, k) + lambda_c(i, j - 1, k)) / 6.0 +
              Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  2.0 +
              Jacobian_c(i, j + 1, k) *
                  (2.0 * mu_c(i, j + 1, k) + lambda_c(i, j + 1, k)) / 2.0 +
              Jacobian_c(i, j + 2, k) *
                  (2.0 * mu_c(i, j + 2, k) + lambda_c(i, j + 2, k)) / 6.0) *
                 u_c_t(2, i, j + 1, k) +
             (-Jacobian_c(i, j, k) * (2.0 * mu_c(i, j, k) + lambda_c(i, j, k)) /
                  8.0 +
              Jacobian_c(i, j + 1, k) *
                  (2.0 * mu_c(i, j + 1, k) + lambda_c(i, j + 1, k)) / 6.0 -
              Jacobian_c(i, j + 2, k) *
                  (2.0 * mu_c(i, j + 2, k) + lambda_c(i, j + 2, k)) / 8.0) *
                 u_c_t(2, i, j + 2, k)) /
                pow(h2_c, 2) / pow(l2, 2) +
            (Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) *
                 (u_c_t(1, i - 2, j - 2, k) / 12.0 -
                  u_c_t(1, i - 2, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i - 2, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i - 2, j + 2, k) / 12.0) /
                 12.0 -
             Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) *
                 (u_c_t(1, i - 1, j - 2, k) / 12.0 -
                  u_c_t(1, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i - 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i - 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) *
                 (u_c_t(1, i + 1, j - 2, k) / 12.0 -
                  u_c_t(1, i + 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 1, j + 2, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) *
                 (u_c_t(1, i + 2, j - 2, k) / 12.0 -
                  u_c_t(1, i + 2, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 2, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_c / h2_c +
            (Jacobian_c(i, j - 2, k) * lambda_c(i, j - 2, k) *
                 (u_c_t(1, i - 2, j - 2, k) / 12.0 -
                  u_c_t(1, i - 1, j - 2, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 1, j - 2, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 2, j - 2, k) / 12.0) /
                 12.0 -
             Jacobian_c(i, j - 1, k) * lambda_c(i, j - 1, k) *
                 (u_c_t(1, i - 2, j - 1, k) / 12.0 -
                  u_c_t(1, i - 1, j - 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 1, j - 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 2, j - 1, k) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(i, j + 1, k) * lambda_c(i, j + 1, k) *
                 (u_c_t(1, i - 2, j + 1, k) / 12.0 -
                  u_c_t(1, i - 1, j + 1, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 1, j + 1, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 2, j + 1, k) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(i, j + 2, k) * lambda_c(i, j + 2, k) *
                 (u_c_t(1, i - 2, j + 2, k) / 12.0 -
                  u_c_t(1, i - 1, j + 2, k) * 2.0 / 3.0 +
                  u_c_t(1, i + 1, j + 2, k) * 2.0 / 3.0 -
                  u_c_t(1, i + 2, j + 2, k) / 12.0) /
                 12.0) /
                l1 / l2 / h1_c / h2_c;
        // third set
        lh_c(i, j, k, 3) =
            lh_c(i, j, k, 3) +
            ((-Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 8.0 +
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0) *
                 u_c_t(3, i - 2, j, k) +
             (Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 6.0 +
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 2.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 6.0) *
                 u_c_t(3, i - 1, j, k) +
             (-Jacobian_c(i - 2, j, k) * mu_c(i - 2, j, k) / 24.0 -
              Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) * 5.0 / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) * 3.0 / 4.0 -
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) * 5.0 / 6.0 -
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 24.0) *
                 u_c_t(3, i - 0, j, k) +
             (Jacobian_c(i - 1, j, k) * mu_c(i - 1, j, k) / 6.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 2.0 +
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 6.0) *
                 u_c_t(3, i + 1, j, k) +
             (-Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0 +
              Jacobian_c(i + 1, j, k) * mu_c(i + 1, j, k) / 6.0 -
              Jacobian_c(i + 2, j, k) * mu_c(i + 2, j, k) / 8.0) *
                 u_c_t(3, i + 2, j, k)) /
                pow(h1_c, 2) / pow(l1, 2) +
            ((-Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 8.0 +
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0) *
                 u_c_t(3, i, j - 2, k) +
             (Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 6.0 +
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 2.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 6.0) *
                 u_c_t(3, i, j - 1, k) +
             (-Jacobian_c(i, j - 2, k) * mu_c(i, j - 2, k) / 24.0 -
              Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) * 5.0 / 6.0 -
              Jacobian_c(i, j, k) * mu_c(i, j, k) * 3.0 / 4.0 -
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) * 5.0 / 6.0 -
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 24.0) *
                 u_c_t(3, i, j - 0, k) +
             (Jacobian_c(i, j - 1, k) * mu_c(i, j - 1, k) / 6.0 +
              Jacobian_c(i, j, k) * mu_c(i, j, k) / 2.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 2.0 +
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 6.0) *
                 u_c_t(3, i, j + 1, k) +
             (-Jacobian_c(i, j, k) * mu_c(i, j, k) / 8.0 +
              Jacobian_c(i, j + 1, k) * mu_c(i, j + 1, k) / 6.0 -
              Jacobian_c(i, j + 2, k) * mu_c(i, j + 2, k) / 8.0) *
                 u_c_t(3, i, j + 2, k)) /
                pow(h2_c, 2) / pow(l2, 2);
      }
    }
  }
  //
  for (k = 1; k <= 4; k++) {
    for (i = 1; i <= n2_c; i++) {
      for (j = 1; j <= n1_c; j++) {
        for (k1 = 1; k1 <= 6; k1++) {
          // mixed derivative 13  23
          // first set equation
          lh_c(j, i, k, 1) =
              lh_c(j, i, k, 1) +
              (Jacobian_c(j - 2, i, k) *
                   (2.0 * mu_c(j - 2, i, k) + lambda_c(j - 2, i, k)) *
                   XI13_c(j - 2, i, k) * bof(k, k1) * u_c_t(1, j - 2, i, k1) /
                   12.0 -
               Jacobian_c(j - 1, i, k) *
                   (2.0 * mu_c(j - 1, i, k) + lambda_c(j - 1, i, k)) *
                   XI13_c(j - 1, i, k) * bof(k, k1) * u_c_t(1, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j + 1, i, k) *
                   (2.0 * mu_c(j + 1, i, k) + lambda_c(j + 1, i, k)) *
                   XI13_c(j + 1, i, k) * bof(k, k1) * u_c_t(1, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j + 2, i, k) *
                   (2.0 * mu_c(j + 2, i, k) + lambda_c(j + 2, i, k)) *
                   XI13_c(j + 2, i, k) * bof(k, k1) * u_c_t(1, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_c / h3_c +
              (Jacobian_c(j - 2, i, k) * lambda_c(j - 2, i, k) *
                   XI23_c(j - 2, i, k) * bof(k, k1) * u_c_t(2, j - 2, i, k1) /
                   12.0 -
               Jacobian_c(j - 1, i, k) * lambda_c(j - 1, i, k) *
                   XI23_c(j - 1, i, k) * bof(k, k1) * u_c_t(2, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j + 1, i, k) * lambda_c(j + 1, i, k) *
                   XI23_c(j + 1, i, k) * bof(k, k1) * u_c_t(2, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j + 2, i, k) * lambda_c(j + 2, i, k) *
                   XI23_c(j + 2, i, k) * bof(k, k1) * u_c_t(2, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_c / h3_c +
              (Jacobian_c(j - 2, i, k) * lambda_c(j - 2, i, k) *
                   XI33_c(j - 2, i, k) * bof(k, k1) * u_c_t(3, j - 2, i, k1) /
                   12.0 -
               Jacobian_c(j - 1, i, k) * lambda_c(j - 1, i, k) *
                   XI33_c(j - 1, i, k) * bof(k, k1) * u_c_t(3, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j + 1, i, k) * lambda_c(j + 1, i, k) *
                   XI33_c(j + 1, i, k) * bof(k, k1) * u_c_t(3, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j + 2, i, k) * lambda_c(j + 2, i, k) *
                   XI33_c(j + 2, i, k) * bof(k, k1) * u_c_t(3, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_c / h3_c +
              (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) *
                   XI23_c(j, i - 2, k) * bof(k, k1) * u_c_t(1, j, i - 2, k1) /
                   12.0 -
               Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) *
                   XI23_c(j, i - 1, k) * bof(k, k1) * u_c_t(1, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) *
                   XI23_c(j, i + 1, k) * bof(k, k1) * u_c_t(1, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) *
                   XI23_c(j, i + 2, k) * bof(k, k1) * u_c_t(1, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_c / h3_c +
              (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) *
                   XI13_c(j, i - 2, k) * bof(k, k1) * u_c_t(2, j, i - 2, k1) /
                   12.0 -
               Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) *
                   XI13_c(j, i - 1, k) * bof(k, k1) * u_c_t(2, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) *
                   XI13_c(j, i + 1, k) * bof(k, k1) * u_c_t(2, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) *
                   XI13_c(j, i + 2, k) * bof(k, k1) * u_c_t(2, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_c / h3_c;

          lh_c(j, i, n3_c + 1 - k, 1) =
              lh_c(j, i, n3_c + 1 - k, 1) +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   (2.0 * mu_c(j - 2, i, n3_c + 1 - k) +
                    lambda_c(j - 2, i, n3_c + 1 - k)) *
                   XI13_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   (2.0 * mu_c(j - 1, i, n3_c + 1 - k) +
                    lambda_c(j - 1, i, n3_c + 1 - k)) *
                   XI13_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   (2.0 * mu_c(j + 1, i, n3_c + 1 - k) +
                    lambda_c(j + 1, i, n3_c + 1 - k)) *
                   XI13_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   (2.0 * mu_c(j + 2, i, n3_c + 1 - k) +
                    lambda_c(j + 2, i, n3_c + 1 - k)) *
                   XI13_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   lambda_c(j - 2, i, n3_c + 1 - k) *
                   XI23_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   lambda_c(j - 1, i, n3_c + 1 - k) *
                   XI23_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   lambda_c(j + 1, i, n3_c + 1 - k) *
                   XI23_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   lambda_c(j + 2, i, n3_c + 1 - k) *
                   XI23_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   lambda_c(j - 2, i, n3_c + 1 - k) *
                   XI33_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   lambda_c(j - 1, i, n3_c + 1 - k) *
                   XI33_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   lambda_c(j + 1, i, n3_c + 1 - k) *
                   XI33_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   lambda_c(j + 2, i, n3_c + 1 - k) *
                   XI33_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   mu_c(j, i - 2, n3_c + 1 - k) *
                   XI23_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   mu_c(j, i - 1, n3_c + 1 - k) *
                   XI23_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   mu_c(j, i + 1, n3_c + 1 - k) *
                   XI23_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   mu_c(j, i + 2, n3_c + 1 - k) *
                   XI23_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   mu_c(j, i - 2, n3_c + 1 - k) *
                   XI13_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   mu_c(j, i - 1, n3_c + 1 - k) *
                   XI13_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   mu_c(j, i + 1, n3_c + 1 - k) *
                   XI13_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   mu_c(j, i + 2, n3_c + 1 - k) *
                   XI13_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c;
          // second set equation
          lh_c(j, i, k, 2) =
              lh_c(j, i, k, 2) +
              (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) *
                   XI23_c(j - 2, i, k) * bof(k, k1) * u_c_t(1, j - 2, i, k1) /
                   12.0 -
               Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) *
                   XI23_c(j - 1, i, k) * bof(k, k1) * u_c_t(1, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) *
                   XI23_c(j + 1, i, k) * bof(k, k1) * u_c_t(1, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) *
                   XI23_c(j + 2, i, k) * bof(k, k1) * u_c_t(1, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_c / h3_c +
              (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) *
                   XI13_c(j - 2, i, k) * bof(k, k1) * u_c_t(2, j - 2, i, k1) /
                   12.0 -
               Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) *
                   XI13_c(j - 1, i, k) * bof(k, k1) * u_c_t(2, j - 1, i, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) *
                   XI13_c(j + 1, i, k) * bof(k, k1) * u_c_t(2, j + 1, i, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) *
                   XI13_c(j + 2, i, k) * bof(k, k1) * u_c_t(2, j + 2, i, k1) /
                   12.0) /
                  l1 / h1_c / h3_c +
              (Jacobian_c(j, i - 2, k) * lambda_c(j, i - 2, k) *
                   XI13_c(j, i - 2, k) * bof(k, k1) * u_c_t(1, j, i - 2, k1) /
                   12.0 -
               Jacobian_c(j, i - 1, k) * lambda_c(j, i - 1, k) *
                   XI13_c(j, i - 1, k) * bof(k, k1) * u_c_t(1, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j, i + 1, k) * lambda_c(j, i + 1, k) *
                   XI13_c(j, i + 1, k) * bof(k, k1) * u_c_t(1, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j, i + 2, k) * lambda_c(j, i + 2, k) *
                   XI13_c(j, i + 2, k) * bof(k, k1) * u_c_t(1, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_c / h3_c +
              (Jacobian_c(j, i - 2, k) *
                   (2.0 * mu_c(j, i - 2, k) + lambda_c(j, i - 2, k)) *
                   XI23_c(j, i - 2, k) * bof(k, k1) * u_c_t(2, j, i - 2, k1) /
                   12.0 -
               Jacobian_c(j, i - 1, k) *
                   (2.0 * mu_c(j, i - 1, k) + lambda_c(j, i - 1, k)) *
                   XI23_c(j, i - 1, k) * bof(k, k1) * u_c_t(2, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j, i + 1, k) *
                   (2.0 * mu_c(j, i + 1, k) + lambda_c(j, i + 1, k)) *
                   XI23_c(j, i + 1, k) * bof(k, k1) * u_c_t(2, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j, i + 2, k) *
                   (2.0 * mu_c(j, i + 2, k) + lambda_c(j, i + 2, k)) *
                   XI23_c(j, i + 2, k) * bof(k, k1) * u_c_t(2, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_c / h3_c +
              (Jacobian_c(j, i - 2, k) * lambda_c(j, i - 2, k) *
                   XI33_c(j, i - 2, k) * bof(k, k1) * u_c_t(3, j, i - 2, k1) /
                   12.0 -
               Jacobian_c(j, i - 1, k) * lambda_c(j, i - 1, k) *
                   XI33_c(j, i - 1, k) * bof(k, k1) * u_c_t(3, j, i - 1, k1) *
                   2.0 / 3.0 +
               Jacobian_c(j, i + 1, k) * lambda_c(j, i + 1, k) *
                   XI33_c(j, i + 1, k) * bof(k, k1) * u_c_t(3, j, i + 1, k1) *
                   2.0 / 3.0 -
               Jacobian_c(j, i + 2, k) * lambda_c(j, i + 2, k) *
                   XI33_c(j, i + 2, k) * bof(k, k1) * u_c_t(3, j, i + 2, k1) /
                   12.0) /
                  l2 / h2_c / h3_c;

          lh_c(j, i, n3_c + 1 - k, 2) =
              lh_c(j, i, n3_c + 1 - k, 2) +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   mu_c(j - 2, i, n3_c + 1 - k) *
                   XI23_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   mu_c(j - 1, i, n3_c + 1 - k) *
                   XI23_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   mu_c(j + 1, i, n3_c + 1 - k) *
                   XI23_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   mu_c(j + 2, i, n3_c + 1 - k) *
                   XI23_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   mu_c(j - 2, i, n3_c + 1 - k) *
                   XI13_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   mu_c(j - 1, i, n3_c + 1 - k) *
                   XI13_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   mu_c(j + 1, i, n3_c + 1 - k) *
                   XI13_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   mu_c(j + 2, i, n3_c + 1 - k) *
                   XI13_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   lambda_c(j, i - 2, n3_c + 1 - k) *
                   XI13_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   lambda_c(j, i - 1, n3_c + 1 - k) *
                   XI13_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   lambda_c(j, i + 1, n3_c + 1 - k) *
                   XI13_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   lambda_c(j, i + 2, n3_c + 1 - k) *
                   XI13_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   (2.0 * mu_c(j, i - 2, n3_c + 1 - k) +
                    lambda_c(j, i - 2, n3_c + 1 - k)) *
                   XI23_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   (2.0 * mu_c(j, i - 1, n3_c + 1 - k) +
                    lambda_c(j, i - 1, n3_c + 1 - k)) *
                   XI23_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   (2.0 * mu_c(j, i + 1, n3_c + 1 - k) +
                    lambda_c(j, i + 1, n3_c + 1 - k)) *
                   XI23_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   (2.0 * mu_c(j, i + 2, n3_c + 1 - k) +
                    lambda_c(j, i + 2, n3_c + 1 - k)) *
                   XI23_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   lambda_c(j, i - 2, n3_c + 1 - k) *
                   XI33_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   lambda_c(j, i - 1, n3_c + 1 - k) *
                   XI33_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   lambda_c(j, i + 1, n3_c + 1 - k) *
                   XI33_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   lambda_c(j, i + 2, n3_c + 1 - k) *
                   XI33_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c;
          // third set equation
          lh_c(j, i, k, 3) = lh_c(j, i, k, 3) +
                             (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) *
                                  XI33_c(j - 2, i, k) * bof(k, k1) *
                                  u_c_t(1, j - 2, i, k1) / 12.0 -
                              Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) *
                                  XI33_c(j - 1, i, k) * bof(k, k1) *
                                  u_c_t(1, j - 1, i, k1) * 2.0 / 3.0 +
                              Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) *
                                  XI33_c(j + 1, i, k) * bof(k, k1) *
                                  u_c_t(1, j + 1, i, k1) * 2.0 / 3.0 -
                              Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) *
                                  XI33_c(j + 2, i, k) * bof(k, k1) *
                                  u_c_t(1, j + 2, i, k1) / 12.0) /
                                 l1 / h1_c / h3_c +
                             (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) *
                                  XI13_c(j - 2, i, k) * bof(k, k1) *
                                  u_c_t(3, j - 2, i, k1) / 12.0 -
                              Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) *
                                  XI13_c(j - 1, i, k) * bof(k, k1) *
                                  u_c_t(3, j - 1, i, k1) * 2.0 / 3.0 +
                              Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) *
                                  XI13_c(j + 1, i, k) * bof(k, k1) *
                                  u_c_t(3, j + 1, i, k1) * 2.0 / 3.0 -
                              Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) *
                                  XI13_c(j + 2, i, k) * bof(k, k1) *
                                  u_c_t(3, j + 2, i, k1) / 12.0) /
                                 l1 / h1_c / h3_c +
                             (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) *
                                  XI33_c(j, i - 2, k) * bof(k, k1) *
                                  u_c_t(2, j, i - 2, k1) / 12.0 -
                              Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) *
                                  XI33_c(j, i - 1, k) * bof(k, k1) *
                                  u_c_t(2, j, i - 1, k1) * 2.0 / 3.0 +
                              Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) *
                                  XI33_c(j, i + 1, k) * bof(k, k1) *
                                  u_c_t(2, j, i + 1, k1) * 2.0 / 3.0 -
                              Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) *
                                  XI33_c(j, i + 2, k) * bof(k, k1) *
                                  u_c_t(2, j, i + 2, k1) / 12.0) /
                                 l2 / h2_c / h3_c +
                             (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) *
                                  XI23_c(j, i - 2, k) * bof(k, k1) *
                                  u_c_t(3, j, i - 2, k1) / 12.0 -
                              Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) *
                                  XI23_c(j, i - 1, k) * bof(k, k1) *
                                  u_c_t(3, j, i - 1, k1) * 2.0 / 3.0 +
                              Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) *
                                  XI23_c(j, i + 1, k) * bof(k, k1) *
                                  u_c_t(3, j, i + 1, k1) * 2.0 / 3.0 -
                              Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) *
                                  XI23_c(j, i + 2, k) * bof(k, k1) *
                                  u_c_t(3, j, i + 2, k1) / 12.0) /
                                 l2 / h2_c / h3_c;

          lh_c(j, i, n3_c + 1 - k, 3) =
              lh_c(j, i, n3_c + 1 - k, 3) +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   mu_c(j - 2, i, n3_c + 1 - k) *
                   XI33_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   mu_c(j - 1, i, n3_c + 1 - k) *
                   XI33_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   mu_c(j + 1, i, n3_c + 1 - k) *
                   XI33_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   mu_c(j + 2, i, n3_c + 1 - k) *
                   XI33_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(1, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j - 2, i, n3_c + 1 - k) *
                   mu_c(j - 2, i, n3_c + 1 - k) *
                   XI13_c(j - 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j - 2, i, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j - 1, i, n3_c + 1 - k) *
                   mu_c(j - 1, i, n3_c + 1 - k) *
                   XI13_c(j - 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j - 1, i, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j + 1, i, n3_c + 1 - k) *
                   mu_c(j + 1, i, n3_c + 1 - k) *
                   XI13_c(j + 1, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j + 1, i, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j + 2, i, n3_c + 1 - k) *
                   mu_c(j + 2, i, n3_c + 1 - k) *
                   XI13_c(j + 2, i, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j + 2, i, n3_c + 1 - k1) / 12.0) /
                  l1 / h1_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   mu_c(j, i - 2, n3_c + 1 - k) *
                   XI33_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   mu_c(j, i - 1, n3_c + 1 - k) *
                   XI33_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   mu_c(j, i + 1, n3_c + 1 - k) *
                   XI33_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   mu_c(j, i + 2, n3_c + 1 - k) *
                   XI33_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(2, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c +
              (-Jacobian_c(j, i - 2, n3_c + 1 - k) *
                   mu_c(j, i - 2, n3_c + 1 - k) *
                   XI23_c(j, i - 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i - 2, n3_c + 1 - k1) / 12.0 +
               Jacobian_c(j, i - 1, n3_c + 1 - k) *
                   mu_c(j, i - 1, n3_c + 1 - k) *
                   XI23_c(j, i - 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i - 1, n3_c + 1 - k1) * 2.0 / 3.0 -
               Jacobian_c(j, i + 1, n3_c + 1 - k) *
                   mu_c(j, i + 1, n3_c + 1 - k) *
                   XI23_c(j, i + 1, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i + 1, n3_c + 1 - k1) * 2.0 / 3.0 +
               Jacobian_c(j, i + 2, n3_c + 1 - k) *
                   mu_c(j, i + 2, n3_c + 1 - k) *
                   XI23_c(j, i + 2, n3_c + 1 - k) * bof(k, k1) *
                   u_c_t(3, j, i + 2, n3_c + 1 - k1) / 12.0) /
                  l2 / h2_c / h3_c;
        }
      }
    }
  }
  //
  for (k = 5; k <= n3_c - 4; k++) {
    for (i = 1; i <= n2_c; i++) {
      for (j = 1; j <= n1_c; j++) {
        // mixed derivative 13  23
        // first set equation
        lh_c(j, i, k, 1) =
            lh_c(j, i, k, 1) +
            (Jacobian_c(j - 2, i, k) *
                 (2.0 * mu_c(j - 2, i, k) + lambda_c(j - 2, i, k)) *
                 XI13_c(j - 2, i, k) *
                 (u_c_t(1, j - 2, i, k - 2) / 12.0 -
                  u_c_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) *
                 (2.0 * mu_c(j - 1, i, k) + lambda_c(j - 1, i, k)) *
                 XI13_c(j - 1, i, k) *
                 (u_c_t(1, j - 1, i, k - 2) / 12.0 -
                  u_c_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) *
                 (2.0 * mu_c(j + 1, i, k) + lambda_c(j + 1, i, k)) *
                 XI13_c(j + 1, i, k) *
                 (u_c_t(1, j + 1, i, k - 2) / 12.0 -
                  u_c_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) *
                 (2.0 * mu_c(j + 2, i, k) + lambda_c(j + 2, i, k)) *
                 XI13_c(j + 2, i, k) *
                 (u_c_t(1, j + 2, i, k - 2) / 12.0 -
                  u_c_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j - 2, i, k) * lambda_c(j - 2, i, k) *
                 XI23_c(j - 2, i, k) *
                 (u_c_t(2, j - 2, i, k - 2) / 12.0 -
                  u_c_t(2, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * lambda_c(j - 1, i, k) *
                 XI23_c(j - 1, i, k) *
                 (u_c_t(2, j - 1, i, k - 2) / 12.0 -
                  u_c_t(2, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * lambda_c(j + 1, i, k) *
                 XI23_c(j + 1, i, k) *
                 (u_c_t(2, j + 1, i, k - 2) / 12.0 -
                  u_c_t(2, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * lambda_c(j + 2, i, k) *
                 XI23_c(j + 2, i, k) *
                 (u_c_t(2, j + 2, i, k - 2) / 12.0 -
                  u_c_t(2, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j - 2, i, k) * lambda_c(j - 2, i, k) *
                 XI33_c(j - 2, i, k) *
                 (u_c_t(3, j - 2, i, k - 2) / 12.0 -
                  u_c_t(3, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * lambda_c(j - 1, i, k) *
                 XI33_c(j - 1, i, k) *
                 (u_c_t(3, j - 1, i, k - 2) / 12.0 -
                  u_c_t(3, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * lambda_c(j + 1, i, k) *
                 XI33_c(j + 1, i, k) *
                 (u_c_t(3, j + 1, i, k - 2) / 12.0 -
                  u_c_t(3, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * lambda_c(j + 2, i, k) *
                 XI33_c(j + 2, i, k) *
                 (u_c_t(3, j + 2, i, k - 2) / 12.0 -
                  u_c_t(3, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) * XI23_c(j, i - 2, k) *
                 (u_c_t(1, j, i - 2, k - 2) / 12.0 -
                  u_c_t(1, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) * XI23_c(j, i - 1, k) *
                 (u_c_t(1, j, i - 1, k - 2) / 12.0 -
                  u_c_t(1, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) * XI23_c(j, i + 1, k) *
                 (u_c_t(1, j, i + 1, k - 2) / 12.0 -
                  u_c_t(1, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) * XI23_c(j, i + 2, k) *
                 (u_c_t(1, j, i + 2, k - 2) / 12.0 -
                  u_c_t(1, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c +
            (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) * XI13_c(j, i - 2, k) *
                 (u_c_t(2, j, i - 2, k - 2) / 12.0 -
                  u_c_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) * XI13_c(j, i - 1, k) *
                 (u_c_t(2, j, i - 1, k - 2) / 12.0 -
                  u_c_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) * XI13_c(j, i + 1, k) *
                 (u_c_t(2, j, i + 1, k - 2) / 12.0 -
                  u_c_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) * XI13_c(j, i + 2, k) *
                 (u_c_t(2, j, i + 2, k - 2) / 12.0 -
                  u_c_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c;
        // second set equation
        lh_c(j, i, k, 2) =
            lh_c(j, i, k, 2) +
            (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) * XI23_c(j - 2, i, k) *
                 (u_c_t(1, j - 2, i, k - 2) / 12.0 -
                  u_c_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) * XI23_c(j - 1, i, k) *
                 (u_c_t(1, j - 1, i, k - 2) / 12.0 -
                  u_c_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) * XI23_c(j + 1, i, k) *
                 (u_c_t(1, j + 1, i, k - 2) / 12.0 -
                  u_c_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) * XI23_c(j + 2, i, k) *
                 (u_c_t(1, j + 2, i, k - 2) / 12.0 -
                  u_c_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) * XI13_c(j - 2, i, k) *
                 (u_c_t(2, j - 2, i, k - 2) / 12.0 -
                  u_c_t(2, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) * XI13_c(j - 1, i, k) *
                 (u_c_t(2, j - 1, i, k - 2) / 12.0 -
                  u_c_t(2, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) * XI13_c(j + 1, i, k) *
                 (u_c_t(2, j + 1, i, k - 2) / 12.0 -
                  u_c_t(2, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) * XI13_c(j + 2, i, k) *
                 (u_c_t(2, j + 2, i, k - 2) / 12.0 -
                  u_c_t(2, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j, i - 2, k) * lambda_c(j, i - 2, k) *
                 XI13_c(j, i - 2, k) *
                 (u_c_t(1, j, i - 2, k - 2) / 12.0 -
                  u_c_t(1, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * lambda_c(j, i - 1, k) *
                 XI13_c(j, i - 1, k) *
                 (u_c_t(1, j, i - 1, k - 2) / 12.0 -
                  u_c_t(1, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * lambda_c(j, i + 1, k) *
                 XI13_c(j, i + 1, k) *
                 (u_c_t(1, j, i + 1, k - 2) / 12.0 -
                  u_c_t(1, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * lambda_c(j, i + 2, k) *
                 XI13_c(j, i + 2, k) *
                 (u_c_t(1, j, i + 2, k - 2) / 12.0 -
                  u_c_t(1, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c +
            (Jacobian_c(j, i - 2, k) *
                 (2.0 * mu_c(j, i - 2, k) + lambda_c(j, i - 2, k)) *
                 XI23_c(j, i - 2, k) *
                 (u_c_t(2, j, i - 2, k - 2) / 12.0 -
                  u_c_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) *
                 (2.0 * mu_c(j, i - 1, k) + lambda_c(j, i - 1, k)) *
                 XI23_c(j, i - 1, k) *
                 (u_c_t(2, j, i - 1, k - 2) / 12.0 -
                  u_c_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) *
                 (2.0 * mu_c(j, i + 1, k) + lambda_c(j, i + 1, k)) *
                 XI23_c(j, i + 1, k) *
                 (u_c_t(2, j, i + 1, k - 2) / 12.0 -
                  u_c_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) *
                 (2.0 * mu_c(j, i + 2, k) + lambda_c(j, i + 2, k)) *
                 XI23_c(j, i + 2, k) *
                 (u_c_t(2, j, i + 2, k - 2) / 12.0 -
                  u_c_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c +
            (Jacobian_c(j, i - 2, k) * lambda_c(j, i - 2, k) *
                 XI33_c(j, i - 2, k) *
                 (u_c_t(3, j, i - 2, k - 2) / 12.0 -
                  u_c_t(3, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * lambda_c(j, i - 1, k) *
                 XI33_c(j, i - 1, k) *
                 (u_c_t(3, j, i - 1, k - 2) / 12.0 -
                  u_c_t(3, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * lambda_c(j, i + 1, k) *
                 XI33_c(j, i + 1, k) *
                 (u_c_t(3, j, i + 1, k - 2) / 12.0 -
                  u_c_t(3, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * lambda_c(j, i + 2, k) *
                 XI33_c(j, i + 2, k) *
                 (u_c_t(3, j, i + 2, k - 2) / 12.0 -
                  u_c_t(3, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c;
        // third set equation
        lh_c(j, i, k, 3) =
            lh_c(j, i, k, 3) +
            (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) * XI33_c(j - 2, i, k) *
                 (u_c_t(1, j - 2, i, k - 2) / 12.0 -
                  u_c_t(1, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) * XI33_c(j - 1, i, k) *
                 (u_c_t(1, j - 1, i, k - 2) / 12.0 -
                  u_c_t(1, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) * XI33_c(j + 1, i, k) *
                 (u_c_t(1, j + 1, i, k - 2) / 12.0 -
                  u_c_t(1, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) * XI33_c(j + 2, i, k) *
                 (u_c_t(1, j + 2, i, k - 2) / 12.0 -
                  u_c_t(1, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(1, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(1, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j - 2, i, k) * mu_c(j - 2, i, k) * XI13_c(j - 2, i, k) *
                 (u_c_t(3, j - 2, i, k - 2) / 12.0 -
                  u_c_t(3, j - 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j - 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j - 2, i, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j - 1, i, k) * mu_c(j - 1, i, k) * XI13_c(j - 1, i, k) *
                 (u_c_t(3, j - 1, i, k - 2) / 12.0 -
                  u_c_t(3, j - 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j - 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j - 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j + 1, i, k) * mu_c(j + 1, i, k) * XI13_c(j + 1, i, k) *
                 (u_c_t(3, j + 1, i, k - 2) / 12.0 -
                  u_c_t(3, j + 1, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j + 1, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j + 1, i, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j + 2, i, k) * mu_c(j + 2, i, k) * XI13_c(j + 2, i, k) *
                 (u_c_t(3, j + 2, i, k - 2) / 12.0 -
                  u_c_t(3, j + 2, i, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j + 2, i, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j + 2, i, k + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) * XI33_c(j, i - 2, k) *
                 (u_c_t(2, j, i - 2, k - 2) / 12.0 -
                  u_c_t(2, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) * XI33_c(j, i - 1, k) *
                 (u_c_t(2, j, i - 1, k - 2) / 12.0 -
                  u_c_t(2, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) * XI33_c(j, i + 1, k) *
                 (u_c_t(2, j, i + 1, k - 2) / 12.0 -
                  u_c_t(2, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) * XI33_c(j, i + 2, k) *
                 (u_c_t(2, j, i + 2, k - 2) / 12.0 -
                  u_c_t(2, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(2, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(2, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c +
            (Jacobian_c(j, i - 2, k) * mu_c(j, i - 2, k) * XI23_c(j, i - 2, k) *
                 (u_c_t(3, j, i - 2, k - 2) / 12.0 -
                  u_c_t(3, j, i - 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i - 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i - 2, k + 2) / 12.0) /
                 12.0 -
             Jacobian_c(j, i - 1, k) * mu_c(j, i - 1, k) * XI23_c(j, i - 1, k) *
                 (u_c_t(3, j, i - 1, k - 2) / 12.0 -
                  u_c_t(3, j, i - 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i - 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i - 1, k + 2) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(j, i + 1, k) * mu_c(j, i + 1, k) * XI23_c(j, i + 1, k) *
                 (u_c_t(3, j, i + 1, k - 2) / 12.0 -
                  u_c_t(3, j, i + 1, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i + 1, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i + 1, k + 2) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(j, i + 2, k) * mu_c(j, i + 2, k) * XI23_c(j, i + 2, k) *
                 (u_c_t(3, j, i + 2, k - 2) / 12.0 -
                  u_c_t(3, j, i + 2, k - 1) * 2.0 / 3.0 +
                  u_c_t(3, j, i + 2, k + 1) * 2.0 / 3.0 -
                  u_c_t(3, j, i + 2, k + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c;
      }
    }
  }
  //
  for (i = 5; i <= n3_c - 4; i++) {
    for (j = 1; j <= n2_c; j++) {
      for (k = 1; k <= n1_c; k++) {
        // mixed derivative 31
        // first set equation
        lh_c(k, j, i, 1) =
            lh_c(k, j, i, 1) +
            (Jacobian_c(k, j, i - 2) *
                 (2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                 XI13_c(k, j, i - 2) *
                 (u_c_t(1, k - 2, j, i - 2) / 12.0 -
                  u_c_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) *
                 (2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                 XI13_c(k, j, i - 1) *
                 (u_c_t(1, k - 2, j, i - 1) / 12.0 -
                  u_c_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) *
                 (2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                 XI13_c(k, j, i + 1) *
                 (u_c_t(1, k - 2, j, i + 1) / 12.0 -
                  u_c_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) *
                 (2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                 XI13_c(k, j, i + 2) *
                 (u_c_t(1, k - 2, j, i + 2) / 12.0 -
                  u_c_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_c / h1_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI23_c(k, j, i - 2) *
                 (u_c_t(2, k - 2, j, i - 2) / 12.0 -
                  u_c_t(2, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI23_c(k, j, i - 1) *
                 (u_c_t(2, k - 2, j, i - 1) / 12.0 -
                  u_c_t(2, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI23_c(k, j, i + 1) *
                 (u_c_t(2, k - 2, j, i + 1) / 12.0 -
                  u_c_t(2, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI23_c(k, j, i + 2) *
                 (u_c_t(2, k - 2, j, i + 2) / 12.0 -
                  u_c_t(2, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_c / h1_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI33_c(k, j, i - 2) *
                 (u_c_t(3, k - 2, j, i - 2) / 12.0 -
                  u_c_t(3, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI33_c(k, j, i - 1) *
                 (u_c_t(3, k - 2, j, i - 1) / 12.0 -
                  u_c_t(3, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI33_c(k, j, i + 1) *
                 (u_c_t(3, k - 2, j, i + 1) / 12.0 -
                  u_c_t(3, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI33_c(k, j, i + 2) *
                 (u_c_t(3, k - 2, j, i + 2) / 12.0 -
                  u_c_t(3, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c;
        // second set equation
        lh_c(k, j, i, 2) =
            lh_c(k, j, i, 2) +
            (Jacobian_c(k, j, i - 2) * lambda_c(k, j, i - 2) *
                 XI23_c(k, j, i - 2) *
                 (u_c_t(1, k - 2, j, i - 2) / 12.0 -
                  u_c_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * lambda_c(k, j, i - 1) *
                 XI23_c(k, j, i - 1) *
                 (u_c_t(1, k - 2, j, i - 1) / 12.0 -
                  u_c_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * lambda_c(k, j, i + 1) *
                 XI23_c(k, j, i + 1) *
                 (u_c_t(1, k - 2, j, i + 1) / 12.0 -
                  u_c_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * lambda_c(k, j, i + 2) *
                 XI23_c(k, j, i + 2) *
                 (u_c_t(1, k - 2, j, i + 2) / 12.0 -
                  u_c_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI13_c(k, j, i - 2) *
                 (u_c_t(2, k - 2, j, i - 2) / 12.0 -
                  u_c_t(2, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI13_c(k, j, i - 1) *
                 (u_c_t(2, k - 2, j, i - 1) / 12.0 -
                  u_c_t(2, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI13_c(k, j, i + 1) *
                 (u_c_t(2, k - 2, j, i + 1) / 12.0 -
                  u_c_t(2, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI13_c(k, j, i + 2) *
                 (u_c_t(2, k - 2, j, i + 2) / 12.0 -
                  u_c_t(2, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(2, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(2, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h1_c / h3_c;
        // third set equation
        lh_c(k, j, i, 3) =
            lh_c(k, j, i, 3) +
            (Jacobian_c(k, j, i - 2) * lambda_c(k, j, i - 2) *
                 XI33_c(k, j, i - 2) *
                 (u_c_t(1, k - 2, j, i - 2) / 12.0 -
                  u_c_t(1, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * lambda_c(k, j, i - 1) *
                 XI33_c(k, j, i - 1) *
                 (u_c_t(1, k - 2, j, i - 1) / 12.0 -
                  u_c_t(1, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * lambda_c(k, j, i + 1) *
                 XI33_c(k, j, i + 1) *
                 (u_c_t(1, k - 2, j, i + 1) / 12.0 -
                  u_c_t(1, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * lambda_c(k, j, i + 2) *
                 XI33_c(k, j, i + 2) *
                 (u_c_t(1, k - 2, j, i + 2) / 12.0 -
                  u_c_t(1, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(1, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(1, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_c / h1_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI13_c(k, j, i - 2) *
                 (u_c_t(3, k - 2, j, i - 2) / 12.0 -
                  u_c_t(3, k - 1, j, i - 2) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i - 2) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI13_c(k, j, i - 1) *
                 (u_c_t(3, k - 2, j, i - 1) / 12.0 -
                  u_c_t(3, k - 1, j, i - 1) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i - 1) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI13_c(k, j, i + 1) *
                 (u_c_t(3, k - 2, j, i + 1) / 12.0 -
                  u_c_t(3, k - 1, j, i + 1) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i + 1) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI13_c(k, j, i + 2) *
                 (u_c_t(3, k - 2, j, i + 2) / 12.0 -
                  u_c_t(3, k - 1, j, i + 2) * 2.0 / 3.0 +
                  u_c_t(3, k + 1, j, i + 2) * 2.0 / 3.0 -
                  u_c_t(3, k + 2, j, i + 2) / 12.0) /
                 12.0) /
                l1 / h3_c / h1_c;

        // 32
        lh_c(k, j, i, 1) =
            lh_c(k, j, i, 1) +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI23_c(k, j, i - 2) *
                 (u_c_t(1, k, j - 2, i - 2) / 12.0 -
                  u_c_t(1, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI23_c(k, j, i - 1) *
                 (u_c_t(1, k, j - 2, i - 1) / 12.0 -
                  u_c_t(1, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI23_c(k, j, i + 1) *
                 (u_c_t(1, k, j - 2, i + 1) / 12.0 -
                  u_c_t(1, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI23_c(k, j, i + 2) *
                 (u_c_t(1, k, j - 2, i + 2) / 12.0 -
                  u_c_t(1, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_c / h2_c +
            (Jacobian_c(k, j, i - 2) * lambda_c(k, j, i - 2) *
                 XI13_c(k, j, i - 2) *
                 (u_c_t(2, k, j - 2, i - 2) / 12.0 -
                  u_c_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * lambda_c(k, j, i - 1) *
                 XI13_c(k, j, i - 1) *
                 (u_c_t(2, k, j - 2, i - 1) / 12.0 -
                  u_c_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * lambda_c(k, j, i + 1) *
                 XI13_c(k, j, i + 1) *
                 (u_c_t(2, k, j - 2, i + 1) / 12.0 -
                  u_c_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * lambda_c(k, j, i + 2) *
                 XI13_c(k, j, i + 2) *
                 (u_c_t(2, k, j - 2, i + 2) / 12.0 -
                  u_c_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_c / h2_c;
        // second set equation
        lh_c(k, j, i, 2) =
            lh_c(k, j, i, 2) +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI13_c(k, j, i - 2) *
                 (u_c_t(1, k, j - 2, i - 2) / 12.0 -
                  u_c_t(1, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI13_c(k, j, i - 1) *
                 (u_c_t(1, k, j - 2, i - 1) / 12.0 -
                  u_c_t(1, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI13_c(k, j, i + 1) *
                 (u_c_t(1, k, j - 2, i + 1) / 12.0 -
                  u_c_t(1, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI13_c(k, j, i + 2) *
                 (u_c_t(1, k, j - 2, i + 2) / 12.0 -
                  u_c_t(1, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(1, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(1, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_c / h2_c +
            (Jacobian_c(k, j, i - 2) *
                 (2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                 XI23_c(k, j, i - 2) *
                 (u_c_t(2, k, j - 2, i - 2) / 12.0 -
                  u_c_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) *
                 (2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                 XI23_c(k, j, i - 1) *
                 (u_c_t(2, k, j - 2, i - 1) / 12.0 -
                  u_c_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) *
                 (2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                 XI23_c(k, j, i + 1) *
                 (u_c_t(2, k, j - 2, i + 1) / 12.0 -
                  u_c_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) *
                 (2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                 XI23_c(k, j, i + 2) *
                 (u_c_t(2, k, j - 2, i + 2) / 12.0 -
                  u_c_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_c / h2_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI33_c(k, j, i - 2) *
                 (u_c_t(3, k, j - 2, i - 2) / 12.0 -
                  u_c_t(3, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI33_c(k, j, i - 1) *
                 (u_c_t(3, k, j - 2, i - 1) / 12.0 -
                  u_c_t(3, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI33_c(k, j, i + 1) *
                 (u_c_t(3, k, j - 2, i + 1) / 12.0 -
                  u_c_t(3, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI33_c(k, j, i + 2) *
                 (u_c_t(3, k, j - 2, i + 2) / 12.0 -
                  u_c_t(3, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h3_c / h2_c;
        // third set equation
        lh_c(k, j, i, 3) =
            lh_c(k, j, i, 3) +
            (Jacobian_c(k, j, i - 2) * lambda_c(k, j, i - 2) *
                 XI33_c(k, j, i - 2) *
                 (u_c_t(2, k, j - 2, i - 2) / 12.0 -
                  u_c_t(2, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * lambda_c(k, j, i - 1) *
                 XI33_c(k, j, i - 1) *
                 (u_c_t(2, k, j - 2, i - 1) / 12.0 -
                  u_c_t(2, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * lambda_c(k, j, i + 1) *
                 XI33_c(k, j, i + 1) *
                 (u_c_t(2, k, j - 2, i + 1) / 12.0 -
                  u_c_t(2, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * lambda_c(k, j, i + 2) *
                 XI33_c(k, j, i + 2) *
                 (u_c_t(2, k, j - 2, i + 2) / 12.0 -
                  u_c_t(2, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(2, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(2, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c +
            (Jacobian_c(k, j, i - 2) * mu_c(k, j, i - 2) * XI23_c(k, j, i - 2) *
                 (u_c_t(3, k, j - 2, i - 2) / 12.0 -
                  u_c_t(3, k, j - 1, i - 2) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i - 2) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i - 2) / 12.0) /
                 12.0 -
             Jacobian_c(k, j, i - 1) * mu_c(k, j, i - 1) * XI23_c(k, j, i - 1) *
                 (u_c_t(3, k, j - 2, i - 1) / 12.0 -
                  u_c_t(3, k, j - 1, i - 1) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i - 1) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i - 1) / 12.0) *
                 2.0 / 3.0 +
             Jacobian_c(k, j, i + 1) * mu_c(k, j, i + 1) * XI23_c(k, j, i + 1) *
                 (u_c_t(3, k, j - 2, i + 1) / 12.0 -
                  u_c_t(3, k, j - 1, i + 1) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i + 1) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i + 1) / 12.0) *
                 2.0 / 3.0 -
             Jacobian_c(k, j, i + 2) * mu_c(k, j, i + 2) * XI23_c(k, j, i + 2) *
                 (u_c_t(3, k, j - 2, i + 2) / 12.0 -
                  u_c_t(3, k, j - 1, i + 2) * 2.0 / 3.0 +
                  u_c_t(3, k, j + 1, i + 2) * 2.0 / 3.0 -
                  u_c_t(3, k, j + 2, i + 2) / 12.0) /
                 12.0) /
                l2 / h2_c / h3_c;
      }
    }
  }
  for (i = 1; i <= 4; i++) {
    for (j = 1; j <= n2_c; j++) {
      for (k = 1; k <= n1_c; k++) {
        for (k1 = 1; k1 <= 6; k1++) {
          // mixed derivative 31
          // first set equation
          lh_c(k, j, i, 1) =
              lh_c(k, j, i, 1) +
              (bof(i, k1) * Jacobian_c(k, j, k1) *
               (2.0 * mu_c(k, j, k1) + lambda_c(k, j, k1)) * XI13_c(k, j, k1) *
               (u_c_t(1, k - 2, j, k1) / 12.0 -
                u_c_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                u_c_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                u_c_t(1, k + 2, j, k1) / 12.0)) /
                  l1 / h3_c / h1_c +
              (bof(i, k1) * Jacobian_c(k, j, k1) * mu_c(k, j, k1) *
               XI23_c(k, j, k1) *
               (u_c_t(2, k - 2, j, k1) / 12.0 -
                u_c_t(2, k - 1, j, k1) * 2.0 / 3.0 +
                u_c_t(2, k + 1, j, k1) * 2.0 / 3.0 -
                u_c_t(2, k + 2, j, k1) / 12.0)) /
                  l1 / h3_c / h1_c +
              (bof(i, k1) * Jacobian_c(k, j, k1) * mu_c(k, j, k1) *
               XI33_c(k, j, k1) *
               (u_c_t(3, k - 2, j, k1) / 12.0 -
                u_c_t(3, k - 1, j, k1) * 2.0 / 3.0 +
                u_c_t(3, k + 1, j, k1) * 2.0 / 3.0 -
                u_c_t(3, k + 2, j, k1) / 12.0)) /
                  l1 / h1_c / h3_c;

          lh_c(k, j, n3_c + 1 - i, 1) =
              lh_c(k, j, n3_c + 1 - i, 1) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               (2.0 * mu_c(k, j, n3_c + 1 - k1) +
                lambda_c(k, j, n3_c + 1 - k1)) *
               XI13_c(k, j, n3_c + 1 - k1) *
               (u_c_t(1, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(1, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(1, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(1, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h3_c / h1_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI23_c(k, j, n3_c + 1 - k1) *
               (u_c_t(2, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(2, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(2, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(2, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h3_c / h1_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI33_c(k, j, n3_c + 1 - k1) *
               (u_c_t(3, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(3, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(3, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(3, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h1_c / h3_c;
          // second set equation
          lh_c(k, j, i, 2) = lh_c(k, j, i, 2) +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              lambda_c(k, j, k1) * XI23_c(k, j, k1) *
                              (u_c_t(1, k - 2, j, k1) / 12.0 -
                               u_c_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                               u_c_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                               u_c_t(1, k + 2, j, k1) / 12.0)) /
                                 l1 / h1_c / h3_c +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              mu_c(k, j, k1) * XI13_c(k, j, k1) *
                              (u_c_t(2, k - 2, j, k1) / 12.0 -
                               u_c_t(2, k - 1, j, k1) * 2.0 / 3.0 +
                               u_c_t(2, k + 1, j, k1) * 2.0 / 3.0 -
                               u_c_t(2, k + 2, j, k1) / 12.0)) /
                                 l1 / h1_c / h3_c;

          lh_c(k, j, n3_c + 1 - i, 2) =
              lh_c(k, j, n3_c + 1 - i, 2) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               lambda_c(k, j, n3_c + 1 - k1) * XI23_c(k, j, n3_c + 1 - k1) *
               (u_c_t(1, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(1, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(1, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(1, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h1_c / h3_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI13_c(k, j, n3_c + 1 - k1) *
               (u_c_t(2, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(2, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(2, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(2, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h1_c / h3_c;
          // third set equation
          lh_c(k, j, i, 3) = lh_c(k, j, i, 3) +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              lambda_c(k, j, k1) * XI33_c(k, j, k1) *
                              (u_c_t(1, k - 2, j, k1) / 12.0 -
                               u_c_t(1, k - 1, j, k1) * 2.0 / 3.0 +
                               u_c_t(1, k + 1, j, k1) * 2.0 / 3.0 -
                               u_c_t(1, k + 2, j, k1) / 12.0)) /
                                 l1 / h1_c / h3_c +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              mu_c(k, j, k1) * XI13_c(k, j, k1) *
                              (u_c_t(3, k - 2, j, k1) / 12.0 -
                               u_c_t(3, k - 1, j, k1) * 2.0 / 3.0 +
                               u_c_t(3, k + 1, j, k1) * 2.0 / 3.0 -
                               u_c_t(3, k + 2, j, k1) / 12.0)) /
                                 l1 / h3_c / h1_c;

          lh_c(k, j, n3_c + 1 - i, 3) =
              lh_c(k, j, n3_c + 1 - i, 3) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               lambda_c(k, j, n3_c + 1 - k1) * XI33_c(k, j, n3_c + 1 - k1) *
               (u_c_t(1, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(1, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(1, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(1, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h3_c / h1_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI13_c(k, j, n3_c + 1 - k1) *
               (u_c_t(3, k - 2, j, n3_c + 1 - k1) / 12.0 -
                u_c_t(3, k - 1, j, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(3, k + 1, j, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(3, k + 2, j, n3_c + 1 - k1) / 12.0)) /
                  l1 / h3_c / h1_c;

          // 31
          // first set equation
          lh_c(k, j, i, 1) = lh_c(k, j, i, 1) +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              mu_c(k, j, k1) * XI23_c(k, j, k1) *
                              (u_c_t(1, k, j - 2, k1) / 12.0 -
                               u_c_t(1, k, j - 1, k1) * 2.0 / 3.0 +
                               u_c_t(1, k, j + 1, k1) * 2.0 / 3.0 -
                               u_c_t(1, k, j + 2, k1) / 12.0)) /
                                 l2 / h3_c / h2_c +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              lambda_c(k, j, k1) * XI13_c(k, j, k1) *
                              (u_c_t(2, k, j - 2, k1) / 12.0 -
                               u_c_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                               u_c_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                               u_c_t(2, k, j + 2, k1) / 12.0)) /
                                 l2 / h3_c / h2_c;

          lh_c(k, j, n3_c + 1 - i, 1) =
              lh_c(k, j, n3_c + 1 - i, 1) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI23_c(k, j, n3_c + 1 - k1) *
               (u_c_t(1, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(1, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(1, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(1, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h3_c / h2_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               lambda_c(k, j, n3_c + 1 - k1) * XI13_c(k, j, n3_c + 1 - k1) *
               (u_c_t(2, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(2, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(2, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(2, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h3_c / h2_c;
          // second set equation
          lh_c(k, j, i, 2) =
              lh_c(k, j, i, 2) +
              (bof(i, k1) * Jacobian_c(k, j, k1) * mu_c(k, j, k1) *
               XI13_c(k, j, k1) *
               (u_c_t(1, k, j - 2, k1) / 12.0 -
                u_c_t(1, k, j - 1, k1) * 2.0 / 3.0 +
                u_c_t(1, k, j + 1, k1) * 2.0 / 3.0 -
                u_c_t(1, k, j + 2, k1) / 12.0)) /
                  l2 / h3_c / h2_c +
              (bof(i, k1) * Jacobian_c(k, j, k1) *
               (2.0 * mu_c(k, j, k1) + lambda_c(k, j, k1)) * XI23_c(k, j, k1) *
               (u_c_t(2, k, j - 2, k1) / 12.0 -
                u_c_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                u_c_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                u_c_t(2, k, j + 2, k1) / 12.0)) /
                  l2 / h3_c / h2_c +
              (bof(i, k1) * Jacobian_c(k, j, k1) * mu_c(k, j, k1) *
               XI33_c(k, j, k1) *
               (u_c_t(3, k, j - 2, k1) / 12.0 -
                u_c_t(3, k, j - 1, k1) * 2.0 / 3.0 +
                u_c_t(3, k, j + 1, k1) * 2.0 / 3.0 -
                u_c_t(3, k, j + 2, k1) / 12.0)) /
                  l2 / h3_c / h2_c;

          lh_c(k, j, n3_c + 1 - i, 2) =
              lh_c(k, j, n3_c + 1 - i, 2) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI13_c(k, j, n3_c + 1 - k1) *
               (u_c_t(1, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(1, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(1, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(1, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h3_c / h2_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               (2.0 * mu_c(k, j, n3_c + 1 - k1) +
                lambda_c(k, j, n3_c + 1 - k1)) *
               XI23_c(k, j, n3_c + 1 - k1) *
               (u_c_t(2, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(2, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(2, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(2, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h3_c / h2_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI33_c(k, j, n3_c + 1 - k1) *
               (u_c_t(3, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(3, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(3, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(3, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h3_c / h2_c;
          // third set equation
          lh_c(k, j, i, 3) = lh_c(k, j, i, 3) +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              lambda_c(k, j, k1) * XI33_c(k, j, k1) *
                              (u_c_t(2, k, j - 2, k1) / 12.0 -
                               u_c_t(2, k, j - 1, k1) * 2.0 / 3.0 +
                               u_c_t(2, k, j + 1, k1) * 2.0 / 3.0 -
                               u_c_t(2, k, j + 2, k1) / 12.0)) /
                                 l2 / h2_c / h3_c +
                             (bof(i, k1) * Jacobian_c(k, j, k1) *
                              mu_c(k, j, k1) * XI23_c(k, j, k1) *
                              (u_c_t(3, k, j - 2, k1) / 12.0 -
                               u_c_t(3, k, j - 1, k1) * 2.0 / 3.0 +
                               u_c_t(3, k, j + 1, k1) * 2.0 / 3.0 -
                               u_c_t(3, k, j + 2, k1) / 12.0)) /
                                 l2 / h2_c / h3_c;

          lh_c(k, j, n3_c + 1 - i, 3) =
              lh_c(k, j, n3_c + 1 - i, 3) +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               lambda_c(k, j, n3_c + 1 - k1) * XI33_c(k, j, n3_c + 1 - k1) *
               (u_c_t(2, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(2, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(2, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(2, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h2_c / h3_c +
              (-bof(i, k1) * Jacobian_c(k, j, n3_c + 1 - k1) *
               mu_c(k, j, n3_c + 1 - k1) * XI23_c(k, j, n3_c + 1 - k1) *
               (u_c_t(3, k, j - 2, n3_c + 1 - k1) / 12.0 -
                u_c_t(3, k, j - 1, n3_c + 1 - k1) * 2.0 / 3.0 +
                u_c_t(3, k, j + 1, n3_c + 1 - k1) * 2.0 / 3.0 -
                u_c_t(3, k, j + 2, n3_c + 1 - k1) / 12.0)) /
                  l2 / h2_c / h3_c;
        }
      }
    }
  }
  //
  for (i = 7; i <= n3_c - 6; i++) {
    for (j = 1; j <= n2_c; j++) {
      for (k = 1; k <= n1_c; k++) {
        // second derivative 33
        // first set equation
        lh_c(k, j, i, 1) =
            lh_c(k, j, i, 1) +
            ((-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI13_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI23_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI13_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI23_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI13_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI23_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  8.0) *
                 u_c_t(1, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI13_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI23_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI13_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI23_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI13_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI23_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI13_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI23_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  6.0) *
                 u_c_t(1, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI13_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI23_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI13_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI23_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI13_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI23_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI13_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI23_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI13_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI23_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  24.0) *
                 u_c_t(1, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI13_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI23_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI13_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI23_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI13_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI23_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI13_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI23_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  6.0) *
                 u_c_t(1, k, j, i + 1) +
             (-Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI13_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI23_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  8.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI13_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI23_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI13_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI23_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  8.0) *
                 u_c_t(1, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 8.0) *
                 u_c_t(2, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 6.0) *
                 u_c_t(2, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 24.0) *
                 u_c_t(2, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 6.0) *
                 u_c_t(2, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 8.0) *
                 u_c_t(2, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 8.0) *
                 u_c_t(3, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0) *
                 u_c_t(3, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 24.0) *
                 u_c_t(3, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 6.0) *
                 u_c_t(3, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 8.0) *
                 u_c_t(3, k, j, i + 2)) /
                pow(h3_c, 2);
        // second set of equation
        lh_c(k, j, i, 2) =
            lh_c(k, j, i, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 8.0) *
                 u_c_t(1, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 6.0) *
                 u_c_t(1, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI23_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 24.0) *
                 u_c_t(1, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI23_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 6.0) *
                 u_c_t(1, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI23_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI23_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI23_c(k, j, i + 2) / 8.0) *
                 u_c_t(1, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI23_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI23_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI23_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  8.0) *
                 u_c_t(2, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI23_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI23_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI23_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI23_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  6.0) *
                 u_c_t(2, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI23_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI33_c(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI23_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI23_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI23_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI23_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  24.0) *
                 u_c_t(2, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI23_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI33_c(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI23_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI23_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI23_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  6.0) *
                 u_c_t(2, k, j, i + 1) +
             (-Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI23_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI33_c(k, j, i), 2))) /
                  8.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI23_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI33_c(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI23_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI33_c(k, j, i + 2), 2))) /
                  8.0) *
                 u_c_t(2, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 8.0) *
                 u_c_t(3, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0) *
                 u_c_t(3, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 24.0) *
                 u_c_t(3, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 6.0) *
                 u_c_t(3, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 8.0) *
                 u_c_t(3, k, j, i + 2)) /
                pow(h3_c, 2);
        // third set equation
        lh_c(k, j, i, 3) =
            lh_c(k, j, i, 3) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 8.0) *
                 u_c_t(1, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0) *
                 u_c_t(1, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI13_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 24.0) *
                 u_c_t(1, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI13_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 6.0) *
                 u_c_t(1, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI13_c(k, j, i) * XI33_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI13_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI13_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 8.0) *
                 u_c_t(1, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 8.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 8.0) *
                 u_c_t(2, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 6.0 +
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 2.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0) *
                 u_c_t(2, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  (mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                  XI23_c(k, j, i - 2) * XI33_c(k, j, i - 2) / 24.0 -
              Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) * 3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) * 5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 24.0) *
                 u_c_t(2, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  (mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                  XI23_c(k, j, i - 1) * XI33_c(k, j, i - 1) / 6.0 +
              Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 2.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 2.0 +
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 6.0) *
                 u_c_t(2, k, j, i + 1) +
             (-Jacobian_c(k, j, i) * (mu_c(k, j, i) + lambda_c(k, j, i)) *
                  XI23_c(k, j, i) * XI33_c(k, j, i) / 8.0 +
              Jacobian_c(k, j, i + 1) *
                  (mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                  XI23_c(k, j, i + 1) * XI33_c(k, j, i + 1) / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  (mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                  XI23_c(k, j, i + 2) * XI33_c(k, j, i + 2) / 8.0) *
                 u_c_t(2, k, j, i + 2)) /
                pow(h3_c, 2) +
            ((-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI33_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI23_c(k, j, i - 2), 2))) /
                  8.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI33_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI23_c(k, j, i - 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI33_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI23_c(k, j, i), 2))) /
                  8.0) *
                 u_c_t(3, k, j, i - 2) +
             (Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI33_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI23_c(k, j, i - 2), 2))) /
                  6.0 +
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI33_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI23_c(k, j, i - 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI33_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI23_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI33_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI23_c(k, j, i + 1), 2))) /
                  6.0) *
                 u_c_t(3, k, j, i - 1) +
             (-Jacobian_c(k, j, i - 2) *
                  ((2.0 * mu_c(k, j, i - 2) + lambda_c(k, j, i - 2)) *
                       pow(XI33_c(k, j, i - 2), 2) +
                   mu_c(k, j, i - 2) * (pow(XI13_c(k, j, i - 2), 2) +
                                        pow(XI23_c(k, j, i - 2), 2))) /
                  24.0 -
              Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI33_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI23_c(k, j, i - 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI33_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI23_c(k, j, i), 2))) *
                  3.0 / 4.0 -
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI33_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI23_c(k, j, i + 1), 2))) *
                  5.0 / 6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI33_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI23_c(k, j, i + 2), 2))) /
                  24.0) *
                 u_c_t(3, k, j, i) +
             (Jacobian_c(k, j, i - 1) *
                  ((2.0 * mu_c(k, j, i - 1) + lambda_c(k, j, i - 1)) *
                       pow(XI33_c(k, j, i - 1), 2) +
                   mu_c(k, j, i - 1) * (pow(XI13_c(k, j, i - 1), 2) +
                                        pow(XI23_c(k, j, i - 1), 2))) /
                  6.0 +
              Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI33_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI23_c(k, j, i), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI33_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI23_c(k, j, i + 1), 2))) /
                  2.0 +
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI33_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI23_c(k, j, i + 2), 2))) /
                  6.0) *
                 u_c_t(3, k, j, i + 1) +
             (-Jacobian_c(k, j, i) *
                  ((2.0 * mu_c(k, j, i) + lambda_c(k, j, i)) *
                       pow(XI33_c(k, j, i), 2) +
                   mu_c(k, j, i) *
                       (pow(XI13_c(k, j, i), 2) + pow(XI23_c(k, j, i), 2))) /
                  8.0 +
              Jacobian_c(k, j, i + 1) *
                  ((2.0 * mu_c(k, j, i + 1) + lambda_c(k, j, i + 1)) *
                       pow(XI33_c(k, j, i + 1), 2) +
                   mu_c(k, j, i + 1) * (pow(XI13_c(k, j, i + 1), 2) +
                                        pow(XI23_c(k, j, i + 1), 2))) /
                  6.0 -
              Jacobian_c(k, j, i + 2) *
                  ((2.0 * mu_c(k, j, i + 2) + lambda_c(k, j, i + 2)) *
                       pow(XI33_c(k, j, i + 2), 2) +
                   mu_c(k, j, i + 2) * (pow(XI13_c(k, j, i + 2), 2) +
                                        pow(XI23_c(k, j, i + 2), 2))) /
                  8.0) *
                 u_c_t(3, k, j, i + 2)) /
                pow(h3_c, 2);
      }
    }
  }
  for (j = 1; j <= n2_c; j++) {
    for (k = 1; k <= n1_c; k++) {
      for (i = 1; i <= 6; i++) {
        for (k1 = 1; k1 <= 8; k1++) {
          for (m = 1; m <= 8; m++) {
            // second derivative 33
            // first set equation
            lh_c(k, j, i, 1) =
                lh_c(k, j, i, 1) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 ((2.0 * mu_c(k, j, m) + lambda_c(k, j, m)) *
                      pow(XI13_c(k, j, m), 2) +
                  mu_c(k, j, m) *
                      (pow(XI23_c(k, j, m), 2) + pow(XI33_c(k, j, m), 2))) *
                 u_c_t(1, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI13_c(k, j, m) *
                 XI23_c(k, j, m) * u_c_t(2, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI13_c(k, j, m) *
                 XI33_c(k, j, m) * u_c_t(3, k, j, k1)) /
                    pow(h3_c, 2);

            lh_c(k, j, n3_c + 1 - i, 1) =
                lh_c(k, j, n3_c + 1 - i, 1) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                   lambda_c(k, j, n3_c + 1 - m)) *
                      pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                  mu_c(k, j, n3_c + 1 - m) *
                      (pow(XI23_c(k, j, n3_c + 1 - m), 2) +
                       pow(XI33_c(k, j, n3_c + 1 - m), 2))) *
                 u_c_t(1, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI13_c(k, j, n3_c + 1 - m) * XI23_c(k, j, n3_c + 1 - m) *
                 u_c_t(2, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI13_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
                 u_c_t(3, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2);
            // second set equation
            lh_c(k, j, i, 2) =
                lh_c(k, j, i, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI13_c(k, j, m) *
                 XI23_c(k, j, m) * u_c_t(1, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 ((2.0 * mu_c(k, j, m) + lambda_c(k, j, m)) *
                      pow(XI23_c(k, j, m), 2) +
                  mu_c(k, j, m) *
                      (pow(XI13_c(k, j, m), 2) + pow(XI33_c(k, j, m), 2))) *
                 u_c_t(2, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI23_c(k, j, m) *
                 XI33_c(k, j, m) * u_c_t(3, k, j, k1)) /
                    pow(h3_c, 2);

            lh_c(k, j, n3_c + 1 - i, 2) =
                lh_c(k, j, n3_c + 1 - i, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI13_c(k, j, n3_c + 1 - m) * XI23_c(k, j, n3_c + 1 - m) *
                 u_c_t(1, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                   lambda_c(k, j, n3_c + 1 - m)) *
                      pow(XI23_c(k, j, n3_c + 1 - m), 2) +
                  mu_c(k, j, n3_c + 1 - m) *
                      (pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                       pow(XI33_c(k, j, n3_c + 1 - m), 2))) *
                 u_c_t(2, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI23_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
                 u_c_t(3, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2);
            // third set equation
            lh_c(k, j, i, 3) =
                lh_c(k, j, i, 3) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI13_c(k, j, m) *
                 XI33_c(k, j, m) * u_c_t(1, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 (mu_c(k, j, m) + lambda_c(k, j, m)) * XI23_c(k, j, m) *
                 XI33_c(k, j, m) * u_c_t(2, k, j, k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, m) *
                 ((2.0 * mu_c(k, j, m) + lambda_c(k, j, m)) *
                      pow(XI33_c(k, j, m), 2) +
                  mu_c(k, j, m) *
                      (pow(XI13_c(k, j, m), 2) + pow(XI23_c(k, j, m), 2))) *
                 u_c_t(3, k, j, k1)) /
                    pow(h3_c, 2);

            lh_c(k, j, n3_c + 1 - i, 3) =
                lh_c(k, j, n3_c + 1 - i, 3) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI13_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
                 u_c_t(1, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 (mu_c(k, j, n3_c + 1 - m) + lambda_c(k, j, n3_c + 1 - m)) *
                 XI23_c(k, j, n3_c + 1 - m) * XI33_c(k, j, n3_c + 1 - m) *
                 u_c_t(2, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2) +
                (acof(i, k1, m) * Jacobian_c(k, j, n3_c + 1 - m) *
                 ((2.0 * mu_c(k, j, n3_c + 1 - m) +
                   lambda_c(k, j, n3_c + 1 - m)) *
                      pow(XI33_c(k, j, n3_c + 1 - m), 2) +
                  mu_c(k, j, n3_c + 1 - m) *
                      (pow(XI13_c(k, j, n3_c + 1 - m), 2) +
                       pow(XI23_c(k, j, n3_c + 1 - m), 2))) *
                 u_c_t(3, k, j, n3_c + 1 - k1)) /
                    pow(h3_c, 2);
          }
        }
      }
      // first set equation
      lh_c(k, j, 1, 1) = lh_c(k, j, 1, 1) +
                         (u_c_t(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                               pow(XI13_c(k, j, 1), 2) +
                           mu_c(k, j, 1) * (pow(XI23_c(k, j, 1), 2) +
                                            pow(XI33_c(k, j, 1), 2)))) /
                             pow(h3_c, 2) +
                         (u_c_t(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI13_c(k, j, 1) * XI23_c(k, j, 1)) /
                             pow(h3_c, 2) +
                         (u_c_t(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI13_c(k, j, 1) * XI33_c(k, j, 1)) /
                             pow(h3_c, 2);

      lh_c(k, j, n3_c, 1) =
          lh_c(k, j, n3_c, 1) +
          (u_c_t(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI13_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI23_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c_t(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c_t(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // second set equation
      lh_c(k, j, 1, 2) = lh_c(k, j, 1, 2) +
                         (u_c_t(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI13_c(k, j, 1) * XI23_c(k, j, 1)) /
                             pow(h3_c, 2) +
                         (u_c_t(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                               pow(XI23_c(k, j, 1), 2) +
                           mu_c(k, j, 1) * (pow(XI13_c(k, j, 1), 2) +
                                            pow(XI33_c(k, j, 1), 2)))) /
                             pow(h3_c, 2) +
                         (u_c_t(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI23_c(k, j, 1) * XI33_c(k, j, 1)) /
                             pow(h3_c, 2);

      lh_c(k, j, n3_c, 2) =
          lh_c(k, j, n3_c, 2) +
          (u_c_t(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI23_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c_t(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI23_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI33_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2) +
          (u_c_t(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2);
      // third set equation
      lh_c(k, j, 1, 3) = lh_c(k, j, 1, 3) +
                         (u_c_t(1, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI13_c(k, j, 1) * XI33_c(k, j, 1)) /
                             pow(h3_c, 2) +
                         (u_c_t(2, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          (mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                          XI23_c(k, j, 1) * XI33_c(k, j, 1)) /
                             pow(h3_c, 2) +
                         (u_c_t(3, k, j, 0) * ghcof(1) * Jacobian_c(k, j, 1) *
                          ((2.0 * mu_c(k, j, 1) + lambda_c(k, j, 1)) *
                               pow(XI33_c(k, j, 1), 2) +
                           mu_c(k, j, 1) * (pow(XI13_c(k, j, 1), 2) +
                                            pow(XI23_c(k, j, 1), 2)))) /
                             pow(h3_c, 2);

      lh_c(k, j, n3_c, 3) =
          lh_c(k, j, n3_c, 3) +
          (u_c_t(1, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI13_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c_t(2, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           (mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) * XI23_c(k, j, n3_c) *
           XI33_c(k, j, n3_c)) /
              pow(h3_c, 2) +
          (u_c_t(3, k, j, n3_c + 1) * ghcof(1) * Jacobian_c(k, j, n3_c) *
           ((2.0 * mu_c(k, j, n3_c) + lambda_c(k, j, n3_c)) *
                pow(XI33_c(k, j, n3_c), 2) +
            mu_c(k, j, n3_c) *
                (pow(XI13_c(k, j, n3_c), 2) + pow(XI23_c(k, j, n3_c), 2)))) /
              pow(h3_c, 2);
    }
  }
  //
  // std::cout<<"END \n";
}  // END UPDATE INTERIOR

void forcing(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
             float_sw4 &f1, float_sw4 &f2, float_sw4 &f3) {
  float_sw4 d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3;
  float_sw4 d1d1u1, d2d2u1, d3d3u1, d1d1u2, d2d2u2, d3d3u2, d1d1u3, d2d2u3,
      d3d3u3;
  float_sw4 d1d2u1, d1d3u1, d2d3u1, d1d2u2, d1d3u2, d2d3u2, d1d2u3, d1d3u3,
      d2d3u3;
  float_sw4 d2d1u1, d3d1u1, d3d2u1, d2d1u2, d3d1u2, d3d2u2, d2d1u3, d3d1u3,
      d3d2u3;
  float_sw4 mu, lambda, rho, d1mu, d2mu, d3mu, d1lambda, d2lambda, d3lambda;

  // material property;
  mu = 3.0 + 1.0 * sin(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * sin(x3);
  lambda =
      21.0 + 1.0 * cos(x1 + 0.10) * cos(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  rho = 2.0 + 1.0 * sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 - 0.20);

  // space derivatives;
  d1u1 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2u1 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3u1 = cos(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d1u2 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2u2 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3u2 = sin(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d1u3 = cos(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d2u3 = sin(x1 + 0.20) * cos(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d3u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * sin(x3 + 0.20) * sin(t);

  d1d1u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2d2u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3d3u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d1d1u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2d2u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3d3u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d1d1u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d2d2u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d3d3u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);

  d1d2u1 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2d1u1 = d1d2u1;
  d1d3u1 = -sin(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d3d1u1 = d1d3u1;
  d2d3u1 = cos(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d3d2u1 = d2d3u1;
  d1d2u2 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2d1u2 = d1d2u2;
  d1d3u2 = cos(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d3d1u2 = d1d3u2;
  d2d3u2 = -sin(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d3d2u2 = d2d3u2;
  d1d2u3 = cos(x1 + 0.20) * cos(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d2d1u3 = d1d2u3;
  d1d3u3 = -cos(x1 + 0.20) * sin(x2 + 0.20) * sin(x3 + 0.20) * sin(t);
  d3d1u3 = d1d3u3;
  d2d3u3 = -sin(x1 + 0.20) * cos(x2 + 0.20) * sin(x3 + 0.20) * sin(t);
  d3d2u3 = d2d3u3;

  d1mu = 3.0 * cos(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * sin(x3);
  d2mu = 3.0 * sin(3.0 * x1 + 0.10) * cos(3.0 * x2 + 0.10) * sin(x3);
  d3mu = sin(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * cos(x3);
  d1lambda = -sin(x1 + 0.10) * cos(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  d2lambda = -cos(x1 + 0.10) * sin(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  d3lambda = cos(x1 + 0.10) * cos(x2 + 0.10) * 2.0 * sin(3.0 * x3) *
             cos(3.0 * x3) * 3.0;

  f1 = rho * (-cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
              (cos(t * t) * 4.0 * t * t + sin(t * t) * 2.0)) -
       (2.0 * d1mu * d1u1 + d1lambda * d1u1 + (2.0 * mu + lambda) * d1d1u1 +
        d1lambda * d2u2 + lambda * d1d2u2 + d1lambda * d3u3 + lambda * d1d3u3 +
        d2mu * d1u2 + mu * d2d1u2 + d2mu * d2u1 + mu * d2d2u1 + d3mu * d1u3 +
        mu * d3d1u3 + d3mu * d3u1 + mu * d3d3u1);

  f2 =
      rho * (-sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
             (cos(t * t) * 4.0 * t * t + sin(t * t) * 2.0)) -
      (d1mu * d1u2 + mu * d1d1u2 + d1mu * d2u1 + mu * d1d2u1 + d2lambda * d1u1 +
       lambda * d2d1u1 + 2.0 * d2mu * d2u2 + d2lambda * d2u2 +
       (2.0 * mu + lambda) * d2d2u2 + d2lambda * d3u3 + lambda * d2d3u3 +
       d3mu * d2u3 + mu * d3d2u3 + d3mu * d3u2 + mu * d3d3u2);

  f3 = rho * (-sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t)) -
       (d1mu * d1u3 + mu * d1d1u3 + d1mu * d3u1 + mu * d1d3u1 + d2mu * d2u3 +
        mu * d2d2u3 + d2mu * d3u2 + mu * d2d3u2 + d3lambda * d1u1 +
        lambda * d3d1u1 + d3lambda * d2u2 + lambda * d3d2u2 +
        2.0 * d3mu * d3u3 + d3lambda * d3u3 + (2.0 * mu + lambda) * d3d3u3);
  // print *,"FORCE",x1,x2,x3,f1,f2,f3
}

void forcing_tt(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 t,
                float_sw4 &f1tt, float_sw4 &f2tt, float_sw4 &f3tt) {
  float_sw4 d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3;
  float_sw4 d1d1u1, d2d2u1, d3d3u1, d1d1u2, d2d2u2, d3d3u2, d1d1u3, d2d2u3,
      d3d3u3;
  float_sw4 d1d2u1, d1d3u1, d2d3u1, d1d2u2, d1d3u2, d2d3u2, d1d2u3, d1d3u3,
      d2d3u3;
  float_sw4 d2d1u1, d3d1u1, d3d2u1, d2d1u2, d3d1u2, d3d2u2, d2d1u3, d3d1u3,
      d3d2u3;
  float_sw4 mu, lambda, rho, d1mu, d2mu, d3mu, d1lambda, d2lambda, d3lambda;

  // material property;
  mu = 3.0 + 1.0 * sin(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * sin(x3);
  lambda =
      21.0 + 1.0 * cos(x1 + 0.10) * cos(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  rho = 2.0 + 1.0 * sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 - 0.20);

  //  second time derivatives to space derivatives;
  d1u1 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2u1 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3u1 = cos(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d1u2 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2u2 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3u2 = sin(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) *
         (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d1u3 = cos(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));
  d2u3 = sin(x1 + 0.20) * cos(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));
  d3u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * sin(x3 + 0.20) * (-sin(t));

  d1d1u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2d2u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d3u1 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d1d1u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2d2u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d3u2 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d1d1u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));
  d2d2u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));
  d3d3u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));

  d1d2u1 = -sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2d1u1 = d1d2u1;
  d1d3u1 = -sin(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d1u1 = d1d3u1;
  d2d3u1 = cos(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d2u1 = d2d3u1;
  d1d2u2 = -cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d2d1u2 = d1d2u2;
  d1d3u2 = cos(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d1u2 = d1d3u2;
  d2d3u2 = -sin(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) *
           (-cos(t * t) * 4.0 * t * t - sin(t * t) * 2.0);
  d3d2u2 = d2d3u2;
  d1d2u3 = cos(x1 + 0.20) * cos(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t));
  d2d1u3 = d1d2u3;
  d1d3u3 = -cos(x1 + 0.20) * sin(x2 + 0.20) * sin(x3 + 0.20) * (-sin(t));
  d3d1u3 = d1d3u3;
  d2d3u3 = -sin(x1 + 0.20) * cos(x2 + 0.20) * sin(x3 + 0.20) * (-sin(t));
  d3d2u3 = d2d3u3;

  d1mu = 3.0 * cos(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * sin(x3);
  d2mu = 3.0 * sin(3.0 * x1 + 0.10) * cos(3.0 * x2 + 0.10) * sin(x3);
  d3mu = sin(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * cos(x3);
  d1lambda = -sin(x1 + 0.10) * cos(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  d2lambda = -cos(x1 + 0.10) * sin(x2 + 0.10) * (pow(sin(3.0 * x3), 2));
  d3lambda = cos(x1 + 0.10) * cos(x2 + 0.10) * 2.0 * sin(3.0 * x3) *
             cos(3.0 * x3) * 3.0;

  float_sw4 t4 = t * t * t * t;

  f1tt =
      rho * (-cos(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) *
             (cos(t * t) * (12.0 - 16.0 * (t4)) - sin(t * t) * 48.0 * t * t)) -
      (2.0 * d1mu * d1u1 + d1lambda * d1u1 + (2.0 * mu + lambda) * d1d1u1 +
       d1lambda * d2u2 + lambda * d1d2u2 + d1lambda * d3u3 + lambda * d1d3u3 +
       d2mu * d1u2 + mu * d2d1u2 + d2mu * d2u1 + mu * d2d2u1 + d3mu * d1u3 +
       mu * d3d1u3 + d3mu * d3u1 + mu * d3d3u1);

  f2tt =
      rho * (-sin(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) *
             (cos(t * t) * (12.0 - 16.0 * (t4)) - sin(t * t) * 48.0 * t * t)) -
      (d1mu * d1u2 + mu * d1d1u2 + d1mu * d2u1 + mu * d1d2u1 + d2lambda * d1u1 +
       lambda * d2d1u1 + 2.0 * d2mu * d2u2 + d2lambda * d2u2 +
       (2.0 * mu + lambda) * d2d2u2 + d2lambda * d3u3 + lambda * d2d3u3 +
       d3mu * d2u3 + mu * d3d2u3 + d3mu * d3u2 + mu * d3d3u2);

  f3tt = rho * (-sin(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * (-sin(t))) -
         (d1mu * d1u3 + mu * d1d1u3 + d1mu * d3u1 + mu * d1d3u1 + d2mu * d2u3 +
          mu * d2d2u3 + d2mu * d3u2 + mu * d2d3u2 + d3lambda * d1u1 +
          lambda * d3d1u1 + d3lambda * d2u2 + lambda * d3d2u2 +
          2.0 * d3mu * d3u3 + d3lambda * d3u3 + (2.0 * mu + lambda) * d3d3u3);
}

void update_gp(Farray &Xgrid_c_1, Farray &Xgrid_c_2, Farray &Xgrid_c_3,
               Sarray &u_c, Farray &Xgrid_f_1, Farray &Xgrid_f_2,
               Farray &Xgrid_f_3, Sarray &u_f, float_sw4 tv, PackArgs &a,
               int index) {
  int i, j, k;

  auto n1_c = a.n1_c;
  auto n2_c = a.n2_c;
  auto n3_c = a.n3_c;

  auto n1_f = a.n1_f;
  auto n2_f = a.n2_f;
  auto n3_f = a.n3_f;

  auto nrg = a.nrg;
  // Update ghost point values on the left and right domain
  // fine mesh
  for (i = 1 - nrg; i <= n3_f + nrg; i++) {
    for (j = 1 - nrg; j <= n2_f + nrg; j++) {
      for (k = 1 - nrg; k <= 0; k++) {
        exact_solution(Xgrid_f_1(k), Xgrid_f_2(j), Xgrid_f_3(k, j, i), tv,
                       u_f(1, k, j, i), u_f(2, k, j, i), u_f(3, k, j, i), 1);
      }
      for (k = n1_f + 1; k <= n1_f + nrg; k++) {
        exact_solution(Xgrid_f_1(k), Xgrid_f_2(j), Xgrid_f_3(k, j, i), tv,
                       u_f(1, k, j, i), u_f(2, k, j, i), u_f(3, k, j, i), 1);
      }
    }
  }
  // coarse mesh
  for (i = 1 - nrg; i <= n3_c + nrg; i++) {
    for (j = 1 - nrg; j <= n2_c + nrg; j++) {
      for (k = 1 - nrg; k <= 0; k++) {
        exact_solution(Xgrid_c_1(k), Xgrid_c_2(j), Xgrid_c_3(k, j, i), tv,
                       u_c(1, k, j, i), u_c(2, k, j, i), u_c(3, k, j, i), 0);
      }
      for (k = n1_c + 1; k <= n1_c + nrg; k++) {
        exact_solution(Xgrid_c_1(k), Xgrid_c_2(j), Xgrid_c_3(k, j, i), tv,
                       u_c(1, k, j, i), u_c(2, k, j, i), u_c(3, k, j, i), 0);
      }
    }
  }
  // Update ghost point values on the front and back domain
  // fine mesh
  for (i = 1 - nrg; i <= n3_f + nrg; i++) {
    for (j = 1 - nrg; j <= 0; j++) {
      for (k = 1; k <= n1_f; k++) {
        exact_solution(Xgrid_f_1(k), Xgrid_f_2(j), Xgrid_f_3(k, j, i), tv,
                       u_f(1, k, j, i), u_f(2, k, j, i), u_f(3, k, j, i), 1);
      }
    }
    for (j = n2_f + 1; j <= n2_f + nrg; j++) {
      for (k = 1; k <= n1_f; k++) {
        exact_solution(Xgrid_f_1(k), Xgrid_f_2(j), Xgrid_f_3(k, j, i), tv,
                       u_f(1, k, j, i), u_f(2, k, j, i), u_f(3, k, j, i), 1);
      }
    }
  }
  // coarse mesh
  for (i = 1 - nrg; i <= n3_c + nrg; i++) {
    for (j = 1 - nrg; j <= 0; j++) {
      for (k = 1; k <= n1_c; k++) {
        exact_solution(Xgrid_c_1(k), Xgrid_c_2(j), Xgrid_c_3(k, j, i), tv,
                       u_c(1, k, j, i), u_c(2, k, j, i), u_c(3, k, j, i), 0);
      }
    }
    for (j = n2_c + 1; j <= n2_c + nrg; j++) {
      for (k = 1; k <= n1_c; k++) {
        exact_solution(Xgrid_c_1(k), Xgrid_c_2(j), Xgrid_c_3(k, j, i), tv,
                       u_c(1, k, j, i), u_c(2, k, j, i), u_c(3, k, j, i), 0);
      }
    }
  }
}

void update_traction(Sarray &traction_data, Farray &Xgrid_f_1,
                     Farray &Xgrid_f_2, Farray &Xgrid_f_3, Sarray &u_f,
                     Sarray &mu_f, Sarray &lambda_f, Sarray &Jacobian_f,
                     Sarray &XI13_f, Sarray &XI23_f, Sarray &XI33_f, Farray &Sb,
                     float_sw4 tv, PackArgs &a, int index) {
  // traction B.C. on the top of the fine mesh

  int i, j, k;

  float_sw4 l1 = a.l1;
  float_sw4 l2 = a.l2;
  float_sw4 l3 = a.l3;

  auto n1_f = a.n1_f;
  auto n2_f = a.n2_f;
  auto n3_f = a.n3_f;

  auto nrg = a.nrg;

  auto h1_f = a.h1_f;
  auto h2_f = a.h2_f;
  auto h3_f = a.h3_f;

  float_sw4 traction_rhs[4];
  for (j = 1; j <= n2_f; j++) {
    for (i = 1; i <= n1_f; i++) {
      top_normal_data(Xgrid_f_1(i), Xgrid_f_2(j), Xgrid_f_3(i, j, n3_f), l1, l2,
                      tv, i, j, traction_data, a);
    }
  }
  //
  for (j = 1; j <= n2_f; j++) {
    for (i = 1; i <= n1_f; i++) {
      traction_rhs[0] = traction_rhs[1] = traction_rhs[2] = traction_rhs[3] =
          0.0;
      // 33
      // first set equation
      for (k = 1; k <= 4; k++) {
        // first component
        traction_rhs[1] =
            traction_rhs[1] -
            1.0 / h3_f * Sb(k) * u_f(1, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI23_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) -
            1.0 / h3_f * Sb(k) * u_f(2, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI23_f(i, j, n3_f) -
            1.0 / h3_f * Sb(k) * u_f(3, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI33_f(i, j, n3_f);
      }
      // 31  32
      // first component
      traction_rhs[1] =
          traction_rhs[1] +
          Jacobian_f(i, j, n3_f) *
              (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) / l1 *
              XI13_f(i, j, n3_f) *
              (u_f(1, i - 2, j, n3_f) / 12.0 -
               u_f(1, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(1, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(1, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l2 * XI23_f(i, j, n3_f) *
              (u_f(1, i, j - 2, n3_f) / 12.0 -
               u_f(1, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(1, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(1, i, j + 2, n3_f) / 12.0) /
              h2_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l1 * XI23_f(i, j, n3_f) *
              (u_f(2, i - 2, j, n3_f) / 12.0 -
               u_f(2, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(2, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(2, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) * lambda_f(i, j, n3_f) / l2 *
              XI13_f(i, j, n3_f) *
              (u_f(2, i, j - 2, n3_f) / 12.0 -
               u_f(2, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(2, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(2, i, j + 2, n3_f) / 12.0) /
              h2_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l1 * XI33_f(i, j, n3_f) *
              (u_f(3, i - 2, j, n3_f) / 12.0 -
               u_f(3, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(3, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(3, i + 2, j, n3_f) / 12.0) /
              h1_f;
      // scale
      traction_rhs[1] =
          traction_rhs[1] -
          sqrt(pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2) +
               pow(XI33_f(i, j, n3_f), 2)) *
              Jacobian_f(i, j, n3_f) * traction_data(i, j, 1);
      // 33
      // second set equation
      for (k = 1; k <= 4; k++) {
        // first component
        traction_rhs[2] =
            traction_rhs[2] -
            1.0 / h3_f * Sb(k) * u_f(1, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI23_f(i, j, n3_f) -
            1.0 / h3_f * Sb(k) * u_f(2, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) -
            1.0 / h3_f * Sb(k) * u_f(3, i, j, n3_f + 1 - k) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f);
      }
      // 31  32
      // first component
      traction_rhs[2] =
          traction_rhs[2] +
          Jacobian_f(i, j, n3_f) * lambda_f(i, j, n3_f) / l1 *
              XI23_f(i, j, n3_f) *
              (u_f(1, i - 2, j, n3_f) / 12.0 -
               u_f(1, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(1, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(1, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l2 * XI13_f(i, j, n3_f) *
              (u_f(1, i, j - 2, n3_f) / 12.0 -
               u_f(1, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(1, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(1, i, j + 2, n3_f) / 12.0) /
              h2_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l1 * XI13_f(i, j, n3_f) *
              (u_f(2, i - 2, j, n3_f) / 12.0 -
               u_f(2, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(2, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(2, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) *
              (2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) / l2 *
              XI23_f(i, j, n3_f) *
              (u_f(2, i, j - 2, n3_f) / 12.0 -
               u_f(2, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(2, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(2, i, j + 2, n3_f) / 12.0) /
              h2_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l2 * XI33_f(i, j, n3_f) *
              (u_f(3, i, j - 2, n3_f) / 12.0 -
               u_f(3, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(3, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(3, i, j + 2, n3_f) / 12.0) /
              h2_f;
      // scale
      traction_rhs[2] =
          traction_rhs[2] -
          sqrt(pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2) +
               pow(XI33_f(i, j, n3_f), 2)) *
              Jacobian_f(i, j, n3_f) * traction_data(i, j, 2);
      // 33
      // third set equation
      for (k = 1; k <= 4; k++) {
        // first component
        traction_rhs[3] = traction_rhs[3] -
                          1.0 / h3_f * Sb(k) * u_f(1, i, j, n3_f + 1 - k) *
                              Jacobian_f(i, j, n3_f) *
                              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                              XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) -
                          1.0 / h3_f * Sb(k) * u_f(2, i, j, n3_f + 1 - k) *
                              Jacobian_f(i, j, n3_f) *
                              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                              XI23_f(i, j, n3_f) * XI33_f(i, j, n3_f) -
                          1.0 / h3_f * Sb(k) * u_f(3, i, j, n3_f + 1 - k) *
                              Jacobian_f(i, j, n3_f) *
                              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                                   pow(XI33_f(i, j, n3_f), 2) +
                               mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                                   pow(XI23_f(i, j, n3_f), 2)));
      }
      // 31  32
      // first component
      traction_rhs[3] =
          traction_rhs[3] +
          Jacobian_f(i, j, n3_f) * lambda_f(i, j, n3_f) / l1 *
              XI33_f(i, j, n3_f) *
              (u_f(1, i - 2, j, n3_f) / 12.0 -
               u_f(1, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(1, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(1, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) * lambda_f(i, j, n3_f) / l2 *
              XI33_f(i, j, n3_f) *
              (u_f(2, i, j - 2, n3_f) / 12.0 -
               u_f(2, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(2, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(2, i, j + 2, n3_f) / 12.0) /
              h2_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l1 * XI13_f(i, j, n3_f) *
              (u_f(3, i - 2, j, n3_f) / 12.0 -
               u_f(3, i - 1, j, n3_f) * 2.0 / 3.0 +
               u_f(3, i + 1, j, n3_f) * 2.0 / 3.0 -
               u_f(3, i + 2, j, n3_f) / 12.0) /
              h1_f +
          Jacobian_f(i, j, n3_f) * mu_f(i, j, n3_f) / l2 * XI23_f(i, j, n3_f) *
              (u_f(3, i, j - 2, n3_f) / 12.0 -
               u_f(3, i, j - 1, n3_f) * 2.0 / 3.0 +
               u_f(3, i, j + 1, n3_f) * 2.0 / 3.0 -
               u_f(3, i, j + 2, n3_f) / 12.0) /
              h2_f;
      // scale
      traction_rhs[3] =
          traction_rhs[3] -
          sqrt(pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2) +
               pow(XI33_f(i, j, n3_f), 2)) *
              Jacobian_f(i, j, n3_f) * traction_data(i, j, 3);

      // Update ghost point at the traction boundary
      // compute the determinant of the coefficient matrix
      float_sw4 mat_det =
          Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI13_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI23_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
              Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI23_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI13_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
              Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI33_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2))) +
          Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
              XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) +
          Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
              XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
              XI23_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) -
          Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
              XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI23_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI13_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
              Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) -
          Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI13_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI23_f(i, j, n3_f), 2) + pow(XI33_f(i, j, n3_f), 2))) *
              Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
              XI33_f(i, j, n3_f) -
          Jacobian_f(i, j, n3_f) *
              ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                   pow(XI33_f(i, j, n3_f), 2) +
               mu_f(i, j, n3_f) *
                   (pow(XI13_f(i, j, n3_f), 2) + pow(XI23_f(i, j, n3_f), 2))) *
              Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
              XI23_f(i, j, n3_f) * Jacobian_f(i, j, n3_f) *
              (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
              XI23_f(i, j, n3_f);
      // update the first component
      u_f(1, i, j, n3_f + 1) =
          h3_f / Sb(0) / mat_det *
          ((Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI33_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI23_f(i, j, n3_f), 2))) -
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI23_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f)) *
               traction_rhs[1] +
           (Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI33_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI23_f(i, j, n3_f), 2)))) *
               traction_rhs[2] +
           (Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2)))) *
               traction_rhs[3]);
      // update the second component
      u_f(2, i, j, n3_f + 1) =
          h3_f / Sb(0) / mat_det *
          ((Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI23_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI33_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI33_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI23_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI23_f(i, j, n3_f)) *
               traction_rhs[1] +
           (Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI23_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI33_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI23_f(i, j, n3_f), 2))) -
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI33_f(i, j, n3_f)) *
               traction_rhs[2] +
           (Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI33_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI23_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI23_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f)) *
               traction_rhs[3]);
      // update the third component
      u_f(3, i, j, n3_f + 1) =
          h3_f / Sb(0) / mat_det *
          ((Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI33_f(i, j, n3_f)) *
               traction_rhs[1] +
           (Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI33_f(i, j, n3_f) -
            Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI23_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI23_f(i, j, n3_f) *
                XI33_f(i, j, n3_f)) *
               traction_rhs[2] +
           (Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI13_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI23_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) *
                Jacobian_f(i, j, n3_f) *
                ((2.0 * mu_f(i, j, n3_f) + lambda_f(i, j, n3_f)) *
                     pow(XI23_f(i, j, n3_f), 2) +
                 mu_f(i, j, n3_f) * (pow(XI13_f(i, j, n3_f), 2) +
                                     pow(XI33_f(i, j, n3_f), 2))) -
            Jacobian_f(i, j, n3_f) * (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) *
                XI13_f(i, j, n3_f) * XI23_f(i, j, n3_f) *
                Jacobian_f(i, j, n3_f) *
                (lambda_f(i, j, n3_f) + mu_f(i, j, n3_f)) * XI13_f(i, j, n3_f) *
                XI23_f(i, j, n3_f)) *
               traction_rhs[3]);
    }
  }
}  // TRACTION

void top_normal_data(float_sw4 x1, float_sw4 x2, float_sw4 x3, float_sw4 l1,
                     float_sw4 l2, float_sw4 t, int i, int j, Sarray &g,
                     PackArgs &a) {
  auto dim = a.dim;
  auto amp = a.amp;
  auto peak = a.peak;

  Farray St(1, dim, 1, dim);
  Farray n(1, dim);

  float_sw4 d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3, mu, lambda;

  // material property;
  mu = 3.0 + 1.0 * sin(3.0 * x1 + 0.10) * sin(3.0 * x2 + 0.10) * sin(x3);
  lambda =
      21.0 + 1.0 * cos(x1 + 0.10) * cos(x2 + 0.10) * (pow(sin(3.0 * x3), 2));

  // space derivatives;
  d1u1 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2u1 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3u1 = cos(x1 + 0.30) * sin(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d1u2 = cos(x1 + 0.30) * cos(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d2u2 = -sin(x1 + 0.30) * sin(x2 + 0.30) * sin(x3 + 0.20) * cos(t * t);
  d3u2 = sin(x1 + 0.30) * cos(x2 + 0.30) * cos(x3 + 0.20) * cos(t * t);
  d1u3 = cos(x1 + 0.20) * sin(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d2u3 = sin(x1 + 0.20) * cos(x2 + 0.20) * cos(x3 + 0.20) * sin(t);
  d3u3 = -sin(x1 + 0.20) * sin(x2 + 0.20) * sin(x3 + 0.20) * sin(t);

  St(1, 1) = (2.0 * mu + lambda) * d1u1 + lambda * d2u2 + lambda * d3u3;

  St(1, 2) = mu * d1u2 + mu * d2u1;

  St(1, 3) = mu * d1u3 + mu * d3u1;

  St(2, 1) = mu * d1u2 + mu * d2u1;

  St(2, 2) = lambda * d1u1 + (2.0 * mu + lambda) * d2u2 + lambda * d3u3;

  St(2, 3) = mu * d2u3 + mu * d3u2;

  St(3, 1) = mu * d1u3 + mu * d3u1;

  St(3, 2) = mu * d2u3 + mu * d3u2;

  St(3, 3) = lambda * d1u1 + lambda * d2u2 + (2.0 * mu + lambda) * d3u3;

  n(1) = amp * exp(-pow((x1 / l1 - 0.50), 2) / peak) * (2.0 * x1 / l1 - 1.0) /
         peak / l1;
  n(2) = amp * exp(-pow((x2 / l2 - 0.50), 2) / peak) * (2.0 * x2 / l2 - 1.0) /
         peak / l2;
  n(3) = 1.0;
  float_sw4 l2norm = sqrt(n(1) * n(1) + n(2) * n(2) + n(3) * n(3));
  for (int i = 1; i <= dim; i++) n(i) = n(i) / l2norm;

  // n = n/sqrt(n(1)**2+n(2)**2+n(3)**2);

  g(i, j, 1) = St(1, 1) * n(1) + St(1, 2) * n(2) + St(1, 3) * n(3);
  g(i, j, 2) = St(2, 1) * n(1) + St(2, 2) * n(2) + St(2, 3) * n(3);
  g(i, j, 3) = St(3, 1) * n(1) + St(3, 2) * n(2) + St(3, 3) * n(3);
}

void update_dirichlet_bc(Farray &Xgrid_c_1, Farray &Xgrid_c_2,
                         Farray &Xgrid_c_3, Sarray &u_c, float_sw4 tv,
                         PackArgs &a, int index)

{
  auto nrg = a.nrg;

  int n1_c = a.n1_c;
  int n2_c = a.n2_c;
  int n3_c = a.n3_c;

  int n1_f = a.n1_f;
  int n2_f = a.n2_f;
  int n3_f = a.n3_f;

  int i, j, k;
  // fine mesh
  for (j = 1 - nrg; j <= n2_f + nrg; j++) {
    for (k = 1 - nrg; k <= n1_f + nrg; k++) {
      i = 1;
      // call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(i),tv, &
      // u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
      i = n3_f;
      // call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
      //          u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
    }
  }
  // coarse mesh
  for (j = 1 - nrg; j <= n2_c + nrg; j++) {
    for (k = 1 - nrg; k <= n1_c + nrg; k++) {
      i = 1;
      exact_solution(Xgrid_c_1(k), Xgrid_c_2(j), Xgrid_c_3(k, j, i), tv,
                     u_c(1, k, j, i), u_c(2, k, j, i), u_c(3, k, j, i), 0);
      i = n3_c;
      // call
      // exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,3),tv, &
      // u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
    }
  }
}
