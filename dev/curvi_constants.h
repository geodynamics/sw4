#ifndef __CURVI_CONSTANTS_
#define __CURVI_CONSTANTS_
const int dim = 3;
const float_sw4 tn = 1.0;
const float_sw4 pi = M_PI;
const float_sw4 l1 = 2.0 * pi, l2 = 2.0 * pi, l3 = 2.0 * pi;  // space interval
const float_sw4 int_pos =
    0.50;  // This is the position in r3 where the interface is located.
const int nrg = 5;
const int n1_c = 25;
const int n2_c = 25;
const float_sw4 h1phy_c = l1 / (n1_c - 1),
                h1phy_f = h1phy_c * 0.50;  // mesh size in physical space, x
const float_sw4 h2phy_c = l2 / (n2_c - 1),
                h2phy_f = h2phy_c * 0.50;  // mesh size in physical space, y
const int n3_c = std::ceil(int_pos * l3 / h1phy_c) +
                 1;  //  number of grid points in direction-3
const int n1_f = n1_c * 2 - 1, n2_f = n2_c * 2 - 1,
          n3_f = std::ceil((1.0 - int_pos) * l3 / (h1phy_f)) + 1;
const float_sw4 h1_c = 1.0 / (n1_c - 1), h2_c = 1.0 / (n2_c - 1),
                h3_c = 1.0 / (n3_c - 1);
const float_sw4 h1_f = 1.0 / (n1_f - 1), h2_f = 1.0 / (n2_f - 1),
                h3_f = 1.0 / (n3_f - 1);
#endif
