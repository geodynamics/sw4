#include "EW.h"
#include "sw4.h"
//-----------------------------------------------------------------------
void EW::solerr3_ci(int ib, int ie, int jb, int je, int kb, int ke, float_sw4 h,
                    float_sw4* __restrict__ uex, float_sw4* __restrict__ u,
                    float_sw4& li, float_sw4& l2, float_sw4& xli,
                    float_sw4 zmin, float_sw4 x0, float_sw4 y0, float_sw4 z0,
                    float_sw4 radius, int imin, int imax, int jmin, int jmax,
                    int kmin, int kmax) {
  li = 0;
  l2 = 0;
  xli = 0;
  float_sw4 sradius2 = radius * radius;
  if (radius < 0) sradius2 = -sradius2;
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const float_sw4 h3 = h * h * h;
  const size_t nijk = nij * (ke - kb + 1);

  for (int c = 0; c < 3; c++) {
    float_sw4 liloc = 0, l2loc = 0, xliloc = 0;
#pragma omp parallel for reduction(max : liloc, xliloc) reduction(+ : l2loc)
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          if (((i - 1) * h - x0) * ((i - 1) * h - x0) +
                  ((j - 1) * h - y0) * ((j - 1) * h - y0) +
                  ((k - 1) * h + zmin - z0) * ((k - 1) * h + zmin - z0) >
              sradius2) {
            size_t ind = i - ib + ni * (j - jb) + nij * (k - kb) + nijk * c;
            float_sw4 err = fabs(u[ind] - uex[ind]);
            //		     if( fabs(u[ind]-uex[ind])>liloc )
            //			liloc = fabs(u[ind]-uex[ind]);
            if (err > liloc) liloc = err;
            if (fabs(uex[ind]) > xliloc) xliloc = fabs(uex[ind]);
            // if (uex[ind] > xliloc) xliloc = uex[ind]; // ORG PRE CURVIMR
            //		     l2loc += h3*(u[ind]-uex[ind])*(u[ind]-uex[ind]);
            l2loc += h3 * err * err;
          }
        }
    l2 += l2loc;
    li = liloc > li ? liloc : li;
    xli = xliloc > xli ? xliloc : xli;
    // std::cout<<"SOL2 "<<c<<" "<<xli<<" "<<xliloc<<"\n";
  }
}

//-----------------------------------------------------------
void EW::solerrgp_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4 h, float_sw4* __restrict__ uex,
                     float_sw4* __restrict__ u, float_sw4& li, float_sw4& l2) {
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const float_sw4 h3 = h * h * h;
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);
  li = 0;
  l2 = 0;
  //   int k=kfirst+1;
  for (int c = 0; c < 3; c++) {
    float_sw4 liloc = 0, l2loc = 0;
#pragma omp parallel for reduction(max : liloc) reduction(+ : l2loc)
    for (int j = jfirst + 2; j <= jlast - 2; j++)
      for (int i = ifirst + 2; i <= ilast - 2; i++) {
        // exact solution in array 'uex'
        size_t ind = base + i + ni * j + nij * (kfirst + 1) + nijk * c;
        float_sw4 err = abs(u[ind] - uex[ind]);
        if (liloc < err) {
          liloc = err;
          //	       cout << " num err= " << err << " at " << i << " " << j <<
          //" " << kfirst+1 << ", c= " << c << endl; 	       cout << "u, uex
          //= " << u[ind] << " " << uex[ind] << endl;
        }
        //	    liloc = liloc > err ? liloc : err;
        l2loc += h3 * err * err;
      }
    li = li > liloc ? li : liloc;
    l2 += l2loc;
  }
  //      k=klast-1;
  for (int c = 0; c < 3; c++) {
    float_sw4 liloc = 0, l2loc = 0;
#pragma omp parallel for reduction(max : liloc) reduction(+ : l2loc)
    for (int j = jfirst + 2; j <= jlast - 2; j++)
      for (int i = ifirst + 2; i <= ilast - 2; i++) {
        size_t ind = base + i + ni * j + nij * (klast - 1) + nijk * c;
        float_sw4 err = abs(u[ind] - uex[ind]);
        if (liloc < err) {
          liloc = err;
          //	       cout << " num err= " << err << " at " << i << " " << j <<
          //" " << klast-1 << ", c= " << c <<  endl; 	       cout << "u, uex
          //= " << u[ind] << " " << uex[ind] << endl;
        }
        //	    liloc = liloc > err ? liloc : err;
        l2loc += h3 * err * err;
      }
    li = li > liloc ? li : liloc;
    l2 += l2loc;
  }
}

//-----------------------------------------------------------------------
void EW::solerr3c_ci(int ib, int ie, int jb, int je, int kb, int ke,
                     float_sw4* __restrict__ uex, float_sw4* __restrict__ u,
                     float_sw4* __restrict__ x, float_sw4* __restrict__ y,
                     float_sw4* __restrict__ z, float_sw4* __restrict__ jac,
                     float_sw4& li, float_sw4& l2, float_sw4& xli, float_sw4 x0,
                     float_sw4 y0, float_sw4 z0, float_sw4 radius, int imin,
                     int imax, int jmin, int jmax, int kmin, int kmax,
                     int usesg, float_sw4* __restrict__ strx,
                     float_sw4* __restrict__ stry) {
  li = 0;
  l2 = 0;
  xli = 0;
  float_sw4 sradius2 = radius * radius;
  if (radius < 0) sradius2 = -sradius2;
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const size_t nijk = nij * (ke - kb + 1);
  const size_t base = -(ib + ni * jb + nij * kb);

  for (int c = 0; c < 3; c++) {
    float_sw4 liloc = 0, xliloc = 0, l2loc = 0;
#pragma omp parallel for reduction(max : liloc, xliloc) reduction(+ : l2loc)
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 dist = (x[ind] - x0) * (x[ind] - x0) +
                           (y[ind] - y0) * (y[ind] - y0) +
                           (z[ind] - z0) * (z[ind] - z0);
          if (dist > sradius2) {
            size_t ind3 = ind + nijk * c;
            float_sw4 err = fabs(u[ind3] - uex[ind3]);
            liloc = liloc > err ? liloc : err;
            xliloc = xliloc > fabs(uex[ind3]) ? xliloc : fabs(uex[ind3]);
            // xliloc = xliloc > uex[ind3] ? xliloc : uex[ind3]; // ORG PRE
            // CURVIMR std::cout<<"SOL"<<ind3<<" "<<xliloc<<"
            // "<<uex[ind3]<<"\n";
            if (usesg != 1)
              l2loc += jac[ind] * err * err;
            else
              l2loc += jac[ind] * err * err / (strx[i - ib] * stry[j - jb]);
          }
        }
    li = liloc > li ? liloc : li;
    l2 += l2loc;
    xli = xliloc > xli ? xliloc : xli;
  }
}

#ifdef XL_BUG_152435FIXED
//-----------------------------------------------------------------------
// This version does not compile due to a bug in XL :
//
// LLNL: SW4 ICE in clangtana with -qsmp=omp (152435).
// Use workaround below until this is fixed.
void EW::meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4* __restrict__ met,
                     float_sw4* __restrict__ metex, float_sw4* __restrict__ jac,
                     float_sw4* __restrict__ jacex, float_sw4 li[5],
                     float_sw4 l2[5], int imin, int imax, int jmin, int jmax,
                     int kmin, int kmax, float_sw4 h) {
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

  for (int c = 0; c < 5; c++) li[c] = l2[c] = 0;

  const float_sw4 isqh = 1 / sqrt(h);
  const float_sw4 ih3 = 1 / (h * h * h);
#pragma omp parallel
  for (int c = 0; c < 5; c++)
#pragma omp for reduction(max : li[:5]) reduction(+ : l2[:5])
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 err;
          if (c < 4)
            err = fabs(met[ind + c * nijk] - metex[ind + c * nijk]) * isqh;
          else
            err = fabs(jac[ind] - jacex[ind]) * ih3;
          li[c] = li[c] > err ? li[c] : err;
          l2[c] += jacex[ind] * err * err;
        }
}
#else
//-----------------------------------------------------------------------
void EW::meterr4c_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
                     int klast, float_sw4* __restrict__ met,
                     float_sw4* __restrict__ metex, float_sw4* __restrict__ jac,
                     float_sw4* __restrict__ jacex, float_sw4 li[5],
                     float_sw4 l2[5], int imin, int imax, int jmin, int jmax,
                     int kmin, int kmax, float_sw4 h) {
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);

  for (int c = 0; c < 5; c++) li[c] = l2[c] = 0;

  float_sw4 tmp_li;
  float_sw4 tmp_l2;
  const float_sw4 isqh = 1 / sqrt(h);
  const float_sw4 ih3 = 1 / (h * h * h);
  for (int c = 0; c < 5; c++) {
    tmp_li = li[c];
    tmp_l2 = l2[c];
#pragma omp parallel for reduction(max : tmp_li) reduction(+ : tmp_l2)
    for (size_t k = kmin; k <= kmax; k++)
      for (size_t j = jmin; j <= jmax; j++)
        for (size_t i = imin; i <= imax; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 err;
          if (c < 4)
            err = fabs(met[ind + c * nijk] - metex[ind + c * nijk]) * isqh;
          else
            err = fabs(jac[ind] - jacex[ind]) * ih3;
          li[c] = li[c] > err ? li[c] : err;
          l2[c] += jacex[ind] * err * err;
        }
    li[c] = tmp_li;
    l2[c] = tmp_l2;
  }
}
#endif
