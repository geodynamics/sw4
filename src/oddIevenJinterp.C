#include <cstdio>

#include "Mspace.h"
#include "Sarray.h"
#include "caliper.h"
#include "policies.h"
#include "sw4.h"

//--------------------- Jacobi ---------------------
void oddIevenJinterpJacobi(float_sw4 rmax[6], Sarray &Uf, Sarray &UfNew,
                           Sarray &Uc, Sarray &Morc, Sarray &Mlrc, Sarray &Morf,
                           Sarray &Mlrf, Sarray &Unextf, Sarray &UnextcInterp,
                           sw4_type a_iStart[], sw4_type a_iEnd[], sw4_type a_jStart[],
                           sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
                           sw4_type a_iStartInt[], sw4_type a_iEndInt[],
                           sw4_type a_jStartInt[], sw4_type a_jEndInt[], sw4_type gf, sw4_type gc,
                           sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
                           float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[],
                           float_sw4 a_ghcof[]) {
  // tmp
  //  printf("Inside oddIevenJinterp! ");

  // sw4_type icb = a_iStartInt[gc];
  sw4_type ifb = a_iStartInt[gf];
  if (ifb % 2 == 0) ifb++;  // make sure ifb is odd

  // sw4_type ice = a_iEndInt[gc];
  sw4_type ife = a_iEndInt[gf];

  // sw4_type jcb = a_jStartInt[gc];
  sw4_type jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++;  // make sure jfb is even

  // sw4_type jce = a_jEndInt[gc];
  sw4_type jfe = a_jEndInt[gf];

  float_sw4 nuf =
      a_Dt * a_Dt / (cof * hf * hf);  // cof=12 for the predictor, cof=1 for the
                                      // corrector (argument to this routine)
  float_sw4 nuc = a_Dt * a_Dt / (cof * hc * hc);
  // float_sw4 ihc = 1/hc, ihf=1/hf;

  const float_sw4 i16 = 1.0 / 16;
  // const float_sw4 i256 = 1.0/256;
  //  const float_sw4 i1024 = 1.0/1024;

  // residuals
  float_sw4 rmax1 = 0, rmax2 = 0, rmax3 = 0;

#pragma omp parallel for reduction(max : rmax1, rmax2, rmax3)
  for (sw4_type j = jfb; j <= jfe; j += 2)
#pragma omp simd
    for (sw4_type i = ifb; i <= ife; i += 2) {
      sw4_type ic, jc;
      float_sw4 b1, a11, r3;
      ic = (i + 1) / 2;
      jc = j / 2;

      //      for( sw4_type c=1 ; c <= 2 ; c++ ) // updated components 1,2 of the
      //      ghost point value of Uf
      //      {
      // Un-roll c-loop
      sw4_type c = 1;
      // All Uc terms
      // NOTE: Uc is not changed by this routine, so the Uc-dependence could be
      // pre-computed
      b1 = UnextcInterp(c, i, j, 1) +
           nuc * a_ghcof[0] * i16 *
               (-Uc(c, ic, jc - 1, 0) * Morc(ic, jc - 1, 1) +
                9 * Uc(c, ic, jc, 0) * Morc(ic, jc, 1) +
                9 * Uc(c, ic, jc + 1, 0) * Morc(ic, jc + 1, 1) -
                Uc(c, ic, jc + 2, 0) * Morc(ic, jc + 2, 1)) -
           Unextf(c, i, j, nkf);

      //      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
      a11 = nuf * a_ghcof[0] * Morf(i, j, nkf);

      // update ghost point value Uf(c,i,j,nkf+1)
      //         Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c, i, j, nkf + 1) =
          relax * b1 / a11 + (1 - relax) * Uf(c, i, j, nkf + 1);
      // change in ghost point value
      r3 = UfNew(c, i, j, nkf + 1) - Uf(c, i, j, nkf + 1);
      rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);

      c = 2;
      // All Uc terms
      // NOTE: Uc is not changed by this routine, so the Uc-dependence could be
      // pre-computed
      b1 = UnextcInterp(c, i, j, 1) +
           nuc * a_ghcof[0] * i16 *
               (-Uc(c, ic, jc - 1, 0) * Morc(ic, jc - 1, 1) +
                9 * Uc(c, ic, jc, 0) * Morc(ic, jc, 1) +
                9 * Uc(c, ic, jc + 1, 0) * Morc(ic, jc + 1, 1) -
                Uc(c, ic, jc + 2, 0) * Morc(ic, jc + 2, 1)) -
           Unextf(c, i, j, nkf);

      //      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
      a11 = nuf * a_ghcof[0] * Morf(i, j, nkf);

      // update ghost point value Uf(c,i,j,nkf+1)
      //         Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c, i, j, nkf + 1) =
          relax * b1 / a11 + (1 - relax) * Uf(c, i, j, nkf + 1);
      // change in ghost point value
      r3 = UfNew(c, i, j, nkf + 1) - Uf(c, i, j, nkf + 1);
      rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
      //      } // end for c=1,2

      // work on component 3 of the ghost point value of Uf
      // ic = (i+1)/2;
      // jc = j/2;
      // All Uc terms
      // right hand side is mismatch in displacement
      // NOTE: Uc is not changed by this routine, so the Uc-dependence could be
      // pre-computed
      b1 = UnextcInterp(3, i, j, 1) +
           nuc * a_ghcof[0] * i16 *
               (-Uc(3, ic, jc - 1, 0) * Mlrc(ic, jc - 1, 1) +
                9 * Uc(3, ic, jc, 0) * Mlrc(ic, jc, 1) +
                9 * Uc(3, ic, jc + 1, 0) * Mlrc(ic, jc + 1, 1) -
                Uc(3, ic, jc + 2, 0) * Mlrc(ic, jc + 2, 1)) -
           Unextf(3, i, j, nkf);

      //    a11 =
      //    nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf)); //
      //    no str
      a11 = nuf * a_ghcof[0] * Mlrf(i, j, nkf);  // no str

      //    r3 = Uf(3,i,j,nkf+1); // save previous value for relaxation below
      // solve for the ghost point value Uf(3,i,j,nkf+1)
      //    Uf(3,i,j,nkf+1) = b1/a11;
      UfNew(3, i, j, nkf + 1) =
          relax * b1 / a11 + (1 - relax) * Uf(3, i, j, nkf + 1);
      r3 = UfNew(3, i, j, nkf + 1) - Uf(3, i, j, nkf + 1);
      rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);

    }  // end for i odd, j even

// update Uf
#pragma omp parallel
  for (sw4_type c = 1; c <= 3; c++)
#pragma omp for
    for (sw4_type j = jfb; j <= jfe; j += 2)
#pragma omp simd
      for (sw4_type i = ifb; i <= ife; i += 2) {
        Uf(c, i, j, nkf + 1) = UfNew(c, i, j, nkf + 1);
      }

  rmax[3] = rmax1;
  rmax[4] = rmax2;
  rmax[5] = rmax3;
}  // end oddIevenJinterpJacobi

void oddIevenJinterpJacobiOpt(
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> &rmax1,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> &rmax2,
    RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> &rmax3,
    float_sw4 *__restrict__ a_uf, float_sw4 *__restrict__ a_ufnew,
    float_sw4 *__restrict__ a_uc, float_sw4 *__restrict__ a_morc,
    float_sw4 *__restrict__ a_mlrc, float_sw4 *__restrict__ a_morf,
    float_sw4 *__restrict__ a_mlrf, float_sw4 *__restrict__ a_unextf,
    float_sw4 *__restrict__ a_uncsw4_type, sw4_type a_iStart[], sw4_type a_iEnd[],
    sw4_type a_jStart[], sw4_type a_jEnd[], sw4_type a_kStart[], sw4_type a_kEnd[],
    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
    sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
    float_sw4 cof, float_sw4 relax, float_sw4 a_sbop[], float_sw4 a_ghcof[]) {
  SW4_MARK_FUNCTION;
  const sw4_type iStartC = a_iStart[gc];
  const sw4_type jStartC = a_jStart[gc];
  const sw4_type kStartC = a_kStart[gc];

  const sw4_type iEndC = a_iEnd[gc];
  const sw4_type jEndC = a_jEnd[gc];
  const sw4_type kEndC = a_kEnd[gc];

  const sw4_type iStartF = a_iStart[gf];
  const sw4_type jStartF = a_jStart[gf];
  const sw4_type kStartF = a_kStart[gf];
  const sw4_type iEndF = a_iEnd[gf];
  const sw4_type jEndF = a_jEnd[gf];
  const sw4_type kEndF = a_kEnd[gf];

  // Bf indexing
  const sw4_type niF = iEndF - iStartF + 1;
  const sw4_type nijF = niF * (jEndF - jStartF + 1);
  const sw4_type nijk_bf = nijF * (1);  // only one k-plane
  const sw4_type base3_bf =
      (iStartF + niF * jStartF + nijF * nkf + nijk_bf);  // only one k=nkf
#define Unextf(c, i, j, k) \
  a_unextf[-base3_bf + i + niF * (j) + nijF * (k) + nijk_bf * (c)]

  const sw4_type nijk_uncsw4_type = nijF * (1);  // only one k-plane
  const sw4_type base3_uncsw4_type =
      (iStartF + niF * jStartF + nijF * 1 + nijk_uncsw4_type);  // only k=1
#define UnextcInterp(c, i, j, k) \
  a_uncsw4_type[-base3_uncsw4_type + i + niF * (j) + nijF * (k) + nijk_uncsw4_type * (c)]

  const sw4_type base_mufs =
      (iStartF + niF * jStartF + nijF * nkf);  // only one k=nkf
#define Morf(i, j, k) \
  a_morf[-base_mufs + i + niF * (j) + nijF * (k)]  // same size as Mufs
#define Mlrf(i, j, k) \
  a_mlrf[-base_mufs + i + niF * (j) + nijF * (k)]  // same size as Mufs

  const sw4_type niC = iEndC - iStartC + 1;
  const sw4_type nijC = niC * (jEndC - jStartC + 1);
  const sw4_type base_morc = (iStartC + niC * jStartC + nijC * 1);  // only one k=1
#define Morc(i, j, k) a_morc[-base_morc + i + niC * (j) + nijC * (k)]
#define Mlrc(i, j, k) \
  a_mlrc[-base_morc + i + niC * (j) + nijC * (k)]  // same size as Morc

  const sw4_type nijk_uc = nijC * (kEndC - kStartC + 1);
  const sw4_type base3_uc = (iStartC + niC * jStartC + nijC * kStartC +
                        nijk_uc * 1);  // c-index has base=1
#define Uc(c, i, j, k) \
  a_uc[-base3_uc + i + niC * (j) + nijC * (k) + nijk_uc * (c)]

  const sw4_type nijk_uf = nijF * (kEndF - kStartF + 1);
  const sw4_type base3_uf = (iStartF + niF * jStartF + nijF * kStartF +
                        nijk_uf * 1);  // c-index has base=1
#define Uf(c, i, j, k) \
  a_uf[-base3_uf + i + niF * (j) + nijF * (k) + nijk_uf * (c)]

  const sw4_type nijk_ufnew = nijF * (1);  // only one k-plane
  const sw4_type base3_ufnew =
      (iStartF + niF * jStartF + nijF * (nkf + 1) +
       nijk_ufnew * 1);  // only one k=nkf+1; c-index has base=1
#define UfNew(c, i, j, k) \
  a_ufnew[-base3_ufnew + i + niF * (j) + nijF * (k) + nijk_ufnew * (c)]

  // previous stuff
  // sw4_type icb = a_iStartInt[gc];
  sw4_type ifb = a_iStartInt[gf];
  if (ifb % 2 == 0) ifb++;  // make sure ifb is odd

  // sw4_type ice = a_iEndInt[gc];
  sw4_type ife = a_iEndInt[gf];

  //  sw4_type jcb = a_jStartInt[gc];
  sw4_type jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++;  // make sure jfb is even

  //  sw4_type jce = a_jEndInt[gc];
  sw4_type jfe = a_jEndInt[gf];

  float_sw4 nuf =
      a_Dt * a_Dt / (cof * hf * hf);  // cof=12 for the predictor, cof=1 for the
                                      // corrector (argument to this routine)
  float_sw4 nuc = a_Dt * a_Dt / (cof * hc * hc);
  //  float_sw4 ihc = 1/hc, ihf=1/hf;

  const float_sw4 i16 = 1.0 / 16;
  //  const float_sw4 i256 = 1.0/256;
  // const float_sw4 i1024 = 1.0/1024;

  // residuals
  //  RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> rrmax[3]={};
  // RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> rmax1(0);
  // RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> rmax2(0);
  // RAJA::ReduceMax<REDUCTION_POLICY, float_sw4> rmax3(0);

  // Lines below are not correct , incorrect stride
  RAJA::RangeSegment j_range(jfb, jfe + 1);
  RAJA::RangeSegment i_range(ifb, ife + 1);
  RAJA::TypedRangeStrideSegment<long> j_srange(jfb, jfe + 1, 2);
  RAJA::TypedRangeStrideSegment<long> i_srange(ifb, ife + 1, 2);
  RAJA::kernel<ODDIEVENJ_EXEC_POL1_ASYNC>(
      RAJA::make_tuple(j_srange, i_srange),
      [=] RAJA_DEVICE(sw4_type j, sw4_type i) {
        // float_sw4 rmax1=0, rmax2=0, rmax3=0;
        // #pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
        //   for( sw4_type j=jfb ; j <= jfe ; j+=2 )
        // #pragma omp simd
        //     for( sw4_type i=ifb ; i <= ife ; i+=2 )
        //     {
        sw4_type ic, jc;
        float_sw4 b1, a11, r3;
        ic = (i + 1) / 2;
        jc = j / 2;

        //      for( sw4_type c=1 ; c <= 2 ; c++ ) // updated components 1,2 of the
        //      ghost point value of Uf
        //      {
        // Un-roll c-loop
        sw4_type c = 1;
        // All Uc terms
        // NOTE: Uc is not changed by this routine, so the Uc-dependence could
        // be pre-computed
        b1 = UnextcInterp(c, i, j, 1) +
             nuc * a_ghcof[0] * i16 *
                 (-Uc(c, ic, jc - 1, 0) * Morc(ic, jc - 1, 1) +
                  9 * Uc(c, ic, jc, 0) * Morc(ic, jc, 1) +
                  9 * Uc(c, ic, jc + 1, 0) * Morc(ic, jc + 1, 1) -
                  Uc(c, ic, jc + 2, 0) * Morc(ic, jc + 2, 1)) -
             Unextf(c, i, j, nkf);

        //      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
        a11 = nuf * a_ghcof[0] * Morf(i, j, nkf);

        // update ghost point value Uf(c,i,j,nkf+1)
        //         Uf(c,i,j,nkf+1) = b1/a11;
        UfNew(c, i, j, nkf + 1) =
            relax * b1 / a11 + (1 - relax) * Uf(c, i, j, nkf + 1);
        // change in ghost point value
        r3 = UfNew(c, i, j, nkf + 1) - Uf(c, i, j, nkf + 1);
        Uf(c, i, j, nkf + 1) = UfNew(c, i, j, nkf + 1);
        // rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
        rmax1.max(fabs(r3));
        c = 2;
        // All Uc terms
        // NOTE: Uc is not changed by this routine, so the Uc-dependence could
        // be pre-computed
        b1 = UnextcInterp(c, i, j, 1) +
             nuc * a_ghcof[0] * i16 *
                 (-Uc(c, ic, jc - 1, 0) * Morc(ic, jc - 1, 1) +
                  9 * Uc(c, ic, jc, 0) * Morc(ic, jc, 1) +
                  9 * Uc(c, ic, jc + 1, 0) * Morc(ic, jc + 1, 1) -
                  Uc(c, ic, jc + 2, 0) * Morc(ic, jc + 2, 1)) -
             Unextf(c, i, j, nkf);

        //      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
        a11 = nuf * a_ghcof[0] * Morf(i, j, nkf);

        // update ghost point value Uf(c,i,j,nkf+1)
        //         Uf(c,i,j,nkf+1) = b1/a11;
        UfNew(c, i, j, nkf + 1) =
            relax * b1 / a11 + (1 - relax) * Uf(c, i, j, nkf + 1);
        // change in ghost point value
        r3 = UfNew(c, i, j, nkf + 1) - Uf(c, i, j, nkf + 1);
        Uf(c, i, j, nkf + 1) = UfNew(c, i, j, nkf + 1);
        // rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
        rmax2.max(fabs(r3));
        //      } // end for c=1,2

        // work on component 3 of the ghost point value of Uf
        // ic = (i+1)/2;
        // jc = j/2;
        // All Uc terms
        // right hand side is mismatch in displacement
        // NOTE: Uc is not changed by this routine, so the Uc-dependence could
        // be pre-computed
        c = 3;
        b1 = UnextcInterp(3, i, j, 1) +
             nuc * a_ghcof[0] * i16 *
                 (-Uc(3, ic, jc - 1, 0) * Mlrc(ic, jc - 1, 1) +
                  9 * Uc(3, ic, jc, 0) * Mlrc(ic, jc, 1) +
                  9 * Uc(3, ic, jc + 1, 0) * Mlrc(ic, jc + 1, 1) -
                  Uc(3, ic, jc + 2, 0) * Mlrc(ic, jc + 2, 1)) -
             Unextf(3, i, j, nkf);

        //    a11 =
        //    nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf));
        //    // no str
        a11 = nuf * a_ghcof[0] * Mlrf(i, j, nkf);  // no str

        //    r3 = Uf(3,i,j,nkf+1); // save previous value for relaxation below
        // solve for the ghost point value Uf(3,i,j,nkf+1)
        //    Uf(3,i,j,nkf+1) = b1/a11;
        UfNew(3, i, j, nkf + 1) =
            relax * b1 / a11 + (1 - relax) * Uf(3, i, j, nkf + 1);
        r3 = UfNew(3, i, j, nkf + 1) - Uf(3, i, j, nkf + 1);
        Uf(c, i, j, nkf + 1) = UfNew(c, i, j, nkf + 1);
        // rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);
        rmax3.max(fabs(r3));
      }  // end for i odd, j even
  );

  // update Uf
  // #pragma omp parallel
  //   for( sw4_type c=1 ; c <= 3 ; c++ )
  // #pragma omp for
  //     for( sw4_type j=jfb ; j <= jfe ; j+=2 )
  // #pragma omp simd
  //       for( sw4_type i=ifb ; i <= ife ; i+=2 )
  //       {
  //   Uf(c,i,j,nkf+1) = UfNew(c,i,j,nkf+1);
  // }

  // RAJA::RangeSegment c_range(1,4);

  // RAJA::kernel<ODDIEVENJ_EXEC_POL2>(
  // 		       RAJA::make_tuple(c_range,j_srange,i_srange),
  // 		       [=]RAJA_DEVICE (sw4_type c,sw4_type j,sw4_type i) {
  // 			 Uf(c,i,j,nkf+1) = UfNew(c,i,j,nkf+1);
  // 		       });

  // SYNC_STREAM;
  // rmax[3] = std::max(rmax[3], static_cast<float_sw4>(rmax1.get()));
  // rmax[4] = std::max(rmax[4], static_cast<float_sw4>(rmax2.get()));
  // rmax[5] = std::max(rmax[5], static_cast<float_sw4>(rmax3.get()));

  // rmax[3]=rmax1;
  // rmax[4]=rmax2;
  // rmax[5]=rmax3;

#undef Unextf
#undef UnextcInterp
#undef Mufs
#undef Mlfs
#undef Morf
#undef Mlrf
#undef Morc
#undef Mlrc
#undef Mucs
#undef Mlcs
#undef Unextc
#undef Bc
#undef BfRestrict
#undef Uc
#undef Uf
#undef UfNew
}  // end oddIevenJinterpJacobiOpt

//--------------------- Reference implementation
void oddIevenJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf,
                     Sarray &Lambdaf, Sarray &Rhof, Sarray &Uc, Sarray &Muc,
                     Sarray &Lambdac, Sarray &Rhoc, Sarray &Morc, Sarray &Mlrc,
                     Sarray &Unextf, Sarray &Bf, Sarray &UnextcInterp,
                     Sarray &Bc, sw4_type a_iStart[], sw4_type a_jStart[],
                     sw4_type a_iStartInt[], sw4_type a_iEndInt[], sw4_type a_jStartInt[],
                     sw4_type a_jEndInt[], sw4_type gf, sw4_type gc, sw4_type nkf, float_sw4 a_Dt,
                     float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
                     float_sw4 *a_strf_x, float_sw4 *a_strf_y,
                     float_sw4 *a_strc_x, float_sw4 *a_strc_y,
                     float_sw4 a_sbop[], float_sw4 a_ghcof[]) {
  SW4_MARK_FUNCTION;
// tmp
//  printf("Inside oddIevenJinterp! ");

// stretching on the coarse side
#define strc_x(i) a_strc_x[(i - a_iStart[gc])]
#define strc_y(j) a_strc_y[(j - a_jStart[gc])]

// stretching on the fine side
#define strf_x(i) a_strf_x[(i - a_iStart[gf])]
#define strf_y(j) a_strf_y[(j - a_jStart[gf])]

  //  sw4_type icb = a_iStartInt[gc];
  sw4_type ifb = a_iStartInt[gf];
  if (ifb % 2 == 0) ifb++;  // make sure ifb is odd

  //  sw4_type ice = a_iEndInt[gc];
  sw4_type ife = a_iEndInt[gf];

  //  sw4_type jcb = a_jStartInt[gc];
  sw4_type jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++;  // make sure jfb is even

  //  sw4_type jce = a_jEndInt[gc];
  sw4_type jfe = a_jEndInt[gf];

  float_sw4 nuf =
      a_Dt * a_Dt / (cof * hf * hf);  // cof=12 for the predictor, cof=1 for the
                                      // corrector (argument to this routine)
  float_sw4 nuc = a_Dt * a_Dt / (cof * hc * hc);
  //  float_sw4 ihc = 1/hc, ihf=1/hf;

  const float_sw4 i16 = 1.0 / 16;
  //  const float_sw4 i256 = 1.0/256;
  //  const float_sw4 i1024 = 1.0/1024;

  // residuals
  float_sw4 rmax1 = 0, rmax2 = 0, rmax3 = 0;

#pragma omp parallel for reduction(max : rmax1, rmax2, rmax3)
  for (sw4_type j = jfb; j <= jfe; j += 2)
#pragma omp simd
    for (sw4_type i = ifb; i <= ife; i += 2) {
      sw4_type ic, jc;
      float_sw4 b1, a11, r3;
      // updated components 1,2 of the ghost point value of Uf
      for (sw4_type c = 1; c <= 2; c++) {
        ic = (i + 1) / 2;
        jc = j / 2;
        // All Unextc terms
        //         b1 =
        //         i16*(-Unextc(c,ic,jc-1,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic,jc+1,1))-Unextc(c,ic,jc+2,1));
        // All Uc terms
        b1 = UnextcInterp(c, i, j, 1) +
             nuc * a_ghcof[0] * i16 *
                 (-Uc(c, ic, jc - 1, 0) * Morc(ic, jc - 1, 1) +
                  9 * Uc(c, ic, jc, 0) * Morc(ic, jc, 1) +
                  9 * Uc(c, ic, jc + 1, 0) * Morc(ic, jc + 1, 1) -
                  Uc(c, ic, jc + 2, 0) * Morc(ic, jc + 2, 1));

        b1 = b1 - Unextf(c, i, j, nkf);
        a11 = nuf * a_ghcof[0] * Muf(i, j, nkf) / (Rhof(i, j, nkf));
        //                  a11 =
        //                  nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j));
        r3 = Uf(c, i, j, nkf + 1);  // save old value for relaxation
                                    // update ghost point value Uf(c,i,j,nkf+1)
        Uf(c, i, j, nkf + 1) = b1 / a11;
        Uf(c, i, j, nkf + 1) = relax * Uf(c, i, j, nkf + 1) + (1 - relax) * r3;
        // change in ghost point value
        r3 = r3 - Uf(c, i, j, nkf + 1);
        if (c == 1)
          rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
        else
          rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
        //		  rmax[c-1+3] = rmax[c-1+3] > fabs(r3) ? rmax[c-1+3] :
        // fabs(r3);
      }  // end for c=1,2

      // work on componet 3 of the ghost point value of Uf
      {
        ic = (i + 1) / 2;
        jc = j / 2;
        // All Unextc terms
        //      b1 =
        //      i16*(-Unextc(3,ic,jc-1,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic,jc+1,1))-Unextc(3,ic,jc+2,1));
        // All Uc terms
        b1 = UnextcInterp(3, i, j, 1) +
             nuc * a_ghcof[0] * i16 *
                 (-Uc(3, ic, jc - 1, 0) * Mlrc(ic, jc - 1, 1) +
                  9 * Uc(3, ic, jc, 0) * Mlrc(ic, jc, 1) +
                  9 * Uc(3, ic, jc + 1, 0) * Mlrc(ic, jc + 1, 1) -
                  Uc(3, ic, jc + 2, 0) * Mlrc(ic, jc + 2, 1));
      }
      // right hand side is mismatch in displacement
      b1 = b1 - Unextf(3, i, j, nkf);
      a11 = nuf * a_ghcof[0] * (2 * Muf(i, j, nkf) + Lambdaf(i, j, nkf)) /
            (Rhof(i, j, nkf));  // no str
                                //               a11 =
      //               nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j));
      //               // with str

      r3 = Uf(3, i, j, nkf + 1);  // save previous value for relaxation below
      // solve for the ghost point value Uf(3,i,j,nkf+1)
      Uf(3, i, j, nkf + 1) = b1 / a11;
      Uf(3, i, j, nkf + 1) = relax * Uf(3, i, j, nkf + 1) + (1 - relax) * r3;
      r3 = r3 - Uf(3, i, j, nkf + 1);
      rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);
      //	       rmax[2+3] = rmax[2+3] > fabs(r3) ? rmax[2+3] : fabs(r3);

    }  // end for i odd, j even
  rmax[3] = rmax1;
  rmax[4] = rmax2;
  rmax[5] = rmax3;
#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
}  // end oddIevenJinterp
