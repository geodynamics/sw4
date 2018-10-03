#include "sw4.h"
#include "Sarray.h"
#include <cstdio>
#include "caliper.h"
//--------------------- Jacobi ------------------------
void evenIevenJinterpJacobi(float_sw4 rmax[6], Sarray &Uf, Sarray &UfNew, Sarray &Uc, 
			    Sarray &Morc, Sarray &Mlrc,
			    Sarray &Morf, Sarray &Mlrf,
			    Sarray &Unextf, Sarray &UnextcInterp,
			    int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
			    int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
			    float_sw4 cof, float_sw4 relax,
			    float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
  SW4_MARK_FUNCTION;
// tmp
//  printf("Inside evenIevenJinterp! ");
  
//  int icb = a_iStartInt[gc];
  int ifb = a_iStartInt[gf];
  if (ifb % 2 == 1) ifb++; // make sure ifb is even

  //  int ice = a_iEndInt[gc];
  int ife = a_iEndInt[gf];
   
  // int jcb = a_jStartInt[gc];
  int jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++; // make sure jfb is even
  

  //  int jce = a_jEndInt[gc];
  int jfe = a_jEndInt[gf];

  float_sw4 nuf = a_Dt*a_Dt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
  float_sw4 nuc = a_Dt*a_Dt/(cof*hc*hc);
  //  float_sw4 ihc = 1/hc, ihf=1/hf;

  // const float_sw4 i16 = 1.0/16;
  const float_sw4 i256 = 1.0/256;
  //const float_sw4 i1024 = 1.0/1024;

// residuals
  float_sw4 rmax1=0, rmax2=0, rmax3=0;

#pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
  for( int j=jfb ; j <= jfe ; j+=2 )
#pragma omp simd
    for( int i=ifb ; i <= ife ; i+=2 )
    {
      int ic, jc;
      ic = i/2;
      jc = j/2;

      float_sw4 b1, a11, r3;

	// for( int c=1 ; c <= 2 ; c++ ) // update components 1,2 of the ghost point value of Uf
	// {
// un-roll the c-loop
      int c=1;
      
// All Uc terms
// NOTE: This routine does not change Uc
      b1 = UnextcInterp(c,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
	Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
	  +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
	  -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
	  -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
	  9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
	  -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
	  )
	+ Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
	-9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
	+Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
	) - Unextf(c,i,j,nkf);

//      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
      a11 = nuf*a_ghcof[0]*Morf(i,j,nkf);

// update ghost point value Uf(c,i,j,nkf+1)
//	  Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c,i,j,nkf+1) = relax* b1/a11 + (1-relax)*Uf(c,i,j,nkf+1);
// change in ghost point value
      r3 = UfNew(c,i,j,nkf+1) - Uf(c,i,j,nkf+1);
      rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
//	} // end for c=1,2
      
      c=2;
// All Uc terms
// NOTE: This routine does not change Uc
      b1 = UnextcInterp(c,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
	Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
	  +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
	  -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
	  -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
	  9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
	  -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
	  )
	+ Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
	-9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
	+Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
	) - Unextf(c,i,j,nkf);

//      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));

// update ghost point value Uf(c,i,j,nkf+1)
//	  Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c,i,j,nkf+1) = relax* b1/a11 + (1-relax)*Uf(c,i,j,nkf+1);
// change in ghost point value
      r3 = UfNew(c,i,j,nkf+1) - Uf(c,i,j,nkf+1);
      rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
               
// work on componet 3 of the ghost point value of Uf
// right hand side is mismatch in displacement                
// All Uc terms
      b1 =  UnextcInterp(3,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(3,ic-1,jc-1,0)*Mlrc(ic-1,jc-1,1)
	-9*(Uc(3,ic,  jc-1,0)*Mlrc(ic,  jc-1,1)+
	    Uc(3,ic+1,jc-1,0)*Mlrc(ic+1,jc-1,1))+
	Uc(3,ic+2,jc-1,0)*Mlrc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+
	  9*( Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+
	      Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1))
	  -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1)  
	  -Uc(3,ic-1,jc+1,0)*Mlrc(ic-1,jc+1,1)+
	  9*( Uc(3,ic,  jc+1,0)*Mlrc(ic,  jc+1,1)+
	      Uc(3,ic+1,jc+1,0)*Mlrc(ic+1,jc+1,1))
	  -Uc(3,ic+2,jc+1,0)*Mlrc(ic+2,jc+1,1)  )
	+ Uc(3,ic-1,jc+2,0)*Mlrc(ic-1,jc+2,1)
	-9*(Uc(3,ic,  jc+2,0)*Mlrc(ic,  jc+2,1) +
	    Uc(3,ic+1,jc+2,0)*Mlrc(ic+1,jc+2,1)) +
	Uc(3,ic+2,jc+2,0)*Mlrc(ic+2,jc+2,1) ) - Unextf(3,i,j,nkf); 

//      a11 = nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf)); // no str
      a11 = nuf*a_ghcof[0]*Mlrf(i,j,nkf);

// solve for the ghost point value Uf(3,i,j,nkf+1)
//	Uf(3,i,j,nkf+1) = b1/a11;
      UfNew(3,i,j,nkf+1) = relax*b1/a11 + (1-relax)*Uf(3,i,j,nkf+1);
      r3 = UfNew(3,i,j,nkf+1) - Uf(3,i,j,nkf+1);
      rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);
    } // end for i even, j even

// update Uf
#pragma omp parallel
  for( int c=1 ; c <= 3 ; c++ ) 
#pragma omp for
    for( int j=jfb ; j <= jfe ; j+=2 )
#pragma omp simd
      for( int i=ifb ; i <= ife ; i+=2 )
      {
	Uf(c,i,j,nkf+1) = UfNew(c,i,j,nkf+1);
      }

  rmax[3] = rmax1;
  rmax[4] = rmax2;
  rmax[5] = rmax3;
} // end evenIevenJinterpJacobi

// ---------- Optimized version -----------
void evenIevenJinterpJacobiOpt(float_sw4 rmax[6], float_sw4* __restrict__ a_uf, 
			       float_sw4* __restrict__ a_ufnew, float_sw4* __restrict__ a_uc,
			       float_sw4* __restrict__ a_morc, float_sw4* __restrict__ a_mlrc,
			       float_sw4* __restrict__ a_morf, float_sw4* __restrict__ a_mlrf,
			       float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_uncint,
			       int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
			       int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
			       int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc,
			       float_sw4 cof, float_sw4 relax,
			       float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
  SW4_MARK_FUNCTION;
  const int iStartC = a_iStart[gc];
  const int jStartC = a_jStart[gc];
  const int kStartC = a_kStart[gc];

  const int iEndC = a_iEnd[gc];
  const int jEndC = a_jEnd[gc];
  const int kEndC = a_kEnd[gc];

  const int iStartF = a_iStart[gf];
  const int jStartF = a_jStart[gf];
  const int kStartF = a_kStart[gf];
  const int iEndF = a_iEnd[gf];
  const int jEndF = a_jEnd[gf];
  const int kEndF = a_kEnd[gf];

// Bf indexing
  const int niF    = iEndF-iStartF+1;
  const int nijF   = niF*(jEndF-jStartF+1);
  const int nijk_bf = nijF*(1); // only one k-plane
  const int base3_bf = (iStartF+niF*jStartF+nijF*nkf+nijk_bf); // only one k=nkf
#define Unextf(c,i,j,k) a_unextf[-base3_bf+i+niF*(j)+nijF*(k)+nijk_bf*(c)] 

  const int nijk_uncint = nijF*(1); // only one k-plane
  const int base3_uncint = (iStartF+niF*jStartF+nijF*1+nijk_uncint); // only k=1
#define UnextcInterp(c,i,j,k) a_uncint[-base3_uncint+i+niF*(j)+nijF*(k)+nijk_uncint*(c)]

  const int base_mufs = (iStartF+niF*jStartF+nijF*nkf); // only one k=nkf
#define Morf(i,j,k)     a_morf[-base_mufs+i+niF*(j)+nijF*(k)] // same size as Mufs
#define Mlrf(i,j,k)     a_mlrf[-base_mufs+i+niF*(j)+nijF*(k)] // same size as Mufs
 
  const int niC    = iEndC-iStartC+1;
  const int nijC   = niC*(jEndC-jStartC+1);
  const int base_morc = (iStartC+niC*jStartC+nijC*1); // only one k=1
#define Morc(i,j,k)     a_morc[-base_morc+i+niC*(j)+nijC*(k)]
#define Mlrc(i,j,k)     a_mlrc[-base_morc+i+niC*(j)+nijC*(k)] // same size as Morc

  const int nijk_uc = nijC*(kEndC-kStartC+1);
  const int base3_uc = (iStartC+niC*jStartC+nijC*kStartC+nijk_uc*1); // c-index has base=1 
#define Uc(c,i,j,k)   a_uc[-base3_uc+i+niC*(j)+nijC*(k)+nijk_uc*(c)]   

  const int nijk_uf = nijF*(kEndF-kStartF+1);
  const int base3_uf = (iStartF+niF*jStartF+nijF*kStartF+nijk_uf*1); // c-index has base=1 
#define Uf(c,i,j,k)   a_uf[-base3_uf+i+niF*(j)+nijF*(k)+nijk_uf*(c)]   

  const int nijk_ufnew = nijF*(1); // only one k-plane
  const int base3_ufnew = (iStartF+niF*jStartF+nijF*(nkf+1)+nijk_ufnew*1); // only one k=nkf+1; c-index has base=1 
#define UfNew(c,i,j,k)   a_ufnew[-base3_ufnew+i+niF*(j)+nijF*(k)+nijk_ufnew*(c)]   

// previous stuff  
//  int icb = a_iStartInt[gc];
  int ifb = a_iStartInt[gf];
  if (ifb % 2 == 1) ifb++; // make sure ifb is even

  //  int ice = a_iEndInt[gc];
  int ife = a_iEndInt[gf];
   
  //  int jcb = a_jStartInt[gc];
  int jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++; // make sure jfb is even
  

  //  int jce = a_jEndInt[gc];
  int jfe = a_jEndInt[gf];

  float_sw4 nuf = a_Dt*a_Dt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
  float_sw4 nuc = a_Dt*a_Dt/(cof*hc*hc);
  // float_sw4 ihc = 1/hc, ihf=1/hf;

  //  const float_sw4 i16 = 1.0/16;
  const float_sw4 i256 = 1.0/256;
  //const float_sw4 i1024 = 1.0/1024;

// residuals
  //float_sw4 rmax1=0, rmax2=0, rmax3=0;
  RAJA::ReduceMax<REDUCTION_POLICY,float_sw4> rmax1(0),rmax2(0),rmax3(0);
  
  RAJA::TypedRangeStrideSegment<long> j_srange(jfb,jfe+1,2);
  RAJA::TypedRangeStrideSegment<long> i_srange(ifb,ife+1,2);
  
// #pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
//   for( int j=jfb ; j <= jfe ; j+=2 )
// #pragma omp simd
//     for( int i=ifb ; i <= ife ; i+=2 )
RAJA::kernel<EVENIEVENJ_EXEC_POL>(
				    RAJA::make_tuple(j_srange,i_srange),
				    [=]RAJA_DEVICE (int j,int i) 
    {
      int ic, jc;
      ic = i/2;
      jc = j/2;

      float_sw4 b1, a11, r3;

	// for( int c=1 ; c <= 2 ; c++ ) // update components 1,2 of the ghost point value of Uf
	// {
// un-roll the c-loop
      int c=1;
      
// All Uc terms
// NOTE: This routine does not change Uc
      b1 = UnextcInterp(c,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
	Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
	  +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
	  -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
	  -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
	  9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
	  -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
	  )
	+ Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
	-9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
	+Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
	) - Unextf(c,i,j,nkf);

//      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
      a11 = nuf*a_ghcof[0]*Morf(i,j,nkf);

// update ghost point value Uf(c,i,j,nkf+1)
//	  Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c,i,j,nkf+1) = relax* b1/a11 + (1-relax)*Uf(c,i,j,nkf+1);
// change in ghost point value
      r3 = UfNew(c,i,j,nkf+1) - Uf(c,i,j,nkf+1);
      //rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
      rmax1.max(fabs(r3));
//	} // end for c=1,2
      
      c=2;
// All Uc terms
// NOTE: This routine does not change Uc
      b1 = UnextcInterp(c,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
	Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
	  +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
	  -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
	  -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
	  9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
	  -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
	  )
	+ Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
	-9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
	+Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
	) - Unextf(c,i,j,nkf);

//      a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));

// update ghost point value Uf(c,i,j,nkf+1)
//	  Uf(c,i,j,nkf+1) = b1/a11;
      UfNew(c,i,j,nkf+1) = relax* b1/a11 + (1-relax)*Uf(c,i,j,nkf+1);
// change in ghost point value
      r3 = UfNew(c,i,j,nkf+1) - Uf(c,i,j,nkf+1);
      //rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
      rmax2.max(fabs(r3));
               
// work on componet 3 of the ghost point value of Uf
// right hand side is mismatch in displacement                
// All Uc terms
      b1 =  UnextcInterp(3,i,j,1) + nuc*a_ghcof[0]*i256*(
	Uc(3,ic-1,jc-1,0)*Mlrc(ic-1,jc-1,1)
	-9*(Uc(3,ic,  jc-1,0)*Mlrc(ic,  jc-1,1)+
	    Uc(3,ic+1,jc-1,0)*Mlrc(ic+1,jc-1,1))+
	Uc(3,ic+2,jc-1,0)*Mlrc(ic+2,jc-1,1)
	+ 9*(
	  -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+
	  9*( Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+
	      Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1))
	  -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1)  
	  -Uc(3,ic-1,jc+1,0)*Mlrc(ic-1,jc+1,1)+
	  9*( Uc(3,ic,  jc+1,0)*Mlrc(ic,  jc+1,1)+
	      Uc(3,ic+1,jc+1,0)*Mlrc(ic+1,jc+1,1))
	  -Uc(3,ic+2,jc+1,0)*Mlrc(ic+2,jc+1,1)  )
	+ Uc(3,ic-1,jc+2,0)*Mlrc(ic-1,jc+2,1)
	-9*(Uc(3,ic,  jc+2,0)*Mlrc(ic,  jc+2,1) +
	    Uc(3,ic+1,jc+2,0)*Mlrc(ic+1,jc+2,1)) +
	Uc(3,ic+2,jc+2,0)*Mlrc(ic+2,jc+2,1) ) - Unextf(3,i,j,nkf); 

//      a11 = nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf)); // no str
      a11 = nuf*a_ghcof[0]*Mlrf(i,j,nkf);

// solve for the ghost point value Uf(3,i,j,nkf+1)
//	Uf(3,i,j,nkf+1) = b1/a11;
      UfNew(3,i,j,nkf+1) = relax*b1/a11 + (1-relax)*Uf(3,i,j,nkf+1);
      r3 = UfNew(3,i,j,nkf+1) - Uf(3,i,j,nkf+1);
      //rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);
      rmax3.max(fabs(r3));
    } ); // end for i even, j even

// update Uf
// #pragma omp parallel
//   for( int c=1 ; c <= 3 ; c++ ) 
// #pragma omp for
//     for( int j=jfb ; j <= jfe ; j+=2 )
// #pragma omp simd
//       for( int i=ifb ; i <= ife ; i+=2 )
      RAJA::RangeSegment c_range(1,4);
      RAJA::kernel<XRHS_POL_ASYNC>(
      RAJA::make_tuple(c_range,j_srange,i_srange),
	[=]RAJA_DEVICE (int c,int j,int i) 
	  {
      Uf(c,i,j,nkf+1) = UfNew(c,i,j,nkf+1);
	  }); //SYNC_STREAM;
      
      rmax[3] = max(rmax[3],static_cast<float_sw4>(rmax1.get()));
      rmax[4] = max(rmax[4],static_cast<float_sw4>(rmax2.get()));
      rmax[5] = max(rmax[5],static_cast<float_sw4>(rmax3.get()));

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
} // end evenIevenJinterpJacobiOpt

//--------------------- Reference ------------------------
void evenIevenJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &UnextcInterp, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
  SW4_MARK_FUNCTION;
// tmp
//  printf("Inside evenIevenJinterp! ");
  
// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-a_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-a_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-a_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-a_jStart[gf])]   

  //  int icb = a_iStartInt[gc];
  int ifb = a_iStartInt[gf];
  if (ifb % 2 == 1) ifb++; // make sure ifb is even

  // int ice = a_iEndInt[gc];
  int ife = a_iEndInt[gf];
   
  // int jcb = a_jStartInt[gc];
  int jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++; // make sure jfb is even
  

  //  int jce = a_jEndInt[gc];
  int jfe = a_jEndInt[gf];

  float_sw4 nuf = a_Dt*a_Dt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
  float_sw4 nuc = a_Dt*a_Dt/(cof*hc*hc);
  //  float_sw4 ihc = 1/hc, ihf=1/hf;

  // const float_sw4 i16 = 1.0/16;
  const float_sw4 i256 = 1.0/256;
  // const float_sw4 i1024 = 1.0/1024;

// residuals
  float_sw4 rmax1=0, rmax2=0, rmax3=0;

#pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
  for( int j=jfb ; j <= jfe ; j+=2 )
#pragma omp simd
    for( int i=ifb ; i <= ife ; i+=2 )
    {
      int ic, jc;
      float_sw4 b1, a11, r3;
// updated components 1,2 of the ghost point value of Uf
	for( int c=1 ; c <= 2 ; c++ )
	{
	  ic = i/2;
	  jc = j/2;
// All Unextc terms
	  // b1 = i256*
	  //   ( Unextc(c,ic-1,jc-1,1)-9*(Unextc(c,ic,jc-1,1)+Unextc(c,ic+1,jc-1,1))+Unextc(c,ic+2,jc-1,1)
	  //     + 9*(-Unextc(c,ic-1,jc,  1)+9*(Unextc(c,ic,jc,  1)+Unextc(c,ic+1,jc,  1))-Unextc(c,ic+2,jc,  1)  
	  //          -Unextc(c,ic-1,jc+1,1)+9*(Unextc(c,ic,jc+1,1)+Unextc(c,ic+1,jc+1,1))-Unextc(c,ic+2,jc+1,1))
	  //     +Unextc(c,ic-1,jc+2,1)-9*(Unextc(c,ic,jc+2,1)+Unextc(c,ic+1,jc+2,1))+Unextc(c,ic+2,jc+2,1) );

// All Uc terms
	  b1 = UnextcInterp(c,i,j,1) + nuc*a_ghcof[0]*i256*(
	    Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1)) +
	    Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
	    + 9*(
	      -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)
	      +9*(Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
	      -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
	      -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
	      9*(Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
	      -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)
	      )
	    + Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
	    -9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1))
	    +Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)
	    );
	  b1 = b1 - Unextf(c,i,j,nkf); 
	  a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
//                  a11 = nuf*a_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j));
	  r3 = Uf(c,i,j,nkf+1); // save old value for relaxation
// update ghost point value Uf(c,i,j,nkf+1)
	  Uf(c,i,j,nkf+1) = b1/a11;
	  Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1)+(1-relax)*r3;
// change in ghost point value
	  r3 = r3 - Uf(c,i,j,nkf+1);
	  if( c == 1 )
	    rmax1 = rmax1 > fabs(r3) ? rmax1 : fabs(r3);
	  else
	    rmax2 = rmax2 > fabs(r3) ? rmax2 : fabs(r3);
	  //		  rmax[c-1+3] = rmax[c-1+3] > fabs(r3) ? rmax[c-1+3] : fabs(r3);
	} // end for c=1,2
               
// work on componet 3 of the ghost point value of Uf
	ic = i/2;
	jc = j/2;
// All Unextc terms
	// b1 = i256*
	//   ( Unextc(3,ic-1,jc-1,1)-9*(Unextc(3,ic,jc-1,1)+Unextc(3,ic+1,jc-1,1))+Unextc(3,ic+2,jc-1,1)
	//     + 9*(-Unextc(3,ic-1,jc,  1)+9*(Unextc(3,ic,jc,  1)+Unextc(3,ic+1,jc,  1))-Unextc(3,ic+2,jc,  1)  
	// 	 -Unextc(3,ic-1,jc+1,1)+9*(Unextc(3,ic,jc+1,1)+Unextc(3,ic+1,jc+1,1))-Unextc(3,ic+2,jc+1,1))
	//     +Unextc(3,ic-1,jc+2,1)-9*(Unextc(3,ic,jc+2,1)+Unextc(3,ic+1,jc+2,1))+Unextc(3,ic+2,jc+2,1) );

// All Uc terms
	b1 =  UnextcInterp(3,i,j,1) + nuc*a_ghcof[0]*i256*(
	  Uc(3,ic-1,jc-1,0)*Mlrc(ic-1,jc-1,1)
	  -9*(Uc(3,ic,  jc-1,0)*Mlrc(ic,  jc-1,1)+
	      Uc(3,ic+1,jc-1,0)*Mlrc(ic+1,jc-1,1))+
	  Uc(3,ic+2,jc-1,0)*Mlrc(ic+2,jc-1,1)
	  + 9*(
	    -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+
	    9*( Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+
		Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1))
	    -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1)  
	    -Uc(3,ic-1,jc+1,0)*Mlrc(ic-1,jc+1,1)+
	    9*( Uc(3,ic,  jc+1,0)*Mlrc(ic,  jc+1,1)+
		Uc(3,ic+1,jc+1,0)*Mlrc(ic+1,jc+1,1))
	    -Uc(3,ic+2,jc+1,0)*Mlrc(ic+2,jc+1,1)  )
	  + Uc(3,ic-1,jc+2,0)*Mlrc(ic-1,jc+2,1)
	  -9*(Uc(3,ic,  jc+2,0)*Mlrc(ic,  jc+2,1) +
	      Uc(3,ic+1,jc+2,0)*Mlrc(ic+1,jc+2,1)) +
	  Uc(3,ic+2,jc+2,0)*Mlrc(ic+2,jc+2,1)   );
// right hand side is mismatch in displacement                
	b1 = b1 - Unextf(3,i,j,nkf); 
	a11 = nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf)); // no str
//               a11 = nuf*a_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf))/(strf_x(i)*strf_y(j)); // with str

	r3 = Uf(3,i,j,nkf+1); // save previous value for relaxation below
// solve for the ghost point value Uf(3,i,j,nkf+1)
	Uf(3,i,j,nkf+1) = b1/a11;
	Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1)+(1-relax)*r3;
	r3 = r3 - Uf(3,i,j,nkf+1);
	rmax3 = rmax3 > fabs(r3) ? rmax3 : fabs(r3);
	//	       rmax[2+3] = rmax[2+3] > fabs(r3) ? rmax[2+3] : fabs(r3);

    } // end for i even, j even
  rmax[3] = rmax1;
  rmax[4] = rmax2;
  rmax[5] = rmax3;
#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
} // end evenIevenJinterp

