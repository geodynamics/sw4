#include "sw4.h"
#include "Sarray.h"

void oddIoddJinterp(float_sw4 rmax[3], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		    Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		    Sarray &Mufs, Sarray &Mlfs,
		    Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		    int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		    int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		    float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		    float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-a_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-a_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-a_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-a_jStart[gf])]   

  int icb = a_iStartInt[gc];
  int ifb = a_iStartInt[gf];

  int ice = a_iEndInt[gc];
  int ife = a_iEndInt[gf];
   
  int jcb = a_jStartInt[gc];
  int jfb = a_jStartInt[gf];

  int jce = a_jEndInt[gc];
  int jfe = a_jEndInt[gf];

  float_sw4 nuf = a_Dt*a_Dt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
  float_sw4 nuc = a_Dt*a_Dt/(cof*hc*hc);
  float_sw4 ihc = 1/hc, ihf=1/hf;

  const float_sw4 i16 = 1.0/16;
  const float_sw4 i256 = 1.0/256;
  const float_sw4 i1024 = 1.0/1024;

// residuals
  float_sw4 rmax1=0, rmax2=0, rmax3=0;

#pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
  for( int jc= jcb ; jc <= jce ; jc++ )
    for( int ic= icb ; ic <= ice ; ic++ )
    {
      float_sw4 a11, a12, a21, a22, b1, b2, r1, r2, r3, deti;
// i odd, j odd
      int i=2*ic-1, j=2*jc-1;
      // setup 2x2 system matrix
      // unknowns: (Uf, Uc)
// eqn 1: continuity of normal stress: NEED stretching
      a11 = 0.25*Mufs(i,j,nkf)*a_sbop[0]*ihf; // ihf = 1/h on the fine grid; Mufs contains stretching
      a12 = Muc(ic,jc,1)*a_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));  // ihc = 1/h on the coarse grid
// eqn 2: continuity of displacement
// nuf = dt^2/(cof * hf^2), nuc = dt^2/(cof * hc^2)
      a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*a_ghcof[0]; 
      a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*a_ghcof[0];
      for( int c=1 ; c <= 2 ; c++ ) //  the 2 tangential components ? 
      {
// apply the restriction operator to the normal stress on the interface (Bf is on the fine grid)
// scale Bf by 1/strf ?
	b1  = i1024*( 
	  Bf(c,i-3,j-3,nkf)-9*Bf(c,i-3,j-1,nkf)-16*Bf(c,i-3,j,nkf)-9*Bf(c,i-3,j+1,nkf)+Bf(c,i-3,j+3,nkf)
	  +9*(-Bf(c,i-1,j-3,nkf)+9*Bf(c,i-1,j-1,nkf)+16*Bf(c,i-1,j,nkf)+9*Bf(c,i-1,j+1,nkf)-Bf(c,i-1,j+3,nkf))
	  +16*(-Bf(c,i,  j-3,nkf)+9*Bf(c,i,  j-1,nkf)+16*Bf(c,i,  j,nkf)+9*Bf(c,i,  j+1,nkf)-Bf(c,i,  j+3,nkf)) // with Bf(i,j)
	  +9*(-Bf(c,i+1,j-3,nkf)+9*Bf(c,i+1,j-1,nkf)+16*Bf(c,i+1,j,nkf)+9*Bf(c,i+1,j+1,nkf)-Bf(c,i+1,j+3,nkf)) +
	  Bf(c,i+3,j-3,nkf)-9*Bf(c,i+3,j-1,nkf)-16*Bf(c,i+3,j,nkf)-9*Bf(c,i+3,j+1,nkf)+Bf(c,i+3,j+3,nkf)
	  );

// scale Muf by 1/strf ?
	b1 = b1 - i1024*a_sbop[0]*ihf*(
	  Mufs(i-3,j-3,nkf)*Uf(c,i-3,j-3,nkf+1) - 9*Mufs(i-3,j-1,nkf)*Uf(c,i-3,j-1,nkf+1)
	  -16*Mufs(i-3,j,  nkf)*Uf(c,i-3,  j,nkf+1) - 9*Mufs(i-3,j+1,nkf)*Uf(c,i-3,j+1,nkf+1)
	  +Mufs(i-3,j+3,nkf)*Uf(c,i-3,j+3,nkf+1) +					    
	  9*(  -Mufs(i-1,j-3,nkf)*Uf(c,i-1,j-3,nkf+1) + 9*Mufs(i-1,j-1,nkf)*Uf(c,i-1,j-1,nkf+1) 
	       +16*Mufs(i-1,j,  nkf)*Uf(c,i-1,j,  nkf+1) + 9*Mufs(i-1,j+1,nkf)*Uf(c,i-1,j+1,nkf+1) 
	       -Mufs(i-1,j+3,nkf)*Uf(c,i-1,j+3,nkf+1) ) +
	  16*(  -Mufs(i,  j-3,nkf)*Uf(c,i,  j-3,nkf+1) + 9*Mufs(i,  j-1,nkf)*Uf(c,i,  j-1,nkf+1) // NOTE: the Uf(i,j) term is in a11
		+ 9*Mufs(i,  j+1,nkf)*Uf(c,i,  j+1,nkf+1)
		-Mufs(i,  j+3,nkf)*Uf(c,i,  j+3,nkf+1) ) + 
	  9*(  -Mufs(i+1,j-3,nkf)*Uf(c,i+1,j-3,nkf+1) + 9*Mufs(i+1,j-1,nkf)*Uf(c,i+1,j-1,nkf+1)
	       +16*Mufs(i+1,j,  nkf)*Uf(c,i+1,j,  nkf+1) + 9*Mufs(i+1,j+1,nkf)*Uf(c,i+1,j+1,nkf+1)
	       -Mufs(i+1,j+3,nkf)*Uf(c,i+1,j+3,nkf+1) ) +
	  Mufs(i+3,j-3,nkf)*Uf(c,i+3,j-3,nkf+1) - 9*Mufs(i+3,j-1,nkf)*Uf(c,i+3,j-1,nkf+1)
	  -16*Mufs(i+3,j,  nkf)*Uf(c,i+3,j,  nkf+1) - 9*Mufs(i+3,j+1,nkf)*Uf(c,i+3,j+1,nkf+1)
	  +Mufs(i+3,j+3,nkf)*Uf(c,i+3,j+3,nkf+1) );

	// NEED stretching term in b1; scale Bc ?
	b1 = b1 - Bc(c,ic,jc,1);

	b2 = Unextc(c,ic,jc,1)-Unextf(c,i,j,nkf); 

	deti=1/(a11*a22-a12*a21);
	r1 = Uf(c,i,j,nkf+1);
	r2 = Uc(c,ic,jc,0);
// solve the linear 2x2 system
	Uf(c,i,j,nkf+1) = deti*( a22*b1-a12*b2);
	Uc(c,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// damp the update of the ghost point values (r1, r2) hold previous values
	Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1) + (1-relax)*r1;
	Uc(c,ic,jc,0)   = relax*Uc(c,ic,jc,0)   + (1-relax)*r2;
// change in solution
	r1 = r1-Uf(c,i,j,nkf+1);
	r2 = r2-Uc(c,ic,jc,0);
	//               rmax[c-1] = rmax[c-1] > fabs(r1) ? rmax[c-1] : fabs(r1);
	//	       rmax[c-1] = rmax[c-1] > fabs(r2) ? rmax[c-1] : fabs(r2);

	if( c==1 )
	{
	  rmax1 = rmax1 > fabs(r1) ? rmax1 : fabs(r1);
	  rmax1 = rmax1 > fabs(r2) ? rmax1 : fabs(r2);
	}
	else
	{
	  rmax2 = rmax2 > fabs(r1) ? rmax2 : fabs(r1);
	  rmax2 = rmax2 > fabs(r2) ? rmax2 : fabs(r2);
	}
	//	       if( c == 2 && ic == 12 && jc == 13 )
	//	       {
	//	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	//		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	//		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	//	       }
      } // end for c=1,2
            
// setup the matrix for the 3rd component of the normal stress (different coefficients)
      // NEED stretching terms in a11 & a12
      a11 = 0.25*Mlfs(i,j,nkf)*a_sbop[0]*ihf; // Mlfs contains stretching
      a12 = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*a_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));

      a21 = nuf/Rhof(i,j,nkf)*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))*a_ghcof[0];
      a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*a_ghcof[0];

// apply the restriction operator to the fine grid normal stress grid function (Bf)
// scale Bf for stretching?
      b1  = i1024*( 
	Bf(3,i-3,j-3,nkf)-9*Bf(3,i-3,j-1,nkf)-16*Bf(3,i-3,j,nkf)-9*Bf(3,i-3,j+1,nkf)+Bf(3,i-3,j+3,nkf) +
	9*(-Bf(3,i-1,j-3,nkf)+9*Bf(3,i-1,j-1,nkf)+16*Bf(3,i-1,j,nkf)+9*Bf(3,i-1,j+1,nkf)-Bf(3,i-1,j+3,nkf)) +
	16*(-Bf(3,i,  j-3,nkf)+9*Bf(3,i,  j-1,nkf)+16*Bf(3,i,  j,nkf)+9*Bf(3,i,  j+1,nkf)-Bf(3,i,  j+3,nkf)) + 
	9*(-Bf(3,i+1,j-3,nkf)+9*Bf(3,i+1,j-1,nkf)+16*Bf(3,i+1,j,nkf)+9*Bf(3,i+1,j+1,nkf)-Bf(3,i+1,j+3,nkf)) +
	Bf(3,i+3,j-3,nkf)-9*Bf(3,i+3,j-1,nkf)-16*Bf(3,i+3,j,nkf)-9*Bf(3,i+3,j+1,nkf)+Bf(3,i+3,j+3,nkf) );

      b1 = b1 - i1024*a_sbop[0]*ihf*(
	Mlfs(i-3,j-3,nkf)*Uf(3,i-3,j-3,nkf+1) - 9*Mlfs(i-3,j-1,nkf)*Uf(3,i-3,j-1,nkf+1)
	-16*Mlfs(i-3,j,  nkf)*Uf(3,i-3,  j,nkf+1) - 9*Mlfs(i-3,j+1,nkf)*Uf(3,i-3,j+1,nkf+1)
	+Mlfs(i-3,j+3,nkf)*Uf(3,i-3,j+3,nkf+1) +					    
	9*(  -Mlfs(i-1,j-3,nkf)*Uf(3,i-1,j-3,nkf+1) + 9*Mlfs(i-1,j-1,nkf)*Uf(3,i-1,j-1,nkf+1) 
	     +16*Mlfs(i-1,j,  nkf)*Uf(3,i-1,j,  nkf+1) + 9*Mlfs(i-1,j+1,nkf)*Uf(3,i-1,j+1,nkf+1) 
	     -Mlfs(i-1,j+3,nkf)*Uf(3,i-1,j+3,nkf+1) ) +
	16*(  -Mlfs(i,  j-3,nkf)*Uf(3,i,  j-3,nkf+1) + 9*Mlfs(i,  j-1,nkf)*Uf(3,i,  j-1,nkf+1) 
	      + 9*Mlfs(i,  j+1,nkf)*Uf(3,i,  j+1,nkf+1)
	      -Mlfs(i,  j+3,nkf)*Uf(3,i,  j+3,nkf+1) ) + 
	9*(  -Mlfs(i+1,j-3,nkf)*Uf(3,i+1,j-3,nkf+1) + 9*Mlfs(i+1,j-1,nkf)*Uf(3,i+1,j-1,nkf+1)
	     +16*Mlfs(i+1,j,  nkf)*Uf(3,i+1,j,  nkf+1) + 9*Mlfs(i+1,j+1,nkf)*Uf(3,i+1,j+1,nkf+1)
	     -Mlfs(i+1,j+3,nkf)*Uf(3,i+1,j+3,nkf+1) ) +
	Mlfs(i+3,j-3,nkf)*Uf(3,i+3,j-3,nkf+1) - 9*Mlfs(i+3,j-1,nkf)*Uf(3,i+3,j-1,nkf+1)
	-16*Mlfs(i+3,j,  nkf)*Uf(3,i+3,j,  nkf+1) - 9*Mlfs(i+3,j+1,nkf)*Uf(3,i+3,j+1,nkf+1)
	+Mlfs(i+3,j+3,nkf)*Uf(3,i+3,j+3,nkf+1) );

// setup the RHS
      b1 = b1 - Bc(3,ic,jc,1); // need stretching terms in Bc
      b2 = Unextc(3,ic,jc,1)-Unextf(3,i,j,nkf); 
      deti=1/(a11*a22-a12*a21);
// previous values
      r1 = Uf(3,i,j,nkf+1);
      r2 = Uc(3,ic,jc,0);
// solve the 2x2 system for component 3 of Uf and Uc
      Uf(3,i,j,nkf+1) = deti*( a22*b1-a12*b2);
      Uc(3,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// relax the updated value
      Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1) + (1-relax)*r1;
      Uc(3,ic,jc,0)   = relax*Uc(3,ic,jc,0)   + (1-relax)*r2;
// change in ghost point values
      r1 = r1-Uf(3,i,j,nkf+1);
      r2 = r2-Uc(3,ic,jc,0);

      //	       rmax[2] = rmax[2] > fabs(r1) ? rmax[2] : fabs(r1);
      //	       rmax[2] = rmax[2] > fabs(r2) ? rmax[2] : fabs(r2);
      rmax3 = rmax3 > fabs(r1) ? rmax3 : fabs(r1);
      rmax3 = rmax3 > fabs(r2) ? rmax3 : fabs(r2);

      //	       int c=3;
      //	       if( c == 3 && ic == 12 && jc == 13 )
      //	       {
      //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
      //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
      //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
      //	       }
    } // end for ic, jc

  rmax[0] = rmax1;
  rmax[1] = rmax2;
  rmax[2] = rmax3;
#undef strc_x
#undef strc_y
#undef strf_x
#undef strf_y
} // end oddIoddJinterp


//--------------------------------- Optimized version
void oddIoddJinterpOpt(float_sw4 rmax[3], float_sw4* __restrict__ a_uf, float_sw4* __restrict__ a_muf, 
		       float_sw4* __restrict__ a_lambdaf, float_sw4* __restrict__ a_rhof, 
		       float_sw4* __restrict__ a_uc, float_sw4* __restrict__ a_muc, 
		       float_sw4* __restrict__ a_lambdac, float_sw4* __restrict__ a_rhoc,
		       float_sw4* __restrict__ a_mufs,   float_sw4* __restrict__ a_mlfs,
		       float_sw4* __restrict__ a_unextf, float_sw4* __restrict__ a_bf, 
		       float_sw4* __restrict__ a_unextc, float_sw4* __restrict__ a_bc,
		       int a_iStart[], int a_iEnd[], int a_jStart[], int a_jEnd[], int a_kStart[], int a_kEnd[], 
		       int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		       int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		       float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		       float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
// stretching on the coarse side
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
  
#define strc_x(i) a_strc_x[(i-iStartC)]   
#define strc_y(j) a_strc_y[(j-jStartC)]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-iStartF)]   
#define strf_y(j) a_strf_y[(j-jStartF)]   

// Bf indexing
  const int niF    = iEndF-iStartF+1;
  const int nijF   = niF*(jEndF-jStartF+1);
  const int nijk_bf = nijF*(1);
  const int base3_bf = (iStartF+niF*jStartF+nijF*nkf+nijk_bf); // only one k=nkf
#define Bf(c,i,j,k)     a_bf[-base3_bf+i+niF*(j)+nijF*(k)+nijk_bf*(c)]   
#define Unextf(c,i,j,k) a_unextf[-base3_bf+i+niF*(j)+nijF*(k)+nijk_bf*(c)] // same size as Bf

  const int base_mufs = (iStartF+niF*jStartF+nijF*nkf); // only one k=nkf
#define Mufs(i,j,k)     a_mufs[-base_mufs+i+niF*(j)+nijF*(k)]
#define Mlfs(i,j,k)     a_mlfs[-base_mufs+i+niF*(j)+nijF*(k)] // same size as Mufs
 
  const int niC    = iEndC-iStartC+1;
  const int nijC   = niC*(jEndC-jStartC+1);
  const int nijk_unextc = nijC*(1);
  const int base3_unextc = (iStartC+niC*jStartC+nijC*1+nijk_unextc); // only one k=1
#define Unextc(c,i,j,k) a_unextc[-base3_unextc+i+niC*(j)+nijC*(k)+nijk_unextc*(c)]   
#define Bc(c,i,j,k)     a_bc[-base3_unextc+i+niC*(j)+nijC*(k)+nijk_unextc*(c)] // same size as Unextc  

  const int nijk_uc = nijC*(kEndC-kStartC+1);
  const int base3_uc = (iStartC+niC*jStartC+nijC*kStartC+nijk_uc*1); // c-index has base=1 
#define Uc(c,i,j,k)   a_uc[-base3_uc+i+niC*(j)+nijC*(k)+nijk_uc*(c)]   

  const int base_rhoc = (iStartC+niC*jStartC+nijC*kStartC);
#define Rhoc(i,j,k)   a_rhoc[-base_rhoc+i+niC*(j)+nijC*(k)]   
#define Muc(i,j,k)    a_muc[-base_rhoc+i+niC*(j)+nijC*(k)]   // same size as Rhoc
#define Lambdac(i,j,k) a_lambdac[-base_rhoc+i+niC*(j)+nijC*(k)] // same size as Rhoc

  const int nijk_uf = nijF*(kEndF-kStartF+1);
  const int base3_uf = (iStartF+niF*jStartF+nijF*kStartF+nijk_uf*1); // c-index has base=1 
#define Uf(c,i,j,k)   a_uf[-base3_uf+i+niF*(j)+nijF*(k)+nijk_uf*(c)]   

  const int base_rhof = (iStartF+niF*jStartF+nijF*kStartF);
#define Rhof(i,j,k)    a_rhof[-base_rhof+i+niF*(j)+nijF*(k)]   
#define Muf(i,j,k)     a_muf[-base_rhof+i+niF*(j)+nijF*(k)]   // same size as Rhof
#define Lambdaf(i,j,k) a_lambdaf[-base_rhof+i+niF*(j)+nijF*(k)]  // same size as Rhof

//
  const int icb = a_iStartInt[gc];
  const int ifb = a_iStartInt[gf];

  const int ice = a_iEndInt[gc];
  const int ife = a_iEndInt[gf];
   
  const int jcb = a_jStartInt[gc];
  const int jfb = a_jStartInt[gf];

  const int jce = a_jEndInt[gc];
  const int jfe = a_jEndInt[gf];


  const float_sw4 nuf = a_Dt*a_Dt/(cof*hf*hf); // cof=12 for the predictor, cof=1 for the corrector (argument to this routine)
  const float_sw4 nuc = a_Dt*a_Dt/(cof*hc*hc);
  const float_sw4 ihc = 1/hc, ihf=1/hf;

  const float_sw4 i16 = 1.0/16;
  const float_sw4 i256 = 1.0/256;
  const float_sw4 i1024 = 1.0/1024;

// residuals
  float_sw4 rmax1=0, rmax2=0, rmax3=0;

#pragma omp parallel for reduction(max:rmax1,rmax2,rmax3)
  for( int jc= jcb ; jc <= jce ; jc++ )
    for( int ic= icb ; ic <= ice ; ic++ )
    {
      float_sw4 a11, a12, a21, a22, b1, b2, r1, r2, r3, deti;
// i odd, j odd
      int i=2*ic-1, j=2*jc-1;
      // setup 2x2 system matrix
      // unknowns: (Uf, Uc)
// eqn 1: continuity of normal stress: NEED stretching
      a11 = 0.25*Mufs(i,j,nkf)*a_sbop[0]*ihf; // ihf = 1/h on the fine grid; Mufs contains stretching
      a12 = Muc(ic,jc,1)*a_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));  // ihc = 1/h on the coarse grid
// eqn 2: continuity of displacement
// nuf = dt^2/(cof * hf^2), nuc = dt^2/(cof * hc^2)
      a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*a_ghcof[0]; 
      a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*a_ghcof[0];
      for( int c=1 ; c <= 2 ; c++ ) //  the 2 tangential components ? 
      {
// apply the restriction operator to the normal stress on the interface (Bf is on the fine grid)
// scale Bf by 1/strf ?
	b1  = i1024*( 
	  Bf(c,i-3,j-3,nkf)-9*Bf(c,i-3,j-1,nkf)-16*Bf(c,i-3,j,nkf)-9*Bf(c,i-3,j+1,nkf)+Bf(c,i-3,j+3,nkf)
	  +9*(-Bf(c,i-1,j-3,nkf)+9*Bf(c,i-1,j-1,nkf)+16*Bf(c,i-1,j,nkf)+9*Bf(c,i-1,j+1,nkf)-Bf(c,i-1,j+3,nkf))
	  +16*(-Bf(c,i,  j-3,nkf)+9*Bf(c,i,  j-1,nkf)+16*Bf(c,i,  j,nkf)+9*Bf(c,i,  j+1,nkf)-Bf(c,i,  j+3,nkf)) // with Bf(i,j)
	  +9*(-Bf(c,i+1,j-3,nkf)+9*Bf(c,i+1,j-1,nkf)+16*Bf(c,i+1,j,nkf)+9*Bf(c,i+1,j+1,nkf)-Bf(c,i+1,j+3,nkf)) +
	  Bf(c,i+3,j-3,nkf)-9*Bf(c,i+3,j-1,nkf)-16*Bf(c,i+3,j,nkf)-9*Bf(c,i+3,j+1,nkf)+Bf(c,i+3,j+3,nkf)
	  );

// scale Muf by 1/strf ?
	b1 = b1 - i1024*a_sbop[0]*ihf*(
	  Mufs(i-3,j-3,nkf)*Uf(c,i-3,j-3,nkf+1) - 9*Mufs(i-3,j-1,nkf)*Uf(c,i-3,j-1,nkf+1)
	  -16*Mufs(i-3,j,  nkf)*Uf(c,i-3,  j,nkf+1) - 9*Mufs(i-3,j+1,nkf)*Uf(c,i-3,j+1,nkf+1)
	  +Mufs(i-3,j+3,nkf)*Uf(c,i-3,j+3,nkf+1) +					    
	  9*(  -Mufs(i-1,j-3,nkf)*Uf(c,i-1,j-3,nkf+1) + 9*Mufs(i-1,j-1,nkf)*Uf(c,i-1,j-1,nkf+1) 
	       +16*Mufs(i-1,j,  nkf)*Uf(c,i-1,j,  nkf+1) + 9*Mufs(i-1,j+1,nkf)*Uf(c,i-1,j+1,nkf+1) 
	       -Mufs(i-1,j+3,nkf)*Uf(c,i-1,j+3,nkf+1) ) +
	  16*(  -Mufs(i,  j-3,nkf)*Uf(c,i,  j-3,nkf+1) + 9*Mufs(i,  j-1,nkf)*Uf(c,i,  j-1,nkf+1) // NOTE: the Uf(i,j) term is in a11
		+ 9*Mufs(i,  j+1,nkf)*Uf(c,i,  j+1,nkf+1)
		-Mufs(i,  j+3,nkf)*Uf(c,i,  j+3,nkf+1) ) + 
	  9*(  -Mufs(i+1,j-3,nkf)*Uf(c,i+1,j-3,nkf+1) + 9*Mufs(i+1,j-1,nkf)*Uf(c,i+1,j-1,nkf+1)
	       +16*Mufs(i+1,j,  nkf)*Uf(c,i+1,j,  nkf+1) + 9*Mufs(i+1,j+1,nkf)*Uf(c,i+1,j+1,nkf+1)
	       -Mufs(i+1,j+3,nkf)*Uf(c,i+1,j+3,nkf+1) ) +
	  Mufs(i+3,j-3,nkf)*Uf(c,i+3,j-3,nkf+1) - 9*Mufs(i+3,j-1,nkf)*Uf(c,i+3,j-1,nkf+1)
	  -16*Mufs(i+3,j,  nkf)*Uf(c,i+3,j,  nkf+1) - 9*Mufs(i+3,j+1,nkf)*Uf(c,i+3,j+1,nkf+1)
	  +Mufs(i+3,j+3,nkf)*Uf(c,i+3,j+3,nkf+1) );

	// NEED stretching term in b1; scale Bc ?
	b1 = b1 - Bc(c,ic,jc,1);

	b2 = Unextc(c,ic,jc,1)-Unextf(c,i,j,nkf); 

	deti=1/(a11*a22-a12*a21);
	r1 = Uf(c,i,j,nkf+1);
	r2 = Uc(c,ic,jc,0);
// solve the linear 2x2 system
	Uf(c,i,j,nkf+1) = deti*( a22*b1-a12*b2);
	Uc(c,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// damp the update of the ghost point values (r1, r2) hold previous values
	Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1) + (1-relax)*r1;
	Uc(c,ic,jc,0)   = relax*Uc(c,ic,jc,0)   + (1-relax)*r2;
// change in solution
	r1 = r1-Uf(c,i,j,nkf+1);
	r2 = r2-Uc(c,ic,jc,0);
	//               rmax[c-1] = rmax[c-1] > fabs(r1) ? rmax[c-1] : fabs(r1);
	//	       rmax[c-1] = rmax[c-1] > fabs(r2) ? rmax[c-1] : fabs(r2);

	if( c==1 )
	{
	  rmax1 = rmax1 > fabs(r1) ? rmax1 : fabs(r1);
	  rmax1 = rmax1 > fabs(r2) ? rmax1 : fabs(r2);
	}
	else
	{
	  rmax2 = rmax2 > fabs(r1) ? rmax2 : fabs(r1);
	  rmax2 = rmax2 > fabs(r2) ? rmax2 : fabs(r2);
	}
	//	       if( c == 2 && ic == 12 && jc == 13 )
	//	       {
	//	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	//		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	//		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	//	       }
      } // end for c=1,2
            
// setup the matrix for the 3rd component of the normal stress (different coefficients)
      // NEED stretching terms in a11 & a12
      a11 = 0.25*Mlfs(i,j,nkf)*a_sbop[0]*ihf; // Mlfs contains stretching
      a12 = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*a_sbop[0]*ihc/(strc_x(ic)*strc_y(jc));

      a21 = nuf/Rhof(i,j,nkf)*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))*a_ghcof[0];
      a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*a_ghcof[0];

// apply the restriction operator to the fine grid normal stress grid function (Bf)
// scale Bf for stretching?
      b1  = i1024*( 
	Bf(3,i-3,j-3,nkf)-9*Bf(3,i-3,j-1,nkf)-16*Bf(3,i-3,j,nkf)-9*Bf(3,i-3,j+1,nkf)+Bf(3,i-3,j+3,nkf) +
	9*(-Bf(3,i-1,j-3,nkf)+9*Bf(3,i-1,j-1,nkf)+16*Bf(3,i-1,j,nkf)+9*Bf(3,i-1,j+1,nkf)-Bf(3,i-1,j+3,nkf)) +
	16*(-Bf(3,i,  j-3,nkf)+9*Bf(3,i,  j-1,nkf)+16*Bf(3,i,  j,nkf)+9*Bf(3,i,  j+1,nkf)-Bf(3,i,  j+3,nkf)) + 
	9*(-Bf(3,i+1,j-3,nkf)+9*Bf(3,i+1,j-1,nkf)+16*Bf(3,i+1,j,nkf)+9*Bf(3,i+1,j+1,nkf)-Bf(3,i+1,j+3,nkf)) +
	Bf(3,i+3,j-3,nkf)-9*Bf(3,i+3,j-1,nkf)-16*Bf(3,i+3,j,nkf)-9*Bf(3,i+3,j+1,nkf)+Bf(3,i+3,j+3,nkf) );

      b1 = b1 - i1024*a_sbop[0]*ihf*(
	Mlfs(i-3,j-3,nkf)*Uf(3,i-3,j-3,nkf+1) - 9*Mlfs(i-3,j-1,nkf)*Uf(3,i-3,j-1,nkf+1)
	-16*Mlfs(i-3,j,  nkf)*Uf(3,i-3,  j,nkf+1) - 9*Mlfs(i-3,j+1,nkf)*Uf(3,i-3,j+1,nkf+1)
	+Mlfs(i-3,j+3,nkf)*Uf(3,i-3,j+3,nkf+1) +					    
	9*(  -Mlfs(i-1,j-3,nkf)*Uf(3,i-1,j-3,nkf+1) + 9*Mlfs(i-1,j-1,nkf)*Uf(3,i-1,j-1,nkf+1) 
	     +16*Mlfs(i-1,j,  nkf)*Uf(3,i-1,j,  nkf+1) + 9*Mlfs(i-1,j+1,nkf)*Uf(3,i-1,j+1,nkf+1) 
	     -Mlfs(i-1,j+3,nkf)*Uf(3,i-1,j+3,nkf+1) ) +
	16*(  -Mlfs(i,  j-3,nkf)*Uf(3,i,  j-3,nkf+1) + 9*Mlfs(i,  j-1,nkf)*Uf(3,i,  j-1,nkf+1) 
	      + 9*Mlfs(i,  j+1,nkf)*Uf(3,i,  j+1,nkf+1)
	      -Mlfs(i,  j+3,nkf)*Uf(3,i,  j+3,nkf+1) ) + 
	9*(  -Mlfs(i+1,j-3,nkf)*Uf(3,i+1,j-3,nkf+1) + 9*Mlfs(i+1,j-1,nkf)*Uf(3,i+1,j-1,nkf+1)
	     +16*Mlfs(i+1,j,  nkf)*Uf(3,i+1,j,  nkf+1) + 9*Mlfs(i+1,j+1,nkf)*Uf(3,i+1,j+1,nkf+1)
	     -Mlfs(i+1,j+3,nkf)*Uf(3,i+1,j+3,nkf+1) ) +
	Mlfs(i+3,j-3,nkf)*Uf(3,i+3,j-3,nkf+1) - 9*Mlfs(i+3,j-1,nkf)*Uf(3,i+3,j-1,nkf+1)
	-16*Mlfs(i+3,j,  nkf)*Uf(3,i+3,j,  nkf+1) - 9*Mlfs(i+3,j+1,nkf)*Uf(3,i+3,j+1,nkf+1)
	+Mlfs(i+3,j+3,nkf)*Uf(3,i+3,j+3,nkf+1) );

// setup the RHS
      b1 = b1 - Bc(3,ic,jc,1); // need stretching terms in Bc
      b2 = Unextc(3,ic,jc,1)-Unextf(3,i,j,nkf); 
      deti=1/(a11*a22-a12*a21);
// previous values
      r1 = Uf(3,i,j,nkf+1);
      r2 = Uc(3,ic,jc,0);
// solve the 2x2 system for component 3 of Uf and Uc
      Uf(3,i,j,nkf+1) = deti*( a22*b1-a12*b2);
      Uc(3,ic,jc,0)   = deti*(-a21*b1+a11*b2);
// relax the updated value
      Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1) + (1-relax)*r1;
      Uc(3,ic,jc,0)   = relax*Uc(3,ic,jc,0)   + (1-relax)*r2;
// change in ghost point values
      r1 = r1-Uf(3,i,j,nkf+1);
      r2 = r2-Uc(3,ic,jc,0);

      //	       rmax[2] = rmax[2] > fabs(r1) ? rmax[2] : fabs(r1);
      //	       rmax[2] = rmax[2] > fabs(r2) ? rmax[2] : fabs(r2);
      rmax3 = rmax3 > fabs(r1) ? rmax3 : fabs(r1);
      rmax3 = rmax3 > fabs(r2) ? rmax3 : fabs(r2);

      //	       int c=3;
      //	       if( c == 3 && ic == 12 && jc == 13 )
      //	       {
      //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
      //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
      //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
      //	       }
    } // end for ic, jc

  rmax[0] = rmax1;
  rmax[1] = rmax2;
  rmax[2] = rmax3;
} // end oddIoddJinterp
