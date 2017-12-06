#include "sw4.h"
#include "Sarray.h"
#include <cstdio>

void evenIevenJinterp(float_sw4 rmax[6], Sarray &Uf, Sarray &Muf, Sarray &Lambdaf, Sarray &Rhof, 
		     Sarray &Uc, Sarray &Muc, Sarray &Lambdac, Sarray &Rhoc,
		     Sarray &Morc, Sarray &Mlrc,
		     Sarray &Unextf, Sarray &Bf, Sarray &Unextc, Sarray &Bc,
		     int a_iStart[], int a_jStart[], int a_iStartInt[], int a_iEndInt[], int a_jStartInt[], int a_jEndInt[],
		     int gf, int gc, int nkf, float_sw4 a_Dt, float_sw4 hf, float_sw4 hc, float_sw4 cof, float_sw4 relax,
		     float_sw4 *a_strf_x, float_sw4 *a_strf_y, float_sw4 *a_strc_x, float_sw4 *a_strc_y, 
		     float_sw4 a_sbop[], float_sw4 a_ghcof[])
{  
// tmp
//  printf("Inside evenIevenJinterp! ");
  
// stretching on the coarse side
#define strc_x(i) a_strc_x[(i-a_iStart[gc])]   
#define strc_y(j) a_strc_y[(j-a_jStart[gc])]   

// stretching on the fine side
#define strf_x(i) a_strf_x[(i-a_iStart[gf])]   
#define strf_y(j) a_strf_y[(j-a_jStart[gf])]   

  int icb = a_iStartInt[gc];
  int ifb = a_iStartInt[gf];
  if (ifb % 2 == 1) ifb++; // make sure ifb is even

  int ice = a_iEndInt[gc];
  int ife = a_iEndInt[gf];
   
  int jcb = a_jStartInt[gc];
  int jfb = a_jStartInt[gf];
  if (jfb % 2 == 1) jfb++; // make sure jfb is even
  

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
  for( int j=jfb ; j <= jfe ; j+=2 )
    for( int i=ifb ; i <= ife ; i+=2 )
    {
      int ic, jc;
      float_sw4 b1, a11, r3;
// updated components 1,2 of the ghost point value of Uf
	for( int c=1 ; c <= 2 ; c++ )
	{
	  {
	    ic = i/2;
	    jc = j/2;
// All Unextc terms
	    b1 = i256*
	      ( Unextc(c,ic-1,jc-1,1)-9*(Unextc(c,ic,jc-1,1)+Unextc(c,ic+1,jc-1,1))+Unextc(c,ic+2,jc-1,1)
		+ 9*(-Unextc(c,ic-1,jc,  1)+9*(Unextc(c,ic,jc,  1)+Unextc(c,ic+1,jc,  1))-Unextc(c,ic+2,jc,  1)  
		     -Unextc(c,ic-1,jc+1,1)+9*(Unextc(c,ic,jc+1,1)+Unextc(c,ic+1,jc+1,1))-Unextc(c,ic+2,jc+1,1))
		+Unextc(c,ic-1,jc+2,1)-9*(Unextc(c,ic,jc+2,1)+Unextc(c,ic+1,jc+2,1))+Unextc(c,ic+2,jc+2,1) );

// All Uc terms
	    b1 = b1 + nuc*a_ghcof[0]*i256*(
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
	  }
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
	{
	  ic = i/2;
	  jc = j/2;
// All Unextc terms
	  b1 = i256*
	    ( Unextc(3,ic-1,jc-1,1)-9*(Unextc(3,ic,jc-1,1)+Unextc(3,ic+1,jc-1,1))+Unextc(3,ic+2,jc-1,1)
	      + 9*(-Unextc(3,ic-1,jc,  1)+9*(Unextc(3,ic,jc,  1)+Unextc(3,ic+1,jc,  1))-Unextc(3,ic+2,jc,  1)  
		   -Unextc(3,ic-1,jc+1,1)+9*(Unextc(3,ic,jc+1,1)+Unextc(3,ic+1,jc+1,1))-Unextc(3,ic+2,jc+1,1))
	      +Unextc(3,ic-1,jc+2,1)-9*(Unextc(3,ic,jc+2,1)+Unextc(3,ic+1,jc+2,1))+Unextc(3,ic+2,jc+2,1) );

// All Uc terms
	  b1 = b1 + nuc*a_ghcof[0]*i256*(
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
	} // end  j even, i even
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
} // end evenIevenJinterp


  
      
