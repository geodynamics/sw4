#include "EW.h"
//-----------------------------------------------------------------------
void EW::consintp( Sarray& Uf, Sarray& Unextf, Sarray& Bf, Sarray& Muf, Sarray& Lambdaf, Sarray& Rhof, double hf,
		   Sarray& Uc, Sarray& Unextc, Sarray& Bc, Sarray& Muc, Sarray& Lambdac, Sarray& Rhoc, double hc,
		   double cof, int gc, int gf, int is_periodic[2] )
{
   // At boundaries to the left and right, at least three ghost points are required
   // e.g., domain in i-direction:   i=-2,-1,0,1,2,...,Ni,Ni+1,Ni+2,Ni+3
   // we solve for ghost points at i=2,..,Ni-1, assuming Dirichlet conditions given on i=1,i=Ni and
   // at the ghost points i=-2,-1,0,Ni+1,Ni+2,Ni+3. 
   // The arrays Bf,Unextf, etc are assumed to be computed with all z-ghost points zero,
   // i.e., Uf(i,Nkf+1) was zero for i=-2,-1,..,Ni+3 when Bf and Unextf were evaluated from Uf 
   // (similarly for Unextc and Bc). Before this routine is called, correct boundary values for
   // Uf(-2,Nkf+1),..Uf(1,Nkf+1) (and similarly at the upper i-boundary, and for Uc) must be imposed.
   // In this way the restriction and prolongation stencils can be computed without any special
   // treatment near the i-boundaries. 

   const double i16 = 1.0/16;
   const double i256 = 1.0/256;
   const double i1024 = 1.0/1024;

   int jcb, jce, icb, ice, jfb, jfe, ifb, ife, nkf;
   double nuf = mDt*mDt/(cof*hf*hf); // cof is an argument to this routine
   double nuc = mDt*mDt/(cof*hc*hc);
   double ihc = 1/hc, ihf=1/hf;
   double jacerr = m_citol+1,jacerr0;
   double a11, a12, a21, a22, b1, b2, r1, r2, r3, deti, relax;
   int it = 0;
   relax = m_cirelfact;
 
   icb = m_iStartInt[gc];
   ifb = m_iStartInt[gf];

   ice = m_iEndInt[gc];
   ife = m_iEndInt[gf];
   
   jcb = m_jStartInt[gc];
   jfb = m_jStartInt[gf];

   jce = m_jEndInt[gc];
   jfe = m_jEndInt[gf];

   nkf = m_global_nz[gf];
   Sarray Mlf(m_iStart[gf],m_iEnd[gf],m_jStart[gf],m_jEnd[gf],nkf,nkf);
   for( int j=m_jStart[gf] ; j<=m_jEnd[gf] ; j++ )
      for( int i=m_iStart[gf] ; i<=m_iEnd[gf] ; i++ )
	 Mlf(i,j,nkf) = 2*Muf(i,j,nkf)+Lambdaf(i,j,nkf); // 2*mu + lambda on the fine grid
   Sarray Morc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   Sarray Mlrc(m_iStart[gc],m_iEnd[gc],m_jStart[gc],m_jEnd[gc],1,1);
   for( int jc=m_jStart[gc] ; jc<=m_jEnd[gc] ; jc++ )
      for( int ic=m_iStart[gc] ; ic<=m_iEnd[gc] ; ic++ )
      {
	 double irho=1/Rhoc(ic,jc,1);
	 Morc(ic,jc,1) = Muc(ic,jc,1)*irho; // mu/rho on the coarse grid
	 Mlrc(ic,jc,1) = (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*irho; // (2*mu+lambda)/rho on the coarse grid
      }

   while( jacerr > m_citol && it < m_cimaxiter )
   {
      double rmax[6]={0,0,0,0,0,0};
// for i=2*ic-1 and j=2*jc-1: Enforce continuity of displacements and normal stresses along the interface
      for( int jc= jcb ; jc <= jce ; jc++ )
	 for( int ic= icb ; ic <= ice ; ic++ )
	 {
// i odd, j odd
	    int i=2*ic-1, j=2*jc-1;
            // setup 2x2 system matrix
	    a11 = 0.25*Muf(i,j,nkf)*m_sbop[0]*ihf; // ihf = 1/h on the fine grid
	    a12 =      Muc(ic,jc,1)*m_sbop[0]*ihc;  // ihc = 1/h on the coarse grid
	    a21 = nuf/Rhof(i,j,nkf)*Muf(i,j,nkf)*m_ghcof[0]; // nuf = dt^2/(cof * h^2) on the fine grid
	    a22 =-nuc/Rhoc(ic,jc,1)*Muc(ic,jc,1)*m_ghcof[0]; // nuc = dt^2/(cof * h^2) on the coarse grid
	    for( int c=1 ; c <= 2 ; c++ ) // AP: what are the 2 components? 
	    {
// apply the restriction operator to the normal stress on the interface (Bf is on the fine grid)
	       b1  = i1024*( 
                    Bf(c,i-3,j-3,nkf)-9*Bf(c,i-3,j-1,nkf)-16*Bf(c,i-3,j,nkf)-9*Bf(c,i-3,j+1,nkf)+Bf(c,i-3,j+3,nkf) +
	        9*(-Bf(c,i-1,j-3,nkf)+9*Bf(c,i-1,j-1,nkf)+16*Bf(c,i-1,j,nkf)+9*Bf(c,i-1,j+1,nkf)-Bf(c,i-1,j+3,nkf)) +
	       16*(-Bf(c,i,  j-3,nkf)+9*Bf(c,i,  j-1,nkf)+16*Bf(c,i,  j,nkf)+9*Bf(c,i,  j+1,nkf)-Bf(c,i,  j+3,nkf)) + 
	        9*(-Bf(c,i+1,j-3,nkf)+9*Bf(c,i+1,j-1,nkf)+16*Bf(c,i+1,j,nkf)+9*Bf(c,i+1,j+1,nkf)-Bf(c,i+1,j+3,nkf)) +
	  	    Bf(c,i+3,j-3,nkf)-9*Bf(c,i+3,j-1,nkf)-16*Bf(c,i+3,j,nkf)-9*Bf(c,i+3,j+1,nkf)+Bf(c,i+3,j+3,nkf) );
	       b1 = b1 - i1024*m_sbop[0]*ihf*(
                Muf(i-3,j-3,nkf)*Uf(c,i-3,j-3,nkf+1) - 9*Muf(i-3,j-1,nkf)*Uf(c,i-3,j-1,nkf+1)
   	    -16*Muf(i-3,j,  nkf)*Uf(c,i-3,  j,nkf+1) - 9*Muf(i-3,j+1,nkf)*Uf(c,i-3,j+1,nkf+1)
               +Muf(i-3,j+3,nkf)*Uf(c,i-3,j+3,nkf+1) +					    
	  9*(  -Muf(i-1,j-3,nkf)*Uf(c,i-1,j-3,nkf+1) + 9*Muf(i-1,j-1,nkf)*Uf(c,i-1,j-1,nkf+1) 
	    +16*Muf(i-1,j,  nkf)*Uf(c,i-1,j,  nkf+1) + 9*Muf(i-1,j+1,nkf)*Uf(c,i-1,j+1,nkf+1) 
	       -Muf(i-1,j+3,nkf)*Uf(c,i-1,j+3,nkf+1) ) +
	 16*(  -Muf(i,  j-3,nkf)*Uf(c,i,  j-3,nkf+1) + 9*Muf(i,  j-1,nkf)*Uf(c,i,  j-1,nkf+1) 
	                                             + 9*Muf(i,  j+1,nkf)*Uf(c,i,  j+1,nkf+1)
               -Muf(i,  j+3,nkf)*Uf(c,i,  j+3,nkf+1) ) + 
	  9*(  -Muf(i+1,j-3,nkf)*Uf(c,i+1,j-3,nkf+1) + 9*Muf(i+1,j-1,nkf)*Uf(c,i+1,j-1,nkf+1)
            +16*Muf(i+1,j,  nkf)*Uf(c,i+1,j,  nkf+1) + 9*Muf(i+1,j+1,nkf)*Uf(c,i+1,j+1,nkf+1)
	       -Muf(i+1,j+3,nkf)*Uf(c,i+1,j+3,nkf+1) ) +
	        Muf(i+3,j-3,nkf)*Uf(c,i+3,j-3,nkf+1) - 9*Muf(i+3,j-1,nkf)*Uf(c,i+3,j-1,nkf+1)
	    -16*Muf(i+3,j,  nkf)*Uf(c,i+3,j,  nkf+1) - 9*Muf(i+3,j+1,nkf)*Uf(c,i+3,j+1,nkf+1)
	       +Muf(i+3,j+3,nkf)*Uf(c,i+3,j+3,nkf+1) );

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
	       r1 = r1-Uf(c,i,j,nkf+1);
	       r2 = r2-Uc(c,ic,jc,0);
	       rmax[c-1] = rmax[c-1] > fabs(r1) ? rmax[c-1] : fabs(r1);
	       rmax[c-1] = rmax[c-1] > fabs(r2) ? rmax[c-1] : fabs(r2);
	       //	       if( c == 2 && ic == 12 && jc == 13 )
	       //	       {
	       //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	       //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	       //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	       //	       }
	    } // end for c=1,2
            
// setup the matrix for the 3rd component of the normal stress (different coefficients)
	    a11 = 0.25*Mlf(i,j,nkf)*m_sbop[0]*ihf;
	    a12 =    (2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_sbop[0]*ihc;
	    a21 = nuf/Rhof(i,j,nkf)*Mlf(i,j,nkf)*m_ghcof[0];
	    a22 =-nuc/Rhoc(ic,jc,1)*(2*Muc(ic,jc,1)+Lambdac(ic,jc,1))*m_ghcof[0];
// apply the restriction operator to the fine grid normal stress grid function (Bf)
	    b1  = i1024*( 
                    Bf(3,i-3,j-3,nkf)-9*Bf(3,i-3,j-1,nkf)-16*Bf(3,i-3,j,nkf)-9*Bf(3,i-3,j+1,nkf)+Bf(3,i-3,j+3,nkf) +
	        9*(-Bf(3,i-1,j-3,nkf)+9*Bf(3,i-1,j-1,nkf)+16*Bf(3,i-1,j,nkf)+9*Bf(3,i-1,j+1,nkf)-Bf(3,i-1,j+3,nkf)) +
	       16*(-Bf(3,i,  j-3,nkf)+9*Bf(3,i,  j-1,nkf)+16*Bf(3,i,  j,nkf)+9*Bf(3,i,  j+1,nkf)-Bf(3,i,  j+3,nkf)) + 
	        9*(-Bf(3,i+1,j-3,nkf)+9*Bf(3,i+1,j-1,nkf)+16*Bf(3,i+1,j,nkf)+9*Bf(3,i+1,j+1,nkf)-Bf(3,i+1,j+3,nkf)) +
	  	    Bf(3,i+3,j-3,nkf)-9*Bf(3,i+3,j-1,nkf)-16*Bf(3,i+3,j,nkf)-9*Bf(3,i+3,j+1,nkf)+Bf(3,i+3,j+3,nkf) );
	    b1 = b1 - i1024*m_sbop[0]*ihf*(
                Mlf(i-3,j-3,nkf)*Uf(3,i-3,j-3,nkf+1) - 9*Mlf(i-3,j-1,nkf)*Uf(3,i-3,j-1,nkf+1)
   	    -16*Mlf(i-3,j,  nkf)*Uf(3,i-3,  j,nkf+1) - 9*Mlf(i-3,j+1,nkf)*Uf(3,i-3,j+1,nkf+1)
               +Mlf(i-3,j+3,nkf)*Uf(3,i-3,j+3,nkf+1) +					    
	  9*(  -Mlf(i-1,j-3,nkf)*Uf(3,i-1,j-3,nkf+1) + 9*Mlf(i-1,j-1,nkf)*Uf(3,i-1,j-1,nkf+1) 
	    +16*Mlf(i-1,j,  nkf)*Uf(3,i-1,j,  nkf+1) + 9*Mlf(i-1,j+1,nkf)*Uf(3,i-1,j+1,nkf+1) 
	       -Mlf(i-1,j+3,nkf)*Uf(3,i-1,j+3,nkf+1) ) +
	 16*(  -Mlf(i,  j-3,nkf)*Uf(3,i,  j-3,nkf+1) + 9*Mlf(i,  j-1,nkf)*Uf(3,i,  j-1,nkf+1) 
	                                             + 9*Mlf(i,  j+1,nkf)*Uf(3,i,  j+1,nkf+1)
               -Mlf(i,  j+3,nkf)*Uf(3,i,  j+3,nkf+1) ) + 
	  9*(  -Mlf(i+1,j-3,nkf)*Uf(3,i+1,j-3,nkf+1) + 9*Mlf(i+1,j-1,nkf)*Uf(3,i+1,j-1,nkf+1)
            +16*Mlf(i+1,j,  nkf)*Uf(3,i+1,j,  nkf+1) + 9*Mlf(i+1,j+1,nkf)*Uf(3,i+1,j+1,nkf+1)
	       -Mlf(i+1,j+3,nkf)*Uf(3,i+1,j+3,nkf+1) ) +
	        Mlf(i+3,j-3,nkf)*Uf(3,i+3,j-3,nkf+1) - 9*Mlf(i+3,j-1,nkf)*Uf(3,i+3,j-1,nkf+1)
	    -16*Mlf(i+3,j,  nkf)*Uf(3,i+3,j,  nkf+1) - 9*Mlf(i+3,j+1,nkf)*Uf(3,i+3,j+1,nkf+1)
	       +Mlf(i+3,j+3,nkf)*Uf(3,i+3,j+3,nkf+1) );

// setup the RHS
	       b1 = b1 - Bc(3,ic,jc,1);
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
	       rmax[2] = rmax[2] > fabs(r1) ? rmax[2] : fabs(r1);
	       rmax[2] = rmax[2] > fabs(r2) ? rmax[2] : fabs(r2);
	       //	       int c=3;
	       //	       if( c == 3 && ic == 12 && jc == 13 )
	       //	       {
	       //	         cout << "i,j " << i << " " << j << " " << b1 << " " << b2 << " " << r1 << " " << r2 << endl;
	       //		 cout << "   " << Uf(c,i,j,nkf+1) << " " << Uc(c,ic,jc,0) << " " << 
	       //		    a21*Uf(c,i,j,nkf+1)+Unextf(c,i,j,nkf) << " " << -a22*Uc(c,ic,jc,0)+Unextc(c,ic,jc,1) << endl;
	       //	       }
	 }
      //      goto skipthis;
      
// Enforce continuity of displacements along the interface (for fine ghost points in between coarse points)
      int ic, jc;
      for( int j=jfb ; j <= jfe ; j++ )
	 for( int i=ifb ; i <= ife ; i++ )
	 {
	    if( !( (i % 2 == 1 && j % 2 == 1 ) ) ) // not both i and j are odd (handled above)
	    {
// updated components 1,2 of the ghost point value of Uf
	       for( int c=1 ; c <= 2 ; c++ )
	       {
		  if( (j % 2 == 0) && (i % 2 == 1) ) // j is even, i is odd
		  {
		     ic = (i+1)/2;
		     jc = j/2;
		     b1 = i16*(-Unextc(c,ic,jc-1,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic,jc+1,1))-Unextc(c,ic,jc+2,1));
		     b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic,jc-1,0)*Morc(ic,jc-1,1) + 
						      9*Uc(c,ic,jc  ,0)*Morc(ic,jc  ,1) + 
						      9*Uc(c,ic,jc+1,0)*Morc(ic,jc+1,1)
						      -Uc(c,ic,jc+2,0)*Morc(ic,jc+2,1) );
		  }
		  if( (j % 2 == 1) && (i % 2 == 0) ) // j is odd, i is even
		  {
		     ic = i/2;
		     jc = (j+1)/2;
		     b1 = i16*(-Unextc(c,ic-1,jc,1)+9*(Unextc(c,ic,jc,1)+Unextc(c,ic+1,jc,1))-Unextc(c,ic+2,jc,1));
		     b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)+ 
						      9*Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+ 
						      9*Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1)
						      -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1));
		  }
		  if( (j % 2 == 0) && (i % 2 == 0) ) // i is even, j is even
		  {
		     ic = i/2;
		     jc = j/2;
		     b1 = i256*
                            ( Unextc(c,ic-1,jc-1,1)-9*(Unextc(c,ic,jc-1,1)+Unextc(c,ic+1,jc-1,1))+Unextc(c,ic+2,jc-1,1)
	 	        + 9*(-Unextc(c,ic-1,jc,  1)+9*(Unextc(c,ic,jc,  1)+Unextc(c,ic+1,jc,  1))-Unextc(c,ic+2,jc,  1)  
			     -Unextc(c,ic-1,jc+1,1)+9*(Unextc(c,ic,jc+1,1)+Unextc(c,ic+1,jc+1,1))-Unextc(c,ic+2,jc+1,1))
		             +Unextc(c,ic-1,jc+2,1)-9*(Unextc(c,ic,jc+2,1)+Unextc(c,ic+1,jc+2,1))+Unextc(c,ic+2,jc+2,1) );

		     b1 = b1 + nuc*m_ghcof[0]*i256*(
						    Uc(c,ic-1,jc-1,0)*Morc(ic-1,jc-1,1)
						-9*(Uc(c,ic,  jc-1,0)*Morc(ic,  jc-1,1)+
						    Uc(c,ic+1,jc-1,0)*Morc(ic+1,jc-1,1))+
						    Uc(c,ic+2,jc-1,0)*Morc(ic+2,jc-1,1)
			      + 9*(
				      -Uc(c,ic-1,jc,0)*Morc(ic-1,jc,1)+
				   9*( Uc(c,ic,  jc,0)*Morc(ic,  jc,1)+
                                       Uc(c,ic+1,jc,0)*Morc(ic+1,jc,1))
                                      -Uc(c,ic+2,jc,0)*Morc(ic+2,jc,1)  
				      -Uc(c,ic-1,jc+1,0)*Morc(ic-1,jc+1,1)+
				   9*( Uc(c,ic,  jc+1,0)*Morc(ic,  jc+1,1)+
                                       Uc(c,ic+1,jc+1,0)*Morc(ic+1,jc+1,1))
				      -Uc(c,ic+2,jc+1,0)*Morc(ic+2,jc+1,1)  )
					      + Uc(c,ic-1,jc+2,0)*Morc(ic-1,jc+2,1)
					    -9*(Uc(c,ic,  jc+2,0)*Morc(ic,  jc+2,1) +
			  		        Uc(c,ic+1,jc+2,0)*Morc(ic+1,jc+2,1)) +
						    Uc(c,ic+2,jc+2,0)*Morc(ic+2,jc+2,1)   );
		  }
		  b1 = b1 - Unextf(c,i,j,nkf);
		  	       a11 = nuf*m_ghcof[0]*Muf(i,j,nkf)/(Rhof(i,j,nkf));
		  r3 = Uf(c,i,j,nkf+1);
	       //	       Uf(c,i,j,nkf+1) = b1/a11;
// update ghost point value Uf(c,i,j,nkf+1)
		  Uf(c,i,j,nkf+1) = b1*Rhof(i,j,nkf)/(nuf*m_ghcof[0]*Muf(i,j,nkf));
		  Uf(c,i,j,nkf+1) = relax*Uf(c,i,j,nkf+1)+(1-relax)*r3;
		  //		  if( i == 4 && j == 7 && c == 1)
		  //		     cout << "in loop " << -a11*Uf(c,i,j,nkf+1) + b1  << endl;
		  r3 = r3 - Uf(c,i,j,nkf+1);
		  rmax[c-1+3] = rmax[c-1+3] > fabs(r3) ? rmax[c-1+3] : fabs(r3);
	       } // end for c=1,2
               
// work on componet 3 of the ghost point value of Uf
	       if( (j % 2 == 0) && (i % 2 == 1) ) // j even, i odd
	       {
		  ic = (i+1)/2;
		  jc = j/2;
		  b1 = i16*(-Unextc(3,ic,jc-1,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic,jc+1,1))-Unextc(3,ic,jc+2,1));
		  b1 = b1 + nuc*m_ghcof[0]*i16*(   - Uc(3,ic,jc-1,0)*Mlrc(ic,jc-1,1) + 
						   9*Uc(3,ic,jc  ,0)*Mlrc(ic,jc  ,1) + 
						   9*Uc(3,ic,jc+1,0)*Mlrc(ic,jc+1,1)
						    -Uc(3,ic,jc+2,0)*Mlrc(ic,jc+2,1) );
	       }
	       if( (j % 2 == 1) && (i % 2 == 0) ) // j odd, i even
	       {
		  ic = i/2;
		  jc = (j+1)/2;
		  b1 = i16*(-Unextc(3,ic-1,jc,1)+9*(Unextc(3,ic,jc,1)+Unextc(3,ic+1,jc,1))-Unextc(3,ic+2,jc,1));
		  b1 = b1 + nuc*m_ghcof[0]*i16*(   -Uc(3,ic-1,jc,0)*Mlrc(ic-1,jc,1)+ 
						  9*Uc(3,ic,  jc,0)*Mlrc(ic,  jc,1)+ 
						  9*Uc(3,ic+1,jc,0)*Mlrc(ic+1,jc,1)
						   -Uc(3,ic+2,jc,0)*Mlrc(ic+2,jc,1));
	       }
	       if( (j % 2 == 0) && (i % 2 == 0) ) // j even, i even
	       {
		  ic = i/2;
		  jc = j/2;
		  b1 = i256*
                            ( Unextc(3,ic-1,jc-1,1)-9*(Unextc(3,ic,jc-1,1)+Unextc(3,ic+1,jc-1,1))+Unextc(3,ic+2,jc-1,1)
	 	        + 9*(-Unextc(3,ic-1,jc,  1)+9*(Unextc(3,ic,jc,  1)+Unextc(3,ic+1,jc,  1))-Unextc(3,ic+2,jc,  1)  
			     -Unextc(3,ic-1,jc+1,1)+9*(Unextc(3,ic,jc+1,1)+Unextc(3,ic+1,jc+1,1))-Unextc(3,ic+2,jc+1,1))
		             +Unextc(3,ic-1,jc+2,1)-9*(Unextc(3,ic,jc+2,1)+Unextc(3,ic+1,jc+2,1))+Unextc(3,ic+2,jc+2,1) );

		  b1 = b1 + nuc*m_ghcof[0]*i256*(
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
	       a11 = nuf*m_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf))/(Rhof(i,j,nkf));
	       r3 = Uf(3,i,j,nkf+1);
               Uf(3,i,j,nkf+1) = b1/a11;
// solve for the ghost point value Uf(3,i,j,nkf+1)
//	       Uf(3,i,j,nkf+1) = b1*Rhof(i,j,nkf)/(nuf*m_ghcof[0]*(2*Muf(i,j,nkf)+Lambdaf(i,j,nkf)));
	       Uf(3,i,j,nkf+1) = relax*Uf(3,i,j,nkf+1)+(1-relax)*r3;
	       r3 = r3 - Uf(3,i,j,nkf+1);
	       rmax[2+3] = rmax[2+3] > fabs(r3) ? rmax[2+3] : fabs(r3);

	    } // end if not ( i%2=1 &&  j%2=1), i.e, not both i and j are odd

            // (i,j) both odd is handled by the first iteration
            
	 } // end for all fine grid points on the interface
      
      //   skipthis:
      communicate_array_2d( Uf, gf, nkf+1 );
      communicate_array_2d( Uc, gc, 0 );
      double jacerrtmp = 0;
      for (int q=0; q<6; q++)
         jacerrtmp += rmax[q];
      
      MPI_Allreduce( &jacerrtmp, &jacerr, 1, MPI_DOUBLE, MPI_MAX, m_cartesian_communicator );
      if( it == 0 )
	 jacerr0 = jacerr;
      if( jacerr0 > 0 )
	 jacerr = jacerr/jacerr0;
      it++;

   } // end while jacerr > eps
   
   if( jacerr > m_citol && proc_zero() )
      cout << "EW::consintp, Warning, no convergence. err = " << jacerr << " tol= " << m_citol << endl;
      
   if( proc_zero() && mVerbose >= 4 )
      cout << "EW::consintp, no of iterations= " << it << " Jac iteration error= " << jacerr << endl;
}
