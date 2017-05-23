!  SW4 LICENSE
! # ----------------------------------------------------------------------
! # SW4 - Seismic Waves, 4th order
! # ----------------------------------------------------------------------
! # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
! # Produced at the Lawrence Livermore National Laboratory. 
! # 
! # Written by:
! # N. Anders Petersson (petersson1@llnl.gov)
! # Bjorn Sjogreen      (sjogreen2@llnl.gov)
! # 
! # LLNL-CODE-643337 
! # 
! # All rights reserved. 
! # 
! # This file is part of SW4, Version: 1.0
! # 
! # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
! # 
! # This program is free software; you can redistribute it and/or modify
! # it under the terms of the GNU General Public License (as published by
! # the Free Software Foundation) version 2, dated June 1991. 
! # 
! # This program is distributed in the hope that it will be useful, but
! # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
! # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
! # conditions of the GNU General Public License for more details. 
! # 
! # You should have received a copy of the GNU General Public License
! # along with this program; if not, write to the Free Software
! # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
c-----------------------------------------------------------------------
c Adds 4th order artificial disssipation for super-grid damping layers
c
c-----------------------------------------------------------------------
	subroutine addsgd4( dt, h, up, u, um, rho, 
     *               dcx, dcy, dcz, strx, stry, strz, cox, coy, coz, 
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** Version with correct density scaling and supergrid stretching.
*** cox, coy, coz are corner factors that reduce the damping near edges and corners
***
***********************************************************************

	implicit none
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 dt, h
	real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8  rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
	real*8 dcz(kfirst:klast), strz(kfirst:klast), coz(kfirst:klast)
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irho

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;

	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c
c There are enough ghost points to always use the interior formula
c
c the corner tapering is applied by replacing
c strx -> strx*coy(j)*coz(k)
c stry -> stry*cox(i)*coz(k)
c strz -> strz*cox(i)*coy(j)
c
c approximately 375 a.o.
c
!$OMP PARALLEL PRIVATE(k,i,j,c,irho)
!$OMP DO
	do k=kfirst+2,klast-2
	  do j=jfirst+2,jlast-2
	    do i=ifirst+2, ilast-2
              irho = 1/rho(i,j,k)
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) - irho*coeff*( 
c x-differences
     + strx(i)*coy(j)*coz(k)*(
     +  rho(i+1,j,k)*dcx(i+1)*
     *              ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
     + -2*rho(i,j,k)*dcx(i)  *
     *              ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
     + +rho(i-1,j,k)*dcx(i-1)*
     *              ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
     + -rho(i+1,j,k)*dcx(i+1)*
     *              (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
     + +2*rho(i,j,k)*dcx(i)  *
     *              (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
     + -rho(i-1,j,k)*dcx(i-1)*
     *              (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
c y-differences
     + stry(j)*cox(i)*coz(k)*(
     + +rho(i,j+1,k)*dcy(j+1)*
     *              ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
     + -2*rho(i,j,k)*dcy(j)  *
     *              ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
     + +rho(i,j-1,k)*dcy(j-1)*
     *              ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
     + -rho(i,j+1,k)*dcy(j+1)*
     *              (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
     + +2*rho(i,j,k)*dcy(j)  *
     *              (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
     + -rho(i,j-1,k)*dcy(j-1)*
     *              (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) +
     +  strz(k)*cox(i)*coy(j)*(
c z-differences
     + +rho(i,j,k+1)*dcz(k+1)* 
     *            ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
     + -2*rho(i,j,k)*dcz(k)  *
     *            ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
     + +rho(i,j,k-1)*dcz(k-1)*
     *            ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
     + -rho(i,j,k+1)*dcz(k+1)*
     *            (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
     + +2*rho(i,j,k)*dcz(k)  *
     *            (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
     + -rho(i,j,k-1)*dcz(k-1)*
     *            (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) 
     + )
*** TOTAL 125 ops for each component = 375 ops per grid point. 
***       3x26  3D array accesses (u,um), 7 rho, gives = 85 elements of 3D arrays per grid point.
	      enddo
	    enddo
	  enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
	end

c-----------------------------------------------------------------------
	subroutine addsgd6( dt, h, up, u, um, rho, 
     *               dcx, dcy, dcz, strx, stry, strz, cox, coy, coz, 
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** 6th order damping, with correct density scaling and supergrid stretching.
***
***
***********************************************************************

	implicit none
	real*8 dt, h
	real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
	real*8 dcz(kfirst:klast), strz(kfirst:klast), coz(kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irhoh

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c	coeff = beta*dt/h
*** Divide by 2 for the averaged variable coefficient rho*dc
	coeff = beta*0.5d0
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c
!$OMP PARALLEL PRIVATE(k,i,j,c,irhoh)
!$OMP DO
	do k=kfirst+3,klast-3
	  do j=jfirst+3,jlast-3
	    do i=ifirst+3, ilast-3
              irhoh = 1/(rho(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) + irhoh*coeff*( 
     +  strx(i)*coy(j)*coz(k)*(
c x-differences
     +    (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*(
     *    u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
     * -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
     + -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*(
     *    u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
     * -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
     + +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*(
     *    u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
     * -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
     +  - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*(
     *    u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
     * -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
     +            ) +  stry(j)*cox(i)*coz(k)*(
c y-differences
     +    (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*(
     *    u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
     * -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     + -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*(
     *    u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
     * -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     + +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*(
     *    u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
     * -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
     +  - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*(
     *    u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
     * -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
     +            ) +  strz(k)*cox(i)*coy(j)*(
c z-differences
     +    ( rho(i,j,k+2)*dcz(k+2) + rho(i,j,k+1)*dcz(k+1) )*(
     *    u(c,i,j,k+3)- 3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) 
     * -(um(c,i,j,k+3)-3*um(c,i,j,k+2)+3*um(c,i,j,k+1)-um(c,i,  j,k)) )
     + -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*(
     *    u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) 
     * -(um(c,i,j,k+2)-3*um(c,i,j,k+1)+3*um(c,i,  j,k)-um(c,i,j,k-1)) )
     + +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*(
     *    u(c,i,j,k+1)- 3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)-u(c,i,j,k-2) 
     * -(um(c,i,j,k+1)-3*um(c,i,  j,k)+3*um(c,i,j,k-1)-um(c,i,j,k-2)) )
     +  - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*(
     *    u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3)
     * -(um(c,i,  j,k)-3*um(c,i,j,k-1)+3*um(c,i,j,k-2)-um(c,i,j,k-3)) )
     *    )  )
	      enddo
	    enddo
	  enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
	end

c-----------------------------------------------------------------------
	subroutine addsgd4c( dt, up, u, um, rho, 
     *               dcx, dcy, strx, stry, jac, cox, coy,
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** Correct density scaling and supergrid stretching.
*** Curvilinear version, assuming uniform grid in x and y.
***
***********************************************************************

	implicit none
	real*8 dt, h
	real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
        real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irhoj

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c the corner tapering is applied by replacing
c strx -> strx*coy(j)
c stry -> stry*cox(i)
c (no SG-stretching in the z-direction for the curvilinear grid)
c
c
c add in the SG damping
c    
!$OMP PARALLEL PRIVATE(k,i,j,c,irhoj)
!$OMP DO
	do k=kfirst+2,klast-2
	  do j=jfirst+2,jlast-2
	    do i=ifirst+2, ilast-2
              irhoj = 1/(rho(i,j,k)*jac(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) - irhoj*coeff*( strx(i)*coy(j)*(
c x-differences
     +  rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*(
     *                 u(c,i+2,j,k)- 2*u(c,i+1,j,k)+ u(c,i,  j,k)
     *              -(um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) )
     + -2*rho(i,j,k)*dcx(i)*jac(i,j,k)*(
     *                 u(c,i+1,j,k)- 2*u(c,i,  j,k)+ u(c,i-1,j,k)
     *              -(um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) )
     + +rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*(
     *                u(c,i,  j,k) - 2*u(c,i-1,j,k)+ u(c,i-2,j,k) 
     *             - (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) ) +
     +   stry(j)*cox(i)*(
c y-differences
     +  rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*(
     *               u(c,i,j+2,k) - 2*u(c,i,j+1,k) + u(c,i,j,  k)  
     *            -(um(c,i,j+2,k) -2*um(c,i,j+1,k) +um(c,i,j,  k)) ) 
     + -2*rho(i,j,k)*dcy(j)*jac(i,j,k)*(
     *               u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k)
     *            -(um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) ) 
     + +rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*(
     *                u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k) 
     *             -(um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) ) )
	      enddo
	    enddo
	  enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
	end

c-----------------------------------------------------------------------
	subroutine addsgd6c( dt, up, u, um, rho, 
     *               dcx, dcy, strx, stry, jac, cox, coy,
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** Correct density scaling and supergrid stretching.
*** Curvilinear version, assuming uniform grid in x and y.
***
***********************************************************************

	implicit none
	real*8 dt, h
	real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
        real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irhoj

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c
c the corner tapering is applied by replacing
c strx -> strx*coy(j)
c stry -> stry*cox(i)
c (no SG-stretching in the z-direction for the curvilinear grid)
c

*** Divide by 2 for the averaged variable coefficient rho*dc*jac
	coeff = beta*0.5d0
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c
!$OMP PARALLEL PRIVATE(k,i,j,c,irhoj)
!$OMP DO
	do k=kfirst+3,klast-3
	  do j=jfirst+3,jlast-3
	    do i=ifirst+3, ilast-3
              irhoj = 1/(rho(i,j,k)*jac(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) + irhoj*coeff*( strx(i)*coy(j)*(
     +    ( rho(i+2,j,k)*dcx(i+2)*jac(i+2,j,k)+
     +                       rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k) )*(
     *     u(c,i+3,j,k)-3*u(c,i+2,j,k)  +3*u(c,i+1,j,k)- u(c,i,  j,k) 
     *  -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i,  j,k)) )
     + -3*( rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)+
     *                            rho(i,j,k)*dcx(i)*jac(i,j,k) )*(
     *    u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i,  j,k)-u(c,i-1,j,k) 
     * -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i,  j,k)-um(c,i-1,j,k)) )
     + +3*( rho(i,j,k)*dcx(i)*jac(i,j,k)+
     *                       rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k) )*(
     *    u(c,i+1,j,k) -3*u(c,i,  j,k)+ 3*u(c,i-1,j,k)- u(c,i-2,j,k) 
     * -(um(c,i+1,j,k)-3*um(c,i,  j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
     +  -( rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)+
     *                      rho(i-2,j,k)*dcx(i-2)*jac(i-2,j,k) )*(
     *    u(c,i,  j,k) -3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
     * -(um(c,i,  j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) ) 
     +            ) +  stry(j)*cox(i)*(
c y-differences
     +    ( rho(i,j+2,k)*dcy(j+2)*jac(i,j+2,k)+
     *                           rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k) )*(
     *    u(c,i,j+3,k)-3*u(c,i,j+2,k) + 3*u(c,i,j+1,k)- u(c,i,  j,k) 
     * -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     + -3*( rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)+
     *                           rho(i,j,k)*dcy(j)*jac(i,j,k) )*(
     *    u(c,i,j+2,k)- 3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)-u(c,i,j-1,k)
     * -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     + +3*( rho(i,j,k)*dcy(j)*jac(i,j,k)+
     *                           rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k) )*(
     *    u(c,i,j+1,k)- 3*u(c,i,  j,k)+ 3*u(c,i,j-1,k)-u(c,i,j-2,k) 
     * -(um(c,i,j+1,k)-3*um(c,i,  j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) ) 
     +  - ( rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)+
     *                           rho(i,j-2,k)*dcy(j-2)*jac(i,j-2,k) )*(
     *    u(c,i,  j,k)-3*u(c,i,j-1,k) + 3*u(c,i,j-2,k)- u(c,i,j-3,k) 
     * -(um(c,i,  j,k)-3*um(c,i,j-1,k)+3*um(c,i,j-2,k)-um(c,i,j-3,k)) ) 
     *                ) )
	      enddo
	    enddo
	  enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
	end
