      subroutine addsgd4( dt, h, up, u, um, rho, 
     *     dcx, dcy, dcz, strx, stry, strz, cox, coy, coz, 
     *     ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

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

	integer i, j, k, c

	if( beta .eq. 0d0 ) return

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
