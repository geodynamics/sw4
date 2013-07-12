c-----------------------------------------------------------------------
c Adds 4th order artificial disssipation for super-grid damping layers
c
c-----------------------------------------------------------------------
	subroutine addsgd4( dt, h, up, u, um, rho, 
     *               dcx, dcy, dcz, strx, stry, strz, 
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** Version with correct density scaling and supergrid stretching.
***
***
***********************************************************************

	implicit none
	real*8 dt, h
	real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8  rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast)
	real*8 dcz(kfirst:klast), strz(kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irho

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c	coeff = beta*dt/h
	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c

c$$$c Operator is applied at the free surface boundary k=kfirst+1
c  Seems this is not needed, there are two ghost points over the free surface,
c  even if only one is used, it was a mistake to add this part
c$$$        k = kfirst+1
c$$$	do j=jfirst+2,jlast-2
c$$$	   do i=ifirst+2, ilast-2
c$$$              irho = 1/rho(i,j,k)
c$$$	      do c=1,3
c$$$		 up(c,i,j,k) = up(c,i,j,k) - irho*coeff*( strx(i)*(
c$$$c x-differences
c$$$     +  rho(i+1,j,k)*dcx(i+1)*
c$$$     *              ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
c$$$     + -2*rho(i,j,k)*dcx(i)  *
c$$$     *              ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
c$$$     + +rho(i-1,j,k)*dcx(i-1)*
c$$$     *              ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
c$$$     + -rho(i+1,j,k)*dcx(i+1)*
c$$$     *              (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
c$$$     + +2*rho(i,j,k)*dcx(i)  *
c$$$     *              (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
c$$$     + -rho(i-1,j,k)*dcx(i-1)*
c$$$     *              (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
c$$$     +   stry(j)*(
c$$$c y-differences
c$$$     + +rho(i,j+1,k)*dcy(j+1)*
c$$$     *              ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
c$$$     + -2*rho(i,j,k)*dcy(j)  *
c$$$     *              ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
c$$$     + +rho(i,j-1,k)*dcy(j-1)*
c$$$     *              ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
c$$$     + -rho(i,j+1,k)*dcy(j+1)*
c$$$     *              (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
c$$$     + +2*rho(i,j,k)*dcy(j)  *
c$$$     *              (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
c$$$     + -rho(i,j-1,k)*dcy(j-1)*
c$$$     *              (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) )
c$$$	      enddo
c$$$	   enddo
c$$$	enddo
c        write(*,*) 'dims ',ifirst,ilast,jfirst,jlast,kfirst,klast
c        do j=jfirst,jlast
c	   write(*,*) j,dcy(j)
c	enddo

c Interior
	do k=kfirst+2,klast-2
	  do j=jfirst+2,jlast-2
	    do i=ifirst+2, ilast-2
              irho = 1/rho(i,j,k)
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) - irho*coeff*( strx(i)*(
c x-differences
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
     +   stry(j)*(
c y-differences
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
     +  strz(k)*(	 
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
     *            (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) 
     + ) )
	      enddo
	    enddo
	  enddo
	enddo
	end

c-----------------------------------------------------------------------
	subroutine addsgd6( dt, h, up, u, um, rho, 
     *               dcx, dcy, dcz, strx, stry, strz, 
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
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast)
	real*8 dcz(kfirst:klast), strz(kfirst:klast)
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
	do k=kfirst+3,klast-3
	  do j=jfirst+3,jlast-3
	    do i=ifirst+3, ilast-3
              irhoh = 1/(rho(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) + irhoh*coeff*( strx(i)*(
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
     +            ) +  stry(j)*(
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
     +            ) +  strz(k)*(
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
	end

c-----------------------------------------------------------------------
	subroutine addsgd4c( dt, up, u, um, rho, 
     *               dcx, dcy, strx, stry, jac,
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
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast)
        real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irhoj

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c	coeff = beta*dt/h
	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c    
	do k=kfirst+2,klast-2
	  do j=jfirst+2,jlast-2
	    do i=ifirst+2, ilast-2
              irhoj = 1/(rho(i,j,k)*jac(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) - irhoj*coeff*( strx(i)*(
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
     +   stry(j)*(
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
	end

c-----------------------------------------------------------------------
	subroutine addsgd6c( dt, up, u, um, rho, 
     *               dcx, dcy, strx, stry, jac,
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
	real*8 dcx(ifirst:ilast), strx(ifirst:ilast)
	real*8 dcy(jfirst:jlast), stry(jfirst:jlast)
        real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff, irhoj

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c	coeff = beta*dt/h

*** Divide by 2 for the averaged variable coefficient rho*dc*jac
	coeff = beta*0.5d0
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c
	do k=kfirst+3,klast-3
	  do j=jfirst+3,jlast-3
	    do i=ifirst+3, ilast-3
              irhoj = 1/(rho(i,j,k)*jac(i,j,k))
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) + irhoj*coeff*( strx(i)*(
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
     +            ) +  stry(j)*(
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
	end
