c-----------------------------------------------------------------------
c Adds 4th order artificial disssipation for super-grid damping layers
c
c-----------------------------------------------------------------------

	subroutine addsgd( dt, h, up, u, um, rho, dcx, dcy,
     +    dcz, ifirst, ilast, jfirst, jlast, kfirst, klast, beta)
	implicit none
	real*8 dt, h
	real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast)
	real*8 dcy(jfirst:jlast)
	real*8 dcz(kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability condition on beta?

c this routine uses un-divided differences in x and t
	real*8 coeff

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c this scaling is used in the 1D test code
	coeff = beta*dt/h
c	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c
c rho u_{tt} = h^3 *(h/dt)*beta*( psi(x) u_{xxt} )_{xx} + (same in y and z) )
c Note: h/dt has the dimension of a velocity!
c Note need to divide by \tilde{rho} = rho/phi
c
c add in the SG damping
c

c Operator is applied at the free surface boundary k=kfirst+1
        k = kfirst+1
	do j=jfirst+2,jlast-2
	   do i=ifirst+2, ilast-2
	      do c=1,3
		 up(c,i,j,k) = up(c,i,j,k) - coeff*(
c x-differences
     +  dcx(i+1) * ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
     + -2*dcx(i) * ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
     + +dcx(i-1) * ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
     + -dcx(i+1) * (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
     + +2*dcx(i) * (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
     + -dcx(i-1) * (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) 
c y-differences
     + +dcy(j+1) * ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
     + -2*dcy(j) * ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
     + +dcy(j-1) * ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
     + -dcy(j+1) * (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
     + +2*dcy(j) * (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
     + -dcy(j-1) * (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) 
     +		     )/rho(i,j,k)
	      enddo
	   enddo
	enddo
c Interior
	do k=kfirst+2,klast-2
	  do j=jfirst+2,jlast-2
	    do i=ifirst+2, ilast-2
	      do c=1,3
		up(c,i,j,k) = up(c,i,j,k) - coeff*(
c x-differences
     + dcx(i+1) * (u(c,i+2,j,k)-2*u(c,i+1,j,k)+u(c,i,j,k))
     + -2*dcx(i) * (u(c,i+1,j,k)-2*u(c,i,j,k)+u(c,i-1,j,k))
     + +dcx(i-1) * (u(c,i,j,k)-2*u(c,i-1,j,k)+u(c,i-2,j,k)) 
     + -dcx(i+1) * (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,j,k)) 
     + +2*dcx(i) * (um(c,i+1,j,k)-2*um(c,i,j,k)+um(c,i-1,j,k)) 
     + -dcx(i-1) * (um(c,i,j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) 
c y-differences
     + +dcy(j+1) * (u(c,i,j+2,k)-2*u(c,i,j+1,k)+u(c,i,j,k)) 
     + -2*dcy(j) * (u(c,i,j+1,k)-2*u(c,i,j,k)+u(c,i,j-1,k))
     + +dcy(j-1) * (u(c,i,j,k)-2*u(c,i,j-1,k)+u(c,i,j-2,k)) 
     + -dcy(j+1) * (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,k)) 
     + +2*dcy(j) * (um(c,i,j+1,k)-2*um(c,i,j,k)+um(c,i,j-1,k)) 
     + -dcy(j-1) * (um(c,i,j,k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) 
c z-differences
     + +dcz(k+1) * ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
     + -2*dcz(k) * ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
     + +dcz(k-1) * ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
     + -dcz(k+1) * (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
     + +2*dcz(k) * (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
     + -dcz(k-1) * (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) 
     + ) / rho(i,j,k)
	      enddo
	    enddo
	  enddo
	enddo
	end

c-----------------------------------------------------------------------
	subroutine addsgd2( dt, h, up, u, um, rho, 
     *               dcx, dcy, dcz, strx, stry, strz, 
     *	             ifirst, ilast, jfirst, jlast, kfirst, klast, beta )

***********************************************************************
*** New version with correct density scaling and supergrid stretching.
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
	real*8 coeff, irho

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
c	coeff = beta*dt/h
	coeff = beta
c beta is the supergrid damping coefficient as entered in the input file
c
c add in the SG damping
c

c Operator is applied at the free surface boundary k=kfirst+1
        k = kfirst+1
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
     *              (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) )
	      enddo
	   enddo
	enddo
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
	subroutine addsgd26( dt, h, up, u, um, rho, 
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
     +    (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*
     *   ( u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i,  j,k) )
     + -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*
     *   ( u(c,i+2,j,k) -3*u(c,i+1,j,k)+ 3*u(c,i,  j,k)- u(c,i-1,j,k) )
     + +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*
     *   ( u(c,i+1,j,k) -3*u(c,i,  j,k)+ 3*u(c,i-1,j,k)- u(c,i-2,j,k) )
     +  - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*
     *   ( u(c,i,  j,k) -3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) )
     +  - (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*
     *  (um(c,i+3,j,k) -3*um(c,i+2,j,k)+ 3*um(c,i+1,j,k)- um(c,i,  j,k))
     + +3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*
     *  (um(c,i+2,j,k) -3*um(c,i+1,j,k)+ 3*um(c,i,  j,k)- um(c,i-1,j,k))
     + -3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*
     *  (um(c,i+1,j,k) -3*um(c,i,  j,k)+ 3*um(c,i-1,j,k)- um(c,i-2,j,k))
     +  + (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*
     *  (um(c,i,  j,k) -3*um(c,i-1,j,k)+ 3*um(c,i-2,j,k)- um(c,i-3,j,k)) 
     +            ) +  stry(j)*(
c y-differences
     +    (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*
     *   ( u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k) )
     + -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*
     *   ( u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) )
     + +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*
     *   ( u(c,i,j+1,k) -3*u(c,i,  j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) )
     +  - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*
     *   ( u(c,i,  j,k) -3*u(c,i,j-1,k)+ 3*u(c,i,j-2,k)- u(c,i,j-3,k) )
     +  - (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*
     * (um(c,i,j+3,k) -3*um(c,i,j+2,k)+ 3*um(c,i,j+1,k)- um(c,i,  j,k))
     + +3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*
     * (um(c,i,j+2,k) -3*um(c,i,j+1,k)+ 3*um(c,i,  j,k)- um(c,i,j-1,k))
     + -3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*
     * (um(c,i,j+1,k) -3*um(c,i,  j,k)+ 3*um(c,i,j-1,k)- um(c,i,j-2,k))
     +  + (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*
     * (um(c,i,  j,k) -3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)- um(c,i,j-3,k))
     +            ) +  strz(k)*(
c z-differences
     +    (rho(i,j,k+2)*dcz(k+2)+rho(i,j,k+1)*dcz(k+1))*
     *   ( u(c,i,j,k+3) -3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) )
     + -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*
     *   ( u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) )
     + +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*
     *   ( u(c,i,j,k+1) -3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)- u(c,i,j,k-2) )
     +  - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*
     *   ( u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3) )
     +  -  (rho(i,j,k+2)*dcz(k+2)+rho(i,j,k+1)*dcz(k+1))*
     * ( um(c,i,j,k+3) -3*um(c,i,j,k+2)+ 3*um(c,i,j,k+1)- um(c,i,  j,k))
     + +3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*
     * ( um(c,i,j,k+2) -3*um(c,i,j,k+1)+ 3*um(c,i,  j,k)- um(c,i,j,k-1))
     + -3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*
     * ( um(c,i,j,k+1) -3*um(c,i,  j,k)+ 3*um(c,i,j,k-1)- um(c,i,j,k-2))
     +  + (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*
     * ( um(c,i,  j,k) -3*um(c,i,j,k-1)+ 3*um(c,i,j,k-2)- um(c,i,j,k-3))
     *    )  )
	      enddo
	    enddo
	  enddo
	enddo
	end
