c-----------------------------------------------------------------------
c Adds 4th order artificial disssipation for super-grid damping layers
c
c-----------------------------------------------------------------------

	subroutine addsgd( dt, h, up, u, um, dcx, dcy,
     +    dcz, ifirst, ilast, jfirst, jlast, kfirst, klast, beta)
	implicit none
	real*8 dt, h
	real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
	real*8 dcx(ifirst:ilast)
	real*8 dcy(jfirst:jlast)
	real*8 dcz(kfirst:klast)
	integer ifirst, ilast, jfirst, jlast, kfirst, klast
	real*8 beta

c time stepping stability gives the restriction
c tilde{beta} * dt/h < 1/6 (without elastic terms)
c this routine uses un-divided differences in x and t
	real*8 coeff

	integer i, j, k, c;

	if( beta .eq. 0d0 ) return;
	coeff = beta*dt/h
c beta is the supergrid damping coefficient as entered in the input file
c
c rho u_{tt} = h^3 rho*beta*( psi(x) u_{xxt} )_{xx} + (same in y and z) )
c
c add in the SG damping
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
     + +dcz(k+1) * (u(c,i,j,k+2)-2*u(c,i,j,k+1)+u(c,i,j,k)) 
     + -2*dcz(k) * (u(c,i,j,k+1)-2*u(c,i,j,k)+u(c,i,j,k-1))
     + +dcz(k-1) * (u(c,i,j,k)-2*u(c,i,j,k-1)+u(c,i,j,k-2)) 
     + -dcz(k+1) * (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k)) 
     + +2*dcz(k) * (um(c,i,j,k+1)-2*um(c,i,j,k)+um(c,i,j,k-1)) 
     + -dcz(k-1) * (um(c,i,j,k)-2*um(c,i,j,k-1)+um(c,i,j,k-2)) 
     + )
	      enddo
	    enddo
	  enddo
	enddo
	end
