      subroutine rayleighfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, u, t, lambda, mu, rho, cr, omega, h, zmin )
      implicit none
c      real*8 CSW, CSWP
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lambda, mu, rho, omega, h, t, zmin
c      real*8 pr, cRel, xi, latilde, csw0, dxi, cr, xi2, amp1, amp2
      real*8 xi, latilde, dxi, cr, xi2, amp1, amp2
      real*8 x, y, z, phase, pi
      pi = 4*ATAN(1d0)
c the limit mu->0 corresponds to Poisson ratio -> 0.5
c      pr = lambda/(2*(lambda+mu))
c$$$      cRel = (0.87d0 + 1.12d0*pr)/(1+pr)
c$$$      xi = cRel
c$$$      latilde = lambda/mu
c$$$      csw0 = csw(xi, latilde)
c$$$      do iter=0,9
c$$$         dxi = -csw0/CSWP(xi, latilde)
c$$$         xi  = xi+dxi
c$$$         csw0 = CSW(xi, latilde)
c$$$      enddo
c$$$c      write(*,101) 'csw = ',csw0, ' xi=',xi
c$$$ 101  format(' Rayleigh', a, g12.5, a, g18.11 )
c$$$      cRel = xi      
c      cr  = cRel*sqrt(mu)
      xi = cr/sqrt(mu/rho)
      xi2 = xi*xi
c      write(*,*),'Rayleigh params:', cRel, mu, cr, xi
      phase = pi/2
      amp1 = 1
      amp2 = -amp1*(2-xi2)/2
c      write(*,*)'Rayleigh:', omega, amp1, amp2, xi2, c, t, phase
      do k=kfirst,klast
        z = (k-1)*h + zmin
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            u(1,i,j,k) = amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           omega/abs(omega) * sqrt(1-xi2) * 
     +           sin(omega*(cr*t + x) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) * 
     *           omega/abs(omega) / sqrt(1-xi2/(2+latilde)) * 
     *           sin(omega*(cr*t + x) + phase)

            u(2,i,j,k) = 0.

            u(3,i,j,k) = amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           cos(omega*(cr*t + x) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) *
     *           cos(omega*(cr*t + x) + phase)

          enddo
        enddo
      enddo
      end
