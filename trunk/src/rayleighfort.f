      subroutine rayleighfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, u, t, lambda, mu, rho, cr, omega, alpha, h, zmin )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lambda, mu, rho, omega, alpha, h, t, zmin
      real*8 xi, latilde, dxi, cr, xi2, amp1, amp2
      real*8 x, y, z, phase, pi, up, xp
      pi = 4*ATAN(1d0)
c the limit mu->0 corresponds to Poisson ratio -> 0.5
c We used to solve the characteristic eqn on every call to this routine
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
      latilde = lambda/mu
      xi = cr/sqrt(mu/rho)
      xi2 = xi*xi
c      write(*,*),'Rayleigh params:', cRel, mu, cr, xi
      phase = 0
      amp1 = 1
      amp2 = -amp1*(2-xi2)/2
c      write(*,*)'Rayleigh:', omega, amp1, amp2, xi2, c, t, phase
      do k=kfirst,klast
        z = (k-1)*h + zmin
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            xp = x*cos(alpha) + y*sin(alpha)
            up = amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           omega/abs(omega) * sqrt(1-xi2) * 
     +           sin(omega*(cr*t + xp) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) * 
     *           omega/abs(omega) / sqrt(1-xi2/(2+latilde)) * 
     *           sin(omega*(cr*t + xp) + phase)

            u(1,i,j,k)= up*cos(alpha)
            u(2,i,j,k)= up*sin(alpha)
            u(3,i,j,k)= amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           cos(omega*(cr*t + xp) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) *
     *           cos(omega*(cr*t + xp) + phase)

          enddo
        enddo
      enddo
      end

c----------------------------------------------------------------------
      subroutine RAYDIRBDRY( bforce, wind, t, lambda, mu, rho, cr, 
     +     omega, alpha, h, zmin )
      implicit none
      integer wind(6)
      real*8 bforce(3,*), h, t, lambda, mu, rho, cr, omega, zmin
      real*8 xi, latilde, xi2, amp1, amp2, alpha
      real*8 x, y, z, phase, pi, up, xp
      integer i, j, k, qq
c
c NOTE: pass in the window for one side, i.e., wind(1,side) in the calling routine
c
      pi = 4*ATAN(1d0)
      latilde = lambda/mu
      xi = cr/sqrt(mu/rho)
      xi2 = xi*xi
c      write(*,*),'Rayleigh params:', cRel, mu, cr, xi
      phase = 0
      amp1 = 1
      amp2 = -amp1*(2-xi2)/2
*** Twilight forced rayleigh wave Dirichlet condition
      qq = 1
      do k=wind(5),wind(6)
c need to add zmin to work in a composite grid setting
        z = (k-1)*h + zmin
        do j=wind(3),wind(4)
          y = (j-1)*h
          do i=wind(1),wind(2)
            x = (i-1)*h
            xp = x*cos(alpha) + y*sin(alpha)
            up = amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           omega/abs(omega) * sqrt(1-xi2) * 
     +           sin(omega*(cr*t + xp) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) * 
     *           omega/abs(omega) / sqrt(1-xi2/(2+latilde)) * 
     *           sin(omega*(cr*t + xp) + phase)

            bforce(1,qq)= up*cos(alpha)
            bforce(2,qq)= up*sin(alpha)
            bforce(3,qq)= amp1 * exp( -abs(omega)*z*sqrt(1-xi2) ) * 
     *           cos(omega*(cr*t + xp) + phase) +
     *           amp2 * exp( -abs(omega)*z*sqrt(1- xi2/(2+latilde)) ) *
     *           cos(omega*(cr*t + xp) + phase)

            qq = qq+1
          enddo
        enddo
      enddo
      end

      function CSW( xi, latilde )
      real*8 xi, latilde
      CSW = sqrt(1-xi*xi)*sqrt(1-xi*xi/(2+latilde)) - (1-0.5d0*xi*xi)**2
      end


      function CSWP( xi, latilde )
      real*8 xi, latilde
      CSWP = -xi*sqrt(1-xi*xi/(2+latilde))/sqrt(1-xi*xi) -
     *      xi*sqrt(1-xi*xi)/(2+latilde)/sqrt(1-xi*xi/(2+latilde))
     *      + 2*xi*(1-0.5d0*xi*xi)
      end

c this routine computes the Rayleigh wave phase velocity, relative to the shear velocity
      real*8 function RVEL( la, mu )
      implicit none
      real*8 la, mu
      real*8 pr, cRel, xi, latilde, csw0, dxi
      real*8 CSW, CSWP
      integer iter
      pr = la/(2*(la+mu))
      cRel = (0.87d0 + 1.12d0*pr)/(1+pr)
      xi = cRel
      latilde = la/mu
      csw0 = CSW(xi, latilde)
      do iter=0,9
        dxi = -csw0/CSWP(xi, latilde)
        xi  = xi+dxi
        csw0 = CSW(xi, latilde)
      enddo
      RVEL = xi      
      end
