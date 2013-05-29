      subroutine GRIDINFO( ifirst, ilast, jfirst, jlast, kfirst, klast,
     *     met, jac, minj, maxj )
      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )
      real*8 fs, ot, ft, os, d3
      parameter( fs= 5d0/6, ot=1d0/12, ft=4d0/3, os=1d0/6, d3=14d0/3 )

      integer i, j, k, ifirst, ilast, jfirst, jlast
      integer kfirst, klast
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 minj, maxj

      maxj = -1d30
      minj =  1d30
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               if( jac(i,j,k).gt.maxj )then
                  maxj = jac(i,j,k)
               endif
               if( jac(i,j,k).lt.minj )then
                  minj = jac(i,j,k)
               endif
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine METRIC( ifirst, ilast, jfirst, jlast, kfirst, klast,
     *                   x, y, z, met, jac, ierr )
      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )
      real*8 fs, ot, ft, os, d3
      parameter( fs= 5d0/6, ot=1d0/12, ft=4d0/3, os=1d0/6, d3=14d0/3 )

      integer ni, nj, nk, i, j, k, ifirst, ilast, jfirst, jlast, ierr
      integer kfirst, klast
      real*8 x(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 y(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 h, zr, zp, zq, sqzr

      ierr=0
      h = x(ifirst+1,jfirst,kfirst)-x(ifirst,jfirst,kfirst)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast

*** k-derivatives
               if( k.ge.kfirst+2 .and. k.le.klast-2 )then
                  zr = c2*(z(i,j,k+2)-z(i,j,k-2)) +
     *                 c1*(z(i,j,k+1)-z(i,j,k-1))
               elseif( k.eq.kfirst )then
c                  zr = -(2+ot)*z(i,j,k)+4*z(i,j,k+1)-3*z(i,j,k+2)+
c     *                               ft*z(i,j,k+3)-0.25d0*z(i,j,k+4)
                  zr=-2.25d0*z(i,j,k)+(4+fs)*z(i,j,k+1)-d3*z(i,j,k+2)+
     *                   3*z(i,j,k+3)-(1+ot)*z(i,j,k+4) +os*z(i,j,k+5)
               elseif( k.eq.kfirst+1 )then
c                  zr = -0.25d0*z(i,j,k-1)-fs*z(i,j,k)+1.5d0*z(i,j,k+1)-
c     *                                   0.5d0*z(i,j,k+2)+ot*z(i,j,k+3)
                  zr = -os*z(i,j,k-1) -1.25d0*z(i,j,k)+(1+ft)*z(i,j,k+1)
     *                - ft*z(i,j,k+2) + 0.5d0*z(i,j,k+3) -ot*z(i,j,k+4)
               elseif( k.eq.klast-1 )then
c                  zr = 0.25d0*z(i,j,klast)+fs*z(i,j,klast-1)
c     *      -1.5d0*z(i,j,klast-2)+0.5d0*z(i,j,klast-3)-ot*z(i,j,klast-4)
                  zr =  os*z(i,j,k+1) +1.25d0*z(i,j,k)-(1+ft)*z(i,j,k-1)
     *                + ft*z(i,j,k-2) - 0.5d0*z(i,j,k-3) + ot*z(i,j,k-4)
               elseif( k.eq.klast )then
c                  zr = (2+ot)*z(i,j,klast)-4*z(i,j,klast-1)
c     *         +3*z(i,j,klast-2)-ft*z(i,j,klast-3)+0.25d0*z(i,j,klast-4)
                  zr= 2.25d0*z(i,j,k)-(4+fs)*z(i,j,k-1)+d3*z(i,j,k-2)-
     *                   3*z(i,j,k-3)+(1+ot)*z(i,j,k-4) -os*z(i,j,k-5)
               endif
               
*** j-derivatives
               if( j.ge.jfirst+2 .and. j.le.jlast-2 )then
                  zq = c2*(z(i,j+2,k)-z(i,j-2,k)) + 
     *                 c1*(z(i,j+1,k)-z(i,j-1,k))
               elseif( j.eq.jfirst )then
c                  zq = -(2+ot)*z(i,j,k)+4*z(i,j+1,k)-3*z(i,j+2,k)+
c     *                               ft*z(i,j+3,k)-0.25d0*z(i,j+4,k)
                  zq=-2.25d0*z(i,j,k)+(4+fs)*z(i,j+1,k)-d3*z(i,j+2,k)+
     *                   3*z(i,j+3,k)-(1+ot)*z(i,j+4,k) +os*z(i,j+5,k)
               elseif( j.eq.jfirst+1 )then
c                  zq = -0.25d0*z(i,j-1,k)-fs*z(i,j,k)+1.5d0*z(i,j+1,k)-
c     *                                   0.5d0*z(i,j+2,k)+ot*z(i,j+3,k)
                  zq = -os*z(i,j-1,k) -1.25d0*z(i,j,k)+(1+ft)*z(i,j+1,k)
     *                 - ft*z(i,j+2,k) + 0.5d0*z(i,j+3,k) -ot*z(i,j+4,k)
               elseif( j.eq.jlast-1 )then
c                  zq = 0.25d0*z(i,jlast,k)+fs*z(i,jlast-1,k)
c     *     -1.5d0*z(i,jlast-2,k)+0.5d0*z(i,jlast-3,k)-ot*z(i,jlast-4,k)
                  zq = os*z(i,j+1,k) +1.25d0*z(i,j,k)-(1+ft)*z(i,j-1,k)
     *               + ft*z(i,j-2,k) - 0.5d0*z(i,j-3,k) + ot*z(i,j-4,k)
               elseif( j.eq.jlast )then
c                  zq = (2+ot)*z(i,jlast,k)-4*z(i,jlast-1,k)
c     *      + 3*z(i,jlast-2,k)- ft*z(i,jlast-3,k)+0.25d0*z(i,jlast-4,k)
                  zq= 2.25d0*z(i,j,k)-(4+fs)*z(i,j-1,k)+d3*z(i,j-2,k)-
     *                   3*z(i,j-3,k)+(1+ot)*z(i,j-4,k) -os*z(i,j-5,k)
               endif

*** i-derivatives
               if( i.ge.ifirst+2 .and. i.le.ilast-2 )then
                  zp= c2*(z(i+2,j,k)-z(i-2,j,k)) + 
     *                c1*(z(i+1,j,k)-z(i-1,j,k))
               elseif( i.eq.ifirst )then
c                  zp = -(2+ot)*z(i,j,k)+4*z(i+1,j,k)-3*z(i+2,j,k)+
c     *                               ft*z(i+3,j,k)-0.25d0*z(i+4,j,k)
                  zp=-2.25d0*z(i,j,k)+(4+fs)*z(i+1,j,k)-d3*z(i+2,j,k)+
     *                   3*z(i+3,j,k)-(1+ot)*z(i+4,j,k) +os*z(i+5,j,k)
               elseif( i.eq.ifirst+1 )then
c                  zp =-0.25d0*z(i-1,j,k)-fs*z(i,j,k)+1.5d0*z(i+1,j,k)-
c     *                                  0.5d0*z(i+2,j,k)+ot*z(i+3,j,k)
                  zp = -os*z(i-1,j,k) -1.25d0*z(i,j,k)+(1+ft)*z(i+1,j,k)
     *                - ft*z(i+2,j,k) + 0.5d0*z(i+3,j,k) - ot*z(i+4,j,k)
               elseif( i.eq.ilast-1)then
c                  zp = 0.25d0*z(ilast,j,k)+fs*z(ilast-1,j,k)
c     *      -1.5d0*z(ilast-2,j,k)+0.5d0*z(ilast-3,j,k)-ot*z(ilast-4,j,k)
                  zp =  os*z(i+1,j,k) +1.25d0*z(i,j,k)-(1+ft)*z(i-1,j,k)
     *                + ft*z(i-2,j,k) - 0.5d0*z(i-3,j,k) + ot*z(i-4,j,k)
               elseif( i.eq.ilast)then
c                  zp = (2+ot)*z(ilast,j,k)-4*z(ilast-1,j,k)
c     *       +3*z(ilast-2,j,k)-  ft*z(ilast-3,j,k)+0.25d0*z(ilast-4,j,k)
                  zp= 2.25d0*z(i,j,k)-(4+fs)*z(i-1,j,k)+d3*z(i-2,j,k)-
     *                   3*z(i-3,j,k)+(1+ot)*z(i-4,j,k) -os*z(i-5,j,k)
               endif

*** Compute the metric
               if( zr.le.0 )then
                  write(*,101) 'Error, zr = ' , zr, ' at ', i, j, k
 101              format(' ', a, g12.5, a, 3(tr1,i3) )
c are you kidding me?                  stop
                  ierr=-1
                  return
               endif
c               if(i.eq.(ifirst+ilast)/2.and.j.eq.(jfirst+jlast)/2 )then
c                  write(*,102) k, zr, z(i,j,k)
c 102              format(' ', i3, tr2, 2(g15.7,tr2))
c               endif

               sqzr = sqrt(zr)
               jac(i,j,k) = h*h*zr
               met(1,i,j,k) = sqzr
               met(2,i,j,k) = -zp/sqzr
               met(3,i,j,k) = -zq/sqzr
               met(4,i,j,k) = h/sqzr
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine METRICEXGH( ifirst, ilast, jfirst, jlast, kfirst, 
     *                       klast, nx, ny, nz, x, y, z, met, jac,
     *                       order, sb, zmax, amp, xc, yc, xl, yl )
*** Exact metric derivatives for the Gaussian hill topography

      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )
      real*8 fs, ot, ft, os, d3
      parameter( fs= 5d0/6, ot=1d0/12, ft=4d0/3, os=1d0/6, d3=14d0/3 )

      integer ni, nj, nk, i, j, k, ifirst, ilast, jfirst, jlast
      integer kfirst, klast, nx, ny, nz, order, l
      real*8 x(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 y(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zmax, amp, xc, yc, xl, yl, ixl2, iyl2, zz
      real*8 h, zr, zp, zq, sqzr, s, tau, taup, tauq, sb, sdb, p1, p2
      h = x(ifirst+1,jfirst,kfirst)-x(ifirst,jfirst,kfirst)
      ixl2 = 1/(xl*xl)
      iyl2 = 1/(yl*yl)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               s = (k-1d0)/(nz-1)
               if( s.lt.sb )then
                  sdb = s/sb
                  tau  = amp*EXP( - (x(i,j,1)-xc)**2*ixl2 -
     *                 (y(i,j,1)-yc)**2*iyl2 )
                  taup = -2*(x(i,j,1)-xc)*ixl2*tau
                  tauq = -2*(y(i,j,1)-yc)*iyl2*tau
                  p1 = 1-sdb
                  p2 = 1
                  do l=2,order-1
                     p1 = p1 + (1-sdb)**l
                     p2 = p2 + l*(1-sdb)**(l-1)
                  enddo
                  zp = taup*( -(1-sdb)+sdb*p1 )
                  zq = tauq*( -(1-sdb)+sdb*p1)
                  zr = (tau+zmax+(zmax+tau-h*sb*(nz-1))*p1 -
     *                       sdb*(zmax+tau-h*sb*(nz-1))*p2 )/sb
                  zz = (1-sdb)*(-tau) + 
     *                     sdb*(zmax+(zmax+tau-h*sb*(nz-1))*p1)
               else
                  zp = 0
                  zq = 0
                  zr = h*(nz-1)
                  zz = zmax + (s-sb)*h*(nz-1)
               endif

c               if(i.eq.(ifirst+ilast)/2.and.j.eq.(jfirst+jlast)/2 )then
c                  write(*,101) k, zr, zr/(nz-1), z(i,j,k), zz
c 101              format(' ', i3, tr2, 4( g15.7,tr2 ) )
c               endif

*** Convert to 'undivided differences'
               zp = zp*h
               zq = zq*h
               zr = zr/(nz-1)
                  
*** Formulas from metric evaluation numerically
               sqzr = sqrt(zr)
               jac(i,j,k)   = h*h*zr
               met(1,i,j,k) = sqzr
               met(2,i,j,k) = -zp/sqzr
               met(3,i,j,k) = -zq/sqzr
               met(4,i,j,k) = h/sqzr
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine FREESURFCURVI( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, nz, side, u, mu, la, met, s, forcing )
      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k, kl, nz, side
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
      real*8 s(0:4), rhs1, rhs2, rhs3, s0i, ac, bc, cc, dc

      if( side.eq.5 )then
         k = 1
         kl= 1
      elseif( side.eq.6 )then
         k = nz
         kl= -1
      endif

      s0i = 1/s(0)
      do j=jfirst+2,jlast-2
         do i=ifirst+2,ilast-2
*** First tangential derivatives
            rhs1 = 
*** pr
     *   (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
     *          c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *          c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )   
*** qr
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
     *  + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )  -
     *                 forcing(1,i,j)

*** (v-eq)
            rhs2 = 
*** pr
     *    la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
*** qr
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
     * + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  ) 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   ) -
     *                  forcing(2,i,j)

*** (w-eq)
            rhs3 = 
*** pr
     *    la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )
*** qr 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   ) 
     *  + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  ) -
     *                  forcing(3,i,j)

*** Normal derivatives
            ac = met(2,i,j,k)**2+met(3,i,j,k)**2+met(4,i,j,k)**2
            bc = 1/(mu(i,j,k)*ac)
            cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac
            dc = cc*(met(2,i,j,k)*rhs1 + met(3,i,j,k)*rhs2 + 
     *           met(4,i,j,k)*rhs3)
            u(1,i,j,k-kl) = -s0i*(  s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+
     *           s(3)*u(1,i,j,k+2*kl)+s(4)*u(1,i,j,k+3*kl) + bc*rhs1 - 
     *                                           dc*met(2,i,j,k) )
            u(2,i,j,k-kl) = -s0i*(  s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+
     *           s(3)*u(2,i,j,k+2*kl)+s(4)*u(2,i,j,k+3*kl) + bc*rhs2 - 
     *                                           dc*met(3,i,j,k) )
            u(3,i,j,k-kl) = -s0i*(  s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+
     *           s(3)*u(3,i,j,k+2*kl)+s(4)*u(3,i,j,k+3*kl) + bc*rhs3 - 
     *                                           dc*met(4,i,j,k) )
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine GETSURFFORCING( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, met, jac, tau, forcing )
***********************************************************************
***
*** Given tau, Cartesian stress tensor on boundary, compute the stress
*** normal to the k=1 boundary of a curvilinear grid.
***
*** tau is ordered as:
***    tau(1) = t_{xx}, tau(2) = t_{xy} tau(3) = t_{xz} 
***    tau(4) = t_{yy}, tau(5) = t_{yz} tau(6) = t_{zz}
***
***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 sqjac
      do j=jfirst,jlast
         do i=ifirst,ilast
            sqjac = SQRT(jac(i,j,k))
            forcing(1,i,j) =  sqjac*( met(2,i,j,k)*tau(1,i,j)+
     *         met(3,i,j,k)*tau(2,i,j)+met(4,i,j,k)*tau(3,i,j) )
            forcing(2,i,j) =  sqjac*( met(2,i,j,k)*tau(2,i,j)+
     *         met(3,i,j,k)*tau(4,i,j)+met(4,i,j,k)*tau(5,i,j) )
            forcing(3,i,j) =  sqjac*( met(2,i,j,k)*tau(3,i,j)+
     *         met(3,i,j,k)*tau(5,i,j)+met(4,i,j,k)*tau(6,i,j) )
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine GETSURFFORCINGGH( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, h, tau, forcing, amp, xc, yc, xl, yl )
***********************************************************************
***
*** Given tau, Cartesian stress tensor on boundary, compute the stress
*** normal to the k=1 boundary of a curvilinear grid.
*** Use analytical normal, as given by the Gaussian Hill test case
***
*** tau is ordered as:
***    tau(1) = t_{xx}, tau(2) = t_{xy} tau(3) = t_{xz} 
***    tau(4) = t_{yy}, tau(5) = t_{yz} tau(6) = t_{zz}
***
***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
      real*8 amp, xc, yc, xl, yl, ixl2, iyl2, efact
      real*8 zp, zq, x, y, h
      ixl2 = 1/(xl*xl)
      iyl2 = 1/(yl*yl)
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x = (i-1)*h
            efact = amp*EXP( - (x-xc)**2*ixl2 - (y-yc)**2*iyl2 )
            zp = 2*(x-xc)*ixl2*efact
            zq = 2*(y-yc)*iyl2*efact
            forcing(1,i,j) =  h*h*( -zp*tau(1,i,j)
     *                              -zq*tau(2,i,j)+tau(3,i,j) )
            forcing(2,i,j) =  h*h*( -zp*tau(2,i,j)
     *                              -zq*tau(4,i,j)+tau(5,i,j) )
            forcing(3,i,j) =  h*h*( -zp*tau(3,i,j)
     *                              -zq*tau(5,i,j)+tau(6,i,j) )
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine CURVILINEAR4( ifirst, ilast, jfirst, jlast, kfirst,
     *                         klast, u, mu, la, met, jac, lu, 
     *                         onesided, acof, bope, ghcof )

      implicit none
      real*8 c1, c2, tf, i6, i144
      parameter( c1=2d0/3, c2=-1d0/12 )
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k, m, q, kstart, onesided(6)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 Lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 cof1, cof2, cof3, cof4, cof5, r1, r2, r3, ijac
      real*8 mux1, mux2, mux3, mux4
      real*8 mucofu2, mucofuv, mucofuw, mucofv2, mucofvw, mucofw2
      real*8 ghcof(6), acof(6,8,8), bope(6,8)
      real*8 dudrp2, dudrp1, dudrm1, dudrm2
      real*8 dvdrp2, dvdrp1, dvdrm1, dvdrm2
      real*8 dwdrp2, dwdrp1, dwdrm1, dwdrm2

*** met(1) is sqrt(J)*px = sqrt(J)*qy
*** met(2) is sqrt(J)*rx
*** met(3) is sqrt(J)*ry
*** met(4) is sqrt(J)*rz

      kstart = kfirst+2
      if( onesided(5).eq.1 )then
         kstart = 7

*** SBP Boundary closure terms
         do k=1,6
            do j=jfirst+2,jlast-2
               do i=ifirst+2,ilast-2
               ijac = 1/jac(i,j,k)
               r1 = 0
               r2 = 0
               r3 = 0

*** pp derivative (u) (u-eq)
          cof1=(2*mu(i-2,j,k)+la(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(2*mu(i-1,j,k)+la(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(2*mu(i+1,j,k)+la(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(2*mu(i+2,j,k)+la(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i+1,j,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i+2,j,k)-u(1,i,j,k))  )

*** qq derivative (u) (u-eq)
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j+2,k)-u(1,i,j,k))  )

*** pp derivative (v) (v-eq)
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r2 = r2 + i6* (
     *               mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i+2,j,k)-u(2,i,j,k))  )

*** qq derivative (v) (v-eq)
          cof1=(2*mu(i,j-2,k)+la(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(2*mu(i,j-1,k)+la(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(2*mu(i,j+1,k)+la(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(2*mu(i,j+2,k)+la(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r2 = r2 + i6* (
     *               mux1*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j+1,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j+2,k)-u(2,i,j,k))  )

*** pp derivative (w) (w-eq)
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r3 = r3 + i6* (
     *               mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i+2,j,k)-u(3,i,j,k))  )

*** qq derivative (w) (w-eq)
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r3 = r3 + i6* (
     *               mux1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j+2,k)-u(3,i,j,k))  )


*** All rr-derivatives at once
*** averaging the coefficient
            do q=1,8
               mucofu2=0
               mucofuv=0
               mucofuw=0
               mucofvw=0
               mucofv2=0
               mucofw2=0
               do m=1,8
                  mucofu2 = mucofu2 +
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*met(2,i,j,m)**2
     *                   + mu(i,j,m)*(met(3,i,j,m)**2+met(4,i,j,m)**2))
                  mucofv2 = mucofv2+
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*met(3,i,j,m)**2
     *                   + mu(i,j,m)*(met(2,i,j,m)**2+met(4,i,j,m)**2))
                  mucofw2 = mucofw2+
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*met(4,i,j,m)**2
     *                   + mu(i,j,m)*(met(2,i,j,m)**2+met(3,i,j,m)**2))
                  mucofuv = mucofuv+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(2,i,j,m)*met(3,i,j,m)
                  mucofuw = mucofuw+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(2,i,j,m)*met(4,i,j,m)
                  mucofvw = mucofvw+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(3,i,j,m)*met(4,i,j,m)
              enddo
*** Computing the second derivative,
              r1 = r1 + mucofu2*u(1,i,j,q) + mucofuv*u(2,i,j,q) + 
     *                                              mucofuw*u(3,i,j,q)
              r2 = r2 + mucofuv*u(1,i,j,q) + mucofv2*u(2,i,j,q) + 
     *                                              mucofvw*u(3,i,j,q)
              r3 = r3 + mucofuw*u(1,i,j,q) + mucofvw*u(2,i,j,q) + 
     *                                              mucofw2*u(3,i,j,q)
            end do

*** Ghost point values, only nonzero for k=1.
            mucofu2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*met(2,i,j,1)**2
     *           + mu(i,j,1)*(met(3,i,j,1)**2+met(4,i,j,1)**2))
            mucofv2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*met(3,i,j,1)**2
     *           + mu(i,j,1)*(met(2,i,j,1)**2+met(4,i,j,1)**2))
            mucofw2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*met(4,i,j,1)**2
     *           + mu(i,j,1)*(met(2,i,j,1)**2+met(3,i,j,1)**2))
            mucofuv = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(2,i,j,1)*met(3,i,j,1)
            mucofuw = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(2,i,j,1)*met(4,i,j,1)
            mucofvw = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(3,i,j,1)*met(4,i,j,1)
            r1 = r1 + mucofu2*u(1,i,j,0) + mucofuv*u(2,i,j,0) + 
     *                                              mucofuw*u(3,i,j,0)
            r2 = r2 + mucofuv*u(1,i,j,0) + mucofv2*u(2,i,j,0) + 
     *                                              mucofvw*u(3,i,j,0)
            r3 = r3 + mucofuw*u(1,i,j,0) + mucofvw*u(2,i,j,0) + 
     *                                              mucofw2*u(3,i,j,0)


*** pq-derivatives (u-eq)
      r1 = r1 + 
     *   c2*(  mu(i,j+2,k)*met(1,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k)) +
     *        c1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))    )
     *     - mu(i,j-2,k)*met(1,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
     *        c1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))     )
     *     ) +
     *   c1*(  mu(i,j+1,k)*met(1,i,j+1,k)*met(1,i,j+1,k)*(
     *          c2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k)) +
     *          c1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))  )
     *      - mu(i,j-1,k)*met(1,i,j-1,k)*met(1,i,j-1,k)*(
     *          c2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k)) + 
     *          c1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))))

*** qp-derivatives (u-eq)
      r1 = r1 + 
     *   c2*(  la(i+2,j,k)*met(1,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k)) +
     *        c1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))    )
     *     - la(i-2,j,k)*met(1,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
     *        c1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))     )
     *     ) +
     *   c1*(  la(i+1,j,k)*met(1,i+1,j,k)*met(1,i+1,j,k)*(
     *          c2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k)) +
     *          c1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))  )
     *      - la(i-1,j,k)*met(1,i-1,j,k)*met(1,i-1,j,k)*(
     *          c2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k)) + 
     *          c1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))))

*** pq-derivatives (v-eq)
      r2 = r2 + 
     *   c2*(  la(i,j+2,k)*met(1,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k)) +
     *        c1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))    )
     *     - la(i,j-2,k)*met(1,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
     *        c1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))     )
     *     ) +
     *   c1*(  la(i,j+1,k)*met(1,i,j+1,k)*met(1,i,j+1,k)*(
     *          c2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k)) +
     *          c1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))  )
     *      - la(i,j-1,k)*met(1,i,j-1,k)*met(1,i,j-1,k)*(
     *          c2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k)) + 
     *          c1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))))

*** qp-derivatives (v-eq)
      r2 = r2 + 
     *   c2*(  mu(i+2,j,k)*met(1,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k)) +
     *        c1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))    )
     *     - mu(i-2,j,k)*met(1,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
     *        c1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))     )
     *     ) +
     *   c1*(  mu(i+1,j,k)*met(1,i+1,j,k)*met(1,i+1,j,k)*(
     *          c2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k)) +
     *          c1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))  )
     *      - mu(i-1,j,k)*met(1,i-1,j,k)*met(1,i-1,j,k)*(
     *          c2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k)) + 
     *          c1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))))

*** rp - derivatives
         dudrm2 = 0
         dudrm1 = 0
         dudrp1 = 0
         dudrp2 = 0
         dvdrm2 = 0
         dvdrm1 = 0
         dvdrp1 = 0
         dvdrp2 = 0
         dwdrm2 = 0
         dwdrm1 = 0
         dwdrp1 = 0
         dwdrp2 = 0
         do q=1,8
            dudrm2 = dudrm2 + bope(k,q)*u(1,i-2,j,q)
            dvdrm2 = dvdrm2 + bope(k,q)*u(2,i-2,j,q)
            dwdrm2 = dwdrm2 + bope(k,q)*u(3,i-2,j,q)
            dudrm1 = dudrm1 + bope(k,q)*u(1,i-1,j,q)
            dvdrm1 = dvdrm1 + bope(k,q)*u(2,i-1,j,q)
            dwdrm1 = dwdrm1 + bope(k,q)*u(3,i-1,j,q)
            dudrp2 = dudrp2 + bope(k,q)*u(1,i+2,j,q)
            dvdrp2 = dvdrp2 + bope(k,q)*u(2,i+2,j,q)
            dwdrp2 = dwdrp2 + bope(k,q)*u(3,i+2,j,q)
            dudrp1 = dudrp1 + bope(k,q)*u(1,i+1,j,q)
            dvdrp1 = dvdrp1 + bope(k,q)*u(2,i+1,j,q)
            dwdrp1 = dwdrp1 + bope(k,q)*u(3,i+1,j,q)
         enddo

*** rp derivatives (u-eq)
      r1 = r1 + c2*(
     *  (2*mu(i+2,j,k)+la(i+2,j,k))*met(2,i+2,j,k)*met(1,i+2,j,k)*dudrp2
     *   + la(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*dvdrp2
     *   + la(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*dwdrp2
     *-((2*mu(i-2,j,k)+la(i-2,j,k))*met(2,i-2,j,k)*met(1,i-2,j,k)*dudrm2
     *   + la(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*dvdrm2
     *   + la(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*dwdrm2 )
     *              ) + c1*(  
     *  (2*mu(i+1,j,k)+la(i+1,j,k))*met(2,i+1,j,k)*met(1,i+1,j,k)*dudrp1
     *   + la(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*dvdrp1
     *   + la(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*dwdrp1 
     *-((2*mu(i-1,j,k)+la(i-1,j,k))*met(2,i-1,j,k)*met(1,i-1,j,k)*dudrm1
     *   + la(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*dvdrm1
     *   + la(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*dwdrm1 ) )

*** rp derivatives (v-eq)
      r2 = r2 + c2*(
     *     mu(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*dudrp2
     *  +  mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*dvdrp2
     *  - (mu(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*dudrm2
     *  +  mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*dvdrm2 )
     *             ) + c1*(  
     *     mu(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*dudrp1
     *  +  mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*dvdrp1
     *  - (mu(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*dudrm1
     *  +  mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*dvdrm1)
     *                     )

*** rp derivatives (w-eq)
      r3 = r3 + c2*(
     *     mu(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*dudrp2
     *  +  mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*dwdrp2
     *  - (mu(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*dudrm2
     *  +  mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*dwdrm2)
     *             ) + c1*(  
     *     mu(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*dudrp1
     *  +  mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*dwdrp1
     *  - (mu(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*dudrm1
     *  +  mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*dwdrm1)
     *                     ) 

*** rq - derivatives
         dudrm2 = 0
         dudrm1 = 0
         dudrp1 = 0
         dudrp2 = 0
         dvdrm2 = 0
         dvdrm1 = 0
         dvdrp1 = 0
         dvdrp2 = 0
         dwdrm2 = 0
         dwdrm1 = 0
         dwdrp1 = 0
         dwdrp2 = 0
         do q=1,8
            dudrm2 = dudrm2 + bope(k,q)*u(1,i,j-2,q)
            dvdrm2 = dvdrm2 + bope(k,q)*u(2,i,j-2,q)
            dwdrm2 = dwdrm2 + bope(k,q)*u(3,i,j-2,q)
            dudrm1 = dudrm1 + bope(k,q)*u(1,i,j-1,q)
            dvdrm1 = dvdrm1 + bope(k,q)*u(2,i,j-1,q)
            dwdrm1 = dwdrm1 + bope(k,q)*u(3,i,j-1,q)
            dudrp2 = dudrp2 + bope(k,q)*u(1,i,j+2,q)
            dvdrp2 = dvdrp2 + bope(k,q)*u(2,i,j+2,q)
            dwdrp2 = dwdrp2 + bope(k,q)*u(3,i,j+2,q)
            dudrp1 = dudrp1 + bope(k,q)*u(1,i,j+1,q)
            dvdrp1 = dvdrp1 + bope(k,q)*u(2,i,j+1,q)
            dwdrp1 = dwdrp1 + bope(k,q)*u(3,i,j+1,q)
         enddo

*** rq derivatives (u-eq)
      r1 = r1 + c2*(
     *      mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*dudrp2
     *   +  mu(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *   - (mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*dudrm2
     *   +  mu(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*dvdrm2)
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*dudrp1
     *   +  mu(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *   - (mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*dudrm1
     *   +  mu(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*dvdrm1)
     *         )

*** rq derivatives (v-eq)
      r2 = r2 + c2*(
     *      la(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*dudrp2
     * +(2*mu(i,j+2,k)+la(i,j+2,k))*met(3,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *    + la(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*dwdrp2
     *  - ( la(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*dudrm2
     * +(2*mu(i,j-2,k)+la(i,j-2,k))*met(3,i,j-2,k)*met(1,i,j-2,k)*dvdrm2
     *    + la(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*dwdrm2 )
     *             ) + c1*(  
     *      la(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*dudrp1
     * +(2*mu(i,j+1,k)+la(i,j+1,k))*met(3,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *    + la(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*dwdrp1
     *  - ( la(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*dudrm1
     * +(2*mu(i,j-1,k)+la(i,j-1,k))*met(3,i,j-1,k)*met(1,i,j-1,k)*dvdrm1
     *    + la(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*dwdrm1 )
     *        )

*** rq derivatives (w-eq)
      r3 = r3 + c2*(
     *     mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*dwdrp2
     *  +  mu(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *  - (mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*dwdrm2
     *  +  mu(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*dvdrm2)
     *             ) + c1*(  
     *     mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*dwdrp1
     *  +  mu(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *  - (mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*dwdrm1
     *  +  mu(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*dvdrm1)
     *               )

*** pr and qr derivatives at once
      do q=1,8
*** (u-eq)
        r1 = r1 + bope(k,q)*( 
*** pr
     *   (2*mu(i,j,q)+la(i,j,q))*met(2,i,j,q)*met(1,i,j,q)*(
     *          c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *          c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   ) 
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i+2,j,q)-u(2,i-2,j,q)) +
     *        c1*(u(2,i+1,j,q)-u(2,i-1,j,q))  ) 
     *  + mu(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i+2,j,q)-u(3,i-2,j,q)) +
     *        c1*(u(3,i+1,j,q)-u(3,i-1,j,q))  ) 
*** qr
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i,j+2,q)-u(1,i,j-2,q)) +
     *        c1*(u(1,i,j+1,q)-u(1,i,j-1,q))   ) 
     *  + la(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  ) )

*** (v-eq)
        r2 = r2 + bope(k,q)*(
*** pr
     *    la(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *        c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   ) 
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i+2,j,q)-u(2,i-2,j,q)) +
     *        c1*(u(2,i+1,j,q)-u(2,i-1,j,q))  ) 
*** qr
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i,j+2,q)-u(1,i,j-2,q)) +
     *        c1*(u(1,i,j+1,q)-u(1,i,j-1,q))   ) 
     * + (2*mu(i,j,q)+la(i,j,q))*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  ) 
     *  + mu(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i,j+2,q)-u(3,i,j-2,q)) +
     *        c1*(u(3,i,j+1,q)-u(3,i,j-1,q))   )  )

*** (w-eq)
        r3 = r3 + bope(k,q)*(
*** pr
     *    la(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *        c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   ) 
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i+2,j,q)-u(3,i-2,j,q)) +
     *        c1*(u(3,i+1,j,q)-u(3,i-1,j,q))  )
*** qr 
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i,j+2,q)-u(3,i,j-2,q)) +
     *        c1*(u(3,i,j+1,q)-u(3,i,j-1,q))   ) 
     *  + la(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  ) )

      enddo
          lu(1,i,j,k) = r1*ijac
          lu(2,i,j,k) = r2*ijac
          lu(3,i,j,k) = r3*ijac
               enddo
            enddo
         enddo
      endif

      do k=kstart,klast-2
         do j=jfirst+2,jlast-2
            do i=ifirst+2,ilast-2
               ijac = 1/jac(i,j,k)

               r1 = 0
*** pp derivative (u)
          cof1=(2*mu(i-2,j,k)+la(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(2*mu(i-1,j,k)+la(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(2*mu(i+1,j,k)+la(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(2*mu(i+2,j,k)+la(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i+1,j,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i+2,j,k)-u(1,i,j,k))  )
*** qq derivative (u)
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j+2,k)-u(1,i,j,k))  )
*** rr derivative (u)
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)**2 +
     *         mu(i,j,k-2)*(met(3,i,j,k-2)**2+met(4,i,j,k-2)**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)**2 +
     *         mu(i,j,k-1)*(met(3,i,j,k-1)**2+met(4,i,j,k-1)**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)**2 +
     *         mu(i,j,k)*(met(3,i,j,k)**2+met(4,i,j,k)**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)**2 +
     *         mu(i,j,k+1)*(met(3,i,j,k+1)**2+met(4,i,j,k+1)**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)**2 +
     *         mu(i,j,k+2)*(met(3,i,j,k+2)**2+met(4,i,j,k+2)**2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j,k+2)-u(1,i,j,k))  )

*** rr derivative (v)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(3,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(3,i,j,k-1)
          cof3=(mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(3,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(3,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(3,i,j,k+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j,k+2)-u(2,i,j,k))  )

*** rr derivative (w)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(4,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(4,i,j,k-1)
          cof3=(mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(4,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(4,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(4,i,j,k+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j,k+1)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )

*** pq-derivatives
      r1 = r1 + 
     *   c2*(  mu(i,j+2,k)*met(1,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k)) +
     *        c1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))    )
     *     - mu(i,j-2,k)*met(1,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
     *        c1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))     )
     *     ) +
     *   c1*(  mu(i,j+1,k)*met(1,i,j+1,k)*met(1,i,j+1,k)*(
     *          c2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k)) +
     *          c1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))  )
     *      - mu(i,j-1,k)*met(1,i,j-1,k)*met(1,i,j-1,k)*(
     *          c2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k)) + 
     *          c1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))))

*** qp-derivatives
      r1 = r1 + 
     *   c2*(  la(i+2,j,k)*met(1,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k)) +
     *        c1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))    )
     *     - la(i-2,j,k)*met(1,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
     *        c1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))     )
     *     ) +
     *   c1*(  la(i+1,j,k)*met(1,i+1,j,k)*met(1,i+1,j,k)*(
     *          c2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k)) +
     *          c1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))  )
     *      - la(i-1,j,k)*met(1,i-1,j,k)*met(1,i-1,j,k)*(
     *          c2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k)) + 
     *          c1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))))

*** pr-derivatives
      r1 = r1 + c2*(
     *  (2*mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   ) 
     *   + mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2)) +
     *        c1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))  ) 
     *   + mu(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2)) +
     *        c1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))  )
     *  - ((2*mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  ) 
     *     + mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2))   ) 
     *     + mu(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2))   ) )
     *             ) + c1*(  
     *     (2*mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1))  )
     *  - ((2*mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1)) +
     *        c1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1)) +
     *        c1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1))   )  ) )

*** rp derivatives
      r1 = r1 + c2*(
     *  (2*mu(i+2,j,k)+la(i+2,j,k))*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   ) 
     *   + la(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2)) +
     *        c1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))  ) 
     *   + la(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2)) +
     *        c1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1))  )
     *  - ((2*mu(i-2,j,k)+la(i-2,j,k))*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  ) 
     *     + la(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))   ) 
     *     + la(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1))   ) )
     *             ) + c1*(  
     *     (2*mu(i+1,j,k)+la(i+1,j,k))*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) ) 
     *     + la(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1)) ) 
     *     + la(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1))  )
     *  - ((2*mu(i-1,j,k)+la(i-1,j,k))*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) ) 
     *     + la(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2)) +
     *        c1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1)) ) 
     *     + la(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2)) +
     *        c1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1))   )  ) )

*** qr derivatives
      r1 = r1 + c2*(
     *    mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2)) +
     *        c1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))   ) 
     *   + la(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  ) 
     *  - ( mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))  ) 
     *     + la(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   ) ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1)) ) 
     *     + la(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )  
     *  - ( mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1)) +
     *        c1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1)) ) 
     *     + la(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) ) ) )

*** rq derivatives
      r1 = r1 + c2*(
     *    mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2)) +
     *        c1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))   ) 
     *   + mu(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  ) 
     *  - ( mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))  ) 
     *     + mu(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   ) ) 
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1)) ) 
     *     + mu(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) )  
     *  - ( mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2)) +
     *        c1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1)) ) 
     *     + mu(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) ) ) )

c          lu(1,i,j,k) = r1/jac(i,j,k)
          lu(1,i,j,k) = r1*ijac
*** v-equation

          r1 = 0
*** pp derivative (v)
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i+2,j,k)-u(2,i,j,k))  )
*** qq derivative (v)
          cof1=(2*mu(i,j-2,k)+la(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(2*mu(i,j-1,k)+la(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(2*mu(i,j+1,k)+la(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(2*mu(i,j+2,k)+la(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j+1,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j+2,k)-u(2,i,j,k))  )
*** rr derivative (u)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(3,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(3,i,j,k-1)
          cof3=(mu(i,j,k)+  la(i,j,k)  )*met(2,i,j,k)*  met(3,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(3,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(3,i,j,k+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j,k+2)-u(1,i,j,k))  )

*** rr derivative (v)
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*met(3,i,j,k-2)**2 +
     *         mu(i,j,k-2)*(met(2,i,j,k-2)**2+met(4,i,j,k-2)**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*met(3,i,j,k-1)**2 +
     *         mu(i,j,k-1)*(met(2,i,j,k-1)**2+met(4,i,j,k-1)**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)**2 +
     *         mu(i,j,k)*(met(2,i,j,k)**2+met(4,i,j,k)**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*met(3,i,j,k+1)**2 +
     *         mu(i,j,k+1)*(met(2,i,j,k+1)**2+met(4,i,j,k+1)**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*met(3,i,j,k+2)**2 +
     *         mu(i,j,k+2)*(met(2,i,j,k+2)**2+met(4,i,j,k+2)**2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j,k+2)-u(2,i,j,k))  )

*** rr derivative (w)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(3,i,j,k-2)*met(4,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(3,i,j,k-1)*met(4,i,j,k-1)
          cof3=(mu(i,j,k)  +la(i,j,k)  )*met(3,i,j,k)*  met(4,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(3,i,j,k+1)*met(4,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(3,i,j,k+2)*met(4,i,j,k+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j,k+1)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )

*** pq-derivatives
      r1 = r1 + 
     *   c2*(  la(i,j+2,k)*met(1,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k)) +
     *        c1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))    )
     *     - la(i,j-2,k)*met(1,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
     *        c1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))     )
     *     ) +
     *   c1*(  la(i,j+1,k)*met(1,i,j+1,k)*met(1,i,j+1,k)*(
     *          c2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k)) +
     *          c1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))  )
     *      - la(i,j-1,k)*met(1,i,j-1,k)*met(1,i,j-1,k)*(
     *          c2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k)) + 
     *          c1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))))

*** qp-derivatives
      r1 = r1 + 
     *   c2*(  mu(i+2,j,k)*met(1,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k)) +
     *        c1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))    )
     *     - mu(i-2,j,k)*met(1,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
     *        c1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))     )
     *     ) +
     *   c1*(  mu(i+1,j,k)*met(1,i+1,j,k)*met(1,i+1,j,k)*(
     *          c2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k)) +
     *          c1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))  )
     *      - mu(i-1,j,k)*met(1,i-1,j,k)*met(1,i-1,j,k)*(
     *          c2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k)) + 
     *          c1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))))

*** pr-derivatives
      r1 = r1 + c2*(
     *  (la(i,j,k+2))*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   ) 
     *   + mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2)) +
     *        c1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))  ) 
     *  - ((la(i,j,k-2))*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  ) 
     *     + mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2)) )  ) 
     *             ) + c1*(  
     *     (la(i,j,k+1))*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1)) ) 
     *  - (la(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1)) +
     *        c1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1)) ) ) )

*** rp derivatives
      r1 = r1 + c2*(
     *  (mu(i+2,j,k))*met(3,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   ) 
     *   + mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2)) +
     *        c1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))  ) 
     *  - (mu(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  )
     *     + mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))   ) )
     *             ) + c1*(  
     *     (mu(i+1,j,k))*met(3,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) ) 
     *     + mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1)) ) 
     *  - (mu(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2)) +
     *        c1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1)) )  ) )

*** qr derivatives
      r1 = r1 + c2*(
     *    mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2)) +
     *        c1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))   ) 
     *   + (2*mu(i,j,k+2)+la(i,j,k+2))*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  ) 
     *   +mu(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2)) +
     *        c1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))   ) 
     *  - ( mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))  ) 
     *    +(2*mu(i,j,k-2)+ la(i,j,k-2))*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   ) +
     *       mu(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))  ) ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1)) ) 
     *    + (2*mu(i,j,k+1)+la(i,j,k+1))*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )
     *    + mu(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1)) )   
     *  - ( mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1)) +
     *        c1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1)) ) 
     *    + (2*mu(i,j,k-1)+la(i,j,k-1))*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) )
     *    +  mu(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1)) +
     *        c1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1)) )  ) )

*** rq derivatives
      r1 = r1 + c2*(
     *    la(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2)) +
     *        c1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))   ) 
     *   +(2*mu(i,j+2,k)+la(i,j+2,k))*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  ) 
     *    + la(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2)) +
     *        c1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))   ) 
     *  - ( la(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))  ) 
     *     +(2*mu(i,j-2,k)+la(i,j-2,k))*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   ) 
     *    + la(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))  ) ) 
     *             ) + c1*(  
     *      la(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1)) ) 
     *     + (2*mu(i,j+1,k)+la(i,j+1,k))*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) ) 
     *     +la(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1)) )  
     *  - ( la(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2)) +
     *        c1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1)) ) 
     *     + (2*mu(i,j-1,k)+la(i,j-1,k))*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) )
     *     + la(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2)) +
     *        c1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1)) )  ) )

c          lu(2,i,j,k) = r1/jac(i,j,k)
          lu(2,i,j,k) = r1*ijac
*** w-equation

          r1 = 0
*** pp derivative (w)
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i+2,j,k)-u(3,i,j,k))  )

*** qq derivative (w)
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j+2,k)-u(3,i,j,k))  )
*** rr derivative (u)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(4,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(4,i,j,k-1)
          cof3=(mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(4,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(4,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(4,i,j,k+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j,k+2)-u(1,i,j,k))  )

*** rr derivative (v)
          cof1=(mu(i,j,k-2)+la(i,j,k-2))*met(3,i,j,k-2)*met(4,i,j,k-2)
          cof2=(mu(i,j,k-1)+la(i,j,k-1))*met(3,i,j,k-1)*met(4,i,j,k-1)
          cof3=(mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(4,i,j,k)
          cof4=(mu(i,j,k+1)+la(i,j,k+1))*met(3,i,j,k+1)*met(4,i,j,k+1)
          cof5=(mu(i,j,k+2)+la(i,j,k+2))*met(3,i,j,k+2)*met(4,i,j,k+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j,k+2)-u(2,i,j,k))  )

*** rr derivative (w)
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*met(4,i,j,k-2)**2 +
     *         mu(i,j,k-2)*(met(2,i,j,k-2)**2+met(3,i,j,k-2)**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*met(4,i,j,k-1)**2 +
     *         mu(i,j,k-1)*(met(2,i,j,k-1)**2+met(3,i,j,k-1)**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*met(4,i,j,k)**2 +
     *         mu(i,j,k)*(met(2,i,j,k)**2+met(3,i,j,k)**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*met(4,i,j,k+1)**2 +
     *         mu(i,j,k+1)*(met(2,i,j,k+1)**2+met(3,i,j,k+1)**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*met(4,i,j,k+2)**2 +
     *         mu(i,j,k+2)*(met(2,i,j,k+2)**2+met(3,i,j,k+2)**2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j,k+1)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )

*** pr-derivatives
c      r1 = r1 
     *
     *     + c2*(
     *  (la(i,j,k+2))*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   ) 
     *   + mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2)) +
     *        c1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))  ) 
     *  - ((la(i,j,k-2))*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  ) 
     *     + mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2)) )  ) 
     *             ) + c1*(  
     *     (la(i,j,k+1))*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1)) ) 
     *  - (la(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1)) +
     *        c1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1)) ) ) )

*** rp derivatives
c      r1 = r1 
     *
     *    + c2*(
     *  (mu(i+2,j,k))*met(4,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   ) 
     *   + mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2)) +
     *        c1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1))  ) 
     *  - (mu(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  )
     *     + mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1))   ) )
     *             ) + c1*(  
     *     (mu(i+1,j,k))*met(4,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) ) 
     *     + mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1)) ) 
     *  - (mu(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2)) +
     *        c1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1)) )  ) )

*** qr derivatives
     *  
c      r1 = r1
     *    + c2*(
     *    mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2)) +
     *        c1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))   ) 
     *   + la(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  ) 
     *  - ( mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))  ) 
     *     + la(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   ) ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1)) ) 
     *     + la(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )  
     *  - ( mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1)) +
     *        c1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1)) ) 
     *     + la(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) ) ) )

*** rq derivatives
     *
c      r1 = r1 
     *     + c2*(
     *    mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2)) +
     *        c1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))   ) 
     *   + mu(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  ) 
     *  - ( mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))  ) 
     *     + mu(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   ) ) 
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1)) ) 
     *     + mu(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) )  
     *  - ( mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2)) +
     *        c1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1)) ) 
     *     + mu(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) ) ) )


c          lu(3,i,j,k) = r1/jac(i,j,k)
          lu(3,i,j,k) = r1*ijac

      enddo
      enddo
      enddo
      end



