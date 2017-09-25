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
      subroutine GRIDINFO( ifirst, ilast, jfirst, jlast, kfirst, klast,
     *     met, jac, minj, maxj ) bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,k) 
!$OMP DO REDUCTION(max:maxj) REDUCTION(min:minj)
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine METRIC( ifirst, ilast, jfirst, jlast, kfirst, klast,
     *                   x, y, z, met, jac, ierr ) bind(c)
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
      ierr = 0
!$OMP PARALLEL PRIVATE(i,j,k,zr,zq,zp,sqzr) 
!$OMP DO REDUCTION(min:ierr)
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
                  ierr=-1
c                  return
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine METRICEXGH( ifirst, ilast, jfirst, jlast, kfirst, 
     *                       klast, nx, ny, nz, x, y, z, met, jac,
     *                       order, sb, zmax, amp, xc, yc, xl, yl ) 
     *   bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,k,l,s,sdb,tau,taup,tauq,p1,p2,zr,zq,zp,
!$OMP*                 zz,sqzr) 
!$OMP DO
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine FREESURFCURVI( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, nz, side, u, mu, la, met, s, forcing ) bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,rhs1,rhs2,rhs3,ac,bc,cc,dc)
!$OMP DO
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine GETSURFFORCING( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, met, jac, tau, forcing ) bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,sqjac)
!$OMP DO
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine GETSURFFORCINGGH( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, h, tau, forcing, amp, xc, yc, xl, yl ) bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,x,y,efact,zp,zq)
!$OMP DO
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
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine SUBSURFFORCING( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, met, jac, tau, forcing ) bind(c)
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
!$OMP PARALLEL PRIVATE(i,j,sqjac)
!$OMP DO
      do j=jfirst,jlast
         do i=ifirst,ilast
            sqjac = SQRT(jac(i,j,k))
            forcing(1,i,j) =  forcing(1,i,j) - 
     *         sqjac*( met(2,i,j,k)*tau(1,i,j)+ 
     *         met(3,i,j,k)*tau(2,i,j)+met(4,i,j,k)*tau(3,i,j) ) 
            forcing(2,i,j) =  forcing(2,i,j) - 
     *         sqjac*( met(2,i,j,k)*tau(2,i,j)+
     *         met(3,i,j,k)*tau(4,i,j)+met(4,i,j,k)*tau(5,i,j) )
            forcing(3,i,j) =  forcing(3,i,j) - 
     *          sqjac*( met(2,i,j,k)*tau(3,i,j)+
     *         met(3,i,j,k)*tau(5,i,j)+met(4,i,j,k)*tau(6,i,j) )
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine ADDBSTRESSC( ifirst, ilast, jfirst, jlast, kfirst, 
     *              klast, nz, u, mu, la, bs, met, side, s, op, ghterm,
     *              usesg, sgstrx, sgstry ) bind(c)
      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k, kl, nz, side, a1, a2, ghterm, usesg
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 bs(3,ifirst:ilast,jfirst:jlast)
      real*8 sgstrx(ifirst:ilast), sgstry(jfirst:jlast), sgx, sgy
      real*8 s(0:4), rhs1, rhs2, rhs3, ac, rtu, un1, vn1, wn1
      real*8 m2sg, m3sg, m4sg, isgx, isgy
      character*1 op

      if( side.eq.5 )then
         k = 1
         kl= 1
      elseif( side.eq.6 )then
         k = nz
         kl= -1
      endif
      if( op.eq.'=' )then
         a1 = 0
         a2 = 1
      elseif( op.eq.'+' )then
         a1 = 1
         a2 = 1
      elseif( op.eq.'-')then
         a1 = 1
         a2 =-1
      endif

!$OMP PARALLEL PRIVATE(i,j,sgx,sgy,isgy,isgx,rhs1,rhs2,rhs3,un1,vn1,
!$OMP*                 wn1,m2sg,m3sg,m4sg,rtu,ac)
      sgx = 1
      sgy = 1
      isgx = 1
      isgy = 1
!$OMP DO
      do j=jfirst+2,jlast-2
         do i=ifirst+2,ilast-2
            if( usesg.eq.1 )then
               sgx = sgstrx(i)
               sgy = sgstry(j)
               isgy = 1/sgy
               isgx = 1/sgx
            endif

*** First, tangential derivatives
            rhs1 = 
*** pr
     *      (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
     *          c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *          c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*sgx*isgy 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*isgy   
*** qr
     *  +    mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*isgx*sgy 
     *  + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )  
c     *               - forcing(1,i,j)

*** (v-eq)
            rhs2 = 
*** pr
     *     la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*sgx*isgy 
*** qr
     *  +    mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )
     * + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*sgy*isgx 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*isgx 
c     *               - forcing(2,i,j)

*** (w-eq)
            rhs3 = 
*** pr
     *      la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*isgy 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*sgx*isgy 
*** qr 
     *  +    mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*sgy*isgx 
     *  + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*isgx
c     *                - forcing(3,i,j)
      
*** then, normal derivatives
           un1 = s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+s(3)*u(1,i,j,k+2*kl)
     *          +s(4)*u(1,i,j,k+3*kl)
           vn1 = s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+s(3)*u(2,i,j,k+2*kl)
     *          +s(4)*u(2,i,j,k+3*kl)
           wn1 = s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+s(3)*u(3,i,j,k+2*kl)
     *          +s(4)*u(3,i,j,k+3*kl)
           if( ghterm .eq. 1 )then
              un1 = un1 + s(0)*u(1,i,j,k-kl)
              vn1 = vn1 + s(0)*u(2,i,j,k-kl)
              wn1 = wn1 + s(0)*u(3,i,j,k-kl)
           endif
           m2sg = SQRT(sgx*isgy)
           m3sg = 1/m2sg
           m4sg = isgx*m2sg

           rtu = un1*m2sg*met(2,i,j,k) + vn1*m3sg*met(3,i,j,k) +
     *           wn1*m4sg*met(4,i,j,k)
           ac  = sgx*isgy*met(2,i,j,k)**2 + sgy*isgx*met(3,i,j,k)**2 
     *                                   + isgx*isgy*met(4,i,j,k)**2
           rhs1 = rhs1 + (mu(i,j,k)+la(i,j,k))*rtu*m2sg*met(2,i,j,k) +
     *                                            mu(i,j,k)*ac*un1
           rhs2 = rhs2 + (mu(i,j,k)+la(i,j,k))*rtu*m3sg*met(3,i,j,k) +
     *                                            mu(i,j,k)*ac*vn1
           rhs3 = rhs3 + (mu(i,j,k)+la(i,j,k))*rtu*m4sg*met(4,i,j,k) +
     *                                            mu(i,j,k)*ac*wn1
           bs(1,i,j) = a1*bs(1,i,j) + a2*rhs1
           bs(2,i,j) = a1*bs(2,i,j) + a2*rhs2
           bs(3,i,j) = a1*bs(3,i,j) + a2*rhs3
        enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end
