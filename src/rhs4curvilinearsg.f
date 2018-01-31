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
      subroutine FREESURFCURVISG( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, nz, side, u, mu, la, met, s, forcing, strx, stry ) 
     *     bind(c)
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
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 s(0:4), rhs1, rhs2, rhs3, s0i, ac, bc, cc, dc
      real*8 istry, istrx, xoysqrt, yoxsqrt, isqrtxy

      if( side.eq.5 )then
         k = 1
         kl= 1
      elseif( side.eq.6 )then
         k = nz
         kl= -1
      endif

      s0i = 1/s(0)
!$OMP PARALLEL PRIVATE(i,j,istry,istrx,rhs1,rhs2,rhs3,ac,bc,cc,dc,
!$OMP*                 xoysqrt,yoxsqrt,isqrtxy)
!$OMP DO      
      do j=jfirst+2,jlast-2
         istry = 1/stry(j)
         do i=ifirst+2,ilast-2
            istrx = 1/strx(i)

*** First tangential derivatives
            rhs1 = 
*** pr
     *   (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
     *          c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *          c1*(u(1,i+1,j,k)-u(1,i-1,j,k))  )*strx(i)*istry 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*istry   
*** qr
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*istrx*stry(j) 
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
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*strx(i)*istry 
*** qr
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
     * + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*stry(j)*istrx 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*istrx -
     *                  forcing(2,i,j)

*** (w-eq)
            rhs3 = 
*** pr
     *    la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*istry 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*strx(i)*istry
*** qr 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*stry(j)*istrx
     *  + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*istrx -
     *                  forcing(3,i,j)

*** Normal derivatives
            ac = strx(i)*istry*met(2,i,j,k)**2+
     *         stry(j)*istrx*met(3,i,j,k)**2+met(4,i,j,k)**2*istry*istrx
            bc = 1/(mu(i,j,k)*ac)
            cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac

            xoysqrt = SQRT(strx(i)*istry)
            yoxsqrt = 1/xoysqrt
            isqrtxy = istrx*xoysqrt
            dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
     *           yoxsqrt*met(3,i,j,k)*rhs2 + isqrtxy*met(4,i,j,k)*rhs3)

            u(1,i,j,k-kl) = -s0i*(  s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+
     *           s(3)*u(1,i,j,k+2*kl)+s(4)*u(1,i,j,k+3*kl) + bc*rhs1 - 
     *                                       dc*met(2,i,j,k)*xoysqrt )
            u(2,i,j,k-kl) = -s0i*(  s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+
     *           s(3)*u(2,i,j,k+2*kl)+s(4)*u(2,i,j,k+3*kl) + bc*rhs2 - 
     *                                       dc*met(3,i,j,k)*yoxsqrt )
            u(3,i,j,k-kl) = -s0i*(  s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+
     *           s(3)*u(3,i,j,k+2*kl)+s(4)*u(3,i,j,k+3*kl) + bc*rhs3 - 
     *                                       dc*met(4,i,j,k)*isqrtxy )
c            ac = strx(i)*istry*met(2,i,j,k)**2+
c     *                stry(j)*istrx*met(3,i,j,k)**2+met(4,i,j,k)**2
c            bc = 1/(mu(i,j,k)*ac)
c            cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac
c
c            xoysqrt = SQRT(strx(i)*istry)
c            yoxsqrt = 1/xoysqrt
c
c            dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
c     *                yoxsqrt*met(3,i,j,k)*rhs2 + met(4,i,j,k)*rhs3)
c
c            u(1,i,j,k-kl) = -s0i*(  s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+
c     *           s(3)*u(1,i,j,k+2*kl)+s(4)*u(1,i,j,k+3*kl) + bc*rhs1 - 
c     *                                       dc*met(2,i,j,k)*xoysqrt )
c            u(2,i,j,k-kl) = -s0i*(  s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+
c     *           s(3)*u(2,i,j,k+2*kl)+s(4)*u(2,i,j,k+3*kl) + bc*rhs2 - 
c     *                                       dc*met(3,i,j,k)*yoxsqrt )
c            u(3,i,j,k-kl) = -s0i*(  s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+
c     *           s(3)*u(3,i,j,k+2*kl)+s(4)*u(3,i,j,k+3*kl) + bc*rhs3 - 
c     *                                       dc*met(4,i,j,k) )
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine GETSURFFORCINGSG( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, met, jac, tau, strx, stry, forcing ) bind(c)
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
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 sqjac, istrx, istry
!$OMP PARALLEL PRIVATE(i,j,istry,istrx,sqjac)
!$OMP DO      
      do j=jfirst,jlast
         istry = 1/stry(j)
         do i=ifirst,ilast
            istrx = 1/strx(i)
            sqjac = SQRT(jac(i,j,k))
            forcing(1,i,j) =  sqjac*( istry*met(2,i,j,k)*tau(1,i,j)+
     *                                istrx*met(3,i,j,k)*tau(2,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(3,i,j) )
            forcing(2,i,j) =  sqjac*( istry*met(2,i,j,k)*tau(2,i,j)+
     *                                istrx*met(3,i,j,k)*tau(4,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(5,i,j) )
            forcing(3,i,j) =  sqjac*( istry*met(2,i,j,k)*tau(3,i,j)+
     *                                istrx*met(3,i,j,k)*tau(5,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(6,i,j) )
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine SUBSURFFORCINGSG( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, k, met, jac, tau, strx, stry, forcing ) bind(c)
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
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 sqjac, istrx, istry
!$OMP PARALLEL PRIVATE(i,j,istry,istrx,sqjac)
!$OMP DO      
      do j=jfirst,jlast
         istry = 1/stry(j)
         do i=ifirst,ilast
            istrx = 1/strx(i)
            sqjac = SQRT(jac(i,j,k))
            forcing(1,i,j) =  forcing(1,i,j) -
     *                        sqjac*( istry*met(2,i,j,k)*tau(1,i,j)+
     *                                istrx*met(3,i,j,k)*tau(2,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(3,i,j) )
            forcing(2,i,j) =  forcing(2,i,j) -
     *                        sqjac*( istry*met(2,i,j,k)*tau(2,i,j)+
     *                                istrx*met(3,i,j,k)*tau(4,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(5,i,j) )
            forcing(3,i,j) =  forcing(3,i,j) -
     *                        sqjac*( istry*met(2,i,j,k)*tau(3,i,j)+
     *                                istrx*met(3,i,j,k)*tau(5,i,j)+
     *                          istrx*istry*met(4,i,j,k)*tau(6,i,j) )
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end


