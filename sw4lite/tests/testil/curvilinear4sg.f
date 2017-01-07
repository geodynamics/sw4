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
      subroutine CURVILINEAR4SG( ifirst, ilast, jfirst, jlast, kfirst,
     *                         klast, u, mu, la, met, jac, lu, 
     *                         onesided, acof, bope, ghcof, strx, stry,
     *                         op )


*** Routine with supergrid stretchings strx and stry. No stretching
*** in z, since top is always topography, and bottom always interface
*** to a deeper Cartesian grid.
*** opcount: 
***      Interior (k>6), 2126 arithmetic ops.
***      Boundary discretization (1<=k<=6 ), 6049 arithmetic ops.

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
      real*8 mux1, mux2, mux3, mux4, a1, sgn
      real*8 mucofu2, mucofuv, mucofuw, mucofv2, mucofvw, mucofw2
      real*8 ghcof(6), acof(6,8,8), bope(6,8)
      real*8 dudrp2, dudrp1, dudrm1, dudrm2
      real*8 dvdrp2, dvdrp1, dvdrm1, dvdrm2
      real*8 dwdrp2, dwdrp1, dwdrm1, dwdrm2
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 istrx, istry, istrxy
      character*1 op

*** met(1) is sqrt(J)*px = sqrt(J)*qy
*** met(2) is sqrt(J)*rx
*** met(3) is sqrt(J)*ry
*** met(4) is sqrt(J)*rz

      if( op.eq.'=' )then
         a1 = 0
         sgn= 1
      elseif( op.eq.'+')then
         a1 = 1
         sgn= 1
      elseif( op.eq.'-')then
         a1 = 1
         sgn=-1
      endif
      kstart = kfirst+2
!$OMP PARALLEL PRIVATE(k,i,j,q,m,mux1,mux2,mux3,mux4,r1,r2,r3,ijac,
!$OMP*   istry,istrx,istrxy,cof1,cof2,cof3,cof4,cof5,mucofu2,mucofuv,
!$OMP*   mucofuw,mucofvw,mucofv2,mucofw2,dudrm2,dudrm1,dudrp1,dudrp2,
!$OMP*       dvdrm2,dvdrm1,dvdrp1,dvdrp2,dwdrm2,dwdrm1,dwdrp1,dwdrp2)
      if( onesided(5).eq.1 )then
         kstart = 7
*** SBP Boundary closure terms
!$OMP DO
         do k=1,6
            do j=jfirst+2,jlast-2
               do i=ifirst+2,ilast-2
*** 5 ops                  
               ijac   = strx(i)*stry(j)/jac(i,j,k)
               istry  = 1/(stry(j))
               istrx  = 1/(strx(i))
               istrxy = istry*istrx

               r1 = 0
               r2 = 0
               r3 = 0

*** pp derivative (u) (u-eq)
*** 53 ops, tot=58
          cof1=(2*mu(i-2,j,k)+la(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
     *                                              *strx(i-2)
          cof2=(2*mu(i-1,j,k)+la(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
     *                                              *strx(i-1)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*strx(i)
          cof4=(2*mu(i+1,j,k)+la(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
     *                                              *strx(i+1)
          cof5=(2*mu(i+2,j,k)+la(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)
     *                                              *strx(i+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i+1,j,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i+2,j,k)-u(1,i,j,k))  )*istry


*** qq derivative (u) (u-eq)
*** 43 ops, tot=101
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)*stry(j-2)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)*stry(j-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*stry(j)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)*stry(j+1)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)*stry(j+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j+2,k)-u(1,i,j,k))  )*istrx

*** pp derivative (v) (v-eq)
*** 43 ops, tot=144
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)*strx(i-2)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)*strx(i-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*strx(i)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)*strx(i+1)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)*strx(i+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r2 = r2 + i6* (
     *               mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i+2,j,k)-u(2,i,j,k))  )*istry
*** qq derivative (v) (v-eq)
*** 53 ops, tot=197
          cof1=(2*mu(i,j-2,k)+la(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
     *                                                     *stry(j-2)
          cof2=(2*mu(i,j-1,k)+la(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
     *                                                     *stry(j-1)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
     *                                                     *stry(j)
          cof4=(2*mu(i,j+1,k)+la(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
     *                                                     *stry(j+1)
          cof5=(2*mu(i,j+2,k)+la(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
     *                                                     *stry(j+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r2 = r2 + i6* (
     *               mux1*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j+1,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j+2,k)-u(2,i,j,k))  )*istrx


*** pp derivative (w) (w-eq)
*** 43 ops, tot=240
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)*strx(i-2)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)*strx(i-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*strx(i)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)*strx(i+1)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)*strx(i+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r3 = r3 + i6* (
     *               mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i+2,j,k)-u(3,i,j,k))  )*istry

*** qq derivative (w) (w-eq)
*** 43 ops, tot=283
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)*stry(j-2)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)*stry(j-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*stry(j)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)*stry(j+1)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)*stry(j+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r3 = r3 + i6* (
     *               mux1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j+2,k)-u(3,i,j,k))  )*istrx


*** All rr-derivatives at once
*** averaging the coefficient
*** 54*8*8+25*8 = 3656 ops, tot=3939
            do q=1,8
               mucofu2=0
               mucofuv=0
               mucofuw=0
               mucofvw=0
               mucofv2=0
               mucofw2=0
               do m=1,8
                  mucofu2 = mucofu2 +
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*
     *                                  (met(2,i,j,m)*strx(i))**2
     *                           + mu(i,j,m)*
     *                   ((met(3,i,j,m)*stry(j))**2+met(4,i,j,m)**2))
                  mucofv2 = mucofv2+
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*
     *                                  (met(3,i,j,m)*stry(j))**2
     *                              + mu(i,j,m)*
     *                    ((met(2,i,j,m)*strx(i))**2+met(4,i,j,m)**2))
                  mucofw2 = mucofw2+
     *             acof(k,q,m)*((2*mu(i,j,m)+la(i,j,m))*met(4,i,j,m)**2
     *                       + mu(i,j,m)*
     *           ((met(2,i,j,m)*strx(i))**2+(met(3,i,j,m)*stry(j))**2))
                  mucofuv = mucofuv+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(2,i,j,m)*met(3,i,j,m)
                  mucofuw = mucofuw+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(2,i,j,m)*met(4,i,j,m)
                  mucofvw = mucofvw+acof(k,q,m)*(mu(i,j,m)+la(i,j,m))*
     *                 met(3,i,j,m)*met(4,i,j,m)
              enddo
*** Computing the second derivative,
              r1 = r1 + istrxy*mucofu2*u(1,i,j,q) + mucofuv*u(2,i,j,q) + 
     *                                         istry*mucofuw*u(3,i,j,q)
              r2 = r2 + mucofuv*u(1,i,j,q) + istrxy*mucofv2*u(2,i,j,q) + 
     *                                         istrx*mucofvw*u(3,i,j,q)
              r3 = r3 + istry*mucofuw*u(1,i,j,q) + 
     *              istrx*mucofvw*u(2,i,j,q) + istrxy*mucofw2*u(3,i,j,q)
            end do

*** Ghost point values, only nonzero for k=1.
*** 72 ops., tot=4011
            mucofu2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*
     *                           (met(2,i,j,1)*strx(i))**2
     *        + mu(i,j,1)*((met(3,i,j,1)*stry(j))**2+met(4,i,j,1)**2))
            mucofv2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*
     *                           (met(3,i,j,1)*stry(j))**2
     *         + mu(i,j,1)*((met(2,i,j,1)*strx(i))**2+met(4,i,j,1)**2))
            mucofw2 = ghcof(k)*((2*mu(i,j,1)+la(i,j,1))*met(4,i,j,1)**2
     *                             + mu(i,j,1)*
     *           ((met(2,i,j,1)*strx(i))**2+(met(3,i,j,1)*stry(j))**2))
            mucofuv = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(2,i,j,1)*met(3,i,j,1)
            mucofuw = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(2,i,j,1)*met(4,i,j,1)
            mucofvw = ghcof(k)*(mu(i,j,1)+la(i,j,1))*
     *                 met(3,i,j,1)*met(4,i,j,1)
            r1 = r1 + istrxy*mucofu2*u(1,i,j,0) + mucofuv*u(2,i,j,0) + 
     *                                        istry*mucofuw*u(3,i,j,0)
            r2 = r2 + mucofuv*u(1,i,j,0) + istrxy*mucofv2*u(2,i,j,0) + 
     *                                        istrx*mucofvw*u(3,i,j,0)
            r3 = r3 + istry*mucofuw*u(1,i,j,0) + 
     *            istrx*mucofvw*u(2,i,j,0) + istrxy*mucofw2*u(3,i,j,0)

*** pq-derivatives (u-eq)
*** 38 ops., tot=4049
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
*** 38 ops. tot=4087
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
*** 38 ops. , tot=4125
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
*** 38 ops., tot=4163
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
*** 24*8 = 192 ops, tot=4355
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
*** 67 ops, tot=4422
      r1 = r1 + ( c2*(
     *  (2*mu(i+2,j,k)+la(i+2,j,k))*met(2,i+2,j,k)*met(1,i+2,j,k)*
     *                                                 strx(i+2)*dudrp2
     *   + la(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*dvdrp2*stry(j)
     *   + la(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*dwdrp2
     *-((2*mu(i-2,j,k)+la(i-2,j,k))*met(2,i-2,j,k)*met(1,i-2,j,k)*
     *                                                 strx(i-2)*dudrm2
     *   + la(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*dvdrm2*stry(j)
     *   + la(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*dwdrm2 )
     *              ) + c1*(  
     *  (2*mu(i+1,j,k)+la(i+1,j,k))*met(2,i+1,j,k)*met(1,i+1,j,k)*
     *                                                 strx(i+1)*dudrp1
     *   + la(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*dvdrp1*stry(j)
     *   + la(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*dwdrp1 
     *-((2*mu(i-1,j,k)+la(i-1,j,k))*met(2,i-1,j,k)*met(1,i-1,j,k)*
     *                                                 strx(i-1)*dudrm1
     *   + la(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*dvdrm1*stry(j)
     *   + la(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*dwdrm1 ) ) )*istry

*** rp derivatives (v-eq)
*** 42 ops, tot=4464
      r2 = r2 + c2*(
     *     mu(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*dudrp2
     *  +  mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*dvdrp2*
     *                                           strx(i+2)*istry
     *  - (mu(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*dudrm2
     *  +  mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*dvdrm2*
     *                                           strx(i-2)*istry )
     *             ) + c1*(  
     *     mu(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*dudrp1
     *  +  mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*dvdrp1*
     *                                           strx(i+1)*istry
     *  - (mu(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*dudrm1
     *  +  mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*dvdrm1*
     *                                           strx(i-1)*istry )
     *                     )

*** rp derivatives (w-eq)
*** 38 ops, tot=4502
      r3 = r3 + istry*(c2*(
     *     mu(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*dudrp2
     *  +  mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*dwdrp2*strx(i+2)
     *  - (mu(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*dudrm2
     *  +  mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*dwdrm2*strx(i-2))
     *             ) + c1*(  
     *     mu(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*dudrp1
     *  +  mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*dwdrp1*strx(i+1)
     *  - (mu(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*dudrm1
     *  +  mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*dwdrm1*strx(i-1))
     *                     ) )

*** rq - derivatives
*** 24*8 = 192 ops , tot=4694
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
*** 42 ops, tot=4736
      r1 = r1 + c2*(
     *      mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*dudrp2*
     *                                           stry(j+2)*istrx
     *   +  mu(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *   - (mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*dudrm2*
     *                                           stry(j-2)*istrx
     *   +  mu(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*dvdrm2)
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*dudrp1*
     *                                           stry(j+1)*istrx
     *   +  mu(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *   - (mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*dudrm1*
     *                                           stry(j-1)*istrx
     *   +  mu(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*dvdrm1)
     *         )

*** rq derivatives (v-eq)
*** 70 ops, tot=4806
      r2 = r2 + c2*(
     *      la(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*dudrp2
     * +(2*mu(i,j+2,k)+la(i,j+2,k))*met(3,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *                                              *stry(j+2)*istrx
     *    + la(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*dwdrp2*istrx
     *  - ( la(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*dudrm2
     * +(2*mu(i,j-2,k)+la(i,j-2,k))*met(3,i,j-2,k)*met(1,i,j-2,k)*dvdrm2
     *                                              *stry(j-2)*istrx
     *    + la(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*dwdrm2*istrx )
     *             ) + c1*(  
     *      la(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*dudrp1
     * +(2*mu(i,j+1,k)+la(i,j+1,k))*met(3,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *                                              *stry(j+1)*istrx
     *    + la(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*dwdrp1*istrx
     *  - ( la(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*dudrm1
     * +(2*mu(i,j-1,k)+la(i,j-1,k))*met(3,i,j-1,k)*met(1,i,j-1,k)*dvdrm1
     *                                              *stry(j-1)*istrx
     *    + la(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*dwdrm1*istrx )
     *        )


*** rq derivatives (w-eq)
*** 39 ops, tot=4845
      r3 = r3 + ( c2*(
     *     mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*dwdrp2*stry(j+2)
     *  +  mu(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*dvdrp2
     *  - (mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*dwdrm2*stry(j-2)
     *  +  mu(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*dvdrm2)
     *             ) + c1*(  
     *     mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*dwdrp1*stry(j+1)
     *  +  mu(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*dvdrp1
     *  - (mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*dwdrm1*stry(j-1)
     *  +  mu(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*dvdrm1)
     *               ) )*istrx

*** pr and qr derivatives at once
*** in loop: 8*(53+53+43) = 1192 ops, tot=6037
      do q=1,8
*** (u-eq)
*** 53 ops
        r1 = r1 + bope(k,q)*( 
*** pr
     *   (2*mu(i,j,q)+la(i,j,q))*met(2,i,j,q)*met(1,i,j,q)*(
     *          c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *          c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   )*strx(i)*istry
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i+2,j,q)-u(2,i-2,j,q)) +
     *        c1*(u(2,i+1,j,q)-u(2,i-1,j,q))  ) 
     *  + mu(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i+2,j,q)-u(3,i-2,j,q)) +
     *        c1*(u(3,i+1,j,q)-u(3,i-1,j,q))  )*istry 
*** qr
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i,j+2,q)-u(1,i,j-2,q)) +
     *        c1*(u(1,i,j+1,q)-u(1,i,j-1,q))   )*stry(j)*istrx
     *  + la(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  ) )

*** (v-eq)
*** 53 ops
        r2 = r2 + bope(k,q)*(
*** pr
     *    la(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *        c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   ) 
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i+2,j,q)-u(2,i-2,j,q)) +
     *        c1*(u(2,i+1,j,q)-u(2,i-1,j,q))  )*strx(i)*istry 
*** qr
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i,j+2,q)-u(1,i,j-2,q)) +
     *        c1*(u(1,i,j+1,q)-u(1,i,j-1,q))   ) 
     * + (2*mu(i,j,q)+la(i,j,q))*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  )*stry(j)*istrx 
     *  + mu(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i,j+2,q)-u(3,i,j-2,q)) +
     *        c1*(u(3,i,j+1,q)-u(3,i,j-1,q))   )*istrx  )

*** (w-eq)
*** 43 ops
        r3 = r3 + bope(k,q)*(
*** pr
     *    la(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(1,i+2,j,q)-u(1,i-2,j,q)) +
     *        c1*(u(1,i+1,j,q)-u(1,i-1,j,q))   )*istry 
     *  + mu(i,j,q)*met(2,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i+2,j,q)-u(3,i-2,j,q)) +
     *        c1*(u(3,i+1,j,q)-u(3,i-1,j,q))  )*strx(i)*istry
*** qr 
     *  + mu(i,j,q)*met(3,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(3,i,j+2,q)-u(3,i,j-2,q)) +
     *        c1*(u(3,i,j+1,q)-u(3,i,j-1,q))   )*stry(j)*istrx 
     *  + la(i,j,q)*met(4,i,j,q)*met(1,i,j,q)*(
     *        c2*(u(2,i,j+2,q)-u(2,i,j-2,q)) +
     *        c1*(u(2,i,j+1,q)-u(2,i,j-1,q))  )*istrx )

      enddo

*** 12 ops, tot=6049
c          lu(1,i,j,k) = r1*ijac
          lu(1,i,j,k) = a1*lu(1,i,j,k) + sgn*r1*ijac
c          lu(2,i,j,k) = r2*ijac
          lu(2,i,j,k) = a1*lu(2,i,j,k) + sgn*r2*ijac
c          lu(3,i,j,k) = r3*ijac
          lu(3,i,j,k) = a1*lu(3,i,j,k) + sgn*r3*ijac
               enddo
            enddo
         enddo
!$OMP END DO
      endif

!$OMP DO
      do k=kstart,klast-2
         do j=jfirst+2,jlast-2
            do i=ifirst+2,ilast-2
*** 5 ops
               ijac = strx(i)*stry(j)/jac(i,j,k)
               istry = 1/(stry(j))
               istrx = 1/(strx(i))
               istrxy = istry*istrx

               r1 = 0

*** pp derivative (u)
*** 53 ops, tot=58
          cof1=(2*mu(i-2,j,k)+la(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)
     *                                              *strx(i-2)
          cof2=(2*mu(i-1,j,k)+la(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)
     *                                              *strx(i-1)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
     *                                              *strx(i)
          cof4=(2*mu(i+1,j,k)+la(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)
     *                                              *strx(i+1)
          cof5=(2*mu(i+2,j,k)+la(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)
     *                                              *strx(i+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i+1,j,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i+2,j,k)-u(1,i,j,k))  )*istry
*** qq derivative (u)
*** 43 ops, tot=101
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)*stry(j-2)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)*stry(j-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*stry(j)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)*stry(j+1)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)*stry(j+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j+2,k)-u(1,i,j,k))  )*istrx
*** rr derivative (u)
*** 5*11+14+14=83 ops, tot=184
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*(met(2,i,j,k-2)*strx(i))**2
     *   +   mu(i,j,k-2)*((met(3,i,j,k-2)*stry(j))**2+met(4,i,j,k-2)**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*(met(2,i,j,k-1)*strx(i))**2
     *   +   mu(i,j,k-1)*((met(3,i,j,k-1)*stry(j))**2+met(4,i,j,k-1)**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*(met(2,i,j,k)*strx(i))**2 +
     *         mu(i,j,k)*((met(3,i,j,k)*stry(j))**2+met(4,i,j,k)**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*(met(2,i,j,k+1)*strx(i))**2
     *   +   mu(i,j,k+1)*((met(3,i,j,k+1)*stry(j))**2+met(4,i,j,k+1)**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*(met(2,i,j,k+2)*strx(i))**2
     *   +   mu(i,j,k+2)*((met(3,i,j,k+2)*stry(j))**2+met(4,i,j,k+2)**2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *               mux2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *               mux3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *               mux4*(u(1,i,j,k+2)-u(1,i,j,k))  )*istrxy

*** rr derivative (v)
*** 42 ops, tot=226
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
*** 43 ops, tot=269
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
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )*istry

*** pq-derivatives
*** 38 ops, tot=307
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
*** 38 ops, tot=345
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
*** 130 ops., tot=475
      r1 = r1 + c2*(
     *  (2*mu(i,j,k+2)+la(i,j,k+2))*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   )*strx(i)*istry 
     *   + mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2)) +
     *        c1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))  ) 
     *   + mu(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2)) +
     *        c1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))  )*istry
     *  - ((2*mu(i,j,k-2)+la(i,j,k-2))*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  )*strx(i)*istry  
     *     + mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2))   ) 
     *     + mu(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2))   )*istry )
     *             ) + c1*(  
     *     (2*mu(i,j,k+1)+la(i,j,k+1))*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) )*strx(i)*istry 
     *     + mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1))  )*istry
     *  - ((2*mu(i,j,k-1)+la(i,j,k-1))*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) )*strx(i)*istry  
     *     + mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1)) +
     *        c1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1)) +
     *        c1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1))   )*istry  ) )

*** rp derivatives
*** 130 ops, tot=605
      r1 = r1 +( c2*(
     *  (2*mu(i+2,j,k)+la(i+2,j,k))*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   )*strx(i+2) 
     *   + la(i+2,j,k)*met(3,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2)) +
     *        c1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))  )*stry(j) 
     *   + la(i+2,j,k)*met(4,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2)) +
     *        c1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1))  )
     *  - ((2*mu(i-2,j,k)+la(i-2,j,k))*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  )*strx(i-2) 
     *     + la(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))   )*stry(j) 
     *     + la(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1))   ) )
     *             ) + c1*(  
     *     (2*mu(i+1,j,k)+la(i+1,j,k))*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) )*strx(i+1) 
     *     + la(i+1,j,k)*met(3,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1)) )*stry(j) 
     *     + la(i+1,j,k)*met(4,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1))  )
     *  - ((2*mu(i-1,j,k)+la(i-1,j,k))*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) )*strx(i-1) 
     *     + la(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2)) +
     *        c1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1)) )*stry(j) 
     *     + la(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2)) +
     *        c1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1))   )  ) ) )*istry

*** qr derivatives
*** 82 ops, tot=687
      r1 = r1 + c2*(
     *    mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2)) +
     *        c1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))   )*stry(j)*istrx 
     *   + la(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  ) 
     *  - ( mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))  )*stry(j)*istrx  
     *     + la(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   ) ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1)) )*stry(j)*istrx  
     *     + la(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )  
     *  - ( mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1)) +
     *        c1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1)) )*stry(j)*istrx  
     *     + la(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) ) ) )

*** rq derivatives
*** 82 ops, tot=769
      r1 = r1 + c2*(
     *    mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2)) +
     *        c1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))   )*stry(j+2)*istrx 
     *   + mu(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  ) 
     *  - ( mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))  )*stry(j-2)*istrx  
     *     + mu(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   ) ) 
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1)) )*stry(j+1)*istrx
     *     + mu(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) )
     *  - ( mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2)) +
     *        c1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1)) )*stry(j-1)*istrx    
     *     + mu(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) ) ) )

c          lu(1,i,j,k) = r1/jac(i,j,k)
c          lu(1,i,j,k) = r1*ijac
*** 4 ops, tot=773
          lu(1,i,j,k) = a1*lu(1,i,j,k) + sgn*r1*ijac
*** v-equation

          r1 = 0
*** pp derivative (v)
*** 43 ops, tot=816
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)*strx(i-2)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)*strx(i-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*strx(i)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)*strx(i+1)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)*strx(i+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i+2,j,k)-u(2,i,j,k))  )*istry

*** qq derivative (v)
*** 53 ops, tot=869
          cof1=(2*mu(i,j-2,k)+la(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)
     *                                                     *stry(j-2)
          cof2=(2*mu(i,j-1,k)+la(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)
     *                                                     *stry(j-1)
          cof3=(2*mu(i,j,k)+la(i,j,k))*met(1,i,j,k)*met(1,i,j,k)
     *                                                     *stry(j)
          cof4=(2*mu(i,j+1,k)+la(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)
     *                                                     *stry(j+1)
          cof5=(2*mu(i,j+2,k)+la(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)
     *                                                     *stry(j+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j+1,k)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j+2,k)-u(2,i,j,k))  )*istrx

*** rr derivative (u)
*** 42 ops, tot=911
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
*** 83 ops, tot=994
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*(met(3,i,j,k-2)*stry(j))**2
     *  +    mu(i,j,k-2)*((met(2,i,j,k-2)*strx(i))**2+met(4,i,j,k-2)**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*(met(3,i,j,k-1)*stry(j))**2
     *  +    mu(i,j,k-1)*((met(2,i,j,k-1)*strx(i))**2+met(4,i,j,k-1)**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*(met(3,i,j,k)*stry(j))**2 +
     *       mu(i,j,k)*((met(2,i,j,k)*strx(i))**2+met(4,i,j,k)**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*(met(3,i,j,k+1)*stry(j))**2
     *  +    mu(i,j,k+1)*((met(2,i,j,k+1)*strx(i))**2+met(4,i,j,k+1)**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*(met(3,i,j,k+2)*stry(j))**2
     *  +    mu(i,j,k+2)*((met(2,i,j,k+2)*strx(i))**2+met(4,i,j,k+2)**2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *               mux2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *               mux3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *               mux4*(u(2,i,j,k+2)-u(2,i,j,k))  )*istrxy

*** rr derivative (w)
*** 43 ops, tot=1037
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
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )*istrx

*** pq-derivatives
*** 38 ops, tot=1075
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
*** 38 ops, tot=1113
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
*** 82 ops, tot=1195
      r1 = r1 + c2*(
     *  (la(i,j,k+2))*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   ) 
     *   + mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2)) +
     *        c1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))  )*strx(i)*istry 
     *  - ((la(i,j,k-2))*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  ) 
     *     + mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2)) )*strx(i)*istry ) 
     *             ) + c1*(  
     *     (la(i,j,k+1))*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) ) 
     *     + mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1)) )*strx(i)*istry  
     *  - (la(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1)) +
     *        c1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1)) )*strx(i)*istry  ) )

*** rp derivatives
*** 82 ops, tot=1277
      r1 = r1 + c2*(
     *  (mu(i+2,j,k))*met(3,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   ) 
     *   + mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2)) +
     *        c1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))  )*strx(i+2)*istry 
     *  - (mu(i-2,j,k)*met(3,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  )
     *     + mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2)) +
     *        c1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))   )*strx(i-2)*istry )
     *             ) + c1*(  
     *     (mu(i+1,j,k))*met(3,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) ) 
     *     + mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2)) +
     *        c1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1)) )*strx(i+1)*istry 
     *  - (mu(i-1,j,k)*met(3,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2)) +
     *        c1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1)) )*strx(i-1)*istry  ) )

*** qr derivatives
*** 130 ops, tot=1407
      r1 = r1 + c2*(
     *    mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2)) +
     *        c1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))   ) 
     *   + (2*mu(i,j,k+2)+la(i,j,k+2))*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  )*stry(j)*istrx 
     *   +mu(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2)) +
     *        c1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))   )*istrx 
     *  - ( mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))  ) 
     *    +(2*mu(i,j,k-2)+ la(i,j,k-2))*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   )*stry(j)*istrx +
     *       mu(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))  )*istrx ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1)) ) 
     *    + (2*mu(i,j,k+1)+la(i,j,k+1))*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )*stry(j)*istrx
     *    + mu(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1)) )*istrx   
     *  - ( mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1)) +
     *        c1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1)) ) 
     *    + (2*mu(i,j,k-1)+la(i,j,k-1))*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) )*stry(j)*istrx
     *    +  mu(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1)) +
     *        c1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1)) )*istrx  ) )


*** rq derivatives
*** 130 ops, tot=1537
      r1 = r1 + c2*(
     *    la(i,j+2,k)*met(2,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2)) +
     *        c1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))   ) 
     *   +(2*mu(i,j+2,k)+la(i,j+2,k))*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  )*stry(j+2)*istrx 
     *    + la(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2)) +
     *        c1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))   )*istrx 
     *  - ( la(i,j-2,k)*met(2,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2)) +
     *        c1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))  ) 
     *     +(2*mu(i,j-2,k)+la(i,j-2,k))*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   )*stry(j-2)*istrx 
     *    + la(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))  )*istrx  ) 
     *             ) + c1*(  
     *      la(i,j+1,k)*met(2,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2)) +
     *        c1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1)) ) 
     *     + (2*mu(i,j+1,k)+la(i,j+1,k))*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) )*stry(j+1)*istrx 
     *     +la(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1)) )*istrx   
     *  - ( la(i,j-1,k)*met(2,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2)) +
     *        c1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1)) ) 
     *     + (2*mu(i,j-1,k)+la(i,j-1,k))*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) )*stry(j-1)*istrx
     *     + la(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2)) +
     *        c1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1)) )*istrx   ) )


c          lu(2,i,j,k) = r1/jac(i,j,k)
c          lu(2,i,j,k) = r1*ijac
*** 4 ops, tot=1541
          lu(2,i,j,k) = a1*lu(2,i,j,k) + sgn*r1*ijac
*** w-equation

          r1 = 0
*** pp derivative (w)
*** 43 ops, tot=1580
          cof1=(mu(i-2,j,k))*met(1,i-2,j,k)*met(1,i-2,j,k)*strx(i-2)
          cof2=(mu(i-1,j,k))*met(1,i-1,j,k)*met(1,i-1,j,k)*strx(i-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*strx(i)
          cof4=(mu(i+1,j,k))*met(1,i+1,j,k)*met(1,i+1,j,k)*strx(i+1)
          cof5=(mu(i+2,j,k))*met(1,i+2,j,k)*met(1,i+2,j,k)*strx(i+2)

            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i+2,j,k)-u(3,i,j,k))  )*istry

*** qq derivative (w)
*** 43 ops, tot=1623
          cof1=(mu(i,j-2,k))*met(1,i,j-2,k)*met(1,i,j-2,k)*stry(j-2)
          cof2=(mu(i,j-1,k))*met(1,i,j-1,k)*met(1,i,j-1,k)*stry(j-1)
          cof3=(mu(i,j,k))*met(1,i,j,k)*met(1,i,j,k)*stry(j)
          cof4=(mu(i,j+1,k))*met(1,i,j+1,k)*met(1,i,j+1,k)*stry(j+1)
          cof5=(mu(i,j+2,k))*met(1,i,j+2,k)*met(1,i,j+2,k)*stry(j+2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j+2,k)-u(3,i,j,k))  )*istrx

*** rr derivative (u)
*** 43 ops, tot=1666
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
     *               mux4*(u(1,i,j,k+2)-u(1,i,j,k))  )*istry

*** rr derivative (v)
*** 43 ops, tot=1709
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
     *               mux4*(u(2,i,j,k+2)-u(2,i,j,k))  )*istrx

*** rr derivative (w)
*** 83 ops, tot=1792
          cof1 = (2*mu(i,j,k-2)+la(i,j,k-2))*met(4,i,j,k-2)**2 +
     *         mu(i,j,k-2)*((met(2,i,j,k-2)*strx(i))**2+
     *                                 (met(3,i,j,k-2)*stry(j))**2)
          cof2 = (2*mu(i,j,k-1)+la(i,j,k-1))*met(4,i,j,k-1)**2 +
     *         mu(i,j,k-1)*((met(2,i,j,k-1)*strx(i))**2+
     *                                 (met(3,i,j,k-1)*stry(j))**2)
          cof3 = (2*mu(i,j,k)+la(i,j,k))*met(4,i,j,k)**2 +
     *         mu(i,j,k)*((met(2,i,j,k)*strx(i))**2+
     *                                 (met(3,i,j,k)*stry(j))**2)
          cof4 = (2*mu(i,j,k+1)+la(i,j,k+1))*met(4,i,j,k+1)**2 +
     *         mu(i,j,k+1)*((met(2,i,j,k+1)*strx(i))**2+
     *                                 (met(3,i,j,k+1)*stry(j))**2)
          cof5 = (2*mu(i,j,k+2)+la(i,j,k+2))*met(4,i,j,k+2)**2 +
     *         mu(i,j,k+2)*((met(2,i,j,k+2)*strx(i))**2+
     *                                 (met(3,i,j,k+2)*stry(j))**2)
            mux1 = cof2 -tf*(cof3+cof1)
            mux2 = cof1 + cof4+3*(cof3+cof2)
            mux3 = cof2 + cof5+3*(cof4+cof3)
            mux4 = cof4-tf*(cof3+cof5)

            r1 = r1 + i6* (
     *               mux1*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     *               mux2*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     *               mux3*(u(3,i,j,k+1)-u(3,i,j,k)) +
     *               mux4*(u(3,i,j,k+2)-u(3,i,j,k))  )*istrxy
*** pr-derivatives
*** 86 ops, tot=1878
c      r1 = r1 
     *
     *     + c2*(
     *  (la(i,j,k+2))*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2)) +
     *        c1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))   )*istry 
     *   + mu(i,j,k+2)*met(2,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2)) +
     *        c1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))  )*strx(i)*istry 
     *  - ((la(i,j,k-2))*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))  )*istry  
     *     + mu(i,j,k-2)*met(2,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2)) )*strx(i)*istry ) 
     *             ) + c1*(  
     *     (la(i,j,k+1))*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1)) )*istry  
     *     + mu(i,j,k+1)*met(2,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1)) )*strx(i)*istry  
     *  - (la(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1)) +
     *        c1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1)) )*istry  
     *     + mu(i,j,k-1)*met(2,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1)) +
     *        c1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1)) )*strx(i)*istry  ) )
*** rp derivatives
*** 79 ops, tot=1957
c      r1 = r1 
     *
     *    + istry*(c2*(
     *  (mu(i+2,j,k))*met(4,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2)) +
     *        c1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))   ) 
     *   + mu(i+2,j,k)*met(2,i+2,j,k)*met(1,i+2,j,k)*(
     *        c2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2)) +
     *        c1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1)) )*strx(i+2) 
     *  - (mu(i-2,j,k)*met(4,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2)) +
     *        c1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))  )
     *     + mu(i-2,j,k)*met(2,i-2,j,k)*met(1,i-2,j,k)*(
     *        c2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2)) +
     *        c1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1)) )*strx(i-2)  )
     *             ) + c1*(  
     *     (mu(i+1,j,k))*met(4,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2)) +
     *        c1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1)) ) 
     *     + mu(i+1,j,k)*met(2,i+1,j,k)*met(1,i+1,j,k)*(
     *        c2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2)) +
     *        c1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1)) )*strx(i+1)  
     *  - (mu(i-1,j,k)*met(4,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2)) +
     *        c1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1)) ) 
     *     + mu(i-1,j,k)*met(2,i-1,j,k)*met(1,i-1,j,k)*(
     *        c2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2)) +
     *        c1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1)) )*strx(i-1)  ) ) )
*** qr derivatives
*** 86 ops, tot=2043
     *  
c      r1 = r1
     *    + c2*(
     *    mu(i,j,k+2)*met(3,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2)) +
     *        c1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))   )*stry(j)*istrx 
     *   + la(i,j,k+2)*met(4,i,j,k+2)*met(1,i,j,k+2)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2)) +
     *        c1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))  )*istrx 
     *  - ( mu(i,j,k-2)*met(3,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))  )*stry(j)*istrx  
     *     + la(i,j,k-2)*met(4,i,j,k-2)*met(1,i,j,k-2)*(
     *        c2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))   )*istrx  ) 
     *             ) + c1*(  
     *      mu(i,j,k+1)*met(3,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1)) )*stry(j)*istrx  
     *     + la(i,j,k+1)*met(4,i,j,k+1)*met(1,i,j,k+1)*(
     *        c2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1)) )*istrx   
     *  - ( mu(i,j,k-1)*met(3,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1)) +
     *        c1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1)) )*stry(j)*istrx  
     *     + la(i,j,k-1)*met(4,i,j,k-1)*met(1,i,j,k-1)*(
     *        c2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1)) +
     *        c1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1)) )*istrx  ) )
*** rq derivatives
***  79 ops, tot=2122
     *
c      r1 = r1 
     *     + istrx*(c2*(
     *    mu(i,j+2,k)*met(3,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2)) +
     *        c1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))   )*stry(j+2) 
     *   + mu(i,j+2,k)*met(4,i,j+2,k)*met(1,i,j+2,k)*(
     *        c2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2)) +
     *        c1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))  ) 
     *  - ( mu(i,j-2,k)*met(3,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2)) +
     *        c1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))  )*stry(j-2) 
     *     + mu(i,j-2,k)*met(4,i,j-2,k)*met(1,i,j-2,k)*(
     *        c2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2)) +
     *        c1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))   ) ) 
     *             ) + c1*(  
     *      mu(i,j+1,k)*met(3,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2)) +
     *        c1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1)) )*stry(j+1) 
     *     + mu(i,j+1,k)*met(4,i,j+1,k)*met(1,i,j+1,k)*(
     *        c2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2)) +
     *        c1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1)) )  
     *  - ( mu(i,j-1,k)*met(3,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2)) +
     *        c1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1)) )*stry(j-1) 
     *     + mu(i,j-1,k)*met(4,i,j-1,k)*met(1,i,j-1,k)*(
     *        c2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2)) +
     *        c1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1)) ) ) ) )

c          lu(3,i,j,k) = r1/jac(i,j,k)
c          lu(3,i,j,k) = r1*ijac
*** 4 ops, tot=2126
            lu(3,i,j,k) = a1*lu(3,i,j,k) + sgn*r1*ijac

      enddo
      enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end
