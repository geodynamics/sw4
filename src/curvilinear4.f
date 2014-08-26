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
      subroutine CURVILINEAR4( ifirst, ilast, jfirst, jlast, kfirst,
     *                         klast, u, mu, la, met, jac, lu, 
     *                         onesided, acof, bope, ghcof, op )

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
      real*8 mux1, mux2, mux3, mux4, sgn, a1
      real*8 mucofu2, mucofuv, mucofuw, mucofv2, mucofvw, mucofw2
      real*8 ghcof(6), acof(6,8,8), bope(6,8)
      real*8 dudrp2, dudrp1, dudrm1, dudrm2
      real*8 dvdrp2, dvdrp1, dvdrm1, dvdrm2
      real*8 dwdrp2, dwdrp1, dwdrm1, dwdrm2
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
          lu(1,i,j,k) = a1*lu(1,i,j,k) + sgn*r1*ijac
          lu(2,i,j,k) = a1*lu(2,i,j,k) + sgn*r2*ijac
          lu(3,i,j,k) = a1*lu(3,i,j,k) + sgn*r3*ijac
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
c          lu(1,i,j,k) = r1*ijac
          lu(1,i,j,k) = a1*lu(1,i,j,k) + sgn*r1*ijac
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
c          lu(2,i,j,k) = r1*ijac
          lu(2,i,j,k) = a1*lu(2,i,j,k) + sgn*r1*ijac
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
c          lu(3,i,j,k) = r1*ijac
          lu(3,i,j,k) = a1*lu(3,i,j,k) + sgn*r1*ijac

      enddo
      enddo
      enddo
      end

