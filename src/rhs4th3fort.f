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
      subroutine rhs4th3fort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, nz, onesided, acof, bope, ghcof,
     +     uacc, u, mu, la, h, op ) bind(c)

*** in the interior: centered approximation of the spatial operator in the elastic wave equation
*** near physical boundaries: one-sided approximation of the spatial operator in the elastic wave equation

      implicit none

      real*8 tf, i6, i144, i12
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144, i12=1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer nz, onesided(6), m, q, kb, mb, qb, k1, k2
      real*8 acof(6,8,8), bope(6,8), ghcof(6)
      real*8 uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 muz1, muz2, muz3, muz4
      real*8 mucof(8), mu1zz, mu2zz, lau2yz
      real*8 lap2mu(8), mu3zz, mu3xz, mu3yz, lau1xz
      real*8 lau3zx, u3zim2, u3zim1, u3zip1, u3zip2
      real*8 lau3zy, u3zjm2, u3zjm1, u3zjp1, u3zjp2
      real*8 mu1zx, u1zim2, u1zim1, u1zip1, u1zip2
      real*8 mu2zy, u2zjm2, u2zjm1, u2zjp1, u2zjp2
      real*8 r1, r2, r3, h, cof, d4a, d4b, a1
      character*1 op
      parameter( d4a=2d0/3, d4b=-1d0/12 )

      cof = 1d0/(h*h)

      if( op .eq.'=' )then
         a1  = 0
      elseif( op.eq.'+' )then
         a1 = 1
      elseif( op.eq.'-' )then
         a1 = 1
         cof = -cof
      endif


      k1 = kfirst+2
      if (onesided(5).eq.1) k1 = 7;
      k2 = klast-2
      if (onesided(6).eq.1) k2 = nz-6;
!$OMP PARALLEL PRIVATE(i,j,k,kb,m,mb,q,qb,mux1,mux2,mux3,mux4,r1,r2,r3,
!$OMP*  muy1,muy2,muy3,muy4,muz1,muz2,muz3,muz4,
!$OMP*  mucof, mu1zz, mu2zz, lau2yz,
!$OMP* lap2mu, mu3zz, mu3xz, mu3yz, lau1xz,
!$OMP* lau3zx, u3zim2, u3zim1, u3zip1, u3zip2,
!$OMP* lau3zy, u3zjm2, u3zjm1, u3zjp1, u3zjp2,
!$OMP* mu1zx, u1zim2, u1zim1, u1zip1, u1zip2,
!$OMP* mu2zy, u2zjm2, u2zjm1, u2zjp1, u2zjp2 )

c the centered stencil can be evaluated 2 points away from the boundary
!$OMP DO
      do k=k1,k2
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))
c
            muz1 = mu(i,j,k-1)-tf*(mu(i,j,k)+mu(i,j,k-2))
            muz2 = mu(i,j,k-2)+mu(i,j,k+1)+3*(mu(i,j,k)+mu(i,j,k-1))
            muz3 = mu(i,j,k-1)+mu(i,j,k+2)+3*(mu(i,j,k+1)+mu(i,j,k))
            muz4 = mu(i,j,k+1)-tf*(mu(i,j,k)+mu(i,j,k+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) 
     *              + muz1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *                muz2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *                muz3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *                muz4*(u(1,i,j,k+2)-u(1,i,j,k)) )

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k))
     *              + muz1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *                muz2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *                muz3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *                muz4*(u(2,i,j,k+2)-u(2,i,j,k)) ) 

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) +
     *             (2*muz1+la(i,j,k-1)-tf*(la(i,j,k)+la(i,j,k-2)))*
     *                     (u(3,i,j,k-2)-u(3,i,j,k))+
     *      (2*muz2+la(i,j,k-2)+la(i,j,k+1)+3*(la(i,j,k)+la(i,j,k-1)))*
     *                     (u(3,i,j,k-1)-u(3,i,j,k))+ 
     *      (2*muz3+la(i,j,k-1)+la(i,j,k+2)+3*(la(i,j,k+1)+la(i,j,k)))*
     *                     (u(3,i,j,k+1)-u(3,i,j,k))+
     *             (2*muz4+la(i,j,k+1)-tf*(la(i,j,k)+la(i,j,k+2)))*
     *                     (u(3,i,j,k+2)-u(3,i,j,k)) )


*** Mixed derivatives:
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (la*w_z)_x
     *          + i144*( la(i-2,j,k)*(u(3,i-2,j,k-2)-u(3,i-2,j,k+2)+
     *                        8*(-u(3,i-2,j,k-1)+u(3,i-2,j,k+1))) - 8*(
     *                   la(i-1,j,k)*(u(3,i-1,j,k-2)-u(3,i-1,j,k+2)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i-1,j,k+1))) )+8*(
     *                   la(i+1,j,k)*(u(3,i+1,j,k-2)-u(3,i+1,j,k+2)+
     *                        8*(-u(3,i+1,j,k-1)+u(3,i+1,j,k+1))) ) - (
     *                   la(i+2,j,k)*(u(3,i+2,j,k-2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i+2,j,k-1)+u(3,i+2,j,k+1))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) 
***   (mu*w_x)_z
     *          + i144*( mu(i,j,k-2)*(u(3,i-2,j,k-2)-u(3,i+2,j,k-2)+
     *                        8*(-u(3,i-1,j,k-2)+u(3,i+1,j,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i-2,j,k-1)-u(3,i+2,j,k-1)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i+1,j,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i-2,j,k+1)-u(3,i+2,j,k+1)+
     *                        8*(-u(3,i-1,j,k+1)+u(3,i+1,j,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i-2,j,k+2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i-1,j,k+2)+u(3,i+1,j,k+2))) )) 

***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) 
*** (la*w_z)_y
     *          + i144*( la(i,j-2,k)*(u(3,i,j-2,k-2)-u(3,i,j-2,k+2)+
     *                        8*(-u(3,i,j-2,k-1)+u(3,i,j-2,k+1))) - 8*(
     *                   la(i,j-1,k)*(u(3,i,j-1,k-2)-u(3,i,j-1,k+2)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j-1,k+1))) )+8*(
     *                   la(i,j+1,k)*(u(3,i,j+1,k-2)-u(3,i,j+1,k+2)+
     *                        8*(-u(3,i,j+1,k-1)+u(3,i,j+1,k+1))) ) - (
     *                   la(i,j+2,k)*(u(3,i,j+2,k-2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j+2,k-1)+u(3,i,j+2,k+1))) ))
*** (mu*w_y)_z
     *          + i144*( mu(i,j,k-2)*(u(3,i,j-2,k-2)-u(3,i,j+2,k-2)+
     *                        8*(-u(3,i,j-1,k-2)+u(3,i,j+1,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i,j-2,k-1)-u(3,i,j+2,k-1)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j+1,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i,j-2,k+1)-u(3,i,j+2,k+1)+
     *                        8*(-u(3,i,j-1,k+1)+u(3,i,j+1,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i,j-2,k+2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j-1,k+2)+u(3,i,j+1,k+2))) )) 
***  (mu*u_z)_x
            r3 = r3 + i144*( mu(i-2,j,k)*(u(1,i-2,j,k-2)-u(1,i-2,j,k+2)+
     *                        8*(-u(1,i-2,j,k-1)+u(1,i-2,j,k+1))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j,k-2)-u(1,i-1,j,k+2)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i-1,j,k+1))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j,k-2)-u(1,i+1,j,k+2)+
     *                        8*(-u(1,i+1,j,k-1)+u(1,i+1,j,k+1))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j,k-2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i+2,j,k-1)+u(1,i+2,j,k+1))) )) 
*** (mu*v_z)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i,j-2,k-2)-u(2,i,j-2,k+2)+
     *                        8*(-u(2,i,j-2,k-1)+u(2,i,j-2,k+1))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i,j-1,k-2)-u(2,i,j-1,k+2)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j-1,k+1))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i,j+1,k-2)-u(2,i,j+1,k+2)+
     *                        8*(-u(2,i,j+1,k-1)+u(2,i,j+1,k+1))) ) - (
     *                   mu(i,j+2,k)*(u(2,i,j+2,k-2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j+2,k-1)+u(2,i,j+2,k+1))) ))
***   (la*u_x)_z
     *          + i144*( la(i,j,k-2)*(u(1,i-2,j,k-2)-u(1,i+2,j,k-2)+
     *                        8*(-u(1,i-1,j,k-2)+u(1,i+1,j,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(1,i-2,j,k-1)-u(1,i+2,j,k-1)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i+1,j,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(1,i-2,j,k+1)-u(1,i+2,j,k+1)+
     *                        8*(-u(1,i-1,j,k+1)+u(1,i+1,j,k+1))) ) - (
     *                   la(i,j,k+2)*(u(1,i-2,j,k+2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i-1,j,k+2)+u(1,i+1,j,k+2))) )) 
*** (la*v_y)_z
     *          + i144*( la(i,j,k-2)*(u(2,i,j-2,k-2)-u(2,i,j+2,k-2)+
     *                        8*(-u(2,i,j-1,k-2)+u(2,i,j+1,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(2,i,j-2,k-1)-u(2,i,j+2,k-1)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j+1,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(2,i,j-2,k+1)-u(2,i,j+2,k+1)+
     *                        8*(-u(2,i,j-1,k+1)+u(2,i,j+1,k+1))) ) - (
     *                   la(i,j,k+2)*(u(2,i,j-2,k+2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j-1,k+2)+u(2,i,j+1,k+2))) )) 

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP ENDDO

c low-k boundary modified stencils
      if (onesided(5).eq.1) then
!$OMP DO
      do k=1,6
c the centered stencil can be used in the x- and y-directions
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) )
c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            do q=1,8
              mucof(q)=0
              do m=1,8
                mucof(q) = mucof(q)+acof(k,q,m)*mu(i,j,m)
              enddo
            end do
c computing the second derivative
            mu1zz = 0
            do q=1,8
              mu1zz = mu1zz + mucof(q)*u(1,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(k)*mu(i,j,1)*u(1,i,j,0))

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) )
c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
            mu2zz = 0
            do q=1,8
              mu2zz = mu2zz + mucof(q)*u(2,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(k)*mu(i,j,1)*u(2,i,j,0))

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
            do q=1,8
              lap2mu(q)=0
              do m=1,8
                lap2mu(q) = lap2mu(q)+acof(k,q,m)*
     +                                (la(i,j,m)+2*mu(i,j,m))
              enddo
            end do
c computing the second derivative
            mu3zz = 0
            do q=1,8
              mu3zz = mu3zz + lap2mu(q)*u(3,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(k)*(la(i,j,1)+2*mu(i,j,1))*
     +           u(3,i,j,0))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) ))
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do q=1,8
              u3zip2 = u3zip2 + bope(k,q)*u(3,i+2,j,q)
              u3zip1 = u3zip1 + bope(k,q)*u(3,i+1,j,q)
              u3zim1 = u3zim1 + bope(k,q)*u(3,i-1,j,q)
              u3zim2 = u3zim2 + bope(k,q)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do q=1,8
              mu3xz = mu3xz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) ))
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do q=1,8
              u3zjp2 = u3zjp2 + bope(k,q)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 + bope(k,q)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 + bope(k,q)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 + bope(k,q)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do q=1,8
              mu3yz = mu3yz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do q=1,8
              u1zip2 = u1zip2 + bope(k,q)*u(1,i+2,j,q)
              u1zip1 = u1zip1 + bope(k,q)*u(1,i+1,j,q)
              u1zim1 = u1zim1 + bope(k,q)*u(1,i-1,j,q)
              u1zim2 = u1zim2 + bope(k,q)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do q=1,8
              u2zjp2 = u2zjp2 + bope(k,q)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 + bope(k,q)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 + bope(k,q)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 + bope(k,q)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do q=1,8
              lau1xz = lau1xz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do q=1,8
              lau2yz = lau2yz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + lau2yz

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP ENDDO
      endif

c high-k boundary
      if (onesided(6).eq.1) then
!$OMP DO
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) )
c all indices ending with 'b' are indices relative to the boundary, going into the domain (1,2,3,...)
            kb = nz-k+1
c all coefficient arrays (acof, bope, ghcof) should be indexed with these indices
c all solution and material property arrays should be indexed with (i,j,k)

c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            do qb=1,8
              mucof(qb)=0
              do mb=1,8
                m = nz-mb+1
                mucof(qb) = mucof(qb)+acof(kb,qb,mb)*mu(i,j,m)
              enddo
            end do
c computing the second derivative
            mu1zz = 0
            do qb=1,8
              q = nz-qb+1
              mu1zz = mu1zz + mucof(qb)*u(1,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(kb)*mu(i,j,nz)*u(1,i,j,nz+1))

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) )
c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
            mu2zz = 0
            do qb=1,8
              q = nz-qb+1
              mu2zz = mu2zz + mucof(qb)*u(2,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(kb)*mu(i,j,nz)*u(2,i,j,nz+1))

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
            do qb=1,8
              lap2mu(qb)=0
              do mb=1,8
                m = nz-mb+1
                lap2mu(qb) = lap2mu(qb)+acof(kb,qb,mb)*
     +                                (la(i,j,m)+2*mu(i,j,m))
              enddo
            end do
c computing the second derivative
            mu3zz = 0
            do qb=1,8
              q = nz-qb+1
              mu3zz = mu3zz + lap2mu(qb)*u(3,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(kb)*(la(i,j,nz)+2*mu(i,j,nz))*
     +           u(3,i,j,nz+1))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) ))
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do qb=1,8
              q = nz-qb+1
              u3zip2 = u3zip2 - bope(kb,qb)*u(3,i+2,j,q)
              u3zip1 = u3zip1 - bope(kb,qb)*u(3,i+1,j,q)
              u3zim1 = u3zim1 - bope(kb,qb)*u(3,i-1,j,q)
              u3zim2 = u3zim2 - bope(kb,qb)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do qb=1,8
              q = nz-qb+1
              mu3xz = mu3xz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) ))
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do qb=1,8
              q = nz-qb+1
              u3zjp2 = u3zjp2 - bope(kb,qb)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 - bope(kb,qb)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 - bope(kb,qb)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 - bope(kb,qb)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do qb=1,8
              q = nz-qb+1
              mu3yz = mu3yz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do qb=1,8
              q = nz-qb+1
              u1zip2 = u1zip2 - bope(kb,qb)*u(1,i+2,j,q)
              u1zip1 = u1zip1 - bope(kb,qb)*u(1,i+1,j,q)
              u1zim1 = u1zim1 - bope(kb,qb)*u(1,i-1,j,q)
              u1zim2 = u1zim2 - bope(kb,qb)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do qb=1,8
              q = nz-qb+1
              u2zjp2 = u2zjp2 - bope(kb,qb)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 - bope(kb,qb)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 - bope(kb,qb)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 - bope(kb,qb)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do qb=1,8
              q = nz-qb+1
              lau1xz = lau1xz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do qb=1,8
              q = nz-qb+1
              lau2yz = lau2yz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + lau2yz

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP ENDDO
      endif
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine rhs4th3fortsgstr( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, nz, onesided, acof, bope, ghcof,
     +     uacc, u, mu, la, h, strx, stry, strz, op )  bind(c)

*** Routine with supergrid stretchings, strx, stry, and strz.
***
*** In the interior: centered approximation of the spatial operator in the elastic wave equation
*** near physical boundaries: one-sided approximation of the spatial operator in the elastic wave equation

      implicit none

      real*8 tf, i6, i144, i12
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144, i12=1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer nz, onesided(6), m, q, kb, mb, qb, k1, k2
      real*8 acof(6,8,8), bope(6,8), ghcof(6)
      real*8 uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 muz1, muz2, muz3, muz4
      real*8 mucof(8), mu1zz, mu2zz, lau2yz
      real*8 lap2mu(8), mu3zz, mu3xz, mu3yz, lau1xz
      real*8 lau3zx, u3zim2, u3zim1, u3zip1, u3zip2
      real*8 lau3zy, u3zjm2, u3zjm1, u3zjp1, u3zjp2
      real*8 mu1zx, u1zim2, u1zim1, u1zip1, u1zip2
      real*8 mu2zy, u2zjm2, u2zjm1, u2zjp1, u2zjp2
      real*8 r1, r2, r3, h, cof, d4a, d4b, a1
      real*8 mucofs, lap2mus, lacofs
      real*8 strx(ifirst:ilast), stry(jfirst:jlast), strz(kfirst:klast)
      character*1 op
      parameter( d4a=2d0/3, d4b=-1d0/12 )

      cof = 1d0/(h*h)
      if( op.eq.'=' )then
         a1 = 0
      elseif( op.eq.'+')then
         a1 = 1
      elseif( op.eq.'-')then
         a1 = 1
         cof = -cof
      endif

      k1 = kfirst+2
      if (onesided(5).eq.1) k1 = 7
      k2 = klast-2
      if (onesided(6).eq.1) k2 = nz-6
c the centered stencil can be evaluated 2 points away from the boundary
!$OMP PARALLEL PRIVATE(k,i,j,mux1,mux2,mux3,mux4,muy1,muy2,muy3,muy4,
!$OMP*   muz1,muz2,muz3,muz4,
!$OMP*   r1,r2,r3,mucof,mu1zz,mu2zz,mu3zz,lap2mu,q,m,u3zip2,u3zip1,
!$OMP*   u3zim1,u3zim2,lau3zx,mu3xz,u3zjp2,u3zjp1,u3zjm1,u3zjm2,lau3zy,
!$OMP*   mu3yz,mu1zx,u1zip2,u1zip1,u1zim1,u1zim2,lap2mus,mucofs,lacofs,
!$OMP*   u2zjp2,u2zjp1,u2zjm1,u2zjm2,mu2zy,lau1xz,lau2yz,kb,qb,mb)
!$OMP DO
      do k=k1,k2
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a, 28x3 = 84 ops
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                     tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                      3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                      3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                     tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                    tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                     3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                     3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                    tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))
            muz1 = mu(i,j,k-1)*strz(k-1)-
     *                    tf*(mu(i,j,k)*strz(k)+mu(i,j,k-2)*strz(k-2))
            muz2 = mu(i,j,k-2)*strz(k-2)+mu(i,j,k+1)*strz(k+1)+
     *                     3*(mu(i,j,k)*strz(k)+mu(i,j,k-1)*strz(k-1))
            muz3 = mu(i,j,k-1)*strz(k-1)+mu(i,j,k+2)*strz(k+2)+
     *                     3*(mu(i,j,k+1)*strz(k+1)+mu(i,j,k)*strz(k))
            muz4 = mu(i,j,k+1)*strz(k+1)-
     *                    tf*(mu(i,j,k)*strz(k)+mu(i,j,k+2)*strz(k+2))
c            write(*,*) 'F muz 1-4 ',muz1, muz2, muz3, muz4
*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda in the same way as we did for mu
*** 75 ops
            r1 = i6*( strx(i)*( (2*mux1+la(i-1,j,k)*strx(i-1)-
     *          tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *           3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *           3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *          tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *                muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) + strz(k)*(
     *                muz1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *                muz2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *                muz3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *                muz4*(u(1,i,j,k+2)-u(1,i,j,k)) ) )

*** 75 ops
            r2 = i6*( strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) ) + stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                 tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *               tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) + strz(k)*(
     *                muz1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *                muz2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *                muz3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *                muz4*(u(2,i,j,k+2)-u(2,i,j,k)) ) )

*** 75 ops
            r3 = i6*( strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) + strz(k)*(
     *             (2*muz1+la(i,j,k-1)*strz(k-1)-
     *                 tf*(la(i,j,k)*strz(k)+la(i,j,k-2)*strz(k-2)))*
     *                     (u(3,i,j,k-2)-u(3,i,j,k))+
     *      (2*muz2+la(i,j,k-2)*strz(k-2)+la(i,j,k+1)*strz(k+1)+
     *                 3*(la(i,j,k)*strz(k)+la(i,j,k-1)*strz(k-1)))*
     *                     (u(3,i,j,k-1)-u(3,i,j,k))+ 
     *      (2*muz3+la(i,j,k-1)*strz(k-1)+la(i,j,k+2)*strz(k+2)+
     *                 3*(la(i,j,k+1)*strz(k+1)+la(i,j,k)*strz(k)))*
     *                     (u(3,i,j,k+1)-u(3,i,j,k))+
     *             (2*muz4+la(i,j,k+1)*strz(k+1)-
     *               tf*(la(i,j,k)*strz(k)+la(i,j,k+2)*strz(k+2)))*
     *                     (u(3,i,j,k+2)-u(3,i,j,k)) ) )


*** Mixed derivatives:
*** 29ops /mixed derivative
*** 116 ops for r1
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (la*w_z)_x
     *          + strx(i)*strz(k)*       
     *            i144*( la(i-2,j,k)*(u(3,i-2,j,k-2)-u(3,i-2,j,k+2)+
     *                        8*(-u(3,i-2,j,k-1)+u(3,i-2,j,k+1))) - 8*(
     *                   la(i-1,j,k)*(u(3,i-1,j,k-2)-u(3,i-1,j,k+2)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i-1,j,k+1))) )+8*(
     *                   la(i+1,j,k)*(u(3,i+1,j,k-2)-u(3,i+1,j,k+2)+
     *                        8*(-u(3,i+1,j,k-1)+u(3,i+1,j,k+1))) ) - (
     *                   la(i+2,j,k)*(u(3,i+2,j,k-2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i+2,j,k-1)+u(3,i+2,j,k+1))) )) 
***   (mu*v_x)_y
     *          + strx(i)*stry(j)*       
     *            i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) 
***   (mu*w_x)_z
     *          + strx(i)*strz(k)*       
     *            i144*( mu(i,j,k-2)*(u(3,i-2,j,k-2)-u(3,i+2,j,k-2)+
     *                        8*(-u(3,i-1,j,k-2)+u(3,i+1,j,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i-2,j,k-1)-u(3,i+2,j,k-1)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i+1,j,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i-2,j,k+1)-u(3,i+2,j,k+1)+
     *                        8*(-u(3,i-1,j,k+1)+u(3,i+1,j,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i-2,j,k+2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i-1,j,k+2)+u(3,i+1,j,k+2))) )) 

*** 116 ops for r2
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *         + strx(i)*stry(j)*
     *            i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) 
*** (la*w_z)_y
     *          + stry(j)*strz(k)*
     *            i144*( la(i,j-2,k)*(u(3,i,j-2,k-2)-u(3,i,j-2,k+2)+
     *                        8*(-u(3,i,j-2,k-1)+u(3,i,j-2,k+1))) - 8*(
     *                   la(i,j-1,k)*(u(3,i,j-1,k-2)-u(3,i,j-1,k+2)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j-1,k+1))) )+8*(
     *                   la(i,j+1,k)*(u(3,i,j+1,k-2)-u(3,i,j+1,k+2)+
     *                        8*(-u(3,i,j+1,k-1)+u(3,i,j+1,k+1))) ) - (
     *                   la(i,j+2,k)*(u(3,i,j+2,k-2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j+2,k-1)+u(3,i,j+2,k+1))) ))
*** (mu*w_y)_z
     *          + stry(j)*strz(k)*
     *            i144*( mu(i,j,k-2)*(u(3,i,j-2,k-2)-u(3,i,j+2,k-2)+
     *                        8*(-u(3,i,j-1,k-2)+u(3,i,j+1,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i,j-2,k-1)-u(3,i,j+2,k-1)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j+1,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i,j-2,k+1)-u(3,i,j+2,k+1)+
     *                        8*(-u(3,i,j-1,k+1)+u(3,i,j+1,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i,j-2,k+2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j-1,k+2)+u(3,i,j+1,k+2))) )) 
*** 116 ops for r3
***  (mu*u_z)_x
            r3 = r3 + strx(i)*strz(k)*
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j,k-2)-u(1,i-2,j,k+2)+
     *                        8*(-u(1,i-2,j,k-1)+u(1,i-2,j,k+1))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j,k-2)-u(1,i-1,j,k+2)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i-1,j,k+1))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j,k-2)-u(1,i+1,j,k+2)+
     *                        8*(-u(1,i+1,j,k-1)+u(1,i+1,j,k+1))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j,k-2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i+2,j,k-1)+u(1,i+2,j,k+1))) )) 
*** (mu*v_z)_y
     *         + stry(j)*strz(k)*
     *            i144*( mu(i,j-2,k)*(u(2,i,j-2,k-2)-u(2,i,j-2,k+2)+
     *                        8*(-u(2,i,j-2,k-1)+u(2,i,j-2,k+1))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i,j-1,k-2)-u(2,i,j-1,k+2)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j-1,k+1))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i,j+1,k-2)-u(2,i,j+1,k+2)+
     *                        8*(-u(2,i,j+1,k-1)+u(2,i,j+1,k+1))) ) - (
     *                   mu(i,j+2,k)*(u(2,i,j+2,k-2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j+2,k-1)+u(2,i,j+2,k+1))) ))
***   (la*u_x)_z
     *         + strx(i)*strz(k)*
     *            i144*( la(i,j,k-2)*(u(1,i-2,j,k-2)-u(1,i+2,j,k-2)+
     *                        8*(-u(1,i-1,j,k-2)+u(1,i+1,j,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(1,i-2,j,k-1)-u(1,i+2,j,k-1)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i+1,j,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(1,i-2,j,k+1)-u(1,i+2,j,k+1)+
     *                        8*(-u(1,i-1,j,k+1)+u(1,i+1,j,k+1))) ) - (
     *                   la(i,j,k+2)*(u(1,i-2,j,k+2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i-1,j,k+2)+u(1,i+1,j,k+2))) )) 
*** (la*v_y)_z
     *         + stry(j)*strz(k)*
     *            i144*( la(i,j,k-2)*(u(2,i,j-2,k-2)-u(2,i,j+2,k-2)+
     *                        8*(-u(2,i,j-1,k-2)+u(2,i,j+1,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(2,i,j-2,k-1)-u(2,i,j+2,k-1)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j+1,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(2,i,j-2,k+1)-u(2,i,j+2,k+1)+
     *                        8*(-u(2,i,j-1,k+1)+u(2,i,j+1,k+1))) ) - (
     *                   la(i,j,k+2)*(u(2,i,j-2,k+2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j-1,k+2)+u(2,i,j+1,k+2))) )) 

*** 9 ops
            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP END DO

*** Total number of ops away from boundary loop, per grid point computed:
*** 84+75*3+3*116+9 = 666 ops
*** Memory access mu,lambda  2*(5+5+5-2) = 26  (stencil is a cross)
***                          u = 3*( 3*5*5-3*5+1) = 183 (stencil is three intersecting planes)
***                          uacc = 3
***             strx,stry,strz = 5*3 = 15
***  Memory access read, 227 vars. = 1816 bytes (double prec.)
***                write, 3 vars. = 24 bytes (double prec.)
***
*** Loop below for the upper boundary (1 <= k <= 6), per grid point:
***    1244 arithmetic ops
***   Memory access, mu, lambda (5 pts in i,j, 8pts in k:) 2*(5+5+8-2)=32      
***                  u       3*(5*5+5*9+5*9-5-5-9+1) : 291 ( u uses ghost point+8 pts in k)
***                  uacc    3
***            strx,stry,strz:  5*3 = 15
***            acof:   6*8 = 48
***            bope:     8
***      Memory access read 397 vars.= 3176 bytes (double)
***                   write, 3 vars = 24 bytes

c low-k boundary modified stencils
      if (onesided(5).eq.1) then
!$OMP DO
      do k=1,6
c the centered stencil can be used in the x- and y-directions
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
*** 56 ops
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                  3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                  3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
c
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                  3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                  3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
*** 62 ops, tot=118
            r1 = i6*(strx(i)*((2*mux1+la(i-1,j,k)*strx(i-1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *                   3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *                   3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *                muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) )
c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient, 
c leave out the z- supergrid stretching strz, since it will
c never be used together with the sbp-boundary operator
*** 8*(19) = 152 ops, tot = 270 (+another 8*33 which is added 73 lines down).
*** Now brought down to 8*(19+19)=304 from 152+8*33=416, newtot=422
            mu1zz = 0
            mu2zz = 0
            mu3zz = 0
            do q=1,8
               mucofs = acof(k,q,1)*mu(i,j,1)+acof(k,q,2)*mu(i,j,2)+
     *acof(k,q,3)*mu(i,j,3)+acof(k,q,4)*mu(i,j,4)+acof(k,q,5)*mu(i,j,5)+
     *acof(k,q,6)*mu(i,j,6)+acof(k,q,7)*mu(i,j,7)+acof(k,q,8)*mu(i,j,8)
               mu1zz = mu1zz + mucofs*u(1,i,j,q)
               mu2zz = mu2zz + mucofs*u(2,i,j,q)
               lacofs = acof(k,q,1)*la(i,j,1)+acof(k,q,2)*la(i,j,2)+
     *acof(k,q,3)*la(i,j,3)+acof(k,q,4)*la(i,j,4)+acof(k,q,5)*la(i,j,5)+
     *acof(k,q,6)*la(i,j,6)+acof(k,q,7)*la(i,j,7)+acof(k,q,8)*la(i,j,8)
               mu3zz = mu3zz +(lacofs+2*mucofs)*u(3,i,j,q)
c               lap2mus = acof(k,q,1)*(la(i,j,1)+2*mu(i,j,1)) +
c     *        acof(k,q,2)*(la(i,j,2)+2*mu(i,j,2)) +
c     *        acof(k,q,3)*(la(i,j,3)+2*mu(i,j,3)) +
c     *        acof(k,q,4)*(la(i,j,4)+2*mu(i,j,4)) +
c     *        acof(k,q,5)*(la(i,j,5)+2*mu(i,j,5)) +
c     *        acof(k,q,6)*(la(i,j,6)+2*mu(i,j,6)) +
c     *        acof(k,q,7)*(la(i,j,7)+2*mu(i,j,7)) +
c     *        acof(k,q,8)*(la(i,j,8)+2*mu(i,j,8))
c               mu3zz = mu3zz + lap2mus*u(3,i,j,q)
            enddo
c            do q=1,8
c              mucof(q)=0
c              do m=1,8
c                mucof(q) = mucof(q)+acof(k,q,m)*mu(i,j,m)
c              enddo
c            end do
c computing the second derivative
c            mu1zz = 0
c            do q=1,8
c              mu1zz = mu1zz + mucof(q)*u(1,i,j,q)
c            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
*** 66 ops, tot=488
            r1 = r1 + (mu1zz + ghcof(k)*mu(i,j,1)*u(1,i,j,0))

            r2 = i6*(strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) )+ stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                   tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                   3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                   3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *                  tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) )
c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
c            mu2zz = 0
c            do q=1,8
c              mu2zz = mu2zz + mucof(q)*u(2,i,j,q)
c            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
*** 30 ops, tot=518
            r2 = r2 + (mu2zz + ghcof(k)*mu(i,j,1)*u(2,i,j,0))

            r3 = i6*(strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
*** 8*33 = 264 ops, tot=630
c            mu3zz = 0
c            do q=1,8
c               lap2mus = acof(k,q,1)*(la(i,j,1)+2*mu(i,j,1)) +
c     *        acof(k,q,2)*(la(i,j,2)+2*mu(i,j,2)) +
c     *        acof(k,q,3)*(la(i,j,3)+2*mu(i,j,3)) +
c     *        acof(k,q,4)*(la(i,j,4)+2*mu(i,j,4)) +
c     *        acof(k,q,5)*(la(i,j,5)+2*mu(i,j,5)) +
c     *        acof(k,q,6)*(la(i,j,6)+2*mu(i,j,6)) +
c     *        acof(k,q,7)*(la(i,j,7)+2*mu(i,j,7)) +
c     *        acof(k,q,8)*(la(i,j,8)+2*mu(i,j,8))
c               mu3zz = mu3zz + lap2mus*u(3,i,j,q)
c            enddo
c            do q=1,8
c              lap2mu(q)=0
c              do m=1,8
c                lap2mu(q) = lap2mu(q)+acof(k,q,m)*
c     +                                (la(i,j,m)+2*mu(i,j,m))
c              enddo
c            end do
cc computing the second derivative
c            mu3zz = 0
c            do q=1,8
c              mu3zz = mu3zz + lap2mu(q)*u(3,i,j,q)
c            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
*** 5 ops, tot=523
            r3 = r3 + (mu3zz + ghcof(k)*(la(i,j,1)+2*mu(i,j,1))*
     +           u(3,i,j,0))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
*** 56 ops, tot=579
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*(
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) )
***   (la*w_z)_x: NOT CENTERED
*** 8*8 + 12= 76 ops, tot=655
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do q=1,8
              u3zip2 = u3zip2 + bope(k,q)*u(3,i+2,j,q)
              u3zip1 = u3zip1 + bope(k,q)*u(3,i+1,j,q)
              u3zim1 = u3zim1 + bope(k,q)*u(3,i-1,j,q)
              u3zim2 = u3zim2 + bope(k,q)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + strx(i)*lau3zx

***   (mu*w_x)_z: NOT CENTERED
*** 8*9+2 = 74 ops, tot=729
            mu3xz=0
            do q=1,8
              mu3xz = mu3xz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + strx(i)*mu3xz

c cross-terms in second component of rhs
*** 56 ops , tot=785
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*(
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) )
*** (la*w_z)_y : NOT CENTERED
*** 8*8+12=76 ops, tot=861
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do q=1,8
              u3zjp2 = u3zjp2 + bope(k,q)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 + bope(k,q)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 + bope(k,q)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 + bope(k,q)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + stry(j)*lau3zy

*** (mu*w_y)_z: NOT CENTERED
*** 8*9+2 = 74 ops, tot=935
            mu3yz=0
            do q=1,8
              mu3yz = mu3yz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + stry(j)*mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
*** 76 ops, tot=1011
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do q=1,8
              u1zip2 = u1zip2 + bope(k,q)*u(1,i+2,j,q)
              u1zip1 = u1zip1 + bope(k,q)*u(1,i+1,j,q)
              u1zim1 = u1zim1 + bope(k,q)*u(1,i-1,j,q)
              u1zim2 = u1zim2 + bope(k,q)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + strx(i)*mu1zx

*** (mu*v_z)_y: NOT CENTERED
*** 76 ops, tot=1087
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do q=1,8
              u2zjp2 = u2zjp2 + bope(k,q)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 + bope(k,q)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 + bope(k,q)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 + bope(k,q)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + stry(j)*mu2zy

***   (la*u_x)_z: NOT CENTERED
*** 74 ops, tot=1161
            lau1xz=0
            do q=1,8
              lau1xz = lau1xz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + strx(i)*lau1xz

*** (la*v_y)_z: NOT CENTERED
*** 74 ops, tot=1235
            lau2yz=0
            do q=1,8
              lau2yz = lau2yz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + stry(j)*lau2yz
*** 9 ops, tot=1244
            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP END DO
      endif

c high-k boundary
      if (onesided(6).eq.1) then
!$OMP DO
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                  3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                  3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
c
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                  3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                  3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*(strx(i)*((2*mux1+la(i-1,j,k)*strx(i-1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *                   3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *                   3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) )

c all indices ending with 'b' are indices relative to the boundary, going into the domain (1,2,3,...)
            kb = nz-k+1
c all coefficient arrays (acof, bope, ghcof) should be indexed with these indices
c all solution and material property arrays should be indexed with (i,j,k)

c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            do qb=1,8
              mucof(qb)=0
              do mb=1,8
                m = nz-mb+1
                mucof(qb) = mucof(qb)+acof(kb,qb,mb)*mu(i,j,m)
              enddo
            end do
c computing the second derivative
            mu1zz = 0
            do qb=1,8
              q = nz-qb+1
              mu1zz = mu1zz + mucof(qb)*u(1,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(kb)*mu(i,j,nz)*u(1,i,j,nz+1))

            r2 = i6*(strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) )+ stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                   tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                   3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                   3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *                  tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) )

c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
            mu2zz = 0
            do qb=1,8
              q = nz-qb+1
              mu2zz = mu2zz + mucof(qb)*u(2,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(kb)*mu(i,j,nz)*u(2,i,j,nz+1))

            r3 = i6*(strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
            do qb=1,8
              lap2mu(qb)=0
              do mb=1,8
                m = nz-mb+1
                lap2mu(qb) = lap2mu(qb)+acof(kb,qb,mb)*
     +                                (la(i,j,m)+2*mu(i,j,m))
              enddo
            end do
c computing the second derivative
            mu3zz = 0
            do qb=1,8
              q = nz-qb+1
              mu3zz = mu3zz + lap2mu(qb)*u(3,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(kb)*(la(i,j,nz)+2*mu(i,j,nz))*
     +           u(3,i,j,nz+1))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*(
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) )
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do qb=1,8
              q = nz-qb+1
              u3zip2 = u3zip2 - bope(kb,qb)*u(3,i+2,j,q)
              u3zip1 = u3zip1 - bope(kb,qb)*u(3,i+1,j,q)
              u3zim1 = u3zim1 - bope(kb,qb)*u(3,i-1,j,q)
              u3zim2 = u3zim2 - bope(kb,qb)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + strx(i)*lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do qb=1,8
              q = nz-qb+1
              mu3xz = mu3xz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + strx(i)*mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*(
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) )
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do qb=1,8
              q = nz-qb+1
              u3zjp2 = u3zjp2 - bope(kb,qb)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 - bope(kb,qb)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 - bope(kb,qb)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 - bope(kb,qb)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + stry(j)*lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do qb=1,8
              q = nz-qb+1
              mu3yz = mu3yz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + stry(j)*mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do qb=1,8
              q = nz-qb+1
              u1zip2 = u1zip2 - bope(kb,qb)*u(1,i+2,j,q)
              u1zip1 = u1zip1 - bope(kb,qb)*u(1,i+1,j,q)
              u1zim1 = u1zim1 - bope(kb,qb)*u(1,i-1,j,q)
              u1zim2 = u1zim2 - bope(kb,qb)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + strx(i)*mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do qb=1,8
              q = nz-qb+1
              u2zjp2 = u2zjp2 - bope(kb,qb)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 - bope(kb,qb)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 - bope(kb,qb)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 - bope(kb,qb)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + stry(j)*mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do qb=1,8
              q = nz-qb+1
              lau1xz = lau1xz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + strx(i)*lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do qb=1,8
              q = nz-qb+1
              lau2yz = lau2yz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + stry(j)*lau2yz

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
!$OMP END DO
      endif
!$OMP END PARALLEL

      end

c----------------------------------------------------------
      subroutine rhserrfort(ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nz, h, fo, u2, lowZ, interZ, highZ)  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nz
      real*8 h, lowZ(3), interZ(3), highZ(3)
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li(3), l2(3), err(3)
c low-z points, error in rhs
      do c=1,3
        lowZ(c) = 0
      enddo
c this test only includes low-k boundary points
!$OMP PARALLEL PRIVATE(k,i,j,c,err)
!$OMP DO REDUCTION(max:lowZ)
      do k=1,6
         do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( lowZ(c).lt.err(c) )then
                lowZ(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO
c interior points, error in rhs
      do c=1,3
        interZ(c) = 0
      enddo
c this test only includes interior points
!$OMP DO REDUCTION(max:interZ)
      do k=7,nz-6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( interZ(c).lt.err(c) )then
                interZ(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO
c high-z points error in rhs
      do c=1,3
        highZ(c) = 0
      enddo
c this test only includes high-k boundary points
!$OMP DO REDUCTION(max:highZ)
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( highZ(c).lt.err(c) )then
                highZ(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      return
      end

c-------------------------------------------------
      subroutine rhouttlumf(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     nz, uacc, lu, fo, rho, lowZ, interZ, highZ)  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, kp, nz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lowZ(3), interZ(3), highZ(3), err(3)
      integer c, i, j, k

c  low-k boundary points
      do c=1,3
        lowZ(c)=0
      enddo
!$OMP PARALLEL PRIVATE(k,i,j,c,err)
!$OMP DO REDUCTION(max:lowZ)
      do k=1,6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) =rho(i,j,k)*uacc(c,i,j,k) - lu(c,i,j,k)
     +             - fo(c,i,j,k)
              if (abs(err(c)).gt.lowZ(c)) lowZ(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO

c evaluate error in the interior
      do c=1,3
        interZ(c)=0
      enddo
c interior points
!$OMP DO REDUCTION(max:interZ)
      do k=7,nz-6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) =rho(i,j,k)*uacc(c,i,j,k) - lu(c,i,j,k)
     +             - fo(c,i,j,k)
              if (abs(err(c)).gt.interZ(c)) interZ(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO

c this test only includes high-k boundary points
      do c=1,3
        highZ(c)=0
      enddo
!$OMP DO REDUCTION(max:highZ)
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) =rho(i,j,k)*uacc(c,i,j,k) - lu(c,i,j,k)
     +             - fo(c,i,j,k)
              if (abs(err(c)).gt.highZ(c)) highZ(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      end

c---------------------------------------------------------------------
      subroutine predfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, lu, fo, rho, dt2 )  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt2
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
c 2nd order accurate predictor of solution at t+dt
!$OMP PARALLEL PRIVATE(k,i,j,c)
!$OMP DO
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              up(c,i,j,k) = 2*u(c,i,j,k) - um(c,i,j,k) + 
     +             dt2*(lu(c,i,j,k) + fo(c,i,j,k))/rho(i,j,k)
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return 
      end

c-----------------------------------------------------------------------
      subroutine corrfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, lu, fo, rho, dt4 )  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt4, i12, dt4i12
      parameter( i12=1d0/12 )
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
c correct solution to 4th order accuracy
      dt4i12 = dt4*i12
!$OMP PARALLEL PRIVATE(k,i,j,c)
!$OMP DO
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              up(c,i,j,k) = up(c,i,j,k) + dt4i12*
     +             (lu(c,i,j,k) + fo(c,i,j,k))/rho(i,j,k)
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return
      end

c-----------------------------------------------------------------------
      subroutine dpdmtfort(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, u2, dt2i)  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt2i
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c evaluate 2nd divided time difference D+D-(u)
!$OMP PARALLEL PRIVATE(k,i,j,c)
!$OMP DO
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              u2(c,i,j,k) = dt2i*
     +             (up(c,i,j,k) - 2*u(c,i,j,k) + um(c,i,j,k))
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return
      end

c-----------------------------------------------------------------------
      subroutine updatememvar( ifirst, ilast, jfirst, jlast, kfirst,
     * klast, alp, alm, up, u, um, omega, dt, domain )  bind(c)

***********************************************************************
*** 
*** domain = 0 --> Entire domain
*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
***
***********************************************************************      

      implicit none
      real*8 i6
      parameter( i6=1d0/6 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      integer domain, k1, k2
      real*8 omega, dt, dto, icp, cm, cof3
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

      dto = dt*omega
      icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
      cm  =     1d0/2 - 1/(2*dto) - dto/4 + dto*dto/12 
      if( domain.eq.0 )then
         k1 = kfirst
         k2 = klast
      elseif( domain.eq.1 )then
         k1 = kfirst
         k2 = kfirst+2
      elseif( domain.eq.2 )then
         k1 = klast-2
         k2 = klast
      endif
!$OMP PARALLEL PRIVATE(k,i,j,c)
!$OMP DO
      do k=k1,k2
         do j=jfirst,jlast
            do i=ifirst,ilast
               do c=1,3
                  alp(c,i,j,k) = icp*(-cm*alm(c,i,j,k) + u(c,i,j,k) +
     *         i6*( dto*dto*u(c,i,j,k) + dto*(up(c,i,j,k)-um(c,i,j,k)) +
     *              (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k))  )
     *                                                  )
               enddo
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine dpdmtfortatt(ifirst, ilast, jfirst, jlast, kfirst, 
     +    klast, up, u, um, dt2i)  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt2i
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c evaluate 2nd divided time difference D+D-(u), and return in um
!$OMP PARALLEL PRIVATE(k,i,j,c)
!$OMP DO
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              um(c,i,j,k) = dt2i*
     +             (up(c,i,j,k) - 2*u(c,i,j,k) + um(c,i,j,k))
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return
      end

c-----------------------------------------------------------------------
      subroutine satt(up, qs, dt, cfreq, 
     +     ifirst, ilast, jfirst, jlast, kfirst, klast)  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt, cfreq, pi, fact
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 qs(ifirst:ilast,jfirst:jlast,kfirst:klast)
c simple attenuation model
      pi = 4*atan(1.d0)
!$OMP PARALLEL PRIVATE(k,i,j,c,fact)
!$OMP DO
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            fact = exp(-pi*cfreq*dt/qs(i,j,k))
c tmp
c$$$            if (i.eq.ifirst .and. j.eq.jfirst .and. k.eq.kfirst) then
c$$$              write(*,*) "dt = ", dt, " qs = ", qs(i,j,k), 
c$$$     +             " cfreq = ", cfreq, " fact = ", fact
c$$$            endif
            do c=1,3
              up(c,i,j,k) = fact*up(c,i,j,k)
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      return
      end

c-----------------------------------------------------------------------
      subroutine SOLVEATTFREEAC( ifirst, ilast, jfirst, jlast, kfirst, 
     *                           klast, alpha, cof, up )  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      real*8 alpha(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 cof
      k = 0
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
      do j=jfirst+2,jlast-2
         do i=ifirst+2,ilast-2
            alpha(1,i,j,k) = alpha(1,i,j,k) + cof*up(1,i,j,k)
            alpha(2,i,j,k) = alpha(2,i,j,k) + cof*up(2,i,j,k)
            alpha(3,i,j,k) = alpha(3,i,j,k) + cof*up(3,i,j,k)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine SOLVEATTFREEC( ifirst, ilast, jfirst, jlast, kfirst, 
     *       klast, u, mu, la, muve, lambdave, bforcerhs, met, s, 
     *       usesg, sgstrx, sgstry )  bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, kl
      integer usesg
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 bforcerhs(3,ifirst:ilast,jfirst:jlast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 muve(ifirst:ilast,jfirst:jlast), s(0:4)
      real*8 lambdave(ifirst:ilast,jfirst:jlast)
      real*8 mupt, lapt, ac, bc, cc, dc, s0i
      real*8 sgstrx(ifirst:ilast), sgstry(jfirst:jlast)
      real*8 sgx, sgy, isgx, isgy, m2sg, m3sg, m4sg
*** Hardcoded for the k=1 surface
      k  = 1
      kl = 1
      sgx = 1
      sgy = 1
      isgx = 1
      isgy = 1
      s0i= 1/s(0)
!$OMP PARALLEL PRIVATE(i,j,mupt,lapt,sgx,sgy,isgx,isgy,m2sg,m3sg,m4sg,
!$OMP*                 ac,bc,cc,dc)
!$OMP DO
      do j=jfirst+2,jlast-2
         do i=ifirst+2,ilast-2
            mupt = mu(i,j,k)-muve(i,j)
            lapt = la(i,j,k)-lambdave(i,j)
            if( usesg.eq.1 )then
               sgx = sgstrx(i)
               sgy = sgstry(j)
               isgy = 1/sgy
               isgx = 1/sgx
            endif
            m2sg = SQRT(sgx*isgy)
            m3sg = 1/m2sg
            m4sg = isgx*m2sg

            ac = sgx*isgy*met(2,i,j,k)**2+isgx*sgy*met(3,i,j,k)**2+
     *                                     isgx*isgy*met(4,i,j,k)**2
            bc = 1/(mupt*ac)
            cc = (mupt+lapt)/(2*mupt+lapt)*bc/ac
            dc = cc*(  m2sg*met(2,i,j,k)*bforcerhs(1,i,j)  
     *               + m3sg*met(3,i,j,k)*bforcerhs(2,i,j)  
     *               + m4sg*met(4,i,j,k)*bforcerhs(3,i,j) )
            u(1,i,j,k-kl) = s0i*(    
     *                bc*bforcerhs(1,i,j) - dc*met(2,i,j,k)*m2sg )
            u(2,i,j,k-kl) = s0i*(
     *                bc*bforcerhs(2,i,j) - dc*met(3,i,j,k)*m3sg )
            u(3,i,j,k-kl) = s0i*( 
     *                bc*bforcerhs(3,i,j) - dc*met(4,i,j,k)*m4sg )
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end

c-----------------------------------------------------------------------
      subroutine ADDBSTRESSWRESC( ifirst, ilast, jfirst, jlast, kfirst,
     *      klast, nz, alphap, alpham, muve, lave, bforcerhs, 
     *      u, um, met, side, dt, omegave, memforce, muvebnd,
     *      lambdavebnd, s, cof, usesg, sgstrx, sgstry )  bind(c)
      implicit none
      real*8 i6, c1, c2
      parameter( i6=1d0/6, c1=2d0/3, c2=-1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, kl
      integer nz, side, usesg
      real*8 alphap(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 alpham(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 muve(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lave(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 muvebnd(ifirst:ilast,jfirst:jlast)
      real*8 lambdavebnd(ifirst:ilast,jfirst:jlast)
      real*8 bforcerhs(3,ifirst:ilast,jfirst:jlast)
      real*8 memforce(3,ifirst:ilast,jfirst:jlast)
      real*8 dt, omegave, s(0:4), cof, omdt, cp, cm, r1, r2, r3
      real*8 rhs1, rhs2, rhs3, un1, vn1, wn1, rtu, ac
      real*8 sgstrx(ifirst:ilast), sgstry(jfirst:jlast)
      real*8 sgx, sgy, isgx, isgy, m2sg, m3sg, m4sg
      if( side.eq.5 )then
         k  = 1
         kl = 1
      elseif( side.eq.6 )then
         k  = nz
         kl = -1
      endif
      omdt = omegave*dt
      cp = 0.5d0 + 1/(2*omdt) + omdt/4 + omdt*omdt/12
      cm = 0.5d0 - 1/(2*omdt) - omdt/4 + omdt*omdt/12
      cof = (omdt+1)/(6*cp)
      sgx = 1
      sgy = 1
      isgx = 1
      isgy = 1
!$OMP PARALLEL PRIVATE(i,j,r1,r2,r3,sgx,sgy,isgx,isgy,rhs1,rhs2,
!$OMP*                 rhs3,un1,vn1,wn1,m2sg,m3sg,m4sg,rtu,ac)
!$OMP DO
      do j=jfirst+2,jlast-2
         do i=ifirst+2,ilast-2
            r1 = (-cm*alpham(1,i,j,k-kl)+(4+omdt*omdt)*i6*u(1,i,j,k-kl)+
     *                i6*(1-omdt)*um(1,i,j,k-kl)+ memforce(1,i,j) )/cp
            r2 = (-cm*alpham(2,i,j,k-kl)+(4+omdt*omdt)*i6*u(2,i,j,k-kl)+
     *                i6*(1-omdt)*um(2,i,j,k-kl)+ memforce(2,i,j) )/cp
            r3 = (-cm*alpham(3,i,j,k-kl)+(4+omdt*omdt)*i6*u(3,i,j,k-kl)+
     *                i6*(1-omdt)*um(3,i,j,k-kl)+ memforce(3,i,j) )/cp
            alphap(1,i,j,k-kl) = r1
            alphap(2,i,j,k-kl) = r2
            alphap(3,i,j,k-kl) = r3
            muvebnd(i,j)     = muvebnd(i,j)     + cof*muve(i,j,k)
            lambdavebnd(i,j) = lambdavebnd(i,j) + cof*lave(i,j,k)

            if( usesg.eq.1 )then
               sgx = sgstrx(i)
               sgy = sgstry(j)
               isgy = 1/sgy
               isgx = 1/sgx
            endif

***   tangential derivatives
            rhs1 = 
*** pr
     *   (2*muve(i,j,k)+lave(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
     *          c2*(alphap(1,i+2,j,k)-alphap(1,i-2,j,k)) +
     *          c1*(alphap(1,i+1,j,k)-alphap(1,i-1,j,k))   )*sgx*isgy 
     *  + muve(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(2,i+2,j,k)-alphap(2,i-2,j,k)) +
     *        c1*(alphap(2,i+1,j,k)-alphap(2,i-1,j,k))  ) 
     *  + muve(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(3,i+2,j,k)-alphap(3,i-2,j,k)) +
     *        c1*(alphap(3,i+1,j,k)-alphap(3,i-1,j,k))  )*isgy   
*** qr
     *  +  muve(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(1,i,j+2,k)-alphap(1,i,j-2,k)) +
     *        c1*(alphap(1,i,j+1,k)-alphap(1,i,j-1,k))   )*isgx*sgy
     *  + lave(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(2,i,j+2,k)-alphap(2,i,j-2,k)) +
     *        c1*(alphap(2,i,j+1,k)-alphap(2,i,j-1,k))  )  

*** (v-eq)
            rhs2 = 
*** pr
     *     lave(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(1,i+2,j,k)-alphap(1,i-2,j,k)) +
     *        c1*(alphap(1,i+1,j,k)-alphap(1,i-1,j,k))   ) 
     *  + muve(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(2,i+2,j,k)-alphap(2,i-2,j,k)) +
     *        c1*(alphap(2,i+1,j,k)-alphap(2,i-1,j,k))  )*sgx*isgy
*** qr
     *  +    muve(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(1,i,j+2,k)-alphap(1,i,j-2,k)) +
     *        c1*(alphap(1,i,j+1,k)-alphap(1,i,j-1,k))   ) 
     * + (2*muve(i,j,k)+lave(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(2,i,j+2,k)-alphap(2,i,j-2,k)) +
     *        c1*(alphap(2,i,j+1,k)-alphap(2,i,j-1,k))  )*sgy*isgx 
     *  + muve(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(3,i,j+2,k)-alphap(3,i,j-2,k)) +
     *        c1*(alphap(3,i,j+1,k)-alphap(3,i,j-1,k))   )*isgx 

*** (w-eq)
            rhs3 = 
*** pr
     *      lave(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(1,i+2,j,k)-alphap(1,i-2,j,k)) +
     *        c1*(alphap(1,i+1,j,k)-alphap(1,i-1,j,k))   )*isgy
     *  + muve(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(3,i+2,j,k)-alphap(3,i-2,j,k)) +
     *        c1*(alphap(3,i+1,j,k)-alphap(3,i-1,j,k))  )*sgx*isgy
*** qr 
     *  +   muve(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(3,i,j+2,k)-alphap(3,i,j-2,k)) +
     *        c1*(alphap(3,i,j+1,k)-alphap(3,i,j-1,k))   )*sgy*isgx
     *  + lave(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(alphap(2,i,j+2,k)-alphap(2,i,j-2,k)) +
     *        c1*(alphap(2,i,j+1,k)-alphap(2,i,j-1,k))  )*isgx  
      
*** normal derivatives
           un1 = s(1)*alphap(1,i,j,k)+s(2)*alphap(1,i,j,k+kl)
     *          +s(3)*alphap(1,i,j,k+2*kl)
     *          +s(4)*alphap(1,i,j,k+3*kl) + s(0)*r1
           vn1 = s(1)*alphap(2,i,j,k)+s(2)*alphap(2,i,j,k+kl)
     *          +s(3)*alphap(2,i,j,k+2*kl)
     *          +s(4)*alphap(2,i,j,k+3*kl) + s(0)*r2
           wn1 = s(1)*alphap(3,i,j,k)+s(2)*alphap(3,i,j,k+kl)
     *          +s(3)*alphap(3,i,j,k+2*kl)
     *          +s(4)*alphap(3,i,j,k+3*kl) + s(0)*r3

           m2sg = SQRT(sgx*isgy)
           m3sg = 1/m2sg
           m4sg = isgx*m2sg

           rtu = un1*m2sg*met(2,i,j,k)+vn1*m3sg*met(3,i,j,k)+
     *           wn1*m4sg*met(4,i,j,k)
           ac  = sgx*isgy*met(2,i,j,k)**2+sgy*isgx*met(3,i,j,k)**2+
     *           isgx*isgy*met(4,i,j,k)**2
           rhs1 = rhs1 +(muve(i,j,k)+lave(i,j,k))*rtu*m2sg*met(2,i,j,k)+
     *                                            muve(i,j,k)*ac*un1
           rhs2 = rhs2 +(muve(i,j,k)+lave(i,j,k))*rtu*m3sg*met(3,i,j,k)+
     *                                            muve(i,j,k)*ac*vn1
           rhs3 = rhs3 +(muve(i,j,k)+lave(i,j,k))*rtu*m4sg*met(4,i,j,k)+
     *                                            muve(i,j,k)*ac*wn1
           bforcerhs(1,i,j) = bforcerhs(1,i,j) + rhs1
           bforcerhs(2,i,j) = bforcerhs(2,i,j) + rhs2
           bforcerhs(3,i,j) = bforcerhs(3,i,j) + rhs3
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end      
