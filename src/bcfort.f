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
      subroutine bcfort( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     wind, nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, t,
     *     bforce1, bforce2, bforce3, bforce4, bforce5, bforce6, 
     +     om, ph, cv, curvilinear )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6,6), i, j, k, bccnd(6), w, kl, qq, curvilinear
      real*8 x, y, z, h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 bforce(3,ifirst:ilast,jfirst:jlast)
c note that the numbering of bforce adds one from C (side goes from 1 in Fortran)
      real*8 bforce1(3,*),  bforce2(3,*)
      real*8 bforce3(3,*),  bforce4(3,*)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 ux, vy, wx, wy, uz, vz, wz, t, om, ph, cv

c the boundary window 'wind' is now an input argument

c loop over all sides of the 3-D domain
      do s=1,6
*** dirichlet condition, bccnd=1
*** supergrid condition, bccnd=2
c now assigning the forcing arrays outside of this routine!
        if( bccnd(s).eq.1 .or. bccnd(s).eq.2)then

            qq=1
            if (s.eq.1) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce1 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce1(1,qq)
                    u(2,i,j,k) = bforce1(2,qq)
                    u(3,i,j,k) = bforce1(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.2) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce2 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce2(1,qq)
                    u(2,i,j,k) = bforce2(2,qq)
                    u(3,i,j,k) = bforce2(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.3) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce3 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce3(1,qq)
                    u(2,i,j,k) = bforce3(2,qq)
                    u(3,i,j,k) = bforce3(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.4) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce4 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce4(1,qq)
                    u(2,i,j,k) = bforce4(2,qq)
                    u(3,i,j,k) = bforce4(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.5) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce5 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce5(1,qq)
                    u(2,i,j,k) = bforce5(2,qq)
                    u(3,i,j,k) = bforce5(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.6) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce6 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce6(1,qq)
                    u(2,i,j,k) = bforce6(2,qq)
                    u(3,i,j,k) = bforce6(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            endif

          elseif( bccnd(s).eq.3 )then
*** Periodic condition, bccnd=3
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i+nx,j,k)
                    u(2,i,j,k) = u(2,i+nx,j,k)
                    u(3,i,j,k) = u(3,i+nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i-nx,j,k)
                    u(2,i,j,k) = u(2,i-nx,j,k)
                    u(3,i,j,k) = u(3,i-nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j+ny,k)
                    u(2,i,j,k) = u(2,i,j+ny,k)
                    u(3,i,j,k) = u(3,i,j+ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j-ny,k)
                    u(2,i,j,k) = u(2,i,j-ny,k)
                    u(3,i,j,k) = u(3,i,j-ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k+nz)
                    u(2,i,j,k) = u(2,i,j,k+nz)
                    u(3,i,j,k) = u(3,i,j,k+nz)
                  enddo
                enddo
              enddo
            elseif (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k-nz)
                    u(2,i,j,k) = u(2,i,j,k-nz)
                    u(3,i,j,k) = u(3,i,j,k-nz)
                  enddo
                enddo
              enddo
            endif
         elseif( bccnd(s).eq.0 )then
*** Free surface condition, bccnd=0
            if( s.ne.5 .and. s.ne.6 )then
              write(*,*) 
     *             'ERROR: Free surface condition ',
     +             'not implemented for side ', s
              stop
            endif

c moved the assignment of bforce5/6 into its own routine

**** Do the free surface condition (kl is the direction)
          if( s.eq.5 .and. curvilinear.eq.0 )then
             k = 1
             kl= 1
             do j=jfirst+2,jlast-2
                do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                   qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                   wx = 
     *                  d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *                  d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) 
                   ux = 
     *                  d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *                  d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) 

                   wy = 
     *                  d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *                  d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) 
                   vy = 
     *                  d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *                  d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) 

                   uz = 0
                   vz = 0
                   wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                   do w=1,4
                      uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                      vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                      wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                   enddo
                   u(1,i,j,k-kl)=(-uz-kl*wx+
     *                    kl*h*bforce5(1,qq)/mu(i,j,k))/sbop(0)
                   u(2,i,j,k-kl)=(-vz-kl*wy+
     *                    kl*h*bforce5(2,qq)/mu(i,j,k))/sbop(0)
                   u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *                    kl*h*bforce5(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *                    /sbop(0)
                enddo
             enddo
          elseif( s.eq.6 .and. curvilinear.eq.0 )then
c s=6
             k = nz
             kl= -1
             do j=jfirst+2,jlast-2
                do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                   qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                   wx = 
     *                  d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *                  d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) 
                   ux = 
     *                  d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *                  d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) 

                   wy = 
     *                  d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *                  d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) 
                   vy =
     *                  d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *                  d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) 

                   uz = 0
                   vz = 0
                   wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                   do w=1,4
                      uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                      vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                      wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                   enddo
                   u(1,i,j,k-kl)=(-uz-kl*wx+
     *                  kl*h*bforce6(1,qq)/mu(i,j,k))/sbop(0)
                   u(2,i,j,k-kl)=(-vz-kl*wy+
     *                  kl*h*bforce6(2,qq)/mu(i,j,k))/sbop(0)
                   u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *                  kl*h*bforce6(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *                  /sbop(0)
                enddo
             enddo
          endif
       endif
      enddo
      end


c-----------------------------------------------------------------------
      subroutine bcfortsg( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     wind, nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, t,
     *     bforce1, bforce2, bforce3, bforce4, bforce5, bforce6, 
     +     om, ph, cv, strx, stry )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6,6), i, j, k, bccnd(6), w, kl, qq
      real*8 x, y, z, h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 bforce(3,ifirst:ilast,jfirst:jlast)
c note that the numbering of bforce adds one from C (side goes from 1 in Fortran)
      real*8 bforce1(3,*),  bforce2(3,*)
      real*8 bforce3(3,*),  bforce4(3,*)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 ux, vy, wx, wy, uz, vz, wz, t, om, ph, cv

c the boundary window 'wind' is now an input argument

c loop over all sides of the 3-D domain
      do s=1,6
*** dirichlet condition, bccnd=1
*** supergrid condition, bccnd=2
c now assigning the forcing arrays outside of this routine!
        if( bccnd(s).eq.1 .or. bccnd(s).eq.2)then

            qq=1
            if (s.eq.1) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce1 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce1(1,qq)
                    u(2,i,j,k) = bforce1(2,qq)
                    u(3,i,j,k) = bforce1(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.2) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce2 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce2(1,qq)
                    u(2,i,j,k) = bforce2(2,qq)
                    u(3,i,j,k) = bforce2(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.3) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce3 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce3(1,qq)
                    u(2,i,j,k) = bforce3(2,qq)
                    u(3,i,j,k) = bforce3(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.4) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce4 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce4(1,qq)
                    u(2,i,j,k) = bforce4(2,qq)
                    u(3,i,j,k) = bforce4(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.5) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce5 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce5(1,qq)
                    u(2,i,j,k) = bforce5(2,qq)
                    u(3,i,j,k) = bforce5(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.6) then
c              call TWDIRBDRY( wind(1,s), h, t, om, cv, ph, bforce6 )
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce6(1,qq)
                    u(2,i,j,k) = bforce6(2,qq)
                    u(3,i,j,k) = bforce6(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            endif

          elseif( bccnd(s).eq.3 )then
*** Periodic condition, bccnd=3
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i+nx,j,k)
                    u(2,i,j,k) = u(2,i+nx,j,k)
                    u(3,i,j,k) = u(3,i+nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i-nx,j,k)
                    u(2,i,j,k) = u(2,i-nx,j,k)
                    u(3,i,j,k) = u(3,i-nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j+ny,k)
                    u(2,i,j,k) = u(2,i,j+ny,k)
                    u(3,i,j,k) = u(3,i,j+ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j-ny,k)
                    u(2,i,j,k) = u(2,i,j-ny,k)
                    u(3,i,j,k) = u(3,i,j-ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k+nz)
                    u(2,i,j,k) = u(2,i,j,k+nz)
                    u(3,i,j,k) = u(3,i,j,k+nz)
                  enddo
                enddo
              enddo
            elseif (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k-nz)
                    u(2,i,j,k) = u(2,i,j,k-nz)
                    u(3,i,j,k) = u(3,i,j,k-nz)
                  enddo
                enddo
              enddo
            endif
          elseif( bccnd(s).eq.0 )then
*** Free surface condition, bccnd=0
            if( s.ne.5 .and. s.ne.6 )then
              write(*,*) 
     *             'ERROR: Free surface condition ',
     +             'not implemented for side ', s
              stop
            endif

c moved the assignment of bforce5/6 into its own routine

**** Do the free surface condition (kl is the direction)
          if( s.eq.5 )then
            k = 1
            kl= 1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = strx(i)*(
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) )
                ux = strx(i)*(
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) )

                wy = stry(j)*(
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) )
                vy = stry(j)*(
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) )

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce5(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce5(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce5(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          else
c s=6
            k = nz
            kl= -1

            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
*** Compute wx, wy, ux, vy by centered differences along boundary
                wx = strx(i)*(
     *               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *               d4b*(u(3,i+2,j,k)-u(3,i-2,j,k)) )
                ux = strx(i)*(
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *               d4b*(u(1,i+2,j,k)-u(1,i-2,j,k)) )

                wy = stry(j)*(
     *               d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *               d4b*(u(3,i,j+2,k)-u(3,i,j-2,k)) )
                vy = stry(j)*(
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *               d4b*(u(2,i,j+2,k)-u(2,i,j-2,k)) )

                uz = 0
                vz = 0
                wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
                do w=1,4
                  uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                  vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                  wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
                enddo
                u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce6(1,qq)/mu(i,j,k))
     *               /sbop(0)
                u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce6(2,qq)/mu(i,j,k))
     *               /sbop(0)
                u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *               kl*h*bforce6(3,qq))/(2*mu(i,j,k)+la(i,j,k)))
     *               /sbop(0)
              enddo
            enddo
          endif

        endif
      enddo
      end


c----------------------------------------------------------------------
      subroutine HDIRICHLET5( ifirst, ilast, jfirst, jlast, kfirst,
     *  klast, iafirst, ialast, jafirst, jalast, kafirst, kalast, u )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer iafirst, ialast, jafirst, jalast, kafirst, kalast
      integer i, j, k
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      do k=kalast+1,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               u(1,i,j,k) = 0
               u(2,i,j,k) = 0
               u(3,i,j,k) = 0
            enddo
         enddo
      enddo
      do k=kfirst,klast
         do j=jfirst,jafirst-1
            do i=ifirst,ilast
               u(1,i,j,k) = 0
               u(2,i,j,k) = 0
               u(3,i,j,k) = 0
            enddo
         enddo
         do j=jalast+1,jlast
            do i=ifirst,ilast
               u(1,i,j,k) = 0
               u(2,i,j,k) = 0
               u(3,i,j,k) = 0
            enddo
         enddo
         do j=jfirst,jlast
            do i=ifirst,iafirst-1
               u(1,i,j,k) = 0
               u(2,i,j,k) = 0
               u(3,i,j,k) = 0
            enddo
            do i=ialast+1,ilast
               u(1,i,j,k) = 0
               u(2,i,j,k) = 0
               u(3,i,j,k) = 0
            enddo
         enddo
      enddo
      end

c----------------------------------------------------------------------
      subroutine TWDIRBDRY( wind, h, t, om, cv, ph, bforce, zmin )
      implicit none
      integer wind(6)
      real*8 bforce(3,*), h, t, om, cv, ph, x, y, z, zmin
      integer i, j, k, qq
c
c NOTE: pass in the window for one side, i.e., wind(1,side) in the calling routine
c
*** Twilight forced dirichlet condition
      qq = 1
      do k=wind(5),wind(6)
c need to add zmin to work in a composite grid setting
        z = (k-1)*h + zmin
        do j=wind(3),wind(4)
          y = (j-1)*h
          do i=wind(1),wind(2)
            x = (i-1)*h
            bforce(1,qq) = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph)
            bforce(2,qq) = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph)
            bforce(3,qq) = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t))
            qq = qq+1
          enddo
        enddo
      enddo
      end

c----------------------------------------------------------------------
      subroutine TWDIRBDRYC(ifirst, ilast, jfirst, jlast, kfirst, klast, 
     *                      wind, t, om, cv, ph, bforce, x, y, z )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, wind(6)
      real*8 bforce(3,*), t, om, cv, ph
      real*8 x(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 y(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer i, j, k, qq
c
c NOTE: pass in the window for one side, i.e., wind(1,side) in the calling routine
c
*** Twilight forced dirichlet condition
      qq = 1
      do k=wind(5),wind(6)
        do j=wind(3),wind(4)
          do i=wind(1),wind(2)

            bforce(1,qq) = sin(om*(x(i,j,k)-cv*t))*sin(om*y(i,j,k)+ph)
     *                                              *sin(om*z(i,j,k)+ph)
            bforce(2,qq) = sin(om*x(i,j,k)+ph)*sin(om*(y(i,j,k)-cv*t))
     *                                              *sin(om*z(i,j,k)+ph)
            bforce(3,qq) = sin(om*x(i,j,k)+ph)*sin(om*y(i,j,k)+ph)
     *                                          *sin(om*(z(i,j,k)-cv*t))
            qq = qq+1
          enddo
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZ( ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, kz, t, omega, c, phase, bforce, mu, lambda, zmin )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, attenuation
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer i, j, kz
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z, zmin
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase

      doubleprecision forces(3)
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t19
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t28
      doubleprecision t29
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t43
      doubleprecision t44
      doubleprecision t49
      doubleprecision t60
      doubleprecision t62
      doubleprecision t65

      z = (kz-1)*h + zmin
c the do loops should span jfirst,jlast and ifirst,ilast
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x=(i-1)*h
        t13 = mu(i,j,kz)
        t15 = omega*x+phase
        t16 = cos(t15)
        t19 = omega*y+phase
        t20 = sin(t19)
        t21 = c*t
        t23 = omega*(z-t21)
        t24 = sin(t23)
        t28 = omega*(x-t21)
        t29 = sin(t28)
        t32 = omega*z+phase
        t33 = cos(t32)
        t34 = t33*omega
        forces(1) = t13*(t16*omega*t20*t24+t29*t20*t34)
        t37 = sin(t15)
        t38 = cos(t19)
        t43 = omega*(y-t21)
        t44 = sin(t43)
        forces(2) = t13*(t37*t38*omega*t24+t37*t44*t34)
        t49 = cos(t23)
        t60 = cos(t28)
        t62 = sin(t32)
        t65 = cos(t43)
        forces(3) = 2*t13*t37*t20*t49*omega+lambda(i,j,kz)*(t6
     #0*omega*t20*t62+t37*t65*omega*t62+t37*t20*t49*omega)
        bforce(1,i,j) = forces(1)
        bforce(2,i,j) = forces(2)
        bforce(3,i,j) = forces(3)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZATT( ifirst, ilast, jfirst, jlast, kfirst,
     +   klast, h, kz, t, omega, c, phase, bforce, mua, lambdaa, zmin )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer i, j, kz
      doubleprecision mua(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambdaa(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z, zmin
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
      doubleprecision t2, t3, t6, t7, t8, t11, t12, t17, t18, t27, t30
      doubleprecision t20, t31, t35, t40, t45, t50, t52, t54
      doubleprecision t16, t23, t24, t34
      doubleprecision forces(3)

      z = (kz-1)*h + zmin
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x=(i-1)*h
           t2 = omega*x+phase
           t3 = sin(t2)
           t6 = omega*y+phase
           t7 = cos(t6)
           t8 = c*t
           t11 = -omega*(z-t8)-phase
           t12 = sin(t11)
           t16 = omega*(x-t8)
           t17 = -t16-phase
           t18 = cos(t17)
           t20 = t12*omega
           forces(1) = mua(i,j,kz)*(t3*omega*t7*t12+t18*t3*t20)
           t23 = cos(t2)
           t24 = sin(t6)
           t27 = sin(t16)
           t30 = -omega*(y-t8)-phase
           t31 = cos(t30)
           t34 = omega*z+phase
           t35 = sin(t34)
           forces(2) = mua(i,j,kz)*(t23*t24*t20-t27*t31*t35*omega)
           t40 = cos(t11)
           t45 = sin(t17)
           t50 = t40*omega
           t52 = sin(t30)
           t54 = cos(t34)
           forces(3) = 2d0*mua(i,j,kz)*t23*t7*t40*omega+lambdaa(i,j,kz)*
     *     (t45*omega*t3*t40+t18*t23*t50+t27*t52*omega*t54+t23*t7*t50)
           bforce(1,i,j) = bforce(1,i,j) - forces(1)
           bforce(2,i,j) = bforce(2,i,j) - forces(2)
           bforce(3,i,j) = bforce(3,i,j) - forces(3)
        enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZSGSTR( ifirst, ilast, jfirst, jlast, kfirst, 
     *                  klast, h, kz, t, om, c, ph, omstrx, omstry,
     *                  bforce, mu, lambda, zmin )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer i, j, kz
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z, zmin
      doubleprecision t
      doubleprecision om
      doubleprecision c
      doubleprecision ph
      doubleprecision omstrx
      doubleprecision omstry

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t11
      doubleprecision t12
      doubleprecision t15
      doubleprecision t17
      doubleprecision t19
      doubleprecision t20
      doubleprecision t22
      doubleprecision t24
      doubleprecision t25
      doubleprecision t29
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t36
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t46
      doubleprecision t51
      doubleprecision t53
      doubleprecision t56
      doubleprecision t6
      doubleprecision t7

      z = (kz-1)*h + zmin
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x=(i-1)*h
            t1 = c*t
            t3 = om*(x-t1)
            t4 = sin(t3)
            t6 = om*y+ph
            t7 = sin(t6)
            t10 = om*z+ph
            t11 = cos(t10)
            t12 = t11*om
            t15 = sin(omstrx*x)
            t17 = 1+t15/2
            t19 = om*x+ph
            t20 = cos(t19)
            t22 = om*t7
            t24 = om*(z-t1)
            t25 = sin(t24)
            forces(1) = mu(i,j,kz)*(t4*t7*t12+t17*t20*t22*t25)
            t29 = sin(t19)
            t31 = om*(y-t1)
            t32 = sin(t31)
            t36 = sin(omstry*y)
            t39 = (1+t36/2)*t29
            t40 = cos(t6)
            forces(2) = mu(i,j,kz)*(t29*t32*t12+t39*t40*om*t25)
            t46 = cos(t24)
            t51 = cos(t3)
            t53 = sin(t10)
            t56 = cos(t31)
            forces(3) = 2*mu(i,j,kz)*t29*t7*t46*om+lambda(i,j,kz)*
     #(t17*t51*t22*t53+t39*t56*om*t53+t29*t7*t46*om)
            bforce(1,i,j) = forces(1)
            bforce(2,i,j) = forces(2)
            bforce(3,i,j) = forces(3)
         enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZSGSTRATT( ifirst, ilast, jfirst, jlast, 
     *       kfirst, klast, h, kz, t, omega, c, phase, omstrx, omstry,
     *       bforce, mua, lambdaa, zmin )

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer i, j, kz
      doubleprecision mua(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambdaa(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z, zmin
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
      doubleprecision omstrx
      doubleprecision omstry
     
      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t12
      doubleprecision t13
      doubleprecision t17
      doubleprecision t19
      doubleprecision t22
      doubleprecision t23
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t35
      doubleprecision t36
      doubleprecision t4
      doubleprecision t40
      doubleprecision t42
      doubleprecision t43
      doubleprecision t45
      doubleprecision t5
      doubleprecision t51
      doubleprecision t56
      doubleprecision t61
      doubleprecision t66
      doubleprecision t68
      doubleprecision t7
      doubleprecision t8
      z = (kz-1)*h + zmin
      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x=(i-1)*h
            t1 = c*t
            t3 = omega*(x-t1)
            t4 = -t3-phase
            t5 = cos(t4)
            t7 = omega*x+phase
            t8 = sin(t7)
            t12 = -omega*(z-t1)-phase
            t13 = sin(t12)
            t17 = sin(omstrx*x)
            t19 = 1+t17/2
            t22 = omega*y+phase
            t23 = cos(t22)
            forces(1) = mua(i,j,kz)*(t5*t8*t13*omega+t19*t8*
     *                          omega*t23*t13)
            t28 = sin(t3)
            t31 = -omega*(y-t1)-phase
            t32 = cos(t31)
            t35 = omega*z+phase
            t36 = sin(t35)
            t40 = sin(omstry*y)
            t42 = 1+t40/2
            t43 = cos(t7)
            t45 = sin(t22)
            forces(2) = mua(i,j,kz)*(-t28*t32*t36*omega+t42*t43*t45*
     *                           omega*t13)
            t51 = cos(t12)
            t56 = sin(t4)
            t61 = omega*t51
            t66 = sin(t31)
            t68 = cos(t35)
            forces(3) = 2*mua(i,j,kz)*t43*t23*t51*omega+lambdaa(i,j,kz)
     #*(t19*(t56*omega*t8*t51+t5*t43*t61)+t42*t28*t66*omega*t68+
     *          t43*t23*t61)
            bforce(1,i,j) = bforce(1,i,j) - forces(1)
            bforce(2,i,j) = bforce(2,i,j) - forces(2)
            bforce(3,i,j) = bforce(3,i,j) - forces(3)
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWSTENSOR( ifirst, ilast, jfirst, jlast, kfirst, 
     *                  klast, kz, t, om, c, ph, xx, yy, zz,
     *                  tau, mu, lambda )
***********************************************************************
***  Stress tensor ordered as tau(1) = t_{xx}, tau(2) = t_{xy}
***  tau(3) = t_{xz}, tau(4) = t_{yy}, tau(5)=t_{yz}, tau(6)=t_{zz}
***
***********************************************************************

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer i, j, kz
      real*8 muu, lambdaa, x, y, z, t, om, c, ph
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision t1
      doubleprecision t11
      doubleprecision t12
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t26
      doubleprecision t3
      doubleprecision t30
      doubleprecision t31
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t4
      doubleprecision t41
      doubleprecision t42
      doubleprecision t46
      doubleprecision t50
      doubleprecision t51
      doubleprecision t54
      doubleprecision t7
      doubleprecision t8

      do j=jfirst,jlast
         do i=ifirst,ilast
            x = xx(i,j,kz)
            y = yy(i,j,kz)
            z = zz(i,j,kz)
            muu = mu(i,j,kz)
            lambdaa = lambda(i,j,kz)
            t1 = c*t
            t3 = om*(x-t1)
            t4 = cos(t3)
            t7 = om*y+ph
            t8 = sin(t7)
            t11 = om*z+ph
            t12 = sin(t11)
            t20 = om*x+ph
            t21 = sin(t20)
            t23 = om*(y-t1)
            t24 = cos(t23)
            t26 = om*t12
            t30 = om*(z-t1)
            t31 = cos(t30)
            t35 = lambdaa*(t4*om*t8*t12+t21*t24*t26+t21*t8*t31*om)
            tau(1,i,j) = 2*muu*t4*om*t8*t12+t35
            t36 = cos(t20)
            t37 = t36*om
            t38 = sin(t23)
            t41 = sin(t3)
            t42 = cos(t7)
            tau(2,i,j) = muu*(t37*t38*t12+t41*t42*t26)
            t46 = sin(t30)
            t50 = cos(t11)
            t51 = t50*om
            tau(3,i,j) = muu*(t37*t8*t46+t41*t8*t51)
            t54 = muu*t21
            tau(4,i,j) = 2*t54*t24*om*t12+t35
            tau(5,i,j) = muu*(t21*t42*om*t46+t21*t38*t51)
            tau(6,i,j) = 2*t54*t8*t31*om+t35
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWSTENSORSG( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, kz, t, om, c, ph, xx, yy, zz,
     *     tau, mu, lambda, omstrx, omstry )

***********************************************************************
***  Stress tensor ordered as tau(1) = t_{xx}, tau(2) = t_{xy}
***  tau(3) = t_{xz}, tau(4) = t_{yy}, tau(5)=t_{yz}, tau(6)=t_{zz}
***
***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer i, j, kz
      real*8 muu, lambdaa, x, y, z, t, om, c, ph
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision omstrx
      doubleprecision omstry

      doubleprecision t12
      doubleprecision t13
      doubleprecision t14
      doubleprecision t16
      doubleprecision t17
      doubleprecision t18
      doubleprecision t2
      doubleprecision t24
      doubleprecision t26
      doubleprecision t28
      doubleprecision t29
      doubleprecision t30
      doubleprecision t32
      doubleprecision t33
      doubleprecision t35
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t44
      doubleprecision t45
      doubleprecision t46
      doubleprecision t47
      doubleprecision t51
      doubleprecision t53
      doubleprecision t54
      doubleprecision t59
      doubleprecision t6
      doubleprecision t60
      doubleprecision t62
      doubleprecision t8
      doubleprecision t9

      do j=jfirst,jlast
         do i=ifirst,ilast
            x = xx(i,j,kz)
            y = yy(i,j,kz)
            z = zz(i,j,kz)
            muu = mu(i,j,kz)
            lambdaa = lambda(i,j,kz)
            t2 = sin(omstrx*x)
            t4 = 1+t2/2
            t6 = c*t
            t8 = om*(x-t6)
            t9 = cos(t8)
            t12 = om*y+ph
            t13 = sin(t12)
            t14 = om*t13
            t16 = om*z+ph
            t17 = sin(t16)
            t18 = t14*t17
            t24 = sin(omstry*y)
            t26 = 1+t24/2
            t28 = om*x+ph
            t29 = sin(t28)
            t30 = t26*t29
            t32 = om*(y-t6)
            t33 = cos(t32)
            t35 = t33*om*t17
            t39 = om*(z-t6)
            t40 = cos(t39)
            t44 = lambdaa*(t4*t9*t18+t30*t35+t29*t13*t40*om)
            tau(1,i,j) = 2*muu*t4*t9*t18+t44
            t45 = cos(t28)
            t46 = t4*t45
            t47 = sin(t32)
            t51 = sin(t8)
            t53 = cos(t12)
            t54 = t53*om
            tau(2,i,j) = muu*(t46*om*t47*t17+t26*t51*t54*t17)
            t59 = cos(t16)
            t60 = t59*om
            t62 = sin(t39)
            tau(3,i,j) = muu*(t51*t13*t60+t46*t14*t62)
            tau(4,i,j) = 2*muu*t26*t29*t35+t44
            tau(5,i,j) = muu*(t29*t47*t60+t30*t54*t62)
            tau(6,i,j) = 2*muu*t29*t13*t40*om+t44
         enddo
      enddo
      end
         
c-----------------------------------------------------------------------
      subroutine TWSTENSORATT( ifirst, ilast, jfirst, jlast, kfirst, 
     *                  klast, kz, t, omega, c, phase, xx, yy, zz,
     *                  tau, mu, lambda )
***********************************************************************
***  Stress tensor ordered as tau(1) = t_{xx}, tau(2) = t_{xy}
***  tau(3) = t_{xz}, tau(4) = t_{yy}, tau(5)=t_{yz}, tau(6)=t_{zz}
***
***********************************************************************

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer i, j, kz
      real*8 muu, lambdaa, x, y, z, t, omega, c, phase
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision stensor(6)
      doubleprecision t1
      doubleprecision t12
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t17
      doubleprecision t19
      doubleprecision t20
      doubleprecision t24
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t36
      doubleprecision t37
      doubleprecision t4
      doubleprecision t41
      doubleprecision t42
      doubleprecision t44
      doubleprecision t48
      doubleprecision t49
      doubleprecision t5
      doubleprecision t61
      doubleprecision t64
      doubleprecision t8
      doubleprecision t9

      do j=jfirst,jlast
         do i=ifirst,ilast
            x = xx(i,j,kz)
            y = yy(i,j,kz)
            z = zz(i,j,kz)
            muu     = mu(i,j,kz)
            lambdaa = lambda(i,j,kz)
        t1 = c*t
        t3 = omega*(x-t1)
        t4 = -t3-phase
        t5 = sin(t4)
        t8 = omega*x+phase
        t9 = sin(t8)
        t12 = -omega*(z-t1)-phase
        t13 = cos(t12)
        t15 = t5*omega*t9*t13
        t16 = cos(t4)
        t17 = cos(t8)
        t19 = omega*t13
        t20 = t16*t17*t19
        t24 = sin(t3)
        t27 = -omega*(y-t1)-phase
        t28 = sin(t27)
        t31 = omega*z+phase
        t32 = cos(t31)
        t36 = omega*y+phase
        t37 = cos(t36)
        t41 = lambdaa*(t15+t20+t24*t28*omega*t32+t17*t37*t19)
        stensor(1) = 2*muu*(t15+t20)+t41
        t42 = cos(t3)
        t44 = cos(t27)
        stensor(2) = muu*t42*omega*t44*t32
        t48 = sin(t12)
        t49 = t48*omega
        stensor(3) = muu*(t16*t9*t49+t9*omega*t37*t48)
        stensor(4) = 2*muu*t24*t28*omega*t32+t41
        t61 = sin(t31)
        t64 = sin(t36)
        stensor(5) = muu*(-t24*t44*t61*omega+t17*t64*t49)
        stensor(6) = 2*muu*t17*t37*t13*omega+t41
        tau(1,i,j) = stensor(1)
        tau(2,i,j) = stensor(2)
        tau(3,i,j) = stensor(3)
        tau(4,i,j) = stensor(4)
        tau(5,i,j) = stensor(5)
        tau(6,i,j) = stensor(6)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine TWSTENSORSGATT( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, kz, t, omega, c, phase, xx, yy, zz,
     *     tau, mu, lambda, omstrx, omstry )

***********************************************************************
***  Stress tensor ordered as tau(1) = t_{xx}, tau(2) = t_{xy}
***  tau(3) = t_{xz}, tau(4) = t_{yy}, tau(5)=t_{yz}, tau(6)=t_{zz}
***
***********************************************************************
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 tau(6,ifirst:ilast,jfirst:jlast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      integer i, j, kz
      real*8 muu, lambdaa, x, y, z, t, omega, c, phase
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision omstrx
      doubleprecision omstry
      doubleprecision stensor(6)
      doubleprecision t10
      doubleprecision t13
      doubleprecision t14
      doubleprecision t17
      doubleprecision t18
      doubleprecision t2
      doubleprecision t21
      doubleprecision t22
      doubleprecision t24
      doubleprecision t26
      doubleprecision t31
      doubleprecision t33
      doubleprecision t34
      doubleprecision t38
      doubleprecision t39
      doubleprecision t4
      doubleprecision t42
      doubleprecision t43
      doubleprecision t44
      doubleprecision t47
      doubleprecision t48
      doubleprecision t5
      doubleprecision t52
      doubleprecision t53
      doubleprecision t55
      doubleprecision t59
      doubleprecision t6
      doubleprecision t72
      doubleprecision t76
      doubleprecision t8
      doubleprecision t9

      do j=jfirst,jlast
         do i=ifirst,ilast
            x = xx(i,j,kz)
            y = yy(i,j,kz)
            z = zz(i,j,kz)
            muu = mu(i,j,kz)
            lambdaa = lambda(i,j,kz)
        t2 = sin(omstrx*x)
        t4 = 1+t2/2
        t5 = muu*t4
        t6 = c*t
        t8 = omega*(x-t6)
        t9 = -t8-phase
        t10 = sin(t9)
        t13 = omega*x+phase
        t14 = sin(t13)
        t17 = -omega*(z-t6)-phase
        t18 = cos(t17)
        t21 = cos(t9)
        t22 = cos(t13)
        t24 = omega*t18
        t26 = t10*omega*t14*t18+t21*t22*t24
        t31 = sin(omstry*y)
        t33 = 1+t31/2
        t34 = sin(t8)
        t38 = -omega*(y-t6)-phase
        t39 = sin(t38)
        t42 = omega*z+phase
        t43 = cos(t42)
        t44 = t39*omega*t43
        t47 = omega*y+phase
        t48 = cos(t47)
        t52 = lambdaa*(t4*t26+t33*t34*t44+t22*t48*t24)
        stensor(1) = 2*t5*t26+t52
        t53 = cos(t8)
        t55 = cos(t38)
        stensor(2) = t5*t53*omega*t55*t43
        t59 = sin(t17)
        stensor(3) = muu*(t21*t14*t59*omega+t4*t14*omega*t48*t59)
        stensor(4) = 2*muu*t33*t34*t44+t52
        t72 = sin(t42)
        t76 = sin(t47)
        stensor(5) = muu*(-t34*t55*t72*omega+t33*t22*t76*omega*t59)
        stensor(6) = 2*muu*t22*t48*t18*omega+t52
        tau(1,i,j) = stensor(1)
        tau(2,i,j) = stensor(2)
        tau(3,i,j) = stensor(3)
        tau(4,i,j) = stensor(4)
        tau(5,i,j) = stensor(5)
        tau(6,i,j) = stensor(6)
      enddo
      enddo
      end
