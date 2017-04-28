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
      subroutine ADDGRADRHO(ifirst, ilast, jfirst, jlast, kfirst, klast,
     *    ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact,
     *    kap, kapacc, um, u, up, uacc, grho, dt, h, onesided )
     * bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer ifirstact, ilastact, jfirstact, jlastact, kfirstact
      integer klastact, onesided(6)
      integer i, j, k
      real*8  kap(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  kapacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  grho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  dt, h, idt, dt2o12, h3, wgh(4), normfact

      idt = 1/dt
      dt2o12 = dt*dt/12
      h3 = h*h*h
      wgh(1) = 17d0/48
      wgh(2) = 59d0/48
      wgh(3) = 43d0/48
      wgh(4) = 49d0/48
      do k=kfirstact,klastact
         do j=jfirstact,jlastact
            do i=ifirstact,ilastact
               normfact = h3
               if( k.le.4 .and. onesided(5).eq.1 )then
                  normfact = h3*wgh(k)
               endif
               grho(i,j,k) = grho(i,j,k) + (
     *     kap(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k))*idt*idt
     *  + dt2o12*kapacc(1,i,j,k)*uacc(1,i,j,k) +
     *     kap(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k))*idt*idt
     *  + dt2o12*kapacc(2,i,j,k)*uacc(2,i,j,k) +
     *     kap(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))*idt*idt
     *  + dt2o12*kapacc(3,i,j,k)*uacc(3,i,j,k) )*normfact
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ADDGRADRHOC(ifirst,ilast, jfirst, jlast, kfirst, klast,
     *    ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact,
     *    kap, kapacc, um, u, up, uacc, grho, dt, jac, onesided )
     * bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer ifirstact, ilastact, jfirstact, jlastact, kfirstact
      integer klastact, onesided(6)
      integer i, j, k
      real*8  kap(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  kapacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  grho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  dt, idt, dt2o12, wgh(4), normfact

      idt = 1/dt
      dt2o12 = dt*dt/12
      wgh(1) = 17d0/48
      wgh(2) = 59d0/48
      wgh(3) = 43d0/48
      wgh(4) = 49d0/48
      do k=kfirstact,klastact
         do j=jfirstact,jlastact
            do i=ifirstact,ilastact
               normfact = jac(i,j,k)
               if( k.le.4 .and. onesided(5).eq.1 )then
                  normfact = normfact*wgh(k)
               endif
               grho(i,j,k) = grho(i,j,k) + (
     *     kap(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k))*idt*idt
     *  + dt2o12*kapacc(1,i,j,k)*uacc(1,i,j,k) +
     *     kap(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k))*idt*idt
     *  + dt2o12*kapacc(2,i,j,k)*uacc(2,i,j,k) +
     *     kap(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))*idt*idt
     *  + dt2o12*kapacc(3,i,j,k)*uacc(3,i,j,k) )*normfact
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine ADDGRADMULA( ifirst, ilast, jfirst, jlast, kfirst,
     *  klast, ifirstact, ilastact, jfirstact, jlastact, kfirstact,
     *   klastact, kap, kapacc, u, uacc, gmu, glambda,
     *    dt, h, onesided, nb, wb, bop )
     * bind(c)
      implicit none
      real*8 d4a, d4b, c6, c8, al1, al2, al3, al4
      parameter( d4a=2d0/3, d4b=-1d0/12, c6=1d0/18, c8=1d0/144 )
      parameter( al1 = 181507d0/1719312, al2 = -1441d0/39984 )
      parameter( al3 = -2593d0/151704,   al4 = 11d0/3528 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer ifirstact, ilastact, jfirstact, jlastact, kfirstact
      integer klastact, kstart, onesided(6), nb, wb
      integer i, j, k, m
      real*8  kap(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  kapacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8  um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8  up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  gmu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  glambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  bop(nb,wb)
      real*8  dt, h, dt2o12, h3, wgh(4), normfact
      real*8  dux, dvy, dwz, dkx, dly, dmz, duax, dvay, dwaz, dkax 
      real*8  dlay, dmaz, stuxy, stkxy, stuaxy, stkaxy
      real*8  stuxz, stkxz, stuaxz, stkaxz, stuyz, stkyz, stuayz, stkayz
      real*8  d3up, d3um, d3kp, d3km, pd, d3uap, d3uam, d3kap, d3kam
      real*8  w8(4), w6m(4), w6p(4)
      real*8  ih2

      h3  = h*h*h
      ih2 = 1/(h*h)
      dt2o12 = dt*dt/12
      wgh(1) = 17d0/48
      wgh(2) = 59d0/48
      wgh(3) = 43d0/48
      wgh(4) = 49d0/48

      kstart = kfirstact
      if( kfirstact .le. 4 .and. onesided(5).eq.1 )then
         kstart = 5

         w8(1) = 0
         w8(2) = 0
         w8(3) = 1
         w8(4) = 1

         w6m(1)= 0
         w6p(1)= 0

         w6m(2)= 0
         w6p(2)= al1

         w6m(3)= al1
         w6p(3)= al1+al2

         w6m(4)= al1+al2
         w6p(4)= al1+al2+al3

         do k=kfirstact,4
            do j=jfirstact,jlastact
               do i=ifirstact,ilastact
                  normfact = h3*wgh(k)
*** Diagonal terms
               dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))
               dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))
c               dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
c     *               d4a*(u(3,i,j,k+1)-u(3,i,j,k-1))
               dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
     *               d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k))
               dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
     *               d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k))
c               dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
c     *               d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1))

               duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
     *               d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k))
               dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
     *               d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k))
c               dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
c     *               d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1))
               dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
     *               d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k))
               dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
     *               d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k))
c               dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
c     *               d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1))
               dwz =0
               dmz =0
               dwaz=0
               dmaz=0
               do m=1,wb
                  dwz = dwz + bop(k,m)*u(3,i,j,m)
                  dmz = dmz + bop(k,m)*kap(3,i,j,m)
                  dwaz= dwaz+ bop(k,m)*uacc(3,i,j,m)
                  dmaz= dmaz+ bop(k,m)*kapacc(3,i,j,m)
               enddo

               glambda(i,j,k) = glambda(i,j,k) + (
     *           (dux+dvy+dwz)*(dkx+dly+dmz) +
     *           dt2o12*( (duax+dvay+dwaz)*(dkx+dly+dmz) +
     *                    (dux+dvy+dwz)*(dkax+dlay+dmaz) )
     *                  )*normfact*ih2
               gmu(i,j,k) = gmu(i,j,k) + 2*(
     *          dux*dkx+dvy*dly+dwz*dmz +
     *           dt2o12*( (duax*dkx+dvay*dly+dwaz*dmz)+
     *                    (dux*dkax+dvy*dlay+dwz*dmaz) )
     *                     )*normfact*ih2
               
*** Off diagonal stresses
***  xy
               stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
     *                 d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
     *                 d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
     *                 d4a*(u(1,i,j+1,k)-u(1,i,j-1,k))
               stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
     *                 d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
     *                 d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
     *                 d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k))
               stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
     *                 d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
     *                 d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
     *                 d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k))
               stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
     *                 d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
     *                 d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
     *                 d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k))
               gmu(i,j,k)= gmu(i,j,k) + (stuxy*stkxy + 
     *            dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*normfact*ih2

***  xz
               stuxz = d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
     *                 d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))
c     *                 d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
c     *                 d4a*(u(1,i,j,k+1)-u(1,i,j,k-1))
               stkxz = d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
     *                 d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k))
c     *                 d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
c     *                 d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1))
               stuaxz = d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
     *                 d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k))
c     *                 d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
c     *                 d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1))
               stkaxz = d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
     *                 d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k))
c     *                 d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
c     *                 d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1))
               do m=1,wb
                  stuxz = stuxz + bop(k,m)*u(1,i,j,m)
                  stkxz = stkxz + bop(k,m)*kap(1,i,j,m)
                  stuaxz= stuaxz+ bop(k,m)*uacc(1,i,j,m)
                  stkaxz= stkaxz+ bop(k,m)*kapacc(1,i,j,m)
               enddo
               gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
     *            dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact*ih2
***  yz
               stuyz = d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
     *                 d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))
c     *                 d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
c     *                 d4a*(u(2,i,j,k+1)-u(2,i,j,k-1))
               stkyz = d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
     *                 d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k))
c     *                 d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
c     *                 d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1))
               stuayz = d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
     *                 d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k))
c     *                 d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
c     *                 d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1))
               stkayz = d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
     *                 d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k))
c     *                 d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
c     *                 d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1))
               do m=1,wb
                  stuyz = stuyz + bop(k,m)*u(2,i,j,m)
                  stkyz = stkyz + bop(k,m)*kap(2,i,j,m)
                  stuayz= stuayz+ bop(k,m)*uacc(2,i,j,m)
                  stkayz= stkayz+ bop(k,m)*kapacc(2,i,j,m)
               enddo
               gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
     *           dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact*ih2

*** Pos. def extra terms
*** x-direction
               d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
     *              3*u(1,i,j,k) -   u(1,i-1,j,k)
               d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
     *             3*u(1,i-1,j,k)-  u(1,i-2,j,k)
               d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
     *              3*kap(1,i,j,k) -   kap(1,i-1,j,k)
               d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
     *              3*kap(1,i-1,j,k) - kap(1,i-2,j,k)
               d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
     *               3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k)
               d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k)
               d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
     *               3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k)
               d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k)
               pd = (c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *      dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2*wgh(k)
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd
               d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
     *               u(2,i-1,j,k)
               d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
     *               u(2,i-2,j,k)
               d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
     *           3*kap(2,i,j,k)-kap(2,i-1,j,k)
               d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
     *           3*kap(2,i-1,j,k)-kap(2,i-2,j,k)
               d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
     *              3*uacc(2,i,j,k)- uacc(2,i-1,j,k)
               d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k)
               d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
     *           3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k)
               d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
     *           3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
     *               u(3,i-1,j,k)
               d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
     *               u(3,i-2,j,k)
               d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
     *           3*kap(3,i,j,k)-kap(3,i-1,j,k)
               d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
     *           3*kap(3,i-1,j,k)-kap(3,i-2,j,k)
               d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
     *            3*uacc(3,i,j,k)- uacc(3,i-1,j,k)
               d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
     *             3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k)
               d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
     *           3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k)
               d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
     *           3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp + d3um*d3km +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

*** y-direction
               d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
     *               u(1,i,j-1,k)
               d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
     *               u(1,i,j-2,k)
               d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
     *              3*kap(1,i,j,k)-kap(1,i,j-1,k)
               d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
     *              3*kap(1,i,j-1,k)-kap(1,i,j-2,k)
               d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
     *               3*uacc(1,i,j,k)- uacc(1,i,j-1,k)
               d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k)
               d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
     *               3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k)
               d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
     *               u(2,i,j-1,k)
               d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
     *               u(2,i,j-2,k)
               d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
     *              3*kap(2,i,j,k)-kap(2,i,j-1,k)
               d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
     *              3*kap(2,i,j-1,k)-kap(2,i,j-2,k)

               d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
     *               3*uacc(2,i,j,k) -   uacc(2,i,j-1,k)
               d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i,j-1,k)-  uacc(2,i,j-2,k)

               d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
     *               3*kapacc(2,i,j,k)-    kapacc(2,i,j-1,k)
               d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
     *               3*kapacc(2,i,j-1,k)-  kapacc(2,i,j-2,k)

               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh(k)
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k) = gmu(i,j,k) + 2*pd

               d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
     *                u(3,i,j-1,k)
               d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
     *               u(3,i,j-2,k)
               d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
     *              3*kap(3,i,j,k)-kap(3,i,j-1,k)
               d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
     *              3*kap(3,i,j-1,k)-kap(3,i,j-2,k)
               d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
     *               3*uacc(3,i,j,k)- uacc(3,i,j-1,k)
               d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
     *               3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k)
               d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
     *               3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k)
               d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                       (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

*** z-direction
               if( k.ge.2 )then
                  d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
     *                u(1,i,j,k-1)
                  d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
     *                u(1,i,j,k-2)
                  d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
     *                 3*kap(1,i,j,k)-kap(1,i,j,k-1)
                  d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
     *                 3*kap(1,i,j,k-1)-kap(1,i,j,k-2)
                  d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
     *                  3*uacc(1,i,j,k)-uacc(1,i,j,k-1)
                  d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
     *                 3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2)
                  d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
     *                  3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1)
                  d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
     *                  3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2)
                  pd = ( ( 0.5d0*( w6p(k)*d3up*d3kp+w6m(k)*d3um*d3km  +
     *        dt2o12*(w6p(k)*d3up*d3kap + w6m(k)*d3um*d3kam+ 
     *                w6p(k)*d3uap*d3kp + w6m(k)*d3uam*d3km)))
     *            + w8(k)*c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3kp+
     *               dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2
                     gmu(i,j,k+1) = gmu(i,j,k+1) + pd
                  endif
                  d3up = u(2,i,j,k+2)-3*u(2,i,j,k+1)+
     *              3*u(2,i,j,k)- u(2,i,j,k-1)
                  d3um =u(2,i,j,k+1)-3*u(2,i,j,k)+
     *             3*u(2,i,j,k-1)- u(2,i,j,k-2)
                  d3kp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
     *              3*kap(2,i,j,k)-kap(2,i,j,k-1)
                  d3km = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
     *              3*kap(2,i,j,k-1)-kap(2,i,j,k-2)
                  d3uap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
     *              3*uacc(2,i,j,k)- uacc(2,i,j,k-1)
                  d3uam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
     *             3*uacc(2,i,j,k-1)- uacc(2,i,j,k-2)
                  d3kap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
     *              3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1)
                  d3kam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
     *              3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2)
                  pd = (( 0.5d0*( w6p(k)*d3up*d3kp+w6m(k)*d3um*d3km  +
     *        dt2o12*(w6p(k)*d3up*d3kap + w6m(k)*d3um*d3kam+
     *                w6p(k)*d3uap*d3kp + w6m(k)*d3uam*d3km)))
     *            + w8(k)*c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3kp+
     *               dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2
                     gmu(i,j,k+1) = gmu(i,j,k+1) + pd
                  endif

                  d3up = u(3,i,j,k+2) - 3*u(3,i,j,k+1) +
     *                 3*u(3,i,j,k)   -   u(3,i,j,k-1)
                  d3um = u(3,i,j,k+1) - 3*u(3,i,j,k) +
     *                 3*u(3,i,j,k-1) -   u(3,i,j,k-2)
                  d3kp = kap(3,i,j,k+2) - 3*kap(3,i,j,k+1) +
     *                 3*kap(3,i,j,k)   -   kap(3,i,j,k-1)
                  d3km = kap(3,i,j,k+1) - 3*kap(3,i,j,k  ) +
     *                 3*kap(3,i,j,k-1) -   kap(3,i,j,k-2)

                  d3uap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
     *                  3*uacc(3,i,j,k)  -  uacc(3,i,j,k-1)
                  d3uam = uacc(3,i,j,k+1)-3*uacc(3,i,j,k) +
     *                  3*uacc(3,i,j,k-1)-  uacc(3,i,j,k-2)

                  d3kap = kapacc(3,i,j,k+2)- 3*kapacc(3,i,j,k+1)+
     *                  3*kapacc(3,i,j,k)  -   kapacc(3,i,j,k-1)
                  d3kam = kapacc(3,i,j,k+1) -3*kapacc(3,i,j,k)+
     *                  3*kapacc(3,i,j,k-1) -  kapacc(3,i,j,k-2)

                  pd = (0.5d0*( w6p(k)*d3up*d3kp  + w6m(k)*d3um*d3km  +
     *              dt2o12*(    w6p(k)*d3up*d3kap + w6m(k)*d3um*d3kam + 
     *                          w6p(k)*d3uap*d3kp + w6m(k)*d3uam*d3km) )
     *            + w8(k)*c8*( (d3up-d3um)*(d3kp-d3km) +
     *                  dt2o12*(  (d3uap-d3uam)*(d3kp-d3km) +
     *                            (d3up-d3um)*(d3kap-d3kam)   )))*h3*ih2
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  gmu(i,j,k) = gmu(i,j,k) + 2*pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3kp+
     *               dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )*h3*ih2
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  2*pd
                  endif
                  endif
               enddo
            enddo
         enddo
      endif
      
      do k=kstart,klastact
         do j=jfirstact,jlastact
            do i=ifirstact,ilastact
               normfact = h3
*** Diagonal terms
               dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))
               dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))
               dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
     *               d4a*(u(3,i,j,k+1)-u(3,i,j,k-1))
               dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
     *               d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k))
               dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
     *               d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k))
               dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
     *               d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1))

               duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
     *                d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k))
               dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
     *                d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k))
               dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
     *                d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1))
               dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
     *                d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k))
               dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
     *                d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k))
               dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
     *                d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1))


               glambda(i,j,k) = glambda(i,j,k) + (
     *           (dux+dvy+dwz)*(dkx+dly+dmz) +
     *           dt2o12*( (duax+dvay+dwaz)*(dkx+dly+dmz) +
     *                    (dux+dvy+dwz)*(dkax+dlay+dmaz) )
     *                  )*normfact*ih2

               gmu(i,j,k) = gmu(i,j,k) + 2*(
     *          dux*dkx+dvy*dly+dwz*dmz +
     *           dt2o12*( (duax*dkx+dvay*dly+dwaz*dmz)+
     *                    (dux*dkax+dvy*dlay+dwz*dmaz) )
     *                     )*normfact*ih2

*** Off diagonal stresses
***  xy
               stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
     *                 d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
     *                 d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
     *                 d4a*(u(1,i,j+1,k)-u(1,i,j-1,k))
               stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
     *                 d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
     *                 d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
     *                 d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k))
               stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
     *                  d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
     *                  d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
     *                  d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k))
               stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
     *                  d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
     *                  d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
     *                  d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k))

               gmu(i,j,k)= gmu(i,j,k) + (stuxy*stkxy + 
     *            dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*normfact*ih2

***  xz
               stuxz = d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
     *                 d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *                 d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
     *                 d4a*(u(1,i,j,k+1)-u(1,i,j,k-1))
               stkxz = d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
     *                 d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k))+
     *                 d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
     *                 d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1))
               stuaxz = d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
     *                  d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k))+
     *                  d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
     *                  d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1))
               stkaxz = d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
     *                  d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k))+
     *                  d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
     *                  d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1))
               gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
     *            dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact*ih2

***  yz
               stuyz = d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
     *                 d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *                 d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
     *                 d4a*(u(2,i,j,k+1)-u(2,i,j,k-1))
               stkyz = d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
     *                 d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k))+
     *                 d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
     *                 d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1))
               stuayz = d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
     *                  d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k))+
     *                  d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
     *                  d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1))
               stkayz = d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
     *                  d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k))+
     *                  d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
     *                  d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1))
               gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
     *           dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact*ih2

*** Pos. def extra terms
*** x-direction
               d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
     *              3*u(1,i,j,k) -   u(1,i-1,j,k)
               d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
     *             3*u(1,i-1,j,k)-  u(1,i-2,j,k)
               d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
     *              3*kap(1,i,j,k) -   kap(1,i-1,j,k)
               d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
     *              3*kap(1,i-1,j,k) - kap(1,i-2,j,k)
               d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
     *               3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k)
               d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k)
               d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
     *               3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k)
               d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k)
               pd = (c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *      dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd

               d3up = u(2,i+2,j,k)- 3*u(2,i+1,j,k) +
     *              3*u(2,i,j,k)  -   u(2,i-1,j,k)
               d3um = u(2,i+1,j,k)- 3*u(2,i,j,k) +
     *              3*u(2,i-1,j,k)-   u(2,i-2,j,k)
               d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
     *              3*kap(2,i,j,k)  -  kap(2,i-1,j,k)
               d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
     *              3*kap(2,i-1,j,k)-  kap(2,i-2,j,k)
               d3uap = uacc(2,i+2,j,k) - 3*uacc(2,i+1,j,k)+
     *               3*uacc(2,i,j,k)   -   uacc(2,i-1,j,k)
               d3uam = uacc(2,i+1,j,k) - 3*uacc(2,i,j,k)+
     *               3*uacc(2,i-1,j,k) -   uacc(2,i-2,j,k)
               d3kap = kapacc(2,i+2,j,k) - 3*kapacc(2,i+1,j,k)+
     *               3*kapacc(2,i,j,k)   -   kapacc(2,i-1,j,k)
               d3kam = kapacc(2,i+1,j,k) - 3*kapacc(2,i,j,k)+
     *               3*kapacc(2,i-1,j,k) -   kapacc(2,i-2,j,k)
               pd = (c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                   (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+
     *              3*u(3,i,j,k)  -  u(3,i-1,j,k)
               d3um = u(3,i+1,j,k)-3*u(3,i,j,k)+
     *              3*u(3,i-1,j,k)-  u(3,i-2,j,k)
               d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
     *              3*kap(3,i,j,k)  -  kap(3,i-1,j,k)
               d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
     *              3*kap(3,i-1,j,k)-  kap(3,i-2,j,k)
               d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
     *               3*uacc(3,i,j,k) -   uacc(3,i-1,j,k)
               d3uam = uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
     *               3*uacc(3,i-1,j,k)-  uacc(3,i-2,j,k)
               d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
     *               3*kapacc(3,i,j,k) -   kapacc(3,i-1,j,k)
               d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i-1,j,k)-  kapacc(3,i-2,j,k)

               pd = ( c6*0.5d0*( d3up*d3kp + d3um*d3km +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

*** y-direction
               d3up = u(1,i,j+2,k) - 3*u(1,i,j+1,k) +
     *              3*u(1,i,j,k)   -   u(1,i,j-1,k)
               d3um = u(1,i,j+1,k) - 3*u(1,i,j,k) +
     *              3*u(1,i,j-1,k) -   u(1,i,j-2,k)
               d3kp = kap(1,i,j+2,k) - 3*kap(1,i,j+1,k) +
     *              3*kap(1,i,j,k)   -   kap(1,i,j-1,k)
               d3km = kap(1,i,j+1,k) - 3*kap(1,i,j,k)+
     *              3*kap(1,i,j-1,k) -   kap(1,i,j-2,k)
               d3uap = uacc(1,i,j+2,k)- 3*uacc(1,i,j+1,k)+
     *               3*uacc(1,i,j,k) -    uacc(1,i,j-1,k)
               d3uam =uacc(1,i,j+1,k) - 3*uacc(1,i,j,k)+
     *              3*uacc(1,i,j-1,k) -   uacc(1,i,j-2,k)
               d3kap = kapacc(1,i,j+2,k) - 3*kapacc(1,i,j+1,k)+
     *               3*kapacc(1,i,j,k)   -   kapacc(1,i,j-1,k)
               d3kam = kapacc(1,i,j+1,k) - 3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i,j-1,k) -   kapacc(1,i,j-2,k)

               pd = ( c6*0.5d0*( d3up*d3kp + d3um*d3km  +
     *      dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam) ) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(2,i,j+2,k) - 3*u(2,i,j+1,k)+
     *              3*u(2,i,j,k)   -   u(2,i,j-1,k)
               d3um =u(2,i,j+1,k)  - 3*u(2,i,j,k)+
     *             3*u(2,i,j-1,k) -    u(2,i,j-2,k)
               d3kp = kap(2,i,j+2,k) - 3*kap(2,i,j+1,k)+
     *              3*kap(2,i,j,k)   -   kap(2,i,j-1,k)
               d3km = kap(2,i,j+1,k) - 3*kap(2,i,j,k)+
     *              3*kap(2,i,j-1,k) -   kap(2,i,j-2,k)
               d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
     *               3*uacc(2,i,j,k)  -  uacc(2,i,j-1,k)
               d3uam = uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
     *               3*uacc(2,i,j-1,k)-  uacc(2,i,j-2,k)
               d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
     *               3*kapacc(2,i,j,k)  -  kapacc(2,i,j-1,k)
               d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
     *               3*kapacc(2,i,j-1,k)-  kapacc(2,i,j-2,k)

               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2

               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd

               d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+
     *              3*u(3,i,j,k)  -  u(3,i,j-1,k)
               d3um = u(3,i,j+1,k)-3*u(3,i,j,k)+
     *              3*u(3,i,j-1,k)-  u(3,i,j-2,k)
               d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
     *              3*kap(3,i,j,k)  -  kap(3,i,j-1,k)
               d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
     *              3*kap(3,i,j-1,k)-  kap(3,i,j-2,k)
               d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
     *               3*uacc(3,i,j,k)  -  uacc(3,i,j-1,k)
               d3uam = uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
     *               3*uacc(3,i,j-1,k)-  uacc(3,i,j-2,k)
               d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
     *               3*kapacc(3,i,j,k)  -  kapacc(3,i,j-1,k)
               d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i,j-1,k)-  kapacc(3,i,j-2,k)
               pd = ( c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam) ) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

*** z-direction
               d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
     *                u(1,i,j,k-1)
               d3um =u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
     *               u(1,i,j,k-2)
               d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
     *              3*kap(1,i,j,k)-kap(1,i,j,k-1)
               d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
     *              3*kap(1,i,j,k-1)-kap(1,i,j,k-2)
               d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
     *               3*uacc(1,i,j,k)-uacc(1,i,j,k-1)
               d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2)
               d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
     *               3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1)
               d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2)
               pd = ( c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(2,i,j,k+2)-3*u(2,i,j,k+1)+
     *              3*u(2,i,j,k)  -  u(2,i,j,k-1)
               d3um =u(2,i,j,k+1)-3*u(2,i,j,k)+
     *             3*u(2,i,j,k-1)-  u(2,i,j,k-2)
               d3kp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
     *              3*kap(2,i,j,k)  -  kap(2,i,j,k-1)
               d3km = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
     *              3*kap(2,i,j,k-1)-  kap(2,i,j,k-2)
               d3uap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
     *               3*uacc(2,i,j,k)  -  uacc(2,i,j,k-1)
               d3uam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i,j,k-1)-  uacc(2,i,j,k-2)
               d3kap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
     *               3*kapacc(2,i,j,k)  -  kapacc(2,i,j,k-1)
               d3kam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
     *               3*kapacc(2,i,j,k-1) - kapacc(2,i,j,k-2)

               pd = ( c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam)) ) )*h3*ih2
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
     *               u(3,i,j,k-1)
               d3um =u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
     *               u(3,i,j,k-2)
               d3kp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
     *           3*kap(3,i,j,k)-kap(3,i,j,k-1)
               d3km = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
     *           3*kap(3,i,j,k-1)-kap(3,i,j,k-2)
               d3uap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
     *               3*uacc(3,i,j,k)-uacc(3,i,j,k-1)
               d3uam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
     *              3*uacc(3,i,j,k-1)- uacc(3,i,j,k-2)
               d3kap = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
     *               3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1)
               d3kam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )*h3*ih2
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k)    + 2*pd
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine ADDGRADMULAC( ifirst, ilast, jfirst, jlast, kfirst,
     *  klast, ifirstact, ilastact, jfirstact, jlastact, kfirstact,
     *   klastact, kap, kapacc, u, uacc, gmu, glambda,
     *    dt, h, met, jac, onesided, nb, wb, bop )
     * bind(c)
      implicit none
      real*8 d4a, d4b, c6, c8, al1, al2, al3, al4
      parameter( d4a=2d0/3, d4b=-1d0/12, c6=1d0/18, c8=1d0/144 )
      parameter( al1 = 181507d0/1719312, al2 = -1441d0/39984 )
      parameter( al3 = -2593d0/151704,   al4 = 11d0/3528 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer ifirstact, ilastact, jfirstact, jlastact, kfirstact
      integer klastact, kstart, onesided(6), nb, wb
      integer i, j, k, m
      real*8  kap(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  kapacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  gmu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  glambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  bop(nb,wb)
      real*8  dt, h, dt2o12, wgh(4), normfact
      real*8  dux, dvy, dwz, dkx, dly, dmz, duax, dvay, dwaz, dkax 
      real*8  dlay, dmaz, stuxy, stkxy, stuaxy, stkaxy
      real*8  stuxz, stkxz, stuaxz, stkaxz, stuyz, stkyz, stuayz, stkayz
      real*8  d3up, d3um, d3kp, d3km, pd, d3uap, d3uam, d3kap, d3kam
      real*8  d3vp, d3vm, d3lp, d3lm, d3vap, d3vam, d3lap, d3lam
      real*8  d3wp, d3wm, d3mp, d3mm, d3wap, d3wam, d3map, d3mam
      real*8  dkz, dlz, duz, dvz, duaz, dvaz, dkaz, dlaz
      real*8  dwzm, dwazm, dmzm, dmazm
      real*8  w8(4), w6m(4), w6p(4), m1sq, mucof, mucofp

      dt2o12 = dt*dt/12
      wgh(1) = 17d0/48
      wgh(2) = 59d0/48
      wgh(3) = 43d0/48
      wgh(4) = 49d0/48

      kstart = kfirstact
      if( kfirstact .le. 4 .and. onesided(5).eq.1 )then
         kstart = 5

         w8(1) = 0
         w8(2) = 0
         w8(3) = 1
         w8(4) = 1

         w6m(1)= 0
         w6p(1)= 0

         w6m(2)= 0
         w6p(2)= al1

         w6m(3)= al1
         w6p(3)= al1+al2

         w6m(4)= al1+al2
         w6p(4)= al1+al2+al3

         do k=kfirstact,4
            do j=jfirstact,jlastact
               do i=ifirstact,ilastact
                  normfact = wgh(k)
*** Diagonal terms
               dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))
               dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))
c               dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
c     *               d4a*(u(3,i,j,k+1)-u(3,i,j,k-1))
               dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
     *               d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k))
               dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
     *               d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k))
c               dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
c     *               d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1))

               duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
     *               d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k))
               dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
     *               d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k))
c               dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
c     *               d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1))
               dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
     *               d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k))
               dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
     *               d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k))
c               dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
c     *               d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1))
               duz = 0
               dvz = 0
               dwz = 0

               dkz = 0
               dlz = 0
               dmz = 0

               duaz=0
               dvaz=0
               dwaz=0

               dkaz=0
               dlaz=0
               dmaz=0
               do m=1,wb
                  duz = duz + bop(k,m)*u(1,i,j,m)
                  dvz = dvz + bop(k,m)*u(2,i,j,m)
                  dwz = dwz + bop(k,m)*u(3,i,j,m)

                  dkz = dkz + bop(k,m)*kap(1,i,j,m)
                  dlz = dlz + bop(k,m)*kap(2,i,j,m)
                  dmz = dmz + bop(k,m)*kap(3,i,j,m)

                  duaz = duaz + bop(k,m)*uacc(1,i,j,m)
                  dvaz = dvaz + bop(k,m)*uacc(2,i,j,m)
                  dwaz = dwaz + bop(k,m)*uacc(3,i,j,m)

                  dkaz= dkaz+ bop(k,m)*kapacc(1,i,j,m)
                  dlaz= dlaz+ bop(k,m)*kapacc(2,i,j,m)
                  dmaz= dmaz+ bop(k,m)*kapacc(3,i,j,m)

               enddo
               dux = met(1,i,j,k)*dux + met(2,i,j,k)*duz
               dvy = met(1,i,j,k)*dvy + met(3,i,j,k)*dvz
               dkx = met(1,i,j,k)*dkx + met(2,i,j,k)*dkz
               dly = met(1,i,j,k)*dly + met(3,i,j,k)*dlz
               duax = met(1,i,j,k)*duax + met(2,i,j,k)*duaz
               dvay = met(1,i,j,k)*dvay + met(3,i,j,k)*dvaz
               dkax = met(1,i,j,k)*dkax + met(2,i,j,k)*dkaz
               dlay = met(1,i,j,k)*dlay + met(3,i,j,k)*dlaz

               dwzm  = dwz*met(4,i,j,k)               
               dmzm  = dmz*met(4,i,j,k)               
               dwazm = dwaz*met(4,i,j,k)               
               dmazm = dmaz*met(4,i,j,k)               

               glambda(i,j,k) = glambda(i,j,k) + (
     *           (dux+dvy+dwzm)*(dkx+dly+dmzm) +
     *           dt2o12*( (duax+dvay+dwazm)*(dkx+dly+dmzm) +
     *                    (dux+dvy+dwzm)*(dkax+dlay+dmazm) )
     *                  )*normfact
               gmu(i,j,k) = gmu(i,j,k) + 2*(
     *          dux*dkx+dvy*dly+dwzm*dmzm +
     *           dt2o12*( (duax*dkx+dvay*dly+dwazm*dmzm)+
     *                    (dux*dkax+dvy*dlay+dwzm*dmazm) )
     *                     )*normfact
               
*** Off diagonal stresses
***  xy
               stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
     *                 d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
     *                 d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
     *                 d4a*(u(1,i,j+1,k)-u(1,i,j-1,k))
               stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
     *                 d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
     *                 d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
     *                 d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k))
               stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
     *                 d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
     *                 d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
     *                 d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k))
               stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
     *                 d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
     *                 d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
     *                 d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k))

               stuxy = met(1,i,j,k)*stuxy+met(2,i,j,k)*dvz +
     *                                    met(3,i,j,k)*duz
               stkxy = met(1,i,j,k)*stkxy+met(2,i,j,k)*dlz +
     *                                    met(3,i,j,k)*dkz
               stuaxy= met(1,i,j,k)*stuaxy+met(2,i,j,k)*dvaz +
     *                                     met(3,i,j,k)*duaz
               stkaxy= met(1,i,j,k)*stkaxy+met(2,i,j,k)*dlaz +
     *                                     met(3,i,j,k)*dkaz

               gmu(i,j,k)= gmu(i,j,k) + (
     *      stuxy*stkxy + dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )*
     *                                            normfact

c               dwz=dwz/met(4,i,j,k)
c               dmz=dmz/met(4,i,j,k)
c               dwaz=dwaz/met(4,i,j,k)
c               dmaz=dmaz/met(4,i,j,k)
***  xz
               stuxz = met(1,i,j,k)*(d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
     *                               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k)) )
c     *                 d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
c     *                 d4a*(u(1,i,j,k+1)-u(1,i,j,k-1))
               stkxz = met(1,i,j,k)*(
     *                 d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
     *                 d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k)) )
c     *                 d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
c     *                 d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1))
               stuaxz = met(1,i,j,k)*(
     *                 d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
     *                 d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k)) )
c     *                 d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
c     *                 d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1))
               stkaxz = met(1,i,j,k)*(
     *                 d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
     *                 d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k)) )
c     *                 d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
c     *                 d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1))


               stuxz = stuxz + met(2,i,j,k)*dwz+met(4,i,j,k)*duz
               stkxz = stkxz + met(2,i,j,k)*dmz+met(4,i,j,k)*dkz
               stuaxz= stuaxz+ met(2,i,j,k)*dwaz+met(4,i,j,k)*duaz
               stkaxz= stkaxz+ met(2,i,j,k)*dmaz+met(4,i,j,k)*dkaz

               gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
     *            dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )*normfact
***  yz
               stuyz = met(1,i,j,k)*(d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
     *                 d4a*(u(3,i,j+1,k)-u(3,i,j-1,k)))
c     *                 d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
c     *                 d4a*(u(2,i,j,k+1)-u(2,i,j,k-1))
               stkyz = met(1,i,j,k)*(
     *                 d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
     *                 d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k)) )
c     *                 d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
c     *                 d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1))
               stuayz = met(1,i,j,k)*(
     *                 d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
     *                 d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k)) )
c     *                 d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
c     *                 d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1))
               stkayz = met(1,i,j,k)*(
     *                 d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
     *                 d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k)) )
c     *                 d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
c     *                 d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1))

               stuyz = stuyz + met(3,i,j,k)*dwz+met(4,i,j,k)*dvz
               stkyz = stkyz + met(3,i,j,k)*dmz+met(4,i,j,k)*dlz
               stuayz= stuayz+ met(3,i,j,k)*dwaz+met(4,i,j,k)*dvaz
               stkayz= stkayz+ met(3,i,j,k)*dmaz+met(4,i,j,k)*dlaz

               gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
     *           dt2o12*( stuyz*stkayz+stuayz*stkyz) )*normfact

*** Pos. def extra terms
*** x-direction
***   metric is constant for the x and y directions, just use as a norm weight
               m1sq = met(1,i,j,k)*met(1,i,j,k)

               d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
     *              3*u(1,i,j,k) -   u(1,i-1,j,k)
               d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
     *             3*u(1,i-1,j,k)-  u(1,i-2,j,k)
               d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
     *              3*kap(1,i,j,k) -   kap(1,i-1,j,k)
               d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
     *              3*kap(1,i-1,j,k) - kap(1,i-2,j,k)
               d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
     *               3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k)
               d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k)
               d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
     *               3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k)
               d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k)
               pd = (c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *      dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam)) ) )*m1sq*wgh(k)
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd
               d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
     *               u(2,i-1,j,k)
               d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
     *               u(2,i-2,j,k)
               d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
     *           3*kap(2,i,j,k)-kap(2,i-1,j,k)
               d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
     *           3*kap(2,i-1,j,k)-kap(2,i-2,j,k)
               d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
     *              3*uacc(2,i,j,k)- uacc(2,i-1,j,k)
               d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k)
               d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
     *           3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k)
               d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
     *           3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
     *               u(3,i-1,j,k)
               d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
     *               u(3,i-2,j,k)
               d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
     *           3*kap(3,i,j,k)-kap(3,i-1,j,k)
               d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
     *           3*kap(3,i-1,j,k)-kap(3,i-2,j,k)
               d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
     *            3*uacc(3,i,j,k)- uacc(3,i-1,j,k)
               d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
     *             3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k)
               d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
     *           3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k)
               d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
     *           3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp + d3um*d3km +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

*** y-direction
               d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
     *               u(1,i,j-1,k)
               d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
     *               u(1,i,j-2,k)
               d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
     *           3*kap(1,i,j,k)-kap(1,i,j-1,k)
               d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
     *           3*kap(1,i,j-1,k)-kap(1,i,j-2,k)
               d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
     *             3*uacc(1,i,j,k)- uacc(1,i,j-1,k)
               d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
     *             3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k)
               d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
     *           3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k)
               d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
     *           3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
     *               u(2,i,j-1,k)
               d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
     *               u(2,i,j-2,k)
               d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
     *           3*kap(2,i,j,k)-kap(2,i,j-1,k)
               d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
     *           3*kap(2,i,j-1,k)-kap(2,i,j-2,k)
               d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
     *               3*uacc(2,i,j,k)-uacc(2,i,j-1,k)
               d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i,j-1,k)-uacc(2,i,j-2,k)
               d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
     *           3*kapacc(2,i,j,k)-kapacc(2,i,j-1,k)
               d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
     *           3*kapacc(2,i,j-1,k)-kapacc(2,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh(k)
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k) = gmu(i,j,k) + 2*pd

               d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
     *                u(3,i,j-1,k)
               d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
     *               u(3,i,j-2,k)
               d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
     *              3*kap(3,i,j,k)-kap(3,i,j-1,k)
               d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
     *              3*kap(3,i,j-1,k)-kap(3,i,j-2,k)
               d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
     *               3*uacc(3,i,j,k)- uacc(3,i,j-1,k)
               d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
     *               3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k)
               d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
     *               3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k)
               d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam))) )*m1sq*wgh(k)
               gmu(i,j,k) = gmu(i,j,k) + pd

*** z-direction
               if( k.ge.2 )then
*** All derivatives are needed
                  d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
     *                u(1,i,j,k-1)
                  d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
     *                u(1,i,j,k-2)
                  d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
     *                 3*kap(1,i,j,k)-kap(1,i,j,k-1)
                  d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
     *                 3*kap(1,i,j,k-1)-kap(1,i,j,k-2)
                  d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
     *                  3*uacc(1,i,j,k)-uacc(1,i,j,k-1)
                  d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
     *                 3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2)
                  d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
     *                  3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1)
                  d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
     *                  3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2)

                  d3vp = u(2,i,j,k+2)-3*u(2,i,j,k+1)+3*u(2,i,j,k)-
     *                u(2,i,j,k-1)
                  d3vm = u(2,i,j,k+1)-3*u(2,i,j,k)+3*u(2,i,j,k-1)-
     *                u(2,i,j,k-2)
                  d3lp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
     *                 3*kap(2,i,j,k)-kap(2,i,j,k-1)
                  d3lm = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
     *                 3*kap(2,i,j,k-1)-kap(2,i,j,k-2)
                  d3vap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
     *                  3*uacc(2,i,j,k)-uacc(2,i,j,k-1)
                  d3vam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
     *                 3*uacc(2,i,j,k-1)-uacc(2,i,j,k-2)
                  d3lap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
     *                  3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1)
                  d3lam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
     *                  3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2)

                  d3wp = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
     *                u(3,i,j,k-1)
                  d3wm = u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
     *                u(3,i,j,k-2)
                  d3mp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
     *                 3*kap(3,i,j,k)-kap(3,i,j,k-1)
                  d3mm = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
     *                 3*kap(3,i,j,k-1)-kap(3,i,j,k-2)
                  d3wap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
     *                  3*uacc(3,i,j,k)-uacc(3,i,j,k-1)
                  d3wam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
     *                 3*uacc(3,i,j,k-1)-uacc(3,i,j,k-2)
                  d3map = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
     *                  3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1)
                  d3mam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
     *                  3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2)

                  mucof=met(2,i,j,k)**2+met(3,i,j,k)**2+met(4,i,j,k)**2
                  mucofp=met(2,i,j,k+1)**2+met(3,i,j,k+1)**2+
     *                   met(4,i,j,k+1)**2

*** u-u
                  pd = ( ( 0.5d0*( w6p(k)*d3up*d3kp+w6m(k)*d3um*d3km  +
     *        dt2o12*(w6p(k)*d3up*d3kap + w6m(k)*d3um*d3kam+ 
     *                w6p(k)*d3uap*d3kp + w6m(k)*d3uam*d3km)))
     *            + w8(k)*c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(2,i,j,k)**2)
                  glambda(i,j,k) = glambda(i,j,k) + pd*met(2,i,j,k)**2
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3kp+
     *               dt2o12*(-al4*d3up*d3kap -al4*d3uap*d3kp ) )
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
     *                                    pd*(mucofp+met(2,i,j,k+1)**2)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
     *                                    pd*met(2,i,j,k+1)**2
                  endif

*** u-v
                  pd = ( ( 0.5d0*( w6p(k)*d3vp*d3kp+w6m(k)*d3vm*d3km  +
     *        dt2o12*(w6p(k)*d3vp*d3kap + w6m(k)*d3vm*d3kam+ 
     *                w6p(k)*d3vap*d3kp + w6m(k)*d3vam*d3km)))
     *            + w8(k)*c8*( (d3vp-d3vm)*(d3kp-d3km) +
     *        dt2o12*(  (d3vap-d3vam)*(d3kp-d3km)+
     *                            (d3vp-d3vm)*(d3kap-d3kam))) )
     *                                *met(2,i,j,k)*met(3,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3vp*d3kp+
     *               dt2o12*(-al4*d3vp*d3kap -al4*d3vap*d3kp ) )*
     *                               met(2,i,j,k+1)*met(3,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** u-w 
                  pd = ( ( 0.5d0*( w6p(k)*d3wp*d3kp+w6m(k)*d3wm*d3km  +
     *        dt2o12*(w6p(k)*d3wp*d3kap + w6m(k)*d3wm*d3kam+ 
     *                w6p(k)*d3wap*d3kp + w6m(k)*d3wam*d3km)))
     *            + w8(k)*c8*( (d3wp-d3wm)*(d3kp-d3km) +
     *        dt2o12*(  (d3wap-d3wam)*(d3kp-d3km)+
     *                            (d3wp-d3wm)*(d3kap-d3kam))) )
     *                                *met(2,i,j,k)*met(4,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3wp*d3kp+
     *               dt2o12*(-al4*d3wp*d3kap -al4*d3wap*d3kp ) )*
     *                               met(2,i,j,k+1)*met(4,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** v-u
                  pd = ( ( 0.5d0*( w6p(k)*d3up*d3lp+w6m(k)*d3um*d3lm  +
     *        dt2o12*(w6p(k)*d3up*d3lap + w6m(k)*d3um*d3lam+ 
     *                w6p(k)*d3uap*d3lp + w6m(k)*d3uam*d3lm)))
     *            + w8(k)*c8*( (d3up-d3um)*(d3lp-d3lm) +
     *        dt2o12*(  (d3uap-d3uam)*(d3lp-d3lm)+
     *                            (d3up-d3um)*(d3lap-d3lam))) )
     *                                *met(2,i,j,k)*met(3,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3lp+
     *               dt2o12*(-al4*d3up*d3lap -al4*d3uap*d3lp ) )*
     *                               met(2,i,j,k+1)*met(3,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** v-v
                  pd = ( ( 0.5d0*( w6p(k)*d3vp*d3lp+w6m(k)*d3vm*d3lm  +
     *        dt2o12*(w6p(k)*d3vp*d3lap + w6m(k)*d3vm*d3lam+ 
     *                w6p(k)*d3vap*d3lp + w6m(k)*d3vam*d3lm)))
     *            + w8(k)*c8*( (d3vp-d3vm)*(d3lp-d3lm) +
     *        dt2o12*(  (d3vap-d3vam)*(d3lp-d3lm)+
     *                            (d3vp-d3vm)*(d3lap-d3lam))) )
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(3,i,j,k)**2)
                  glambda(i,j,k) = glambda(i,j,k) + pd*met(3,i,j,k)**2
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3vp*d3lp+
     *               dt2o12*(-al4*d3vp*d3lap -al4*d3vap*d3lp ) )
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
     *                                    pd*met(3,i,j,k+1)**2
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
     *                                    pd*(mucofp+met(3,i,j,k+1)**2)
                  endif

*** v-w
                  pd = ( ( 0.5d0*( w6p(k)*d3wp*d3lp+w6m(k)*d3wm*d3lm  +
     *        dt2o12*(w6p(k)*d3wp*d3lap + w6m(k)*d3wm*d3lam+ 
     *                w6p(k)*d3wap*d3lp + w6m(k)*d3wam*d3lm)))
     *            + w8(k)*c8*( (d3wp-d3wm)*(d3lp-d3lm) +
     *        dt2o12*(  (d3wap-d3wam)*(d3lp-d3lm)+
     *                            (d3wp-d3wm)*(d3lap-d3lam))) )
     *                                *met(3,i,j,k)*met(4,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3wp*d3lp+
     *               dt2o12*(-al4*d3wp*d3lap -al4*d3wap*d3lp ) )*
     *                               met(3,i,j,k+1)*met(4,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** w-u
                  pd = ( ( 0.5d0*( w6p(k)*d3up*d3mp+w6m(k)*d3um*d3mm  +
     *        dt2o12*(w6p(k)*d3up*d3map + w6m(k)*d3um*d3mam+ 
     *                w6p(k)*d3uap*d3mp + w6m(k)*d3uam*d3mm)))
     *            + w8(k)*c8*( (d3up-d3um)*(d3mp-d3mm) +
     *        dt2o12*(  (d3uap-d3uam)*(d3mp-d3mm)+
     *                            (d3up-d3um)*(d3map-d3mam))) )
     *                                *met(2,i,j,k)*met(4,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3up*d3mp+
     *               dt2o12*(-al4*d3up*d3map -al4*d3uap*d3mp ) )*
     *                               met(2,i,j,k+1)*met(4,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** w-v
                  pd = ( ( 0.5d0*( w6p(k)*d3vp*d3mp+w6m(k)*d3vm*d3mm  +
     *        dt2o12*(w6p(k)*d3vp*d3map + w6m(k)*d3vm*d3mam+ 
     *                w6p(k)*d3vap*d3mp + w6m(k)*d3vam*d3mm)))
     *            + w8(k)*c8*( (d3vp-d3vm)*(d3mp-d3mm) +
     *        dt2o12*(  (d3vap-d3vam)*(d3mp-d3mm)+
     *                            (d3vp-d3vm)*(d3map-d3mam))) )
     *                                *met(3,i,j,k)*met(4,i,j,k)
                  gmu(i,j,k) = gmu(i,j,k) + pd
                  glambda(i,j,k) = glambda(i,j,k) + pd
                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3vp*d3mp+
     *               dt2o12*(-al4*d3vp*d3map -al4*d3vap*d3mp ) )*
     *                               met(3,i,j,k+1)*met(4,i,j,k+1)
                     glambda(i,j,k+1) = glambda(i,j,k+1) + pd
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  pd
                  endif

*** w-w
                  pd = ( ( 0.5d0*( w6p(k)*d3wp*d3mp+w6m(k)*d3wm*d3mm  +
     *        dt2o12*(w6p(k)*d3wp*d3map + w6m(k)*d3wm*d3mam+ 
     *                w6p(k)*d3wap*d3mp + w6m(k)*d3wam*d3mm)))
     *            + w8(k)*c8*( (d3wp-d3wm)*(d3mp-d3mm) +
     *        dt2o12*(  (d3wap-d3wam)*(d3mp-d3mm)+
     *                            (d3wp-d3wm)*(d3map-d3mam))) )
                  gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(4,i,j,k)**2)
                  glambda(i,j,k) = glambda(i,j,k) + pd*met(4,i,j,k)**2

                  if( k.eq.4 )then
                     pd = 0.5d0*( -al4*d3wp*d3mp+
     *               dt2o12*(-al4*d3wp*d3map -al4*d3wap*d3mp ) )
                     glambda(i,j,k+1) = glambda(i,j,k+1) + 
     *                                    pd*met(4,i,j,k+1)**2
                     gmu(i,j,k+1)     = gmu(i,j,k+1)   +  
     *                                    pd*(mucofp+met(4,i,j,k+1)**2)
                  endif
                  endif
               enddo
            enddo
         enddo
      endif
      
      do k=kstart,klastact
         do j=jfirstact,jlastact
            do i=ifirstact,ilastact
*** Diagonal terms
               dux = d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))+
     *               d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))
               dvy = d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))+
     *               d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))
               dwz = d4b*(u(3,i,j,k+2)-u(3,i,j,k-2))+
     *               d4a*(u(3,i,j,k+1)-u(3,i,j,k-1))
               dkx = d4b*(kap(1,i+2,j,k)-kap(1,i-2,j,k))+
     *               d4a*(kap(1,i+1,j,k)-kap(1,i-1,j,k))
               dly = d4b*(kap(2,i,j+2,k)-kap(2,i,j-2,k))+
     *               d4a*(kap(2,i,j+1,k)-kap(2,i,j-1,k))
               dmz = d4b*(kap(3,i,j,k+2)-kap(3,i,j,k-2))+
     *               d4a*(kap(3,i,j,k+1)-kap(3,i,j,k-1))

               duax = d4b*(uacc(1,i+2,j,k)-uacc(1,i-2,j,k))+
     *               d4a*(uacc(1,i+1,j,k)-uacc(1,i-1,j,k))
               dvay = d4b*(uacc(2,i,j+2,k)-uacc(2,i,j-2,k))+
     *               d4a*(uacc(2,i,j+1,k)-uacc(2,i,j-1,k))
               dwaz = d4b*(uacc(3,i,j,k+2)-uacc(3,i,j,k-2))+
     *               d4a*(uacc(3,i,j,k+1)-uacc(3,i,j,k-1))
               dkax = d4b*(kapacc(1,i+2,j,k)-kapacc(1,i-2,j,k))+
     *               d4a*(kapacc(1,i+1,j,k)-kapacc(1,i-1,j,k))
               dlay = d4b*(kapacc(2,i,j+2,k)-kapacc(2,i,j-2,k))+
     *               d4a*(kapacc(2,i,j+1,k)-kapacc(2,i,j-1,k))
               dmaz = d4b*(kapacc(3,i,j,k+2)-kapacc(3,i,j,k-2))+
     *               d4a*(kapacc(3,i,j,k+1)-kapacc(3,i,j,k-1))

               duz = d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
     *               d4a*(u(1,i,j,k+1)-u(1,i,j,k-1))
               dvz = d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
     *               d4a*(u(2,i,j,k+1)-u(2,i,j,k-1))

               dkz = d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
     *               d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1))
               dlz = d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
     *               d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1))

               duaz= d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
     *               d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1))
               dvaz= d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
     *               d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1))

               dkaz= d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
     *               d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1))
               dlaz= d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
     *               d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1))

               dux = met(1,i,j,k)*dux + met(2,i,j,k)*duz
               dvy = met(1,i,j,k)*dvy + met(3,i,j,k)*dvz
               dkx = met(1,i,j,k)*dkx + met(2,i,j,k)*dkz
               dly = met(1,i,j,k)*dly + met(3,i,j,k)*dlz
               duax = met(1,i,j,k)*duax + met(2,i,j,k)*duaz
               dvay = met(1,i,j,k)*dvay + met(3,i,j,k)*dvaz
               dkax = met(1,i,j,k)*dkax + met(2,i,j,k)*dkaz
               dlay = met(1,i,j,k)*dlay + met(3,i,j,k)*dlaz

               dwzm  = dwz*met(4,i,j,k)               
               dmzm  = dmz*met(4,i,j,k)               
               dwazm = dwaz*met(4,i,j,k)               
               dmazm = dmaz*met(4,i,j,k)               

               glambda(i,j,k) = glambda(i,j,k) + (
     *           (dux+dvy+dwzm)*(dkx+dly+dmzm) +
     *           dt2o12*( (duax+dvay+dwazm)*(dkx+dly+dmzm) +
     *                    (dux+dvy+dwzm)*(dkax+dlay+dmazm) )
     *                  )
               gmu(i,j,k) = gmu(i,j,k) + 2*(
     *          dux*dkx+dvy*dly+dwzm*dmzm +
     *           dt2o12*( (duax*dkx+dvay*dly+dwazm*dmzm)+
     *                    (dux*dkax+dvy*dlay+dwzm*dmazm) )
     *                     )
               
*** Off diagonal stresses
***  xy
               stuxy = d4b*(u(2,i+2,j,k)-u(2,i-2,j,k))+
     *                 d4a*(u(2,i+1,j,k)-u(2,i-1,j,k))+
     *                 d4b*(u(1,i,j+2,k)-u(1,i,j-2,k))+
     *                 d4a*(u(1,i,j+1,k)-u(1,i,j-1,k))
               stkxy = d4b*(kap(2,i+2,j,k)-kap(2,i-2,j,k))+
     *                 d4a*(kap(2,i+1,j,k)-kap(2,i-1,j,k))+
     *                 d4b*(kap(1,i,j+2,k)-kap(1,i,j-2,k))+
     *                 d4a*(kap(1,i,j+1,k)-kap(1,i,j-1,k))
               stuaxy = d4b*(uacc(2,i+2,j,k)-uacc(2,i-2,j,k))+
     *                 d4a*(uacc(2,i+1,j,k)-uacc(2,i-1,j,k))+
     *                 d4b*(uacc(1,i,j+2,k)-uacc(1,i,j-2,k))+
     *                 d4a*(uacc(1,i,j+1,k)-uacc(1,i,j-1,k))
               stkaxy = d4b*(kapacc(2,i+2,j,k)-kapacc(2,i-2,j,k))+
     *                 d4a*(kapacc(2,i+1,j,k)-kapacc(2,i-1,j,k))+
     *                 d4b*(kapacc(1,i,j+2,k)-kapacc(1,i,j-2,k))+
     *                 d4a*(kapacc(1,i,j+1,k)-kapacc(1,i,j-1,k))

               stuxy = met(1,i,j,k)*stuxy+met(2,i,j,k)*dvz +
     *                                    met(3,i,j,k)*duz
               stkxy = met(1,i,j,k)*stkxy+met(2,i,j,k)*dlz +
     *                                    met(3,i,j,k)*dkz
               stuaxy= met(1,i,j,k)*stuaxy+met(2,i,j,k)*dvaz +
     *                                     met(3,i,j,k)*duaz
               stkaxy= met(1,i,j,k)*stkaxy+met(2,i,j,k)*dlaz +
     *                                     met(3,i,j,k)*dkaz

               gmu(i,j,k)= gmu(i,j,k) + (
     *      stuxy*stkxy + dt2o12*(stuxy*stkaxy+stuaxy*stkxy) )

***  xz
               stuxz = met(1,i,j,k)*(d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))+
     *                               d4a*(u(3,i+1,j,k)-u(3,i-1,j,k)) )
c     *                 d4b*(u(1,i,j,k+2)-u(1,i,j,k-2))+
c     *                 d4a*(u(1,i,j,k+1)-u(1,i,j,k-1))
               stkxz = met(1,i,j,k)*(
     *                 d4b*(kap(3,i+2,j,k)-kap(3,i-2,j,k))+
     *                 d4a*(kap(3,i+1,j,k)-kap(3,i-1,j,k)) )
c     *                 d4b*(kap(1,i,j,k+2)-kap(1,i,j,k-2))+
c     *                 d4a*(kap(1,i,j,k+1)-kap(1,i,j,k-1))
               stuaxz = met(1,i,j,k)*(
     *                 d4b*(uacc(3,i+2,j,k)-uacc(3,i-2,j,k))+
     *                 d4a*(uacc(3,i+1,j,k)-uacc(3,i-1,j,k)) )
c     *                 d4b*(uacc(1,i,j,k+2)-uacc(1,i,j,k-2))+
c     *                 d4a*(uacc(1,i,j,k+1)-uacc(1,i,j,k-1))
               stkaxz = met(1,i,j,k)*(
     *                 d4b*(kapacc(3,i+2,j,k)-kapacc(3,i-2,j,k))+
     *                 d4a*(kapacc(3,i+1,j,k)-kapacc(3,i-1,j,k)) )
c     *                 d4b*(kapacc(1,i,j,k+2)-kapacc(1,i,j,k-2))+
c     *                 d4a*(kapacc(1,i,j,k+1)-kapacc(1,i,j,k-1))

               stuxz = stuxz + met(2,i,j,k)*dwz+met(4,i,j,k)*duz
               stkxz = stkxz + met(2,i,j,k)*dmz+met(4,i,j,k)*dkz
               stuaxz= stuaxz+ met(2,i,j,k)*dwaz+met(4,i,j,k)*duaz
               stkaxz= stkaxz+ met(2,i,j,k)*dmaz+met(4,i,j,k)*dkaz

               gmu(i,j,k)= gmu(i,j,k) + ( stuxz*stkxz +
     *            dt2o12*(stuxz*stkaxz + stuaxz*stkxz) )
***  yz
               stuyz = met(1,i,j,k)*(d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))+
     *                 d4a*(u(3,i,j+1,k)-u(3,i,j-1,k)))
c     *                 d4b*(u(2,i,j,k+2)-u(2,i,j,k-2))+
c     *                 d4a*(u(2,i,j,k+1)-u(2,i,j,k-1))
               stkyz = met(1,i,j,k)*(
     *                 d4b*(kap(3,i,j+2,k)-kap(3,i,j-2,k))+
     *                 d4a*(kap(3,i,j+1,k)-kap(3,i,j-1,k)) )
c     *                 d4b*(kap(2,i,j,k+2)-kap(2,i,j,k-2))+
c     *                 d4a*(kap(2,i,j,k+1)-kap(2,i,j,k-1))
               stuayz = met(1,i,j,k)*(
     *                 d4b*(uacc(3,i,j+2,k)-uacc(3,i,j-2,k))+
     *                 d4a*(uacc(3,i,j+1,k)-uacc(3,i,j-1,k)) )
c     *                 d4b*(uacc(2,i,j,k+2)-uacc(2,i,j,k-2))+
c     *                 d4a*(uacc(2,i,j,k+1)-uacc(2,i,j,k-1))
               stkayz = met(1,i,j,k)*(
     *                 d4b*(kapacc(3,i,j+2,k)-kapacc(3,i,j-2,k))+
     *                 d4a*(kapacc(3,i,j+1,k)-kapacc(3,i,j-1,k)) )
c     *                 d4b*(kapacc(2,i,j,k+2)-kapacc(2,i,j,k-2))+
c     *                 d4a*(kapacc(2,i,j,k+1)-kapacc(2,i,j,k-1))

               stuyz = stuyz + met(3,i,j,k)*dwz+met(4,i,j,k)*dvz
               stkyz = stkyz + met(3,i,j,k)*dmz+met(4,i,j,k)*dlz
               stuayz= stuayz+ met(3,i,j,k)*dwaz+met(4,i,j,k)*dvaz
               stkayz= stkayz+ met(3,i,j,k)*dmaz+met(4,i,j,k)*dlaz

               gmu(i,j,k)= gmu(i,j,k) + ( stuyz*stkyz +
     *           dt2o12*( stuyz*stkayz+stuayz*stkyz) )

*** Pos. def extra terms
*** x-direction
***   metric is constant for the x and y directions, just use as a norm weight
               m1sq = met(1,i,j,k)*met(1,i,j,k)

               d3up = u(1,i+2,j,k)-3*u(1,i+1,j,k) +
     *              3*u(1,i,j,k) -   u(1,i-1,j,k)
               d3um =u(1,i+1,j,k)-3*u(1,i,j,k) +
     *             3*u(1,i-1,j,k)-  u(1,i-2,j,k)
               d3kp = kap(1,i+2,j,k)-3*kap(1,i+1,j,k) +
     *              3*kap(1,i,j,k) -   kap(1,i-1,j,k)
               d3km = kap(1,i+1,j,k)-3*kap(1,i,j,k) +
     *              3*kap(1,i-1,j,k) - kap(1,i-2,j,k)
               d3uap = uacc(1,i+2,j,k)-3*uacc(1,i+1,j,k)+
     *               3*uacc(1,i,j,k)  -  uacc(1,i-1,j,k)
               d3uam =uacc(1,i+1,j,k)-3*uacc(1,i,j,k)+
     *              3*uacc(1,i-1,j,k) - uacc(1,i-2,j,k)
               d3kap = kapacc(1,i+2,j,k)-3*kapacc(1,i+1,j,k)+
     *               3*kapacc(1,i,j,k) -   kapacc(1,i-1,j,k)
               d3kam = kapacc(1,i+1,j,k)-3*kapacc(1,i,j,k)+
     *               3*kapacc(1,i-1,j,k)-  kapacc(1,i-2,j,k)
               pd = (c6*0.5d0*( d3up*d3kp+d3um*d3km  +
     *      dt2o12*(d3up*d3kap + d3um*d3kam + d3uap*d3kp + d3uam*d3km) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam)) ) )*m1sq
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k)     = gmu(i,j,k) + 2*pd
               d3up = u(2,i+2,j,k)-3*u(2,i+1,j,k)+3*u(2,i,j,k)-
     *               u(2,i-1,j,k)
               d3um =u(2,i+1,j,k)-3*u(2,i,j,k)+3*u(2,i-1,j,k)-
     *               u(2,i-2,j,k)
               d3kp = kap(2,i+2,j,k)-3*kap(2,i+1,j,k)+
     *           3*kap(2,i,j,k)-kap(2,i-1,j,k)
               d3km = kap(2,i+1,j,k)-3*kap(2,i,j,k)+
     *           3*kap(2,i-1,j,k)-kap(2,i-2,j,k)
               d3uap = uacc(2,i+2,j,k)-3*uacc(2,i+1,j,k)+
     *              3*uacc(2,i,j,k)- uacc(2,i-1,j,k)
               d3uam =uacc(2,i+1,j,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i-1,j,k)- uacc(2,i-2,j,k)
               d3kap = kapacc(2,i+2,j,k)-3*kapacc(2,i+1,j,k)+
     *           3*kapacc(2,i,j,k)-kapacc(2,i-1,j,k)
               d3kam = kapacc(2,i+1,j,k)-3*kapacc(2,i,j,k)+
     *           3*kapacc(2,i-1,j,k)-kapacc(2,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km) ))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(3,i+2,j,k)-3*u(3,i+1,j,k)+3*u(3,i,j,k)-
     *               u(3,i-1,j,k)
               d3um =u(3,i+1,j,k)-3*u(3,i,j,k)+3*u(3,i-1,j,k)-
     *               u(3,i-2,j,k)
               d3kp = kap(3,i+2,j,k)-3*kap(3,i+1,j,k)+
     *           3*kap(3,i,j,k)-kap(3,i-1,j,k)
               d3km = kap(3,i+1,j,k)-3*kap(3,i,j,k)+
     *           3*kap(3,i-1,j,k)-kap(3,i-2,j,k)
               d3uap = uacc(3,i+2,j,k)-3*uacc(3,i+1,j,k)+
     *            3*uacc(3,i,j,k)- uacc(3,i-1,j,k)
               d3uam =uacc(3,i+1,j,k)-3*uacc(3,i,j,k)+
     *             3*uacc(3,i-1,j,k)- uacc(3,i-2,j,k)
               d3kap = kapacc(3,i+2,j,k)-3*kapacc(3,i+1,j,k)+
     *           3*kapacc(3,i,j,k)-kapacc(3,i-1,j,k)
               d3kam = kapacc(3,i+1,j,k)-3*kapacc(3,i,j,k)+
     *           3*kapacc(3,i-1,j,k)-kapacc(3,i-2,j,k)
               pd = (c6*( 0.5d0*( d3up*d3kp + d3um*d3km +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq
               gmu(i,j,k) = gmu(i,j,k) + pd

*** y-direction
               d3up = u(1,i,j+2,k)-3*u(1,i,j+1,k)+3*u(1,i,j,k)-
     *               u(1,i,j-1,k)
               d3um =u(1,i,j+1,k)-3*u(1,i,j,k)+3*u(1,i,j-1,k)-
     *               u(1,i,j-2,k)
               d3kp = kap(1,i,j+2,k)-3*kap(1,i,j+1,k)+
     *           3*kap(1,i,j,k)-kap(1,i,j-1,k)
               d3km = kap(1,i,j+1,k)-3*kap(1,i,j,k)+
     *           3*kap(1,i,j-1,k)-kap(1,i,j-2,k)
               d3uap = uacc(1,i,j+2,k)-3*uacc(1,i,j+1,k)+
     *             3*uacc(1,i,j,k)- uacc(1,i,j-1,k)
               d3uam =uacc(1,i,j+1,k)-3*uacc(1,i,j,k)+
     *             3*uacc(1,i,j-1,k)- uacc(1,i,j-2,k)
               d3kap = kapacc(1,i,j+2,k)-3*kapacc(1,i,j+1,k)+
     *           3*kapacc(1,i,j,k)-kapacc(1,i,j-1,k)
               d3kam = kapacc(1,i,j+1,k)-3*kapacc(1,i,j,k)+
     *           3*kapacc(1,i,j-1,k)-kapacc(1,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)) )
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq
               gmu(i,j,k) = gmu(i,j,k) + pd

               d3up = u(2,i,j+2,k)-3*u(2,i,j+1,k)+3*u(2,i,j,k)-
     *               u(2,i,j-1,k)
               d3um =u(2,i,j+1,k)-3*u(2,i,j,k)+3*u(2,i,j-1,k)-
     *               u(2,i,j-2,k)
               d3kp = kap(2,i,j+2,k)-3*kap(2,i,j+1,k)+
     *           3*kap(2,i,j,k)-kap(2,i,j-1,k)
               d3km = kap(2,i,j+1,k)-3*kap(2,i,j,k)+
     *           3*kap(2,i,j-1,k)-kap(2,i,j-2,k)
               d3uap = uacc(2,i,j+2,k)-3*uacc(2,i,j+1,k)+
     *               3*uacc(2,i,j,k)-uacc(2,i,j-1,k)
               d3uam =uacc(2,i,j+1,k)-3*uacc(2,i,j,k)+
     *              3*uacc(2,i,j-1,k)-uacc(2,i,j-2,k)
               d3kap = kapacc(2,i,j+2,k)-3*kapacc(2,i,j+1,k)+
     *           3*kapacc(2,i,j,k)-kapacc(2,i,j-1,k)
               d3kam = kapacc(2,i,j+1,k)-3*kapacc(2,i,j,k)+
     *           3*kapacc(2,i,j-1,k)-kapacc(2,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                     (d3up-d3um)*(d3kap-d3kam))) )*m1sq
               glambda(i,j,k) = glambda(i,j,k) + pd
               gmu(i,j,k) = gmu(i,j,k) + 2*pd

               d3up = u(3,i,j+2,k)-3*u(3,i,j+1,k)+3*u(3,i,j,k)-
     *                u(3,i,j-1,k)
               d3um =u(3,i,j+1,k)-3*u(3,i,j,k)+3*u(3,i,j-1,k)-
     *               u(3,i,j-2,k)
               d3kp = kap(3,i,j+2,k)-3*kap(3,i,j+1,k)+
     *              3*kap(3,i,j,k)-kap(3,i,j-1,k)
               d3km = kap(3,i,j+1,k)-3*kap(3,i,j,k)+
     *              3*kap(3,i,j-1,k)-kap(3,i,j-2,k)
               d3uap = uacc(3,i,j+2,k)-3*uacc(3,i,j+1,k)+
     *               3*uacc(3,i,j,k)- uacc(3,i,j-1,k)
               d3uam =uacc(3,i,j+1,k)-3*uacc(3,i,j,k)+
     *               3*uacc(3,i,j-1,k)-uacc(3,i,j-2,k)
               d3kap = kapacc(3,i,j+2,k)-3*kapacc(3,i,j+1,k)+
     *               3*kapacc(3,i,j,k)-kapacc(3,i,j-1,k)
               d3kam = kapacc(3,i,j+1,k)-3*kapacc(3,i,j,k)+
     *               3*kapacc(3,i,j-1,k)-kapacc(3,i,j-2,k)
               pd = (c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap+d3um*d3kam+ d3uap*d3kp + d3uam*d3km)))
     *            + c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                      (d3up-d3um)*(d3kap-d3kam))) )*m1sq
               gmu(i,j,k) = gmu(i,j,k) + pd

*** z-direction
*** All derivatives are needed
               d3up = u(1,i,j,k+2)-3*u(1,i,j,k+1)+3*u(1,i,j,k)-
     *                u(1,i,j,k-1)
               d3um = u(1,i,j,k+1)-3*u(1,i,j,k)+3*u(1,i,j,k-1)-
     *                u(1,i,j,k-2)
               d3kp = kap(1,i,j,k+2)-3*kap(1,i,j,k+1)+
     *                 3*kap(1,i,j,k)-kap(1,i,j,k-1)
               d3km = kap(1,i,j,k+1)-3*kap(1,i,j,k)+
     *                 3*kap(1,i,j,k-1)-kap(1,i,j,k-2)
               d3uap = uacc(1,i,j,k+2)-3*uacc(1,i,j,k+1)+
     *                  3*uacc(1,i,j,k)-uacc(1,i,j,k-1)
               d3uam =uacc(1,i,j,k+1)-3*uacc(1,i,j,k)+
     *                 3*uacc(1,i,j,k-1)-uacc(1,i,j,k-2)
               d3kap = kapacc(1,i,j,k+2)-3*kapacc(1,i,j,k+1)+
     *                  3*kapacc(1,i,j,k)-kapacc(1,i,j,k-1)
               d3kam = kapacc(1,i,j,k+1)-3*kapacc(1,i,j,k)+
     *                  3*kapacc(1,i,j,k-1)-kapacc(1,i,j,k-2)

               d3vp = u(2,i,j,k+2)-3*u(2,i,j,k+1)+3*u(2,i,j,k)-
     *                u(2,i,j,k-1)
               d3vm = u(2,i,j,k+1)-3*u(2,i,j,k)+3*u(2,i,j,k-1)-
     *                u(2,i,j,k-2)
               d3lp = kap(2,i,j,k+2)-3*kap(2,i,j,k+1)+
     *                 3*kap(2,i,j,k)-kap(2,i,j,k-1)
               d3lm = kap(2,i,j,k+1)-3*kap(2,i,j,k)+
     *                 3*kap(2,i,j,k-1)-kap(2,i,j,k-2)
               d3vap = uacc(2,i,j,k+2)-3*uacc(2,i,j,k+1)+
     *                  3*uacc(2,i,j,k)-uacc(2,i,j,k-1)
               d3vam =uacc(2,i,j,k+1)-3*uacc(2,i,j,k)+
     *                 3*uacc(2,i,j,k-1)-uacc(2,i,j,k-2)
               d3lap = kapacc(2,i,j,k+2)-3*kapacc(2,i,j,k+1)+
     *                  3*kapacc(2,i,j,k)-kapacc(2,i,j,k-1)
               d3lam = kapacc(2,i,j,k+1)-3*kapacc(2,i,j,k)+
     *                  3*kapacc(2,i,j,k-1)-kapacc(2,i,j,k-2)

               d3wp = u(3,i,j,k+2)-3*u(3,i,j,k+1)+3*u(3,i,j,k)-
     *                u(3,i,j,k-1)
               d3wm = u(3,i,j,k+1)-3*u(3,i,j,k)+3*u(3,i,j,k-1)-
     *                u(3,i,j,k-2)
               d3mp = kap(3,i,j,k+2)-3*kap(3,i,j,k+1)+
     *                 3*kap(3,i,j,k)-kap(3,i,j,k-1)
               d3mm = kap(3,i,j,k+1)-3*kap(3,i,j,k)+
     *                 3*kap(3,i,j,k-1)-kap(3,i,j,k-2)
               d3wap = uacc(3,i,j,k+2)-3*uacc(3,i,j,k+1)+
     *                  3*uacc(3,i,j,k)-uacc(3,i,j,k-1)
               d3wam =uacc(3,i,j,k+1)-3*uacc(3,i,j,k)+
     *                 3*uacc(3,i,j,k-1)-uacc(3,i,j,k-2)
               d3map = kapacc(3,i,j,k+2)-3*kapacc(3,i,j,k+1)+
     *                  3*kapacc(3,i,j,k)-kapacc(3,i,j,k-1)
               d3mam = kapacc(3,i,j,k+1)-3*kapacc(3,i,j,k)+
     *                  3*kapacc(3,i,j,k-1)-kapacc(3,i,j,k-2)

               mucof=met(2,i,j,k)**2+met(3,i,j,k)**2+met(4,i,j,k)**2
*** u-u
               pd = ( c6*( 0.5d0*( d3up*d3kp+d3um*d3km  +
     *        dt2o12*(d3up*d3kap + d3um*d3kam+ 
     *                d3uap*d3kp + d3uam*d3km)))
     *            +  c8*( (d3up-d3um)*(d3kp-d3km) +
     *        dt2o12*(  (d3uap-d3uam)*(d3kp-d3km)+
     *                            (d3up-d3um)*(d3kap-d3kam))) )
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(2,i,j,k)**2)
               glambda(i,j,k) = glambda(i,j,k) + pd*met(2,i,j,k)**2

*** u-v
               pd = ( c6*( 0.5d0*( d3vp*d3kp+d3vm*d3km  +
     *        dt2o12*(d3vp*d3kap + d3vm*d3kam+ 
     *                d3vap*d3kp + d3vam*d3km)))
     *            + c8*( (d3vp-d3vm)*(d3kp-d3km) +
     *        dt2o12*(  (d3vap-d3vam)*(d3kp-d3km)+
     *                            (d3vp-d3vm)*(d3kap-d3kam))) )
     *                                *met(2,i,j,k)*met(3,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** u-w 
               pd = ( c6*( 0.5d0*( d3wp*d3kp+d3wm*d3km  +
     *        dt2o12*(d3wp*d3kap + d3wm*d3kam+ 
     *                d3wap*d3kp + d3wam*d3km)))
     *            + c8*( (d3wp-d3wm)*(d3kp-d3km) +
     *        dt2o12*(  (d3wap-d3wam)*(d3kp-d3km)+
     *                            (d3wp-d3wm)*(d3kap-d3kam))) )
     *                                *met(2,i,j,k)*met(4,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** v-u
               pd = ( c6*( 0.5d0*( d3up*d3lp+d3um*d3lm  +
     *        dt2o12*(d3up*d3lap + d3um*d3lam+ 
     *                d3uap*d3lp + d3uam*d3lm)))
     *            + c8*( (d3up-d3um)*(d3lp-d3lm) +
     *        dt2o12*(  (d3uap-d3uam)*(d3lp-d3lm)+
     *                            (d3up-d3um)*(d3lap-d3lam))) )
     *                                *met(2,i,j,k)*met(3,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** v-v
               pd = ( c6*( 0.5d0*( d3vp*d3lp+d3vm*d3lm  +
     *        dt2o12*(d3vp*d3lap + d3vm*d3lam+ 
     *                d3vap*d3lp + d3vam*d3lm)))
     *            + c8*( (d3vp-d3vm)*(d3lp-d3lm) +
     *        dt2o12*(  (d3vap-d3vam)*(d3lp-d3lm)+
     *                            (d3vp-d3vm)*(d3lap-d3lam))) )
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(3,i,j,k)**2)
               glambda(i,j,k) = glambda(i,j,k) + pd*met(3,i,j,k)**2

*** v-w
               pd = ( c6*( 0.5d0*( d3wp*d3lp+d3wm*d3lm  +
     *        dt2o12*(d3wp*d3lap + d3wm*d3lam+ 
     *                d3wap*d3lp + d3wam*d3lm)))
     *            + c8*( (d3wp-d3wm)*(d3lp-d3lm) +
     *        dt2o12*(  (d3wap-d3wam)*(d3lp-d3lm)+
     *                            (d3wp-d3wm)*(d3lap-d3lam))) )
     *                                *met(3,i,j,k)*met(4,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** w-u
               pd = ( c6*( 0.5d0*( d3up*d3mp+d3um*d3mm  +
     *        dt2o12*(d3up*d3map + d3um*d3mam+ 
     *                d3uap*d3mp + d3uam*d3mm)))
     *            + c8*( (d3up-d3um)*(d3mp-d3mm) +
     *        dt2o12*(  (d3uap-d3uam)*(d3mp-d3mm)+
     *                            (d3up-d3um)*(d3map-d3mam))) )
     *                                *met(2,i,j,k)*met(4,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** w-v
               pd = ( c6*( 0.5d0*( d3vp*d3mp+d3vm*d3mm  +
     *        dt2o12*(d3vp*d3map + d3vm*d3mam+ 
     *                d3vap*d3mp + d3vam*d3mm)))
     *            + c8*( (d3vp-d3vm)*(d3mp-d3mm) +
     *        dt2o12*(  (d3vap-d3vam)*(d3mp-d3mm)+
     *                            (d3vp-d3vm)*(d3map-d3mam))) )
     *                                *met(3,i,j,k)*met(4,i,j,k)
               gmu(i,j,k) = gmu(i,j,k) + pd
               glambda(i,j,k) = glambda(i,j,k) + pd

*** w-w
               pd = ( c6*( 0.5d0*( d3wp*d3mp+d3wm*d3mm  +
     *        dt2o12*(d3wp*d3map + d3wm*d3mam+ 
     *                d3wap*d3mp + d3wam*d3mm)))
     *            + c8*( (d3wp-d3wm)*(d3mp-d3mm) +
     *        dt2o12*(  (d3wap-d3wam)*(d3mp-d3mm)+
     *                            (d3wp-d3wm)*(d3map-d3mam))) )
               gmu(i,j,k) = gmu(i,j,k) + pd*(mucof+met(4,i,j,k)**2)
               glambda(i,j,k) = glambda(i,j,k) + pd*met(4,i,j,k)**2

            enddo
         enddo
      enddo

      end
