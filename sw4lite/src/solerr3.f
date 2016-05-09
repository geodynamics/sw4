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
c------------------------------------------------------------
      subroutine solerr3(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, uex, u, li, l2, xli, zmin, x0, y0, z0, radius,
     +     imin, imax, jmin, jmax, kmin, kmax)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer imin, imax, jmin, jmax, kmin, kmax
      real*8 h, zmin, x0, y0, z0, radius, sradius2, dist
      real*8 uex(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li, l2, err(3), xli

      li = 0
      l2 = 0
c max norm of exact solution
      xli = 0
c tmp
c      write(*,*)'if=', ifirst, 'il=', ilast, 'jf=', jfirst, 'jl=',
c     +     jlast, 'kf=', kfirst, 'kl=', klast
c this test includes interior points (the arrays include an extra ghost point in k)
c
c Possibility to exclude all points within the distance 'radius' from the
c point '(x0,y0,z0)', set radius to a negative value to include all points.
c
      sradius2 = radius*radius
      if( radius.lt.0 )then
         sradius2 = -sradius2
      endif
C       do k=kfirst+2,klast-2
C         do j=jfirst+2,jlast-2
C           do i=ifirst+2,ilast-2
      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            dist = ((i-1)*h-x0)**2+((j-1)*h-y0)**2+((k-1)*h+zmin-z0)**2
            if( dist .gt. sradius2 )then
c exact solution in array 'uex'
               do c=1,3
                  err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
               enddo
               if( li.lt.max(err(1),err(2),err(3)) )then
                  li = max(err(1),err(2),err(3))
               endif
               if( xli.lt.max(uex(1,i,j,k),uex(2,i,j,k),uex(3,i,j,k)) )
     +              xli = max(uex(1,i,j,k),uex(2,i,j,k),uex(3,i,j,k))

               l2 = l2 + 
     +              h*h*h* (err(1)**2 + err(2)**2 + err(3)**2)
            endif
          enddo
        enddo
      enddo
c$$$      write(*,101) 'Solution errors in max- and L2-norm: ', li, l2
c$$$ 101  format(' ', a, 2(g15.7,tr2))
      return
      end

c------------------------------------------------------------
      subroutine solerrgp(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, uex, u, li, l2 )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 h
      real*8 uex(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li, l2, err(3)

      li = 0
      l2 = 0
c tmp
c      write(*,*)'if=', ifirst, 'il=', ilast, 'jf=', jfirst, 'jl=',
c     +     jlast, 'kf=', kfirst, 'kl=', klast
c this test only includes the ghost points
      k=kfirst+1
      do j=jfirst+2,jlast-2
        do i=ifirst+2,ilast-2
c exact solution in array 'uex'
          do c=1,3
            err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
          enddo
          if( li.lt.max(err(1),err(2),err(3)) )then
            li = max(err(1),err(2),err(3))
          endif
          l2 = l2 + 
     +         h*h*h* (err(1)**2 + err(2)**2 + err(3)**2)
        enddo
      enddo

      k=klast-1
      do j=jfirst+2,jlast-2
        do i=ifirst+2,ilast-2
c exact solution in array 'uex'
          do c=1,3
            err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
          enddo
          if( li.lt.max(err(1),err(2),err(3)) )then
            li = max(err(1),err(2),err(3))
          endif
          l2 = l2 + 
     +         h*h*h* (err(1)**2 + err(2)**2 + err(3)**2)
        enddo
      enddo
c$$$      write(*,101) 'Solution errors in max- and L2-norm: ', li, l2
c$$$ 101  format(' ', a, 2(g15.7,tr2))
      return
      end

c-----------------------------------------------------------------------
      subroutine solerr3c(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     uex, u, x, y, z, jac, li, l2, xli, x0, y0, z0, radius,
     +     imin, imax, jmin, jmax, kmin, kmax, usesg, strx, stry )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer imin, imax, jmin, jmax, kmin, kmax, usesg
      real*8 x0, y0, z0, radius, sradius2, dist
      real*8 uex(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 x(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 y(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      integer c, k, j, i
      real*8 li, l2, err(3), xli

      li = 0
      l2 = 0
c max norm of exact solution
      xli = 0
c tmp
c      write(*,*)'if=', ifirst, 'il=', ilast, 'jf=', jfirst, 'jl=',
c     +     jlast, 'kf=', kfirst, 'kl=', klast
c this test includes interior points (the arrays include an extra ghost point in k)
c
c Possibility to exclude all points within the distance 'radius' from the
c point '(x0,y0,z0)', set radius to a negative value to include all points.
c
      sradius2 = radius*radius
      if( radius.lt.0 )then
         sradius2 = -sradius2
      endif
C       do k=kfirst+2,klast-2
C         do j=jfirst+2,jlast-2
C           do i=ifirst+2,ilast-2
      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            dist = (x(i,j,k)-x0)**2+(y(i,j,k)-y0)**2+(z(i,j,k)-z0)**2
            if( dist .gt. sradius2 )then
c exact solution in array 'uex'
               do c=1,3
                  err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
               enddo
               if( li.lt.max(err(1),err(2),err(3)) )then
                  li = max(err(1),err(2),err(3))
               endif
               if( xli.lt.max(uex(1,i,j,k),uex(2,i,j,k),uex(3,i,j,k)) )
     +              xli = max(uex(1,i,j,k),uex(2,i,j,k),uex(3,i,j,k))

               if( usesg.ne.1 )then
                  l2 = l2 + 
     +              jac(i,j,k)*(err(1)**2 + err(2)**2 + err(3)**2)
               else
                  l2 = l2 + 
     +  jac(i,j,k)*(err(1)**2 + err(2)**2 + err(3)**2)/(strx(i)*stry(j))
               endif
            endif
          enddo
        enddo
      enddo
c$$$      write(*,101) 'Solution errors in max- and L2-norm: ', li, l2
c$$$ 101  format(' ', a, 2(g15.7,tr2))
      return
      end

c-----------------------------------------------------------------------
      subroutine meterr4c(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     met, metex, jac, jacex, li, l2, 
     +     imin, imax, jmin, jmax, kmin, kmax, h )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer imin, imax, jmin, jmax, kmin, kmax
      real*8 metex(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jacex(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 h

      integer c, k, j, i
      real*8 li(5), l2(5), err(5)

      do c=1,5
         li(c) = 0
         l2(c) = 0
      enddo
      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
             do c=1,4
                err(c) = ABS( met(c,i,j,k) - metex(c,i,j,k) )/sqrt(h)
             enddo
             err(5) = ABS( jac(i,j,k)-jacex(i,j,k))/(h*h*h)
             do c=1,5
                if( li(c).lt.err(c) )then
                   li(c) = err(c)
                endif
                l2(c) = l2(c) + jacex(i,j,k)*(err(c)**2)
             enddo
          enddo
       enddo
      enddo
      return
      end
