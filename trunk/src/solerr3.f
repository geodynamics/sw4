c------------------------------------------------------------
      subroutine solerr3(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, uex, u, li, l2, zmin, x0, y0, z0, radius )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 h, zmin, x0, y0, z0, radius, sradius2, dist
      real*8 uex(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li, l2, err(3)

      li = 0
      l2 = 0
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
      do k=kfirst+2,klast-2
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            dist = ((i-1)*h-x0)**2+((j-1)*h-y0)**2+((k-1)*h+zmin-z0)**2
            if( dist .gt. sradius2 )then
c exact solution in array 'uex'
               do c=1,3
                  err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
               enddo
               if( li.lt.max(err(1),err(2),err(3)) )then
                  li = max(err(1),err(2),err(3))
               endif
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
     +     h, uex, u, li, l2)
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

