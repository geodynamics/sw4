      subroutine ADDGRADRHO(ifirst, ilast, jfirst, jlast, kfirst, klast,
     *    ifirstact, ilastact, jfirstact, jlastact, kfirstact, klastact,
     *    kap, kapacc, um, u, up, uacc, grho, dt, h, onesided )
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
               grho(i,j,k) = grho(i,j,k) - (
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
      real*8  dt, h, idt, dt2o12, wgh(4), normfact

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
               grho(i,j,k) = grho(i,j,k) - (
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
