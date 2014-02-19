      subroutine FACTORIZEINTERFACEMATRICES( ib, ie, jb, je, nk, 
     *                                       mu,  la, rho,
     *                      ibf, ief, jbf, jef, nkf, muf, laf, rhof,
     *                      h, dt, a1, ipiv1, a2, ipiv2 )

***********************************************************************
***
*** LU Decomposition of 5x5 matrix blocks along the mesh refinement interface.
***
***********************************************************************

      implicit none
      real*8 half
      parameter( half=0.5d0 )

      integer ni, nj, nk, nif, njf, nkf, ipiv1(5,*), ipiv2(5,*), info
      integer ib, ie, jb, je, ibf, ief, jbf, jef
      integer i, j, ind, kk, mm
      real*8   mu(ib:ie,jb:je,0:nk+1)
      real*8   la(ib:ie,jb:je,0:nk+1)
      real*8  rho(ib:ie,jb:je,0:nk+1)
      real*8   muf(ibf:ief,jbf:jef,0:nkf+1)
      real*8   laf(ibf:ief,jbf:jef,0:nkf+1)
      real*8  rhof(ibf:ief,jbf:jef,0:nkf+1)
      real*8  h(2), dt, dt2, nu2, hc, hf, a1(5,5,*), a2(5,5,*)

      hc = h(1)
      hf = h(2)

      dt2 = dt*dt
      nu2 = dt*dt/(hf*hf)

*** Setup 5x5 blocks :
      ni = ie-ib-1
      nj = je-jb-1
*** u,v components
      do i=ib+2,ie-2
         do j=jb+2,je-2
            ind = i-1-ib+(ni-2)*(j-2-jb)
            a1(1,1,ind) = -nu2*half*(mu(i,j,1)+mu(i,j,0))/(4*rho(i,j,1))
            a1(1,2,ind) = nu2*half*(muf(2*i-1,2*j-1,nkf+1)+
     *                 muf(2*i-1,2*j-1,nkf))/(rhof(2*i-1,2*j-1,nkf))
            a1(1,3,ind) = 0
            a1(1,4,ind) = 0
            a1(1,5,ind) = 0

            a1(2,1,ind) = (mu(i,j,1)+mu(i,j,0))/8
            a1(2,2,ind) = (muf(2*i-1,2*j-1,nkf+1)+muf(2*i-1,2*j-1,nkf))
     *                       /16
            a1(2,3,ind) = (muf(2*i,2*j,nkf+1)+muf(2*i,2*j,nkf))
     *                       /64
            a1(2,4,ind) = (muf(2*i-1,2*j,nkf+1)+muf(2*i-1,2*j,nkf))
     *                       /32
            a1(2,5,ind) = (muf(2*i,2*j-1,nkf+1)+muf(2*i,2*j-1,nkf))
     *                       /32

            a1(3,1,ind) = 0
            a1(3,2,ind) =-nu2*(muf(2*i-1,2*j-1,nkf+1)+
     *              muf(2*i-1,2*j-1,nkf))/(8*rhof(2*i-1,2*j-1,nkf))
            a1(3,3,ind) = nu2*(muf(2*i,2*j,nkf+1)+muf(2*i,2*j,nkf))/
     *                                  (2*rhof(2*i,2*j,nkf))
            a1(3,4,ind) = 0
            a1(3,5,ind) = 0

            a1(4,1,ind) = 0
            a1(4,2,ind) = -nu2*(muf(2*i-1,2*j-1,nkf+1)+
     *              muf(2*i-1,2*j-1,nkf))/(4*rhof(2*i-1,2*j-1,nkf))
            a1(4,3,ind) = 0
            a1(4,4,ind) = nu2*(muf(2*i-1,2*j,nkf+1)+muf(2*i-1,2*j,nkf))/
     *                                  (2*rhof(2*i-1,2*j,nkf))
            a1(4,5,ind) = 0

            a1(5,1,ind) = 0
            a1(5,2,ind) = -nu2*(muf(2*i-1,2*j-1,nkf+1)+
     *              muf(2*i-1,2*j-1,nkf))/(4*rhof(2*i-1,2*j-1,nkf))
            a1(5,3,ind) = 0
            a1(5,4,ind) = 0
            a1(5,5,ind) = nu2*(muf(2*i,2*j-1,nkf+1)+muf(2*i,2*j-1,nkf))/
     *                               (2*rhof(2*i,2*j-1,nkf))
            call DGETRF(5,5,a1(1,1,ind),5,ipiv1(1,ind),info)
            if( info.ne.0 )then
               write(*,*) 'DGETRF: A1, info = ',info
               stop
            endif
         enddo
      enddo

*** w component
      do i=ib+2,ie-2
         do j=jb+2,je-2
            ind = i-1-ib+(ni-2)*(j-2-jb)
            a2(1,1,ind) = -nu2*half*(2*mu(i,j,1)+2*mu(i,j,0)+
     *                         la(i,j,1)+la(i,j,0))/(4*rho(i,j,1))
            a2(1,2,ind) = nu2*half*(
     *           2*muf(2*i-1,2*j-1,nkf+1)+laf(2*i-1,2*j-1,nkf+1) +
     *           2*muf(2*i-1,2*j-1,nkf)+laf(2*i-1,2*j-1,nkf))/
     *                       (rhof(2*i-1,2*j-1,nkf))
            a2(1,3,ind) = 0
            a2(1,4,ind) = 0
            a2(1,5,ind) = 0

            a2(2,1,ind) = (2*mu(i,j,1)+2*mu(i,j,0)+
     *                                   la(i,j,1)+la(i,j,0))/8
            a2(2,2,ind)=  (2*muf(2*i-1,2*j-1,nkf+1)+
     *           laf(2*i-1,2*j-1,nkf+1)+ 2*muf(2*i-1,2*j-1,nkf)+
     *           laf(2*i-1,2*j-1,nkf) )/16
            a2(2,3,ind) = (2*muf(2*i,2*j,nkf+1)+2*muf(2*i,2*j,nkf)+
     *                        laf(2*i,2*j,nkf+1)+  laf(2*i,2*j,nkf) )
     *                       /64
            a2(2,4,ind) = (2*muf(2*i-1,2*j,nkf+1)+2*muf(2*i-1,2*j,nkf)
     *                       +laf(2*i-1,2*j,nkf+1)+  laf(2*i-1,2*j,nkf))
     *                       /32
            a2(2,5,ind) = (2*muf(2*i,2*j-1,nkf+1)+2*muf(2*i,2*j-1,nkf)
     *                      + laf(2*i,2*j-1,nkf+1)+  laf(2*i,2*j-1,nkf))
     *                       /32

            a2(3,1,ind) = 0
            a2(3,2,ind) =-nu2*(2*muf(2*i-1,2*j-1,nkf+1)+
     *                laf(2*i-1,2*j-1,nkf+1)+
     *              2*muf(2*i-1,2*j-1,nkf)+
     *                laf(2*i-1,2*j-1,nkf))/(8*rhof(2*i-1,2*j-1,nkf))
            a2(3,3,ind) = nu2*(2*muf(2*i,2*j,nkf+1)+laf(2*i,2*j,nkf+1)+
     *                         2*muf(2*i,2*j,nkf)+laf(2*i,2*j,nkf))/
     *                                  (2*rhof(2*i,2*j,nkf))
            a2(3,4,ind) = 0
            a2(3,5,ind) = 0

            a2(4,1,ind) = 0
            a2(4,2,ind) = -nu2*(2*muf(2*i-1,2*j-1,nkf+1)+
     *           laf(2*i-1,2*j-1,nkf+1)+2*muf(2*i-1,2*j-1,nkf)+
     *           laf(2*i-1,2*j-1,nkf))/(4*rhof(2*i-1,2*j-1,nkf))
            a2(4,3,ind) = 0
            a2(4,4,ind) = nu2*(2*muf(2*i-1,2*j,nkf+1)+
     *  laf(2*i-1,2*j,nkf+1)+ 2*muf(2*i-1,2*j,nkf)+laf(2*i-1,2*j,nkf))/
     *                                (2*rhof(2*i-1,2*j,nkf))
            a2(4,5,ind) = 0

            a2(5,1,ind) = 0
            a2(5,2,ind) = -nu2*(2*muf(2*i-1,2*j-1,nkf+1)+
     *         laf(2*i-1,2*j-1,nkf+1)+ 2*muf(2*i-1,2*j-1,nkf)+
     *         laf(2*i-1,2*j-1,nkf) )/(4*rhof(2*i-1,2*j-1,nkf))
            a2(5,3,ind) = 0
            a2(5,4,ind) = 0
            a2(5,5,ind) = nu2*(2*muf(2*i,2*j-1,nkf+1)+
     *          laf(2*i,2*j-1,nkf+1)+ 2*muf(2*i,2*j-1,nkf)+
     *          laf(2*i,2*j-1,nkf) )/(2*rhof(2*i,2*j-1,nkf))
            call DGETRF(5,5,a2(1,1,ind),5,ipiv2(1,ind),info)
            if( info.ne.0 )then
               write(*,*) 'DGETRF: A2, info = ',info
               stop
            endif
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine CIKPLANE( ni, nj, nk, up, u, um, side, fo, dt, h, 
     *                 mu, la, rho, use_supergrid, beta, dcx, dcy, dcz )

      implicit none

      real*8 half, fourth
      parameter( half = 1d0/2, fourth=1d0/4 )

      integer ni, nj, nk, i, j, k, kp, km, side, use_supergrid, c
      real*8  up(3,0:ni+1,0:nj+1), u(3,0:ni+1,0:nj+1,0:nk+1)
      real*8  um(3,0:ni+1,0:nj+1,0:nk+1)
      real*8  fo(3,0:ni+1,0:nj+1), dt2, la2, cof, dt, h
      real*8  mu(0:ni+1,0:nj+1,0:nk+1), la(0:ni+1,0:nj+1,0:nk+1)
      real*8  rho(0:ni+1,0:nj+1,0:nk+1)
      real*8  beta, dcx(0:ni+1), dcy(0:nj+1), dcz(0:nk+1)
      real*8  vyp, vym, wyp, wym, wxp, wxm, uxp, uxm
      real*8  dymudzv, dyladzw, dxmudzu, dxladzw
      real*8  dyladxu, dymudxv, dxmudyu, dxladyv
      real*8  mupx, mumx, mupy, mumy, mupz, mumz

      dt2 = dt*dt
      la2 = dt*dt/(h*h)

      if( side.eq.5 )then
         k  = 1
         kp = k+1
         km = k
      else
         k  = nk
         kp = k
         km = k-1
      endif

      do j=1,nj
         do i=1,ni

            vyp    = half*(u(2,i,j+1,kp)-u(2,i,j-1,kp))
            vym    = half*(u(2,i,j+1,km)-u(2,i,j-1,km))
            wyp    = half*(u(3,i,j+1,kp)-u(3,i,j-1,kp))
            wym    = half*(u(3,i,j+1,km)-u(3,i,j-1,km))
            wxp    = half*(u(3,i+1,j,kp)-u(3,i-1,j,kp))
            wxm    = half*(u(3,i+1,j,km)-u(3,i-1,j,km))
            uxp    = half*(u(1,i+1,j,kp)-u(1,i-1,j,kp))
            uxm    = half*(u(1,i+1,j,km)-u(1,i-1,j,km))

            dymudzv=  half*(
     *               mu(i,j+1,k)*(u(2,i,j+1,kp)-u(2,i,j+1,km))-
     *               mu(i,j-1,k)*(u(2,i,j-1,kp)-u(2,i,j-1,km)))
            dyladzw = half*(
     *               la(i,j+1,k)*(u(3,i,j+1,kp)-u(3,i,j+1,km))-
     *               la(i,j-1,k)*(u(3,i,j-1,kp)-u(3,i,j-1,km)) )
            dxmudzu = half*(
     *               mu(i+1,j,k)*(u(1,i+1,j,kp)-u(1,i+1,j,km))-
     *               mu(i-1,j,k)*(u(1,i-1,j,kp)-u(1,i-1,j,km)) )
            dxladzw=  half*(
     *               la(i+1,j,k)*(u(3,i+1,j,kp)-u(3,i+1,j,km))-
     *               la(i-1,j,k)*(u(3,i-1,j,kp)-u(3,i-1,j,km)))

            dyladxu = half*(
     *             la(i,j+1,k)*half*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k)) -
     *             la(i,j-1,k)*half*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k)) )
            dymudxv = half*(
     *             mu(i,j+1,k)*half*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k)) -
     *             mu(i,j-1,k)*half*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k)) )
            dxmudyu = half*(
     *             mu(i+1,j,k)*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))*half-
     *             mu(i-1,j,k)*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))*half )
            dxladyv = half*(
     *             la(i+1,j,k)*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))*half -
     *             la(i-1,j,k)*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))*half )

            cof = la2/rho(i,j,k)
            mupx = half*(mu(i,j,k)+mu(i+1,j,k))
            mumx = half*(mu(i,j,k)+mu(i-1,j,k))
            mupy = half*(mu(i,j+1,k)+mu(i,j,k))
            mumy = half*(mu(i,j-1,k)+mu(i,j,k))
            mupz = half*(mu(i,j,k+1)+mu(i,j,k))
            mumz = half*(mu(i,j,k-1)+mu(i,j,k))

            up(1,i,j) = 2*u(1,i,j,k)-um(1,i,j,k) + cof*( 
     * (2*mupx+half*(la(i+1,j,k)+la(i,j,k)))
     *                                *(u(1,i+1,j,k)-u(1,i,j,k)) -
     * (2*mumx+half*(la(i-1,j,k)+la(i,j,k)))
     *                                *(u(1,i,j,k)-u(1,i-1,j,k)) +
     *  dxladyv + dxladzw + dymudxv + mu(i,j,kp)*wxp - mu(i,j,km)*wxm +
     *              mupy*(u(1,i,j+1,k)-u(1,i,j,k)) - 
     *              mumy*(u(1,i,j,k)-u(1,i,j-1,k)) +
     *              mupz*(u(1,i,j,k+1)-u(1,i,j,k)) - 
     *              mumz*(u(1,i,j,k)-u(1,i,j,k-1)) )  + 
     *                      dt2*fo(1,i,j)

            up(2,i,j) = 2*u(2,i,j,k)-um(2,i,j,k) + cof*(
     *              mupx*(u(2,i+1,j,k)-u(2,i,j,k))-
     *              mumx*(u(2,i,j,k)-u(2,i-1,j,k))+
     *  dxmudyu + dyladxu + dyladzw + mu(i,j,kp)*wyp - mu(i,j,km)*wym + 
     *          (2*mupy+ half*(la(i,j,k)+la(i,j+1,k)))
     *                              *(u(2,i,j+1,k)-u(2,i,j,k))-
     *          (2*mumy+ half*(la(i,j,k)+la(i,j-1,k)))
     *                              *(u(2,i,j,k)-u(2,i,j-1,k)) +
     *      mupz*(u(2,i,j,k+1)-u(2,i,j,k)) - 
     *      mumz*(u(2,i,j,k)-u(2,i,j,k-1)) )+
     *            dt2*fo(2,i,j)

            up(3,i,j) = 2*u(3,i,j,k)-um(3,i,j,k)+cof*(
     * mupx*(u(3,i+1,j,k)-u(3,i,j,k))-mumx*(u(3,i,j,k)-u(3,i-1,j,k))+
     *  dxmudzu + dymudzv + la(i,j,kp)*(uxp+vyp) - la(i,j,km)*(uxm+vym)+ 
     *           mupy*(u(3,i,j+1,k)-u(3,i,j,k))-
     *           mumy*(u(3,i,j,k)-u(3,i,j-1,k))+
     *           (2*mupz+half*(la(i,j,k+1)+la(i,j,k)) )*
     *               (u(3,i,j,k+1)-u(3,i,j,k)) - 
     *           (2*mumz+half*(la(i,j,k-1)+la(i,j,k)) )*
     *               (u(3,i,j,k)-u(3,i,j,k-1)) ) +
     *            dt2*fo(3,i,j)
            if( use_supergrid.eq.1 )then
               do c=1,3
                  up(c,i,j) = up(c,i,j) + beta*(
     *         + dcx(i+1)*(u(c,i+1,j,k)-u(c,i,j,k))   
     *         - dcx(i)*  (u(c,i,j,k)-u(c,i-1,j,k))
     *         - dcx(i+1)*(um(c,i+1,j,k)-um(c,i,j,k)) 
     *         + dcx(i)*  (um(c,i,j,k)-um(c,i-1,j,k))
     *         + dcy(j+1)*(u(c,i,j+1,k)-u(c,i,j,k))  
     *         - dcy(j)*  (u(c,i,j,k)-u(c,i,j-1,k))
     *         - dcy(j+1)*(um(c,i,j+1,k)-um(c,i,j,k)) 
     *         + dcy(j)*  (um(c,i,j,k)-um(c,i,j-1,k))
     *         + dcz(k+1)*(u(c,i,j,k+1)-u(c,i,j,k))
     *         - dcz(k)*  (u(c,i,j,k)-u(c,i,j,k-1))
     *         - dcz(k+1)*(um(c,i,j,k+1)-um(c,i,j,k)) 
     *         + dcz(k)*  (um(c,i,j,k)-um(c,i,j,k-1)) )
               enddo
            endif
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine CIBCKPLANE( ni, nj, u, bctype, bd1, bd2, bd3, bd4 )

      implicit none

      integer ni, nj, i, j, bctype(4), c
      real*8 u(3,0:ni+1,0:nj+1)
      real*8 bd1(3,0:nj+1), bd2(3,0:nj+1)
      real*8 bd3(3,0:ni+1), bd4(3,0:ni+1)

      if( bctype(1).eq.1 )then
         i = 1
         do j=1,nj
            do c=1,3
               u(c,i,j) = bd1(c,j)
            enddo
         enddo
      endif

      if( bctype(2).eq.1 )then
         i = ni
         do j=1,nj
            do c=1,3
               u(c,i,j) = bd2(c,j)
            enddo
         enddo
      endif

      if( bctype(3).eq.1 )then
         j = 1
         do i=1,ni
            do c=1,3
               u(c,i,j) = bd3(c,i)
            enddo
         enddo
      endif

      if( bctype(4).eq.1 )then
         j = nj
         do i=1,ni
            do c=1,3
               u(c,i,j) = bd4(c,i)
            enddo
         enddo
      endif
      
      end

c-----------------------------------------------------------------------
      subroutine CIBCRESIDUAL( ni,  nj,  nk,  u,  mu,  la, 
     *                       bres, side, h )

      implicit none

      real*8 half, fourth
      parameter( half = 1d0/2, fourth=1d0/4 )

      integer ni, nj, nk, nif, njf, nkf, i, j, k, side
      real*8  u(3,0:ni+1,0:nj+1,0:nk+1), mu(0:ni+1,0:nj+1,0:nk+1)
      real*8  la(0:ni+1,0:nj+1,0:nk+1)
      real*8  bres(3,0:ni+1,0:nj+1), h, hi

      hi = 1/h
      if( side.eq.5 )then
         k = 1
      elseif( side.eq.6 )then
         k = nk
      endif
      do j=1,nj
         do i=1,ni
            bres(1,i,j) = hi*(fourth*( (mu(i,j,k)+mu(i,j,k+1))*
     *                                    (u(1,i,j,k+1)-u(1,i,j,k)) +
     *                             (mu(i,j,k)+mu(i,j,k-1))* 
     *                                    (u(1,i,j,k)-u(1,i,j,k-1))) + 
     *                   half*mu(i,j,k)*(u(3,i+1,j,k)-u(3,i-1,j,k)))
 
            bres(2,i,j) = hi*(fourth*( (mu(i,j,k)+mu(i,j,k+1))*
     *                                    (u(2,i,j,k+1)-u(2,i,j,k)) +
     *                             (mu(i,j,k)+mu(i,j,k-1))* 
     *                                    (u(2,i,j,k)-u(2,i,j,k-1))) + 
     *                   half*mu(i,j,k)*(u(3,i,j+1,k)-u(3,i,j-1,k)) )

            bres(3,i,j) = hi*(fourth*( (2*mu(i,j,k)+2*mu(i,j,k+1)+
     *                                la(i,j,k)+la(i,j,k+1))*
     *                                    (u(3,i,j,k+1)-u(3,i,j,k)) +
     *                             (2*mu(i,j,k)+2*mu(i,j,k-1)+ 
     *                                la(i,j,k)+la(i,j,k-1))*
     *                                    (u(3,i,j,k)-u(3,i,j,k-1))) + 
     *                half*la(i,j,k)*(
     *      u(1,i+1,j,k)-u(1,i-1,j,k)+u(2,i,j+1,k)-u(2,i,j-1,k) ))
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine CALCNORM( ni, nj, nk, utmp, u )

      implicit none
      integer ni, nj, nk, i, j
      real*8  utmp(3,0:ni+1,0:nj+1,0:nk+1)
      real*8  u(3,0:ni+1,0:nj+1), nrmu

      nrmu = 0
      do j=2,nj-1
         do i=2,ni-1
            nrmu = nrmu + abs(utmp(1,i,j,nk)-u(1,i,j))
         enddo
      enddo
      write(*,*) 'Norm of difference is ', nrmu

      end
c-----------------------------------------------------------------------
      subroutine CIBCKPLANE2( ni, nj, nk, u, k, bctype, bd1, bd2, bd3, 
     *     bd4 )

      implicit none

      integer ni, nj, nk, i, j, bctype(4), c, k
      real*8 u(3,0:ni+1,0:nj+1)
      real*8 bd1(3,0:nj+1,0:nk+1), bd2(3,0:nj+1,0:nk+1)
      real*8 bd3(3,0:ni+1,0:nk+1), bd4(3,0:ni+1,0:nk+1)

      if( bctype(1).eq.1 )then
         i = 1
         do j=1,nj
            do c=1,3
               u(c,i,j) = bd1(c,j,k)
            enddo
         enddo
      endif

      if( bctype(2).eq.1 )then
         i = ni
         do j=1,nj
            do c=1,3
               u(c,i,j) = bd2(c,j,k)
            enddo
         enddo
      endif

      if( bctype(3).eq.1 )then
         j = 1
         do i=1,ni
            do c=1,3
               u(c,i,j) = bd3(c,i,k)
            enddo
         enddo
      endif

      if( bctype(4).eq.1 )then
         j = nj
         do i=1,ni
            do c=1,3
               u(c,i,j) = bd4(c,i,k)
            enddo
         enddo
      endif
      
      end
c-----------------------------------------------------------------------
      subroutine printmatrix( m, n, a )

      implicit none
      integer m, n, i, j
      real*8  a(m,n)

      do i=1,m
         write(*,101) (a(i,j),j=1,n)
 101     format(' ', 12(g15.7,tr2))
      enddo

      end


      
c-----------------------------------------------------------------------
      subroutine CIKPLANEATT( ni, nj, nk, u, lu, side, h, 
     *                 mu, la, rho )

***********************************************************************
***
*** Accumulate operator L(u), routine computes: lu := lu - L(u).
*** lu is defined on a single coordinate plane, for use at grid refinement
*** boundaries.
***
***********************************************************************
      implicit none

      real*8 half, fourth
      parameter( half = 1d0/2, fourth=1d0/4 )

      integer ni, nj, nk, i, j, k, kp, km, side, c
      real*8  lu(3,0:ni+1,0:nj+1), u(3,0:ni+1,0:nj+1,0:nk+1)
      real*8  cof, h
      real*8  mu(0:ni+1,0:nj+1,0:nk+1), la(0:ni+1,0:nj+1,0:nk+1)
      real*8  rho(0:ni+1,0:nj+1,0:nk+1)
      real*8  vyp, vym, wyp, wym, wxp, wxm, uxp, uxm
      real*8  dymudzv, dyladzw, dxmudzu, dxladzw
      real*8  dyladxu, dymudxv, dxmudyu, dxladyv
      real*8  mupx, mumx, mupy, mumy, mupz, mumz

      if( side.eq.5 )then
         k  = 1
         kp = k+1
         km = k
      else
         k  = nk
         kp = k
         km = k-1
      endif

      do j=1,nj
         do i=1,ni

            vyp    = half*(u(2,i,j+1,kp)-u(2,i,j-1,kp))
            vym    = half*(u(2,i,j+1,km)-u(2,i,j-1,km))
            wyp    = half*(u(3,i,j+1,kp)-u(3,i,j-1,kp))
            wym    = half*(u(3,i,j+1,km)-u(3,i,j-1,km))
            wxp    = half*(u(3,i+1,j,kp)-u(3,i-1,j,kp))
            wxm    = half*(u(3,i+1,j,km)-u(3,i-1,j,km))
            uxp    = half*(u(1,i+1,j,kp)-u(1,i-1,j,kp))
            uxm    = half*(u(1,i+1,j,km)-u(1,i-1,j,km))

            dymudzv=  half*(
     *               mu(i,j+1,k)*(u(2,i,j+1,kp)-u(2,i,j+1,km))-
     *               mu(i,j-1,k)*(u(2,i,j-1,kp)-u(2,i,j-1,km)))
            dyladzw = half*(
     *               la(i,j+1,k)*(u(3,i,j+1,kp)-u(3,i,j+1,km))-
     *               la(i,j-1,k)*(u(3,i,j-1,kp)-u(3,i,j-1,km)) )
            dxmudzu = half*(
     *               mu(i+1,j,k)*(u(1,i+1,j,kp)-u(1,i+1,j,km))-
     *               mu(i-1,j,k)*(u(1,i-1,j,kp)-u(1,i-1,j,km)) )
            dxladzw=  half*(
     *               la(i+1,j,k)*(u(3,i+1,j,kp)-u(3,i+1,j,km))-
     *               la(i-1,j,k)*(u(3,i-1,j,kp)-u(3,i-1,j,km)))

            dyladxu = half*(
     *             la(i,j+1,k)*half*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k)) -
     *             la(i,j-1,k)*half*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k)) )
            dymudxv = half*(
     *             mu(i,j+1,k)*half*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k)) -
     *             mu(i,j-1,k)*half*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k)) )
            dxmudyu = half*(
     *             mu(i+1,j,k)*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))*half-
     *             mu(i-1,j,k)*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))*half )
            dxladyv = half*(
     *             la(i+1,j,k)*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))*half -
     *             la(i-1,j,k)*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))*half )

            cof = 1/(h*h*rho(i,j,k))
            mupx = half*(mu(i,j,k)+mu(i+1,j,k))
            mumx = half*(mu(i,j,k)+mu(i-1,j,k))
            mupy = half*(mu(i,j+1,k)+mu(i,j,k))
            mumy = half*(mu(i,j-1,k)+mu(i,j,k))
            mupz = half*(mu(i,j,k+1)+mu(i,j,k))
            mumz = half*(mu(i,j,k-1)+mu(i,j,k))

            lu(1,i,j) = lu(1,i,j) - cof*( 
     * (2*mupx+half*(la(i+1,j,k)+la(i,j,k)))
     *                                *(u(1,i+1,j,k)-u(1,i,j,k)) -
     * (2*mumx+half*(la(i-1,j,k)+la(i,j,k)))
     *                                *(u(1,i,j,k)-u(1,i-1,j,k)) +
     *  dxladyv + dxladzw + dymudxv + mu(i,j,kp)*wxp - mu(i,j,km)*wxm +
     *              mupy*(u(1,i,j+1,k)-u(1,i,j,k)) - 
     *              mumy*(u(1,i,j,k)-u(1,i,j-1,k)) +
     *              mupz*(u(1,i,j,k+1)-u(1,i,j,k)) - 
     *              mumz*(u(1,i,j,k)-u(1,i,j,k-1)) )

            lu(2,i,j) = lu(2,i,j) - cof*(
     *              mupx*(u(2,i+1,j,k)-u(2,i,j,k))-
     *              mumx*(u(2,i,j,k)-u(2,i-1,j,k))+
     *  dxmudyu + dyladxu + dyladzw + mu(i,j,kp)*wyp - mu(i,j,km)*wym + 
     *          (2*mupy+ half*(la(i,j,k)+la(i,j+1,k)))
     *                              *(u(2,i,j+1,k)-u(2,i,j,k))-
     *          (2*mumy+ half*(la(i,j,k)+la(i,j-1,k)))
     *                              *(u(2,i,j,k)-u(2,i,j-1,k)) +
     *      mupz*(u(2,i,j,k+1)-u(2,i,j,k)) - 
     *      mumz*(u(2,i,j,k)-u(2,i,j,k-1)) )

            lu(3,i,j) = lu(3,i,j) - cof*(
     * mupx*(u(3,i+1,j,k)-u(3,i,j,k))-mumx*(u(3,i,j,k)-u(3,i-1,j,k))+
     *  dxmudzu + dymudzv + la(i,j,kp)*(uxp+vyp) - la(i,j,km)*(uxm+vym)+ 
     *           mupy*(u(3,i,j+1,k)-u(3,i,j,k))-
     *           mumy*(u(3,i,j,k)-u(3,i,j-1,k))+
     *           (2*mupz+half*(la(i,j,k+1)+la(i,j,k)) )*
     *               (u(3,i,j,k+1)-u(3,i,j,k)) - 
     *           (2*mumz+half*(la(i,j,k-1)+la(i,j,k)) )*
     *               (u(3,i,j,k)-u(3,i,j,k-1)) ) 

         enddo
      enddo

      end
      
c-----------------------------------------------------------------------
      subroutine CIBCRESIDUALATT( ni,  nj,  nk,  u,  mu,  la, 
     *     bres, side, h )

***********************************************************************
***
*** Computes bres := bres - sigma(u)
***  can be called repeatedly for each memory variable 
***  in an attenuation model.
***
***********************************************************************

      implicit none

      real*8 half, fourth
      parameter( half = 1d0/2, fourth=1d0/4 )

      integer ni, nj, nk, nif, njf, nkf, i, j, k, side
      real*8  u(3,0:ni+1,0:nj+1,0:nk+1), mu(0:ni+1,0:nj+1,0:nk+1)
      real*8  la(0:ni+1,0:nj+1,0:nk+1)
      real*8  bres(3,0:ni+1,0:nj+1), h, hi

      hi = 1/h
      if( side.eq.5 )then
         k = 1
      elseif( side.eq.6 )then
         k = nk
      endif
      do j=1,nj
         do i=1,ni
            bres(1,i,j) = bres(1,i,j) - 
     *                   hi*(fourth*( (mu(i,j,k)+mu(i,j,k+1))*
     *                                    (u(1,i,j,k+1)-u(1,i,j,k)) +
     *                             (mu(i,j,k)+mu(i,j,k-1))* 
     *                                    (u(1,i,j,k)-u(1,i,j,k-1))) + 
     *                   half*mu(i,j,k)*(u(3,i+1,j,k)-u(3,i-1,j,k)))
 
            bres(2,i,j) = bres(2,i,j) -
     *                         hi*(fourth*( (mu(i,j,k)+mu(i,j,k+1))*
     *                                    (u(2,i,j,k+1)-u(2,i,j,k)) +
     *                             (mu(i,j,k)+mu(i,j,k-1))* 
     *                                    (u(2,i,j,k)-u(2,i,j,k-1))) + 
     *                   half*mu(i,j,k)*(u(3,i,j+1,k)-u(3,i,j-1,k)) )

            bres(3,i,j) = bres(3,i,j) - 
     *                        hi*(fourth*( (2*mu(i,j,k)+2*mu(i,j,k+1)+
     *                                        la(i,j,k)+  la(i,j,k+1))*
     *                                    (u(3,i,j,k+1)-u(3,i,j,k)) +
     *                             (2*mu(i,j,k)+2*mu(i,j,k-1)+ 
     *                                la(i,j,k)+  la(i,j,k-1))*
     *                                    (u(3,i,j,k)-u(3,i,j,k-1))) + 
     *                half*la(i,j,k)*(
     *      u(1,i+1,j,k)-u(1,i-1,j,k)+u(2,i,j+1,k)-u(2,i,j-1,k) ) ) 
         enddo
      enddo

      end

