c ----------------------------------------------------------------------
       ! bcsf.f
      subroutine BC_FREE_SURFACE_F( ni, nj, nk, u, rho, mu, la, side,
     *                              curlcoef, h, bforcing )
c This routine will fail if mu=0 and curlcoeff=0, or if mu=0 and lambda=0, at any point along the boundary
c Be aware of setting all velocities to zero in supergrid!
      implicit none

      real*8 half
      parameter( half = 1d0/2 )

      integer ni, nj, nk, i, j, side
      real*8 u(3,ni,nj,nk), rho(ni,nj,nk), mu(ni,nj,nk), la(ni,nj,nk)
      real*8 imy, curlcoef, ga, gap, gam, du, dv, bforcing(3,ni,nj), h
      logical onesided_at_corners

      if( side.ne.4 )then
c         write(*,*) 'Free surface condition only implemented for top side'
      endif
      onesided_at_corners = .false.
      if( onesided_at_corners )then
         do j=2,nj-1
            do i=2,ni-1
               if( mu(i,j,1)+mu(i,j,2).gt.0 )then
                  imy = 2d0/(mu(i,j,1)+mu(i,j,2))
                  if( i.gt.2 .and. i.lt.ni-1 )then
                     u(1,i,j,1) = u(1,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(1,i,j,3)-u(1,i,j,2))+
     *                 imy*mu(i,j,2)*(u(3,i+1,j,2)-u(3,i-1,j,2))
                  elseif( i.eq.2 )then
                     u(1,i,j,1) = u(1,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(1,i,j,3)-u(1,i,j,2))+
     *                 imy*mu(i,j,2)*(u(3,i+1,j,2)-u(3,i,j,2))*2
                  else
                     u(1,i,j,1) = u(1,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(1,i,j,3)-u(1,i,j,2))+
     *                 imy*mu(i,j,2)*(u(3,i,j,2)-u(3,i-1,j,2))*2
                  endif
                  if( j.gt.2 .and. j.lt.nj-1 )then
                     u(2,i,j,1) = u(2,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(2,i,j,3)-u(2,i,j,2))+
     *            imy*mu(i,j,2)*(u(3,i,j+1,2)-u(3,i,j-1,2))
                  elseif( j.eq.2 )then
                     u(2,i,j,1) = u(2,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(2,i,j,3)-u(2,i,j,2))+
     *            imy*mu(i,j,2)*(u(3,i,j+1,2)-u(3,i,j,2))*2
                  else
                     u(2,i,j,1) = u(2,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(2,i,j,3)-u(2,i,j,2))+
     *            imy*mu(i,j,2)*(u(3,i,j,2)-u(3,i,j-1,2))*2
                  endif
               elseif( curlcoef.ne.0 )then
                  ga = 0
                  gam= 0
                  gap= 0
                  if( mu(i,j,2).eq.0 )then
                     ga = curlcoef
                  endif
                  if( mu(i,j,1).eq.0 )then
                     gam = curlcoef
                  endif
                  if( mu(i,j,3).eq.0 )then
                     gap = curlcoef
                  endif
                  imy = 2d0/(gam+ga)
c                  if( i.gt.2 .and. i.lt.ni-1 )then
                     u(1,i,j,1) = u(1,i,j,2) +
     *                    half*imy*(gap+ga)*(u(1,i,j,3)-u(1,i,j,2))-
     *                    imy*ga*(u(3,i+1,j,2)-u(3,i-1,j,2))
c                  elseif( i.eq.2 )then
c                     u(1,i,j,1) = u(1,i,j,2) +
c     *                    half*imy*(gap+ga)*(u(1,i,j,3)-u(1,i,j,2))-
c     *                    imy*ga*(u(3,i+1,j,2)-u(3,i,j,2))*2
c                  else
c                     u(1,i,j,1) = u(1,i,j,2) +
c     *                    half*imy*(gap+ga)*(u(1,i,j,3)-u(1,i,j,2))-
c     *                    imy*ga*(u(3,i,j,2)-u(3,i-1,j,2))*2
c                  endif
c                  if( j.gt.2 .and. j.lt.nj-1 )then
                     u(2,i,j,1) = u(2,i,j,2) +
     *                    half*imy*(gap+ga)*(u(2,i,j,3)-u(2,i,j,2))-
     *                    imy*ga*(u(3,i,j+1,2)-u(3,i,j-1,2))
c                  elseif( j.eq.2 )then
c                     u(2,i,j,1) = u(2,i,j,2) +
c     *           half*imy*(gap+ga)*(u(2,i,j,3)-u(2,i,j,2))-
c     *                imy*ga*(u(3,i,j+1,2)-u(3,i,j,2))*2
c                  else
c                     u(2,i,j,1) = u(2,i,j,2) +
c     *           half*imy*(gap+ga)*(u(2,i,j,3)-u(2,i,j,2))-
c     *                imy*ga*(u(3,i,j,2)-u(3,i,j-1,2))*2
c                  endif
               endif
               u(1,i,j,1) = u(1,i,j,1) - 2*h*bforcing(1,i,j)*imy
               u(2,i,j,1) = u(2,i,j,1) - 2*h*bforcing(2,i,j)*imy
               if( i.gt.2 .and. i.lt.ni-1 )then
                  du = u(1,i+1,j,2)-u(1,i-1,j,2)
               elseif( i.eq.2 )then
                  du = (u(1,i+1,j,2)-u(1,i,j,2))*2
               else
                  du = (u(1,i,j,2)-u(1,i-1,j,2))*2
               endif
               if( j.gt.2 .and. j.lt.nj-1 )then
                  dv = u(2,i,j+1,2)-u(2,i,j-1,2)
               elseif( j.eq.2 )then
                  dv = (u(2,i,j+1,2)-u(2,i,j,2))*2
               else
                  dv = (u(2,i,j,2)-u(2,i,j-1,2))*2
               endif
               if( curlcoef.ne.0 )then
                  du = u(1,i+1,j,2)-u(1,i-1,j,2)
                  dv = u(2,i,j+1,2)-u(2,i,j-1,2)
               endif
               imy =1d0/(mu(i,j,1)+mu(i,j,2)+half*(la(i,j,1)+la(i,j,2)))
               u(3,i,j,1) = u(3,i,j,2) +
     *            imy*(mu(i,j,3)+mu(i,j,2)+half*(la(i,j,3)+la(i,j,2)))*
     *                        (u(3,i,j,3)-u(3,i,j,2)) +
     *             imy*la(i,j,2)*( du + dv ) -2*h*bforcing(3,i,j)*imy
            enddo
         enddo
      else
         do j=2,nj-1
            do i=2,ni-1
               if( mu(i,j,1)+mu(i,j,2).gt.0 )then
                  imy = 2d0/(mu(i,j,1)+mu(i,j,2))
                  u(1,i,j,1) = u(1,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(1,i,j,3)-u(1,i,j,2))+
     *                 imy*mu(i,j,2)*(u(3,i+1,j,2)-u(3,i-1,j,2))
                  u(2,i,j,1) = u(2,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(2,i,j,3)-u(2,i,j,2))+
     *                 imy*mu(i,j,2)*(u(3,i,j+1,2)-u(3,i,j-1,2))
               elseif( curlcoef.ne.0 )then
                  ga = 0
                  gam= 0
                  gap= 0
                  if( mu(i,j,2).eq.0 )then
                     ga = curlcoef
                  endif
                  if( mu(i,j,1).eq.0 )then
                     gam = curlcoef
                  endif
                  if( mu(i,j,3).eq.0 )then
                     gap = curlcoef
                  endif
                  imy = 2d0/(gam+ga)
                  u(1,i,j,1) = u(1,i,j,2) +
     *                 half*imy*(gap+ga)*(u(1,i,j,3)-u(1,i,j,2))-
     *                 imy*ga*(u(3,i+1,j,2)-u(3,i-1,j,2))
                  u(2,i,j,1) = u(2,i,j,2) +
     *                 half*imy*(gap+ga)*(u(2,i,j,3)-u(2,i,j,2))-
     *                 imy*ga*(u(3,i,j+1,2)-u(3,i,j-1,2))
               endif
               imy =1d0/(mu(i,j,1)+mu(i,j,2)+half*(la(i,j,1)+la(i,j,2)))
               u(3,i,j,1) = u(3,i,j,2) +
     *             imy*(mu(i,j,3)+mu(i,j,2)+half*(la(i,j,3)+la(i,j,2)))*
     *              (u(3,i,j,3)-u(3,i,j,2)) +
     *              imy*la(i,j,2)*( u(1,i+1,j,2)-u(1,i-1,j,2)+
     *              u(2,i,j+1,2)-u(2,i,j-1,2))
            enddo
         enddo
      endif

      end

c-----------------------------------------------------------------------
      subroutine BC_HIGDON2_F( ni, nj, nk, up, u, um, rho, mu, la, side,
     *                         dt, h )

      implicit none

      integer ni, nj, nk, side, wind(6), oi, oj, ok
      integer i, j, k, c
      real*8  up(3,ni,nj,nk), u(3,ni,nj,nk), um(3,ni,nj,nk)
      real*8  rho(ni,nj,nk), la(ni,nj,nk), mu(ni,nj,nk), dt, h
      real*8  nu, cp, cs, csp, a, b, dti, hi, forc(3)
      real*8  vm, vp, vpp, vmp, cofi

      a = 0.0
      b = 0.0
      oi=0
      oj=0
      ok=0

      forc(1) = 0
      forc(2) = 0
      forc(3) = 0

      wind(1) = 1
      wind(2) = ni
      wind(3) = 1
      wind(4) = nj
      wind(5) = 1
      wind(6) = nk
      if( side.eq.0 )then
         oi =  1
         wind(2) = 1
      elseif( side.eq.1 )then
         oi = -1
         wind(1) = ni
      elseif( side.eq.2 )then
         oj =  1
         wind(4) = 1
      elseif( side.eq.3 )then
         oj = -1
         wind(3) = nj
      elseif( side.eq.4 )then
         ok =  1
         wind(6) = 1
      elseif( side.eq.5 )then
         ok = -1
         wind(5) = nk
      endif

      dti = 1/dt
      hi  = 1/h

      do k=wind(5),wind(6)
         do j=wind(3),wind(4)
            do i=wind(1),wind(2)
               cp = SQRT(  (2*(mu(i+oi,j+oj,k+ok)+mu(i,j,k))+
     *                         la(i+oi,j+oj,k+ok)+la(i+oi,j+oj,k+ok))/
     *                       (rho(i+oi,j+oj,k+ok)+rho(i,j,k))        )
               cs = SQRT( (mu(i+oi,j+oj,k+ok)+mu(i,j,k))/
     *                    (rho(i+oi,j+oj,k+ok)+rho(i,j,k)) )
               csp=SQRT((mu(i+oi,j+oj,k+ok)+mu(i+2*oi,j+2*oj,k+2*ok))
     *             /(rho(i+oi,j+oj,k+ok)+rho(i+2*oi,j+2*oj,k+2*ok)) )
               do c=1,3
                  vm  = (1-a)*(u(c,i,j,k)-um(c,i,j,k))*dti+
     *                a*(u(c,i+oi,j+oj,k+ok)-um(c,i+oi,j+oj,k+ok))*dti-
     *               cs*( (1-b)*(u(c,i+oi,j+oj,k+ok)-u(c,i,j,k))*hi+ 
     *                     (b)* (um(c,i+oi,j+oj,k+ok)-um(c,i,j,k))*hi )
                  vmp = (1-a)*(u(c,i+oi,j+oj,k+ok)-um(c,i+oi,j+oj,k+ok)
     *                                                            )*dti+
     *                a*(u(c,i+2*oi,j+2*oj,k+2*ok)-
     *                                um(c,i+2*oi,j+2*oj,k+2*ok))*dti-
     *               csp*( (1-b)*(u(c,i+2*oi,j+2*oj,k+2*ok)- 
     *                                       u(c,i+oi,j+oj,k+ok))*hi+ 
     *                     (b)* (um(c,i+2*oi,j+2*oj,k+2*ok)-
     *                                        um(c,i+oi,j+oj,k+ok))*hi )
                  vpp = (1-a)*(up(c,i+oi,j+oj,k+ok) -
     *                          u(c,i+oi,j+oj,k+ok)     )*dti +
     *                    a*(up(c,i+2*oi,j+2*oj,k+2*ok) -
     *                        u(c,i+2*oi,j+2*oj,k+2*ok) )*dti-
     *               csp*( (1-b)*( up(c,i+2*oi,j+2*oj,k+2*ok) - 
     *                             up(c,i+oi,j+oj,k+ok)  )*hi+ 
     *                     (b)* ( u(c,i+2*oi,j+2*oj,k+2*ok)-
     *                            u(c,i+oi,j+oj,k+ok)    )*hi )
                  nu = cs*dt*hi
                  cofi = 1/( (1-a) + nu*(1-b) )
                  vp = cofi*( (-a+(1-b)*nu)*vpp +(a+nu*b)*vmp+
     *             ((1-a)-b*nu)*vm + dt*( forc(c) ) )
                  nu = cp*dt*hi
                  cofi = 1/( (1-a) + nu*(1-b) )
                  up(c,i,j,k) = cofi*((-a+(1-b)*nu)*up(c,i+oi,j+oj,k+ok) 
     *                      +(a+nu*b)*u(c,i+oi,j+oj,k+ok) +
     *                     ((1-a)-b*nu)*u(c,i,j,k) + dt*vp )
               enddo
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine BC_CE1_F( ni, nj, nk, up, u, um, rho, mu, la, 
     *                     side, dt, h, wind, forc )

      implicit none

      integer ni, nj, nk, wind(6), oi, oj, ok
      integer i, j, k, n, t1, t2, s, side, ind
      real*8  up(3,ni,nj,nk), u(3,ni,nj,nk), um(3,ni,nj,nk)
      real*8  rho(ni,nj,nk), la(ni,nj,nk), mu(ni,nj,nk), dt, h
      real*8  nu, cp, cs, forc(3,*), den, rhoiav, di

      oi=0
      oj=0
      ok=0
c      forc(1) = 0
c      forc(2) = 0
c      forc(3) = 0
c      wind(1) = 1
c      wind(2) = ni
c      wind(3) = 1
c      wind(4) = nj
c      wind(5) = 1
c      wind(6) = nk

      if( side.eq.0 )then
         oi =  1
c         wind(2) = 1
         n  = 1
         t1 = 2
         t2 = 3
      elseif( side.eq.1 )then
         oi = -1
c         wind(1) = ni
         n  = 1
         t1 = 2
         t2 = 3
      elseif( side.eq.2 )then
         oj =  1
c         wind(4) = 1
         n  = 2
         t1 = 1
         t2 = 3
      elseif( side.eq.3 )then
         oj = -1
c         wind(3) = nj
         n  = 2
         t1 = 1
         t2 = 3
      elseif( side.eq.4 )then
         ok =  1
c         wind(6) = 1
         n  = 3
         t1 = 1
         t2 = 2
      elseif( side.eq.5 )then
         ok = -1
c         wind(5) = nk
         n  = 3
         t1 = 1
         t2 = 2
      endif

      nu = dt/h
      s = 1-2*MOD(side+1,2)
      ind = 1
c
c BUG: wind is base 0 (when m_ghost_point=1), but all fortran arrays are base 1
c
      do k=wind(5),wind(6)
         do j=wind(3),wind(4)
            do i=wind(1),wind(2)
               rhoiav = 1d0/(rho(i+oi,j+oj,k+ok)+rho(i,j,k))
               cp = SQRT(  (2*(mu(i+oi,j+oj,k+ok)+mu(i,j,k))+
     *                         la(i+oi,j+oj,k+ok)+la(i+oi,j+oj,k+ok))
     *                                        *rhoiav )
               cs = SQRT( (mu(i+oi,j+oj,k+ok)+mu(i,j,k))
     *                    *rhoiav )
               den = 1/(nu*cp+1)
               up(n,i,j,k) = u(n,i+oi,j+oj,k+ok)+
     *            (nu*cp-1)/(nu*cp+1)*(up(n,i+oi,j+oj,k+ok)-u(n,i,j,k)) 
     *       +    s*dt*2*cp*den*forc(n,ind)
               di  = 1d0/(nu*cs+1)
               den = (nu*cs-1)*di
               up(t1,i,j,k) = u(t1,i+oi,j+oj,k+ok)+
     *                den*(up(t1,i+oi,j+oj,k+ok)-u(t1,i,j,k))
     *      +          s*dt*di*2*cs*forc(t1,ind)
               up(t2,i,j,k) = u(t2,i+oi,j+oj,k+ok)+
     *                 den*(up(t2,i+oi,j+oj,k+ok)-u(t2,i,j,k))
     *             +   s*dt*di*2*cs*forc(t2,ind)

               ind = ind + 1
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine BC_FREE_SURFACE_ATT_F( ni, nj, nk, u, mu, la, h, na,
     *     bvars, bforcing, curlcoef )

      implicit none
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, a, na
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, h, bforcing(3,ni,nj), bvars(3,na,ni,nj), a1, a2, a3
      real*8 curlcoef, ga, gam, gap

*** Free surface on top
      do j=2,nj-1
         do i=2,ni-1
            a1 = 0
            a2 = 0
            a3 = 0
            do a=1,na
               a1 = a1 + bvars(1,a,i,j)
               a2 = a2 + bvars(2,a,i,j)
               a3 = a3 + bvars(3,a,i,j)
            enddo
            if( mu(i,j,1)+mu(i,j,2).gt.0 )then

               imy = 2d0/(mu(i,j,1)+mu(i,j,2))
 
               u(1,i,j,1) = u(1,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(1,i,j,3)-u(1,i,j,2))+
     *            imy*mu(i,j,2)*(u(3,i+1,j,2)-u(3,i-1,j,2)) 
     *               - 2*h*(bforcing(1,i,j)+a1)*imy
 
               u(2,i,j,1) = u(2,i,j,2) +
     *           half*imy*(mu(i,j,3)+mu(i,j,2))*(u(2,i,j,3)-u(2,i,j,2))+
     *             imy*mu(i,j,2)*(u(3,i,j+1,2)-u(3,i,j-1,2)) 
     *               - 2*h*(bforcing(2,i,j)+a2)*imy
            elseif( curlcoef.ne.0 )then
               ga = 0
               gam= 0
               gap= 0
               if( mu(i,j,2).eq.0 )then
                  ga = curlcoef
               endif
               if( mu(i,j,1).eq.0 )then
                  gam = curlcoef
               endif
               if( mu(i,j,3).eq.0 )then
                  gap = curlcoef
               endif
               imy = 2d0/(gam+ga)
               u(1,i,j,1) = u(1,i,j,2) +
     *                    half*imy*(gap+ga)*(u(1,i,j,3)-u(1,i,j,2))-
     *                    imy*ga*(u(3,i+1,j,2)-u(3,i-1,j,2))
     *                - 2*h*(bforcing(1,i,j)+a1)*imy
               u(2,i,j,1) = u(2,i,j,2) +
     *                    half*imy*(gap+ga)*(u(2,i,j,3)-u(2,i,j,2))-
     *                    imy*ga*(u(3,i,j+1,2)-u(3,i,j-1,2))
     *                - 2*h*(bforcing(2,i,j)+a2)*imy
            endif
            imy = 1d0/(mu(i,j,1)+mu(i,j,2)+half*(la(i,j,1)+la(i,j,2)))
            u(3,i,j,1) = u(3,i,j,2) +
     *            imy*(mu(i,j,3)+mu(i,j,2)+half*(la(i,j,3)+la(i,j,2)))*
     *                        (u(3,i,j,3)-u(3,i,j,2)) +
     *            imy*la(i,j,2)*( u(1,i+1,j,2)-u(1,i-1,j,2)+
     *                            u(2,i,j+1,2)-u(2,i,j-1,2))
     *                - 2*h*(bforcing(3,i,j)+a3)*imy
         enddo
      enddo
 
      end

c-----------------------------------------------------------------------
      subroutine BC_NEUMANN_F( ni, nj, nk, u, side, wind, h, bforcing )
 
      implicit none

      integer ni, nj, nk, side, wind(6), ind, il, jl, kl, sg, i, j, k
      real*8  u(3,ni,nj,nk), bforcing(3,*), h

      il = 0
      jl = 0
      kl = 0
      if( side.eq.0 )then
         il = 1
         sg = -1
      elseif( side.eq.1 )then
         il = -1
         sg = 1
      elseif( side.eq.2 )then
         jl = 1
         sg = -1
      elseif( side.eq.3 )then
         jl = -1
         sg = 1
      elseif( side.eq.4 )then
         kl = 1
         sg = -1
      else
         kl = -1
         sg = 1
      endif
      ind = 1
      do k=wind(5),wind(6)
         do j=wind(3),wind(4)
            do i=wind(1),wind(2)
               u(1,i,j,k) = u(1,i+il,j+jl,k+kl) + 2*sg*h*bforcing(1,ind)
               u(2,i,j,k) = u(2,i+il,j+jl,k+kl) + 2*sg*h*bforcing(2,ind)
               u(3,i,j,k) = u(3,i+il,j+jl,k+kl) + 2*sg*h*bforcing(3,ind)
               ind = ind + 1
            enddo
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine BC_FREE_I( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), il
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,nj,nk), uy, uz, vy, wz, h, curlcoef
      real*8 gam, ga, gap

      il = 1
      i  = ni-1
      if( side.eq.0 )then
         il = -1
         i = 2
      endif
         do k=2,nk-1
            do j=2,nj-1
               if( j.eq.2 .and. onesided(3).eq.1 )then
                  vy = 2*(u(2,i,j+1,k)-u(2,i,j,k))
                  uy = 2*(u(1,i,j+1,k)-u(1,i,j,k))
               elseif( j.eq.nj-1 .and.onesided(4).eq.1 )then
                  vy = 2*(u(2,i,j,k)-u(2,i,j-1,k))
                  uy = 2*(u(1,i,j,k)-u(1,i,j-1,k))
               else
                  vy = u(2,i,j+1,k)-u(2,i,j-1,k)
                  uy = u(1,i,j+1,k)-u(1,i,j-1,k)
               endif
               if( k.eq.2 .and. onesided(5).eq.1 )then
                  uz = 2*(u(1,i,j,k+1)-u(1,i,j,k))
                  wz = 2*(u(3,i,j,k+1)-u(3,i,j,k))
               elseif( k.eq.nk-1 .and. onesided(6).eq.1 )then
                  uz = 2*(u(1,i,j,k)-u(1,i,j,k-1))
                  wz = 2*(u(3,i,j,k)-u(3,i,j,k-1))
               else
                  uz = u(1,i,j,k+1)-u(1,i,j,k-1)
                  wz = u(3,i,j,k+1)-u(3,i,j,k-1)
               endif
               imy = 2d0/( 2*mu(i,j,k)   +la(i,j,k)+
     *                     2*mu(i+il,j,k)+la(i+il,j,k) )

               u(1,i+il,j,k) = u(1,i,j,k) - imy*( (mu(i,j,k)+
     *          mu(i-il,j,k)+ half*(la(i,j,k)+la(i-il,j,k)) )*
     *                  (u(1,i,j,k)-u(1,i-il,j,k))
     *                 +il*(la(i,j,k)*(vy+wz) -2*h*bcforcing(1,j,k) ) )

               if( mu(i,j,k)+mu(i+il,j,k).gt.0 )then
                  imy = 1d0/(half*(mu(i,j,k)+mu(i+il,j,k)))

                  u(2,i+il,j,k) = u(2,i,j,k) - imy*( half*(mu(i,j,k)+
     *               mu(i-il,j,k))*(u(2,i,j,k)-u(2,i-il,j,k)) 
     *                     +il*(mu(i,j,k)*uy -2*h*bcforcing(2,j,k) ))

                  u(3,i+il,j,k) = u(3,i,j,k) - imy*( half*(mu(i,j,k)+
     *               mu(i-il,j,k))*(u(3,i,j,k)-u(3,i-il,j,k)) 
     *                     +il*(mu(i,j,k)*uz -2*h*bcforcing(3,j,k) ))
               elseif( curlcoef.ne.0 )then
                  ga = 0
                  gam= 0
                  gap= 0
                  if( mu(i,j,k).eq.0 )then
                     ga = curlcoef
                  endif
                  if( mu(i+il,j,k).eq.0 )then
                     gam = curlcoef
                  endif
                  if( mu(i-il,j,k).eq.0 )then
                     gap = curlcoef
                  endif
                  imy = 2d0/(gam+ga)
                  u(2,i+il,j,k) = u(2,i,j,k) +
     *                    half*imy*(gap+ga)*(u(2,i-il,j,k)-u(2,i,j,k))+
     *                    il*imy*ga*2*uy
                  u(3,i+il,j,k) = u(3,i,j,k) +
     *                    half*imy*(gap+ga)*(u(3,i-il,j,k)-u(3,i,j,k))+
     *                    il*imy*ga*2*uz
               else
                  u(2,i+il,j,k) = u(2,i,j,k)
                  u(3,i+il,j,k) = u(3,i,j,k)
               endif
            enddo
         enddo
      end

c-----------------------------------------------------------------------
      subroutine BC_FREE_J( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), jl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,nk,ni), vz, ux, vx, wz, h, curlcoef
      real*8 gam, ga, gap

      jl = 1
      j  = nj-1
      if( side.eq.2 )then
         jl = -1
         j  = 2
      endif
         do k=2,nk-1
            do i=2,ni-1
               if( i.eq.2 .and. onesided(1).eq.1 )then
                  vx = 2*(u(2,i+1,j,k)-u(2,i,j,k))
                  ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
               elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
                  vx = 2*(u(2,i,j,k)-u(2,i-1,j,k))
                  ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
               else
                  vx = u(2,i+1,j,k)-u(2,i-1,j,k)
                  ux = u(1,i+1,j,k)-u(1,i-1,j,k)
               endif
               if( k.eq.2 .and. onesided(5).eq.1 )then
                  vz = 2*(u(2,i,j,k+1)-u(2,i,j,k))
                  wz = 2*(u(3,i,j,k+1)-u(3,i,j,k))
               elseif( k.eq.nk-1 .and. onesided(6).eq.1 )then
                  vz = 2*(u(2,i,j,k)-u(2,i,j,k-1))
                  wz = 2*(u(3,i,j,k)-u(3,i,j,k-1))
               else
                  vz = u(2,i,j,k+1)-u(2,i,j,k-1)
                  wz = u(3,i,j,k+1)-u(3,i,j,k-1)
               endif
               imy = 1d0/(half*(2*mu(i,j,k)+la(i,j,k)+
     *                             2*mu(i,j+jl,k)+la(i,j+jl,k)))

               u(2,i,j+jl,k) = u(2,i,j,k) - imy*( (mu(i,j,k)+
     *              mu(i,j-jl,k)+half*(la(i,j,k)+la(i,j-jl,k)))*
     *                             (u(2,i,j,k)-u(2,i,j-jl,k))
     *               +jl*(la(i,j,k)*(ux+wz) -2*h*bcforcing(2,i,k) ))

               if( mu(i,j,k)+mu(i,j+jl,k).gt.0 )then
                  imy = 1d0/(half*(mu(i,j,k)+mu(i,j+jl,k)))

                  u(1,i,j+jl,k) = u(1,i,j,k) - imy*( half*(mu(i,j,k)+
     *                 mu(i,j-jl,k))*(u(1,i,j,k)-u(1,i,j-jl,k)) 
     *                 +jl*(mu(i,j,k)*vx -2*h*bcforcing(1,i,k)))
                  u(3,i,j+jl,k) = u(3,i,j,k) - imy*( half*(mu(i,j,k)+
     *                 mu(i,j-jl,k))*(u(3,i,j,k)-u(3,i,j-jl,k)) 
     *                 +jl*(mu(i,j,k)*vz -2*h*bcforcing(3,i,k)))
               elseif( curlcoef.ne.0 )then
                  ga = 0
                  gam= 0
                  gap= 0
                  if( mu(i,j,k).eq.0 )then
                     ga = curlcoef
                  endif
                  if( mu(i,j+jl,k).eq.0 )then
                     gam = curlcoef
                  endif
                  if( mu(i,j-jl,k).eq.0 )then
                     gap = curlcoef
                  endif
                  imy = 2d0/(gam+ga)
                  u(1,i,j+jl,k) = u(1,i,j,k) +
     *                    half*imy*(gap+ga)*(u(1,i,j-jl,k)-u(1,i,j,k))+
     *                    jl*imy*ga*2*vx
                  u(3,i,j+jl,k) = u(3,i,j,k) +
     *                    half*imy*(gap+ga)*(u(3,i,j-jl,k)-u(3,i,j,k))+
     *                    jl*imy*ga*2*vz
               else
                  u(1,i,j+jl,k) = u(1,i,j,k)
                  u(3,i,j+jl,k) = u(3,i,j,k)
               endif
            enddo
         enddo

      end

c-----------------------------------------------------------------------
      subroutine BC_FREE_K( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), kl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,ni,nj), wx, ux, vy, wy, h, curlcoef
      real*8 gam, ga, gap

      kl = 1
      k = nk-1
      if( side.eq.4 )then
         kl = -1
         k = 2
      endif
      do j=2,nj-1
         do i=2,ni-1
               if( i.eq.2 .and. onesided(1).eq.1 )then
                  wx = 2*(u(3,i+1,j,k)-u(3,i,j,k))
                  ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
               elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
                  wx = 2*(u(3,i,j,k)-u(3,i-1,j,k))
                  ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
               else
                  wx = u(3,i+1,j,k)-u(3,i-1,j,k)
                  ux = u(1,i+1,j,k)-u(1,i-1,j,k)
               endif
               if( j.eq.2 .and. onesided(3).eq.1 )then
                  vy = 2*(u(2,i,j+1,k)-u(2,i,j,k))
                  wy = 2*(u(3,i,j+1,k)-u(3,i,j,k))
               elseif( j.eq.nj-1 .and. onesided(4).eq.1 )then
                  vy = 2*(u(2,i,j,k)-u(2,i,j-1,k))
                  wy = 2*(u(3,i,j,k)-u(3,i,j-1,k))
               else
                  vy = u(2,i,j+1,k)-u(2,i,j-1,k)
                  wy = u(3,i,j+1,k)-u(3,i,j-1,k)
               endif

               imy = 1d0/(half*(2*mu(i,j,k)+la(i,j,k)+
     *                             2*mu(i,j,k+kl)+la(i,j,k+kl)))

               u(3,i,j,k+kl) = u(3,i,j,k) - imy*( (mu(i,j,k)+
     *           mu(i,j,k-kl) + half*(la(i,j,k)+la(i,j,k-kl)))*
     *               (u(3,i,j,k)-u(3,i,j,k-kl))
     *                +kl*(la(i,j,k)*(ux+vy) -2*h*bcforcing(3,i,j) ))

               if( mu(i,j,k)+mu(i,j,k+kl).gt.0 )then
                  imy = 2d0/(mu(i,j,k)+mu(i,j,k+kl))

                  u(1,i,j,k+kl) = u(1,i,j,k) - imy*( half*(mu(i,j,k)+
     *                 mu(i,j,k-kl))*(u(1,i,j,k)-u(1,i,j,k-kl)) 
     *                      +kl*(mu(i,j,k)*wx - 2*h*bcforcing(1,i,j)))

                  u(2,i,j,k+kl) = u(2,i,j,k) - imy*( half*(mu(i,j,k)+
     *                 mu(i,j,k-kl))*(u(2,i,j,k)-u(2,i,j,k-kl)) 
     *                      +kl*(mu(i,j,k)*wy - 2*h*bcforcing(2,i,j)))
               elseif( curlcoef.ne.0 )then
                  ga = 0
                  gam= 0
                  gap= 0
                  if( mu(i,j,k).eq.0 )then
                     ga = curlcoef
                  endif
                  if( mu(i,j,k+kl).eq.0 )then
                     gam = curlcoef
                  endif
                  if( mu(i,j,k-kl).eq.0 )then
                     gap = curlcoef
                  endif
                  imy = 2d0/(gam+ga)
                  u(1,i,j,k+kl) = u(1,i,j,k) +
     *                    half*imy*(gap+ga)*(u(1,i,j,k-kl)-u(1,i,j,k))+
     *                    kl*imy*ga*2*wx
                  u(2,i,j,k+kl) = u(2,i,j,k) +
     *                    half*imy*(gap+ga)*(u(2,i,j,k-kl)-u(2,i,j,k))+
     *                    kl*imy*ga*2*wy
               else
                  u(1,i,j,k+kl) = u(1,i,j,k)
                  u(2,i,j,k+kl) = u(2,i,j,k)
               endif
            enddo
         enddo

      end


c-----------------------------------------------------------------------
      subroutine BC_FREE_I_2DJ( ni, nj, nk, u, mu, la, h, bcforcing,
     *     side, onesided, curlcoef )

      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), il
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,nj,nk), uy, vy, h, curlcoef
      real*8 gam, ga, gap

      il = 1
      i  = ni-1
      k  = 1
      if( side.eq.0 )then
         il = -1
         i  = 2
      endif
      do j=2,nj-1
         if( j.eq.2 .and. onesided(3).eq.1 )then
            vy = 2*(u(2,i,j+1,k)-u(2,i,j,k))
            uy = 2*(u(1,i,j+1,k)-u(1,i,j,k))
         elseif( j.eq.nj-1 .and.onesided(4).eq.1 )then
            vy = 2*(u(2,i,j,k)-u(2,i,j-1,k))
            uy = 2*(u(1,i,j,k)-u(1,i,j-1,k))
         else
            vy = u(2,i,j+1,k)-u(2,i,j-1,k)
            uy = u(1,i,j+1,k)-u(1,i,j-1,k)
         endif
         imy = 2d0/( 2*mu(i,j,k)+la(i,j,k)+
     *        2*mu(i+il,j,k)+la(i+il,j,k) )

         u(1,i+il,j,k) = u(1,i,j,k) - imy*( (mu(i,j,k)+
     *        mu(i-il,j,k)+ half*(la(i,j,k)+la(i-il,j,k)) )*
     *        (u(1,i,j,k)-u(1,i-il,j,k))
     *        +il*(la(i,j,k)*(vy) -2*h*bcforcing(1,j,k) ) )

         if( mu(i,j,k)+mu(i+il,j,k).gt.0 )then
            imy = 1d0/(half*(mu(i,j,k)+mu(i+il,j,k)))

            u(2,i+il,j,k) = u(2,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i-il,j,k))*(u(2,i,j,k)-u(2,i-il,j,k)) 
     *           +il*(mu(i,j,k)*uy -2*h*bcforcing(2,j,k) ))

            u(3,i+il,j,k) = u(3,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i-il,j,k))*(u(3,i,j,k)-u(3,i-il,j,k)) 
     *           +il*( -2*h*bcforcing(3,j,k) ))
         elseif( curlcoef.ne.0 )then
            ga = 0
            gam= 0
            gap= 0
            if( mu(i,j,k).eq.0 )then
               ga = curlcoef
            endif
            if( mu(i+il,j,k).eq.0 )then
               gam = curlcoef
            endif
            if( mu(i-il,j,k).eq.0 )then
               gap = curlcoef
            endif
            imy = 2d0/(gam+ga)
            u(2,i+il,j,k) = u(2,i,j,k) +
     *           half*imy*(gap+ga)*(u(2,i-il,j,k)-u(2,i,j,k))+
     *           il*imy*ga*2*uy
            u(3,i+il,j,k) = u(3,i,j,k) +
     *           half*imy*(gap+ga)*(u(3,i-il,j,k)-u(3,i,j,k))
         else
            u(2,i+il,j,k) = u(2,i,j,k)
            u(3,i+il,j,k) = u(3,i,j,k)
         endif
      enddo

      end


c-----------------------------------------------------------------------
      subroutine BC_FREE_I_2DK( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), il
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,nj,nk), uz, wz, h, curlcoef
      real*8 gam, ga, gap

      il = 1
      i  = ni-1
      j  = 1
      if( side.eq.0 )then
         il = -1
         i = 2
      endif
      do k=2,nk-1
         if( k.eq.2 .and. onesided(5).eq.1 )then
            uz = 2*(u(1,i,j,k+1)-u(1,i,j,k))
            wz = 2*(u(3,i,j,k+1)-u(3,i,j,k))
         elseif( k.eq.nk-1 .and. onesided(6).eq.1 )then
            uz = 2*(u(1,i,j,k)-u(1,i,j,k-1))
            wz = 2*(u(3,i,j,k)-u(3,i,j,k-1))
         else
            uz = u(1,i,j,k+1)-u(1,i,j,k-1)
            wz = u(3,i,j,k+1)-u(3,i,j,k-1)
         endif
         imy = 2d0/( 2*mu(i,j,k)   +la(i,j,k)+
     *        2*mu(i+il,j,k)+la(i+il,j,k) )

         u(1,i+il,j,k) = u(1,i,j,k) - imy*( (mu(i,j,k)+
     *        mu(i-il,j,k)+ half*(la(i,j,k)+la(i-il,j,k)) )*
     *        (u(1,i,j,k)-u(1,i-il,j,k))
     *        +il*(la(i,j,k)*(wz) -2*h*bcforcing(1,j,k) ) )

         if( mu(i,j,k)+mu(i+il,j,k).gt.0 )then
            imy = 1d0/(half*(mu(i,j,k)+mu(i+il,j,k)))

            u(2,i+il,j,k) = u(2,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i-il,j,k))*(u(2,i,j,k)-u(2,i-il,j,k)) 
     *           +il*( -2*h*bcforcing(2,j,k) ))

            u(3,i+il,j,k) = u(3,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i-il,j,k))*(u(3,i,j,k)-u(3,i-il,j,k)) 
     *           +il*(mu(i,j,k)*uz -2*h*bcforcing(3,j,k) ))
         elseif( curlcoef.ne.0 )then
            ga = 0
            gam= 0
            gap= 0
            if( mu(i,j,k).eq.0 )then
               ga = curlcoef
            endif
            if( mu(i+il,j,k).eq.0 )then
               gam = curlcoef
            endif
            if( mu(i-il,j,k).eq.0 )then
               gap = curlcoef
            endif
            imy = 2d0/(gam+ga)
            u(2,i+il,j,k) = u(2,i,j,k) +
     *           half*imy*(gap+ga)*(u(2,i-il,j,k)-u(2,i,j,k))
            u(3,i+il,j,k) = u(3,i,j,k) +
     *           half*imy*(gap+ga)*(u(3,i-il,j,k)-u(3,i,j,k))+
     *           il*imy*ga*2*uz
         else
            u(2,i+il,j,k) = u(2,i,j,k)
            u(3,i+il,j,k) = u(3,i,j,k)
         endif
      enddo

      end


c-----------------------------------------------------------------------
      subroutine BC_FREE_J_2DI( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )

      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), jl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,nk,ni), ux, vx, h, curlcoef
      real*8 gam, ga, gap

      jl = 1
      j  = nj-1
      k  = 1
      if( side.eq.2 )then
         jl = -1
         j  = 2
      endif
      do i=2,ni-1
         if( i.eq.2 .and. onesided(1).eq.1 )then
            vx = 2*(u(2,i+1,j,k)-u(2,i,j,k))
            ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
         elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
            vx = 2*(u(2,i,j,k)-u(2,i-1,j,k))
            ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
         else
            vx = u(2,i+1,j,k)-u(2,i-1,j,k)
            ux = u(1,i+1,j,k)-u(1,i-1,j,k)
         endif
         imy = 1d0/(half*(2*mu(i,j,k)+la(i,j,k)+
     *        2*mu(i,j+jl,k)+la(i,j+jl,k)))

         u(2,i,j+jl,k) = u(2,i,j,k) - imy*( (mu(i,j,k)+
     *        mu(i,j-jl,k)+half*(la(i,j,k)+la(i,j-jl,k)))*
     *        (u(2,i,j,k)-u(2,i,j-jl,k))
     *        +jl*(la(i,j,k)*(ux) -2*h*bcforcing(2,i,k) ))

         if( mu(i,j,k)+mu(i,j+jl,k).gt.0 )then
            imy = 1d0/(half*(mu(i,j,k)+mu(i,j+jl,k)))

            u(1,i,j+jl,k) = u(1,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i,j-jl,k))*(u(1,i,j,k)-u(1,i,j-jl,k)) 
     *           +jl*(mu(i,j,k)*vx -2*h*bcforcing(1,i,k)))
            u(3,i,j+jl,k) = u(3,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i,j-jl,k))*(u(3,i,j,k)-u(3,i,j-jl,k)) 
     *           +jl*( -2*h*bcforcing(3,i,k)))
         elseif( curlcoef.ne.0 )then
            ga = 0
            gam= 0
            gap= 0
            if( mu(i,j,k).eq.0 )then
               ga = curlcoef
            endif
            if( mu(i,j+jl,k).eq.0 )then
               gam = curlcoef
            endif
            if( mu(i,j-jl,k).eq.0 )then
               gap = curlcoef
            endif
            imy = 2d0/(gam+ga)
            u(1,i,j+jl,k) = u(1,i,j,k) +
     *           half*imy*(gap+ga)*(u(1,i,j-jl,k)-u(1,i,j,k))+
     *           jl*imy*ga*2*vx
            u(3,i,j+jl,k) = u(3,i,j,k) +
     *           half*imy*(gap+ga)*(u(3,i,j-jl,k)-u(3,i,j,k))
         else
            u(1,i,j+jl,k) = u(1,i,j,k)
            u(3,i,j+jl,k) = u(3,i,j,k)
         endif
      enddo

      end


c-----------------------------------------------------------------------
      subroutine BC_FREE_K_2DI( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), kl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8 imy, bcforcing(3,ni,nj), wx, ux, vy, wy, h, curlcoef
      real*8 gam, ga, gap

      kl = 1
      k = nk-1
      j = 1
      if( side.eq.4 )then
         kl = -1
         k = 2
      endif
      do i=2,ni-1
         if( i.eq.2 .and. onesided(1).eq.1 )then
            wx = 2*(u(3,i+1,j,k)-u(3,i,j,k))
            ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
         elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
            wx = 2*(u(3,i,j,k)-u(3,i-1,j,k))
            ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
         else
            wx = u(3,i+1,j,k)-u(3,i-1,j,k)
            ux = u(1,i+1,j,k)-u(1,i-1,j,k)
         endif
         imy = 1d0/(half*(2*mu(i,j,k)+la(i,j,k)+
     *                             2*mu(i,j,k+kl)+la(i,j,k+kl)))

         u(3,i,j,k+kl) = u(3,i,j,k) - imy*( (mu(i,j,k)+
     *        mu(i,j,k-kl) + half*(la(i,j,k)+la(i,j,k-kl)))*
     *        (u(3,i,j,k)-u(3,i,j,k-kl))
     *        +kl*(la(i,j,k)*(ux) -2*h*bcforcing(3,i,j) ))

         if( mu(i,j,k)+mu(i,j,k+kl).gt.0 )then
            imy = 2d0/(mu(i,j,k)+mu(i,j,k+kl))

            u(1,i,j,k+kl) = u(1,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i,j,k-kl))*(u(1,i,j,k)-u(1,i,j,k-kl)) 
     *           +kl*(mu(i,j,k)*wx - 2*h*bcforcing(1,i,j)))

            u(2,i,j,k+kl) = u(2,i,j,k) - imy*( half*(mu(i,j,k)+
     *           mu(i,j,k-kl))*(u(2,i,j,k)-u(2,i,j,k-kl)) 
     *           +kl*( - 2*h*bcforcing(2,i,j)))
         elseif( curlcoef.ne.0 )then
            ga = 0
            gam= 0
            gap= 0
            if( mu(i,j,k).eq.0 )then
               ga = curlcoef
            endif
            if( mu(i,j,k+kl).eq.0 )then
               gam = curlcoef
            endif
            if( mu(i,j,k-kl).eq.0 )then
               gap = curlcoef
            endif
            imy = 2d0/(gam+ga)
            u(1,i,j,k+kl) = u(1,i,j,k) +
     *           half*imy*(gap+ga)*(u(1,i,j,k-kl)-u(1,i,j,k))+
     *           kl*imy*ga*2*wx
            u(2,i,j,k+kl) = u(2,i,j,k) +
     *           half*imy*(gap+ga)*(u(2,i,j,k-kl)-u(2,i,j,k))
         else
            u(1,i,j,k+kl) = u(1,i,j,k)
            u(2,i,j,k+kl) = u(2,i,j,k)
         endif
      enddo

      end

c-----------------------------------------------------------------------
      subroutine BCCE1EDGES( ni, nj, nk, u, mu, la, rho, h, dt,
     *    side1, side2, bcforcingi, bcforcingj, bcforcingk )

      implicit none

      real*8 half
      parameter( half = 0.5d0 )

      integer side1, side2, ni, nj, nk 
      integer inormal, onormal, i, j, k, il, jl, kl
      real*8  u(3,ni,nj,nk)
      real*8  mu(ni,nj,nk), la(ni,nj,nk), rho(ni,nj,nk)
      real*8  bcforcingi(3,*), bcforcingj(3,*), bcforcingk(3,*)
      real*8  h, dt, fq, cpcsq, cp, cs, rhi

      if( (side1.eq.0 .or. side1.eq.1) )then
         inormal = 1
         if( side1.eq.0 )then
            i  = 2
            il =-1
         else
            i  = ni-1
            il = 1
         endif
      elseif( side1.eq.2.or.side1.eq.3 )then
         inormal = 2
         if( side1.eq.2 )then
            j  = 2
            jl =-1
         else
            j  = nj-1
            jl = 1
         endif
      else
         inormal = 3
         if( side1.eq.4 )then
            k  = 2
            kl =-1
         else
            k  = nk-1
            kl = 1
         endif
      endif

      if( (side2.eq.0 .or. side2.eq.1) )then
         onormal = 1
         if( side2.eq.0 )then
            i  = 2
            il =-1
         else
            i  = ni-1
            il = 1
         endif
      elseif( side2.eq.2.or.side2.eq.3 )then
         onormal = 2
         if( side2.eq.2 )then
            j  = 2
            jl =-1
         else
            j  = nj-1
            jl = 1
         endif
      else
         onormal = 3
         if( side2.eq.4 )then
            k  = 2
            kl =-1
         else
            k  = nk-1
            kl = 1
         endif
      endif

      if( inormal+onormal.eq.3 )then
*** IJ-edge
         do k=2,nk-1
            rhi = 1/(rho(i,j,k)+rho(i+il,j,k)+rho(i,j+jl,k)+
     *                                        rho(i+il,j+jl,k))
            cs = SQRT((mu(i,j,k)+mu(i+il,j,k)+mu(i,j+jl,k)+
     *                                    mu(i+il,j+jl,k))*rhi)
            cp = SQRT((2*mu(i,j,k)+2*mu(i+il,j,k)+2*mu(i,j+jl,k)+
     *                                        2*mu(i+il,j+jl,k)+
     *                la(i,j,k)+la(i+il,j,k)+la(i,j+jl,k)+
     *                                        la(i+il,j+jl,k))*rhi )
            cpcsq = (cp+cs)/(cp-cs)
            fq = 2*h/(cp-cs)

            u(1,i+il,j+jl,k) = u(1,i,j,k) + 
     *             cpcsq*(u(1,i,j+jl,k)-u(1,i+il,j,k)) + 
     *               fq*(il*cp*bcforcingi(1,k)-jl*cs*bcforcingj(1,k))
            u(2,i+il,j+jl,k) = u(2,i,j,k) - 
     *             cpcsq*(u(2,i,j+jl,k)-u(2,i+il,j,k)) -
     *               fq*(cs*il*bcforcingi(2,k)-cp*jl*bcforcingj(2,k))
         enddo
      elseif( inormal+onormal.eq.4 )then
*** IK-edge
         do j=2,nj-1
            rhi = 1/(rho(i,j,k)+rho(i+il,j,k)+rho(i,j,k+kl)+
     *                                        rho(i+il,j,k+kl))
            cs = SQRT((mu(i,j,k)+mu(i+il,j,k)+mu(i,j,k+kl)+
     *                                    mu(i+il,j,k+kl))*rhi)
            cp = SQRT((2*mu(i,j,k)+2*mu(i+il,j,k)+2*mu(i,j,k+kl)+
     *                                        2*mu(i+il,j,k+kl)+
     *                la(i,j,k)+la(i+il,j,k)+la(i,j,k+kl)+
     *                                        la(i+il,j,k+kl))*rhi )
            cpcsq = (cp+cs)/(cp-cs)
            fq = 2*h/(cp-cs)
            u(1,i+il,j,k+kl) = u(1,i,j,k) + 
     *             cpcsq*(u(1,i,j,k+kl)-u(1,i+il,j,k)) +
     *               fq*(il*cp*bcforcingi(1,j)-cs*kl*bcforcingk(1,j))
            u(3,i+il,j,k+kl) = u(3,i,j,k) - 
     *             cpcsq*(u(3,i,j,k+kl)-u(3,i+il,j,k)) -
     *               fq*(cs*il*bcforcingi(3,j)-cp*kl*bcforcingk(3,j))
         enddo
      else
*** JK-edge
         do i=2,ni-1
            rhi = 1/(rho(i,j,k)+rho(i,j+jl,k)+rho(i,j,k+kl)+
     *                                        rho(i,j+jl,k+kl))
            cs = SQRT((mu(i,j,k)+mu(i,j+jl,k)+mu(i,j,k+kl)+
     *                                    mu(i,j+jl,k+kl))*rhi)
            cp = SQRT((2*mu(i,j,k)+2*mu(i,j+jl,k)+2*mu(i,j,k+kl)+
     *                                        2*mu(i,j+jl,k+kl)+
     *                la(i,j,k)+la(i,j+jl,k)+la(i,j,k+kl)+
     *                                        la(i,j+jl,k+kl))*rhi )
            cpcsq = (cp+cs)/(cp-cs)
            fq = 2*h/(cp-cs)

            u(2,i,j+jl,k+kl) = u(2,i,j,k) + 
     *             cpcsq*(u(2,i,j,k+kl)-u(2,i,j+jl,k)) + 
     *               fq*(jl*cp*bcforcingj(2,i)-cs*kl*bcforcingk(2,i))
            u(3,i,j+jl,k+kl) = u(3,i,j,k) - 
     *             cpcsq*(u(3,i,j,k+kl)-u(3,i,j+jl,k)) -
     *               fq*(cs*jl*bcforcingj(3,i)-cp*kl*bcforcingk(3,i))
         enddo
      endif

      end

c-----------------------------------------------------------------------
      subroutine BCCE1_2DCORNERS( ni, nj, nk, u, mu, la, rho, h, dt,
     *    side1, side2, bcforcingi, bcforcingj, bcforcingk )

      implicit none

      real*8 half
      parameter( half = 0.5d0 )

      integer side1, side2, ni, nj, nk 
      integer inormal, onormal, i, j, k, il, jl, kl
      real*8  u(3,ni,nj,nk)
      real*8  mu(ni,nj,nk), la(ni,nj,nk), rho(ni,nj,nk)
      real*8  bcforcingi(3,*), bcforcingj(3,*), bcforcingk(3,*)
      real*8  h, dt, fq, cpcsq, cp, cs, rhi

      if( (side1.eq.0 .or. side1.eq.1) )then
         inormal = 1
         if( side1.eq.0 )then
            i  = 2
            il =-1
         else
            i  = ni-1
            il = 1
         endif
      elseif( side1.eq.2.or.side1.eq.3 )then
         inormal = 2
         if( side1.eq.2 )then
            j  = 2
            jl =-1
         else
            j  = nj-1
            jl = 1
         endif
      else
         inormal = 3
         if( side1.eq.4 )then
            k  = 2
            kl =-1
         else
            k  = nk-1
            kl = 1
         endif
      endif

      if( (side2.eq.0 .or. side2.eq.1) )then
         onormal = 1
         if( side2.eq.0 )then
            i  = 2
            il =-1
         else
            i  = ni-1
            il = 1
         endif
      elseif( side2.eq.2.or.side2.eq.3 )then
         onormal = 2
         if( side2.eq.2 )then
            j  = 2
            jl =-1
         else
            j  = nj-1
            jl = 1
         endif
      else
         onormal = 3
         if( side2.eq.4 )then
            k  = 2
            kl =-1
         else
            k  = nk-1
            kl = 1
         endif
      endif

      if( inormal+onormal.eq.3 )then
*** IJ-corner
         k = 1
         rhi = 1/(rho(i,j,k)+rho(i+il,j,k)+rho(i,j+jl,k)+
     *                                        rho(i+il,j+jl,k))
         cs = SQRT((mu(i,j,k)+mu(i+il,j,k)+mu(i,j+jl,k)+
     *                                    mu(i+il,j+jl,k))*rhi)
         cp = SQRT((2*mu(i,j,k)+2*mu(i+il,j,k)+2*mu(i,j+jl,k)+
     *                                        2*mu(i+il,j+jl,k)+
     *                la(i,j,k)+la(i+il,j,k)+la(i,j+jl,k)+
     *                                        la(i+il,j+jl,k))*rhi )
         cpcsq = (cp+cs)/(cp-cs)
         fq = 2*h/(cp-cs)

         u(1,i+il,j+jl,k) = u(1,i,j,k) + 
     *             cpcsq*(u(1,i,j+jl,k)-u(1,i+il,j,k)) + 
     *               fq*(il*cp*bcforcingi(1,k)-jl*cs*bcforcingj(1,k))
         u(2,i+il,j+jl,k) = u(2,i,j,k) - 
     *             cpcsq*(u(2,i,j+jl,k)-u(2,i+il,j,k)) -
     *               fq*(cs*il*bcforcingi(2,k)-cp*jl*bcforcingj(2,k))
      elseif( inormal+onormal.eq.4 )then
*** IK-corner
         j = 1
         rhi = 1/(rho(i,j,k)+rho(i+il,j,k)+rho(i,j,k+kl)+
     *                                        rho(i+il,j,k+kl))
         cs = SQRT((mu(i,j,k)+mu(i+il,j,k)+mu(i,j,k+kl)+
     *                                    mu(i+il,j,k+kl))*rhi)
         cp = SQRT((2*mu(i,j,k)+2*mu(i+il,j,k)+2*mu(i,j,k+kl)+
     *                                        2*mu(i+il,j,k+kl)+
     *                la(i,j,k)+la(i+il,j,k)+la(i,j,k+kl)+
     *                                        la(i+il,j,k+kl))*rhi )
         cpcsq = (cp+cs)/(cp-cs)
         fq = 2*h/(cp-cs)
         u(1,i+il,j,k+kl) = u(1,i,j,k) + 
     *        cpcsq*(u(1,i,j,k+kl)-u(1,i+il,j,k)) +
     *        fq*(il*cp*bcforcingi(1,j)-cs*kl*bcforcingk(1,j))
         u(3,i+il,j,k+kl) = u(3,i,j,k) - 
     *        cpcsq*(u(3,i,j,k+kl)-u(3,i+il,j,k)) -
     *        fq*(cs*il*bcforcingi(3,j)-cp*kl*bcforcingk(3,j))
      endif

      end

c-----------------------------------------------------------------------
      subroutine ACC_BC_FREE_I( ni, nj, nk, u, mu, la, h, bcforcing,
     *     side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), il
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8  bcforcing(3,nj,nk), uy, uz, vy, wz, h, curlcoef
      real*8 hi

      il = 1
      i  = ni-1
      if( side.eq.0 )then
         il = -1
         i = 2
      endif
      hi = 1/(2*h)
      do k=2,nk-1
         do j=2,nj-1
            if( j.eq.2 .and. onesided(3).eq.1 )then
               vy = 2*(u(2,i,j+1,k)-u(2,i,j,k))
               uy = 2*(u(1,i,j+1,k)-u(1,i,j,k))
            elseif( j.eq.nj-1 .and.onesided(4).eq.1 )then
               vy = 2*(u(2,i,j,k)-u(2,i,j-1,k))
               uy = 2*(u(1,i,j,k)-u(1,i,j-1,k))
            else
               vy = u(2,i,j+1,k)-u(2,i,j-1,k)
               uy = u(1,i,j+1,k)-u(1,i,j-1,k)
            endif
            if( k.eq.2 .and. onesided(5).eq.1 )then
               uz = 2*(u(1,i,j,k+1)-u(1,i,j,k))
               wz = 2*(u(3,i,j,k+1)-u(3,i,j,k))
            elseif( k.eq.nk-1 .and. onesided(6).eq.1 )then
               uz = 2*(u(1,i,j,k)-u(1,i,j,k-1))
               wz = 2*(u(3,i,j,k)-u(3,i,j,k-1))
            else
               uz = u(1,i,j,k+1)-u(1,i,j,k-1)
               wz = u(3,i,j,k+1)-u(3,i,j,k-1)
            endif

            bcforcing(1,j,k) = bcforcing(1,j,k) + hi*(
     *      il*(mu(i+il,j,k)+mu(i,j,k)+half*(la(i+il,j,k)+la(i,j,k)))*
     *           (u(1,i+il,j,k) - u(1,i,j,k) ) +  
     *      il*(mu(i,j,k)+mu(i-il,j,k)+ half*(la(i,j,k)+la(i-il,j,k)))*
     *                  (u(1,i,j,k)-u(1,i-il,j,k))
     *                 +  la(i,j,k)*(vy+wz))

            bcforcing(2,j,k) = bcforcing(2,j,k) + hi*(
     *      il*half*(mu(i,j,k)+mu(i+il,j,k))*(u(2,i+il,j,k)-u(2,i,j,k))+
     *      il*half*(mu(i,j,k)+mu(i-il,j,k))*(u(2,i,j,k)-u(2,i-il,j,k)) 
     *                     +mu(i,j,k)*uy)

            bcforcing(3,j,k) = bcforcing(3,j,k) + hi*(
     *      il*half*(mu(i,j,k)+mu(i+il,j,k))*(u(3,i+il,j,k)-u(3,i,j,k))+
     *      il*half*(mu(i,j,k)+mu(i-il,j,k))*(u(3,i,j,k)-u(3,i-il,j,k))
     *                     +mu(i,j,k)*uz)
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine ACC_BC_FREE_J( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), jl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8  bcforcing(3,nk,ni), vz, ux, vx, wz, h, curlcoef
      real*8 hi

      jl = 1
      j  = nj-1
      if( side.eq.2 )then
         jl = -1
         j  = 2
      endif
      hi = 1/(2*h)
      do k=2,nk-1
         do i=2,ni-1
            if( i.eq.2 .and. onesided(1).eq.1 )then
               vx = 2*(u(2,i+1,j,k)-u(2,i,j,k))
               ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
            elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
               vx = 2*(u(2,i,j,k)-u(2,i-1,j,k))
               ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
            else
               vx = u(2,i+1,j,k)-u(2,i-1,j,k)
               ux = u(1,i+1,j,k)-u(1,i-1,j,k)
            endif
            if( k.eq.2 .and. onesided(5).eq.1 )then
               vz = 2*(u(2,i,j,k+1)-u(2,i,j,k))
               wz = 2*(u(3,i,j,k+1)-u(3,i,j,k))
            elseif( k.eq.nk-1 .and. onesided(6).eq.1 )then
               vz = 2*(u(2,i,j,k)-u(2,i,j,k-1))
               wz = 2*(u(3,i,j,k)-u(3,i,j,k-1))
            else
               vz = u(2,i,j,k+1)-u(2,i,j,k-1)
               wz = u(3,i,j,k+1)-u(3,i,j,k-1)
            endif

            bcforcing(2,i,k) = bcforcing(2,i,k) + hi*(
     *    jl*(mu(i,j,k)+mu(i,j+jl,k)+half*(la(i,j,k)+la(i,j+jl,k)))*
     *         (u(2,i,j+jl,k) - u(2,i,j,k)) + 
     *    jl*(mu(i,j,k)+mu(i,j-jl,k)+half*(la(i,j,k)+la(i,j-jl,k)))*
     *         (u(2,i,j,k)-u(2,i,j-jl,k))
     *           +la(i,j,k)*(ux+wz) )

            bcforcing(1,i,k) = bcforcing(1,i,k) + hi*(
     *      jl*half*(mu(i,j,k)+mu(i,j+jl,k))*(u(1,i,j+jl,k)-u(1,i,j,k))+
     *      jl*half*(mu(i,j,k)+mu(i,j-jl,k))*(u(1,i,j,k)-u(1,i,j-jl,k))
     *              + mu(i,j,k)*vx )

            bcforcing(3,i,k) = bcforcing(3,i,k) + hi*( 
     *      jl*half*(mu(i,j,k)+mu(i,j+jl,k))*(u(1,3,j+jl,k)-u(3,i,j,k))+
     *      jl*half*(mu(i,j,k)+mu(i,j-jl,k))*(u(1,3,j,k)-u(3,i,j-jl,k))
     *              +mu(i,j,k)*vz )
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine ACC_BC_FREE_K( ni, nj, nk, u, mu, la, h, bcforcing,
     *                      side, onesided, curlcoef )
      implicit none
 
      real*8 half
      parameter( half = 1d0/2 )
 
      integer i, j, k, ni, nj, nk, side, onesided(6), kl
      real*8  mu(ni,nj,nk), la(ni,nj,nk), u(3,ni,nj,nk)
      real*8  bcforcing(3,ni,nj), wx, ux, vy, wy, h, curlcoef
      real*8 hi

      kl = 1
      k = nk-1
      if( side.eq.4 )then
         kl = -1
         k = 2
      endif
      hi = 1/(2*h)
      do j=2,nj-1
         do i=2,ni-1
            if( i.eq.2 .and. onesided(1).eq.1 )then
               wx = 2*(u(3,i+1,j,k)-u(3,i,j,k))
               ux = 2*(u(1,i+1,j,k)-u(1,i,j,k))
            elseif( i.eq.ni-1 .and.onesided(2).eq.1 )then
               wx = 2*(u(3,i,j,k)-u(3,i-1,j,k))
               ux = 2*(u(1,i,j,k)-u(1,i-1,j,k))
            else
               wx = u(3,i+1,j,k)-u(3,i-1,j,k)
               ux = u(1,i+1,j,k)-u(1,i-1,j,k)
            endif
            if( j.eq.2 .and. onesided(3).eq.1 )then
               vy = 2*(u(2,i,j+1,k)-u(2,i,j,k))
               wy = 2*(u(3,i,j+1,k)-u(3,i,j,k))
            elseif( j.eq.nj-1 .and. onesided(4).eq.1 )then
               vy = 2*(u(2,i,j,k)-u(2,i,j-1,k))
               wy = 2*(u(3,i,j,k)-u(3,i,j-1,k))
            else
               vy = u(2,i,j+1,k)-u(2,i,j-1,k)
               wy = u(3,i,j+1,k)-u(3,i,j-1,k)
            endif

            bcforcing(3,i,j) = bcforcing(3,i,j) + hi*(
     *      kl*(mu(i,j,k+kl)+mu(i,j,k)+half*(la(i,j,k+kl)+la(i,j,k)))*
     *           (u(3,i,j,k+kl)-u(3,i,j,k)) + 
     *      kl*(mu(i,j,k)+mu(i,j,k-kl)+half*(la(i,j,k)+la(i,j,k-kl)))*
     *           (u(3,i,j,k)-u(3,i,j,k-kl))
     *           +la(i,j,k)*(ux+vy) )
            bcforcing(1,i,j) = bcforcing(1,i,j) + hi*(
     *      kl*half*(mu(i,j,k)+mu(i,j,k+kl))*(u(1,i,j,k+kl)-u(1,i,j,k))+
     *      kl*half*(mu(i,j,k)+mu(i,j,k-kl))*(u(1,i,j,k)-u(1,i,j,k-kl))
     *              +mu(i,j,k)*wx )
            bcforcing(2,i,j) = bcforcing(2,i,j) + hi*(
     *      kl*half*(mu(i,j,k)+mu(i,j,k+kl))*(u(2,i,j,k+kl)-u(2,i,j,k))+
     *      kl*half*(mu(i,j,k)+mu(i,j,k-kl))*(u(2,i,j,k)-u(2,i,j,k-kl)) 
     *              + mu(i,j,k)*wy )
         enddo
      enddo

      end

