      subroutine INVTESTMTRL( ib, ie, jb, je, kb, ke, rho, cs, cp, h,
     *                        zmin, nr )
      implicit none
      integer ib, ie, jb, je, kb, ke, nr
      integer i, j, k
      real*8 rho(ib:ie,jb:je,kb:ke), cs(ib:ie,jb:je,kb:ke), zmin, h
      real*8 cp(ib:ie,jb:je,kb:ke), x, y, z, ep, om
      real*8 phip, phis, phir, omx, omy, omz, pi2
      real*8 xminbox, xmaxbox, yminbox, ymaxbox, zminbox, zmaxbox

      ep = 0.01d0
*** 2*pi
      pi2 = ATAN(1d0)*8
      om = pi2
      xminbox = 0.4d0
      xmaxbox = 0.6d0
      yminbox = 0.4d0
      ymaxbox = 0.6d0
      zminbox = 0.1d0
      zmaxbox = 0.3d0
      if( nr.eq.3 )then
         omx = pi2*5d0/30000
         omy = pi2*5d0/30000
         omz = pi2*3d0/17000
         phip = 0.3d0
         phir = 0.17d0
         phis = 0.08d0
      endif
      do k=kb,ke
         z = zmin + (k-1)*h
         do j=jb,je
            y = (j-1)*h
            do i=ib,ie
               x = (i-1)*h
               if( nr.eq.1 )then
                  rho(i,j,k) = 1+ep*sin(om*x+0.13d0)*sin(om*y)*sin(om*z)
                  cs(i,j,k)  = 2+ep*cos(om*x)*sin(om*y)*cos(om*z+0.01d0)
                  cp(i,j,k)  = 4+ep*sin(om*x+0.4d0)*sin(om*y)*
     *                                              cos(om*z+0.1d0)
               elseif( nr.eq.2 )then
                  rho(i,j,k) = 1
                  cs(i,j,k)  = 2
                  cp(i,j,k)  = 4
                  if( x.gt.xminbox .and. x.le.xmaxbox .and. y.ge.yminbox
     *          .and. y.le.ymaxbox .and. z.ge.zminbox .and. z.le.zmaxbox
     *                                                            )then
                     rho(i,j,k) = 1.3d0
                     cs(i,j,k)  = 0.5d0
                     cp(i,j,k)  = 2
                  endif
               elseif( nr.eq.3 )then
                  rho(i,j,k)=2650+
     *                  50*sin(omx*x)*sin(omy*y)*cos(omz*z+phir)
                  cp(i,j,k)=5000+
     *                1000*sin(omx*x+phip)*sin(omy*y)*cos(omz*z)
                  cs(i,j,k) = cp(i,j,k)/( 1.8660d0 + 
     *                  0.1339*cos(omx*x)*sin(omy*y+phis)*cos(omz*z) )
c                  cs(i,j,k) =2732+cos(omx*x)*sin(omy*y+phis)*cos(omz*z)
c     *                 732*
               endif
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine INVTESTMTRLC( ib, ie, jb, je, kb, ke, rho, cs, cp, xx,
     *                         yy, zz, nr )
      implicit none
      integer ib, ie, jb, je, kb, ke, nr
      integer i, j, k
      real*8 rho(ib:ie,jb:je,kb:ke), cs(ib:ie,jb:je,kb:ke)
      real*8 cp(ib:ie,jb:je,kb:ke), xx(ib:ie,jb:je,kb:ke)
      real*8 yy(ib:ie,jb:je,kb:ke), zz(ib:ie,jb:je,kb:ke), x, y, z, om
      real*8 xminbox, xmaxbox, yminbox, ymaxbox, zminbox, zmaxbox, ep
      real*8 phip, phis, phir, omx, omy, omz, pi2

      ep = 0.01d0
*** 2*pi
*** 2*pi
      pi2 = ATAN(1d0)*8
      om = pi2
      xminbox = 0.4d0
      xmaxbox = 0.6d0
      yminbox = 0.4d0
      ymaxbox = 0.6d0
      zminbox = 0.1d0
      zmaxbox = 0.3d0
      if( nr.eq.3 )then
         omx = pi2*5d0/30000
         omy = pi2*5d0/30000
         omz = pi2*3d0/17000
         phip = 0.3d0
         phir = 0.17d0
         phis = 0.08d0
      endif
      do k=kb,ke
         do j=jb,je
            do i=ib,ie
               x = xx(i,j,k)
               y = yy(i,j,k)
               z = zz(i,j,k)
               if( nr.eq.1 )then
                  rho(i,j,k) = 1+ep*sin(om*x+0.13d0)*sin(om*y)*sin(om*z)
                  cs(i,j,k)  = 2+ep*cos(om*x)*sin(om*y)*cos(om*z+0.01d0)
                  cp(i,j,k)  = 4+ep*sin(om*x+0.4d0)*sin(om*y)*
     *                                              cos(om*z+0.1d0)
               elseif( nr.eq.2 )then
                  rho(i,j,k) = 1
                  cs(i,j,k)  = 2
                  cp(i,j,k)  = 4
                  if( x.gt.xminbox .and. x.le.xmaxbox .and. y.ge.yminbox
     *          .and. y.le.ymaxbox .and. z.ge.zminbox .and. z.le.zmaxbox
     *                                                            )then
                     rho(i,j,k) = 1.3d0
                     cs(i,j,k)  = 0.5d0
                     cp(i,j,k)  = 2
                  endif
               elseif( nr.eq.3 )then
                  rho(i,j,k)=2650+
     *                  50*sin(omx*x)*sin(omy*y)*cos(omz*z+phir)
                  cs(i,j,k) =2732+
     *                 732*cos(omx*x)*sin(omy*y+phis)*cos(omz*z)
                  cp(i,j,k)=5000+
     *                1000*sin(omx*x+phip)*sin(omy*y)*cos(omz*z)
               endif
            enddo
         enddo
      enddo
      end
