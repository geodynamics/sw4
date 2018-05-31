      subroutine INTERPOLATEMTRL( nx, ny, nz, xmin, ymin, zmin, hx, hy,
     *     hz, rho, mu, lambda, ib, ie, jb, je, kb, ke,
     *     ibact, ieact, jbact, jeact, kbact, keact, rhogrid, mugrid,
     *     lambdagrid, hf, zmingrid )
*** Trilinear interpolation
      implicit none
      integer nx, ny, nz, ib, ie, jb, je, kb, ke, ibact, ieact
      integer jbact, jeact, kbact, keact, i, j, k, ic, jc, kc
      real*8  xmin, ymin, zmin, hx, hy, hz, hf, zmingrid
      real*8  wgx, wgy, wgz, x, y, z
      real*8  rho(0:nx+1,0:ny+1,0:nz+1), mu(0:nx+1,0:ny+1,0:nz+1)
      real*8  lambda(0:nx+1,0:ny+1,0:nz+1)
      real*8  rhogrid(ib:ie,jb:je,kb:ke), mugrid(ib:ie,jb:je,kb:ke)
      real*8  lambdagrid(ib:ie,jb:je,kb:ke)

      do k=kbact,keact
         z = (k-1)*hf + zmingrid
         kc = (z-zmin)/hz
         if( kc.lt.0 )then
            kc  = 0
            wgz = 0
         elseif( kc.gt.nz )then
            kc  = nz
            wgz = 1
         else
            wgz = (z-zmin)/hz-kc
         endif
         do j=jbact, jeact
            y = (j-1)*hf
            jc = (y-ymin)/hy
            if( jc.lt.0 )then
               jc  = 0
               wgy = 0
            elseif( jc.gt.ny )then
               jc  = ny
               wgy = 1
            else
               wgy = (y-ymin)/hy-jc
            endif
            do i=ibact, ieact
               x = (i-1)*hf
               ic = (x-xmin)/hx
               if( ic.lt.0 )then
                  ic  = 0
                  wgx = 0
               elseif( ic.gt.nx )then
                  ic  = nx
                  wgx = 1
               else
                  wgx = (x-xmin)/hx-ic
               endif
               rhogrid(i,j,k) = rhogrid(i,j,k) +
     *      (1-wgz)*(
     *        (1-wgy)*( (1-wgx)*rho(ic,jc,kc) + wgx*rho(ic+1,jc,kc) ) +
     *         wgy*( (1-wgx)*rho(ic,jc+1,kc) + wgx*rho(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     *    (1-wgy)*( (1-wgx)*rho(ic,jc,kc+1) + wgx*rho(ic+1,jc,kc+1) ) +
     *     wgy*( (1-wgx)*rho(ic,jc+1,kc+1) + wgx*rho(ic+1,jc+1,kc+1) ) )
               mugrid(i,j,k) = mugrid(i,j,k) +
     *      (1-wgz)*(
     *        (1-wgy)*( (1-wgx)*mu(ic,jc,kc) + wgx*mu(ic+1,jc,kc) ) +
     *         wgy*( (1-wgx)*mu(ic,jc+1,kc) + wgx*mu(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     *    (1-wgy)*( (1-wgx)*mu(ic,jc,kc+1) + wgx*mu(ic+1,jc,kc+1) ) +
     *     wgy*( (1-wgx)*mu(ic,jc+1,kc+1) + wgx*mu(ic+1,jc+1,kc+1) ) )
               lambdagrid(i,j,k) = lambdagrid(i,j,k) + 
     *      (1-wgz)*(
     *   (1-wgy)*( (1-wgx)*lambda(ic,jc,kc) + wgx*lambda(ic+1,jc,kc) ) +
     *   wgy*( (1-wgx)*lambda(ic,jc+1,kc) + wgx*lambda(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     * (1-wgy)*( (1-wgx)*lambda(ic,jc,kc+1) + wgx*lambda(ic+1,jc,kc+1))+
     *  wgy*((1-wgx)*lambda(ic,jc+1,kc+1) + wgx*lambda(ic+1,jc+1,kc+1)))
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine INTERPOLATEMTRLC( nx, ny, nz, xmin, ymin, zmin, hx, hy,
     *     hz, rho, mu, lambda, ib, ie, jb, je, kb, ke,
     *     ibact, ieact, jbact, jeact, kbact, keact, rhogrid, mugrid,
     *     lambdagrid, hf, zgrid )
*** Trilinear interpolation
      implicit none
      integer nx, ny, nz, ib, ie, jb, je, kb, ke, ibact, ieact
      integer jbact, jeact, kbact, keact, i, j, k, ic, jc, kc
      real*8  xmin, ymin, zmin, hx, hy, hz, hf, zmingrid
      real*8  wgx, wgy, wgz, x, y, z
      real*8  rho(0:nx+1,0:ny+1,0:nz+1), mu(0:nx+1,0:ny+1,0:nz+1)
      real*8  lambda(0:nx+1,0:ny+1,0:nz+1)
      real*8  rhogrid(ib:ie,jb:je,kb:ke), mugrid(ib:ie,jb:je,kb:ke)
      real*8  lambdagrid(ib:ie,jb:je,kb:ke)
      real*8  zgrid(ib:ie,jb:je,kb:ke)

      do k=kbact,keact
         do j=jbact, jeact
            y = (j-1)*hf
            jc = (y-ymin)/hy
            if( jc.lt.0 )then
               jc  = 0
               wgy = 0
            elseif( jc.gt.ny )then
               jc  = ny
               wgy = 1
            else
               wgy = (y-ymin)/hy-jc
            endif
            do i=ibact, ieact
               kc = (zgrid(i,j,k)-zmin)/hz
               if( kc.lt.0 )then
                  kc  = 0
                  wgz = 0
               elseif( kc.gt.nz )then
                  kc  = nz
                  wgz = 1
               else
                  wgz = (z-zmin)/hz-kc
               endif
               x = (i-1)*hf
               ic = (x-xmin)/hx
               if( ic.lt.0 )then
                  ic  = 0
                  wgx = 0
               elseif( ic.gt.nx )then
                  ic  = nx
                  wgx = 1
               else
                  wgx = (x-xmin)/hx-ic
               endif
               rhogrid(i,j,k) = rhogrid(i,j,k) +
     *      (1-wgz)*(
     *        (1-wgy)*( (1-wgx)*rho(ic,jc,kc) + wgx*rho(ic+1,jc,kc) )+
     *         wgy*( (1-wgx)*rho(ic,jc+1,kc) + wgx*rho(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     *    (1-wgy)*( (1-wgx)*rho(ic,jc,kc+1) + wgx*rho(ic+1,jc,kc+1) )+
     *     wgy*( (1-wgx)*rho(ic,jc+1,kc+1) + wgx*rho(ic+1,jc+1,kc+1) ) )
               mugrid(i,j,k) = mugrid(i,j,k) + 
     *      (1-wgz)*(
     *        (1-wgy)*( (1-wgx)*mu(ic,jc,kc) + wgx*mu(ic+1,jc,kc) )+
     *         wgy*( (1-wgx)*mu(ic,jc+1,kc) + wgx*mu(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     *    (1-wgy)*( (1-wgx)*mu(ic,jc,kc+1) + wgx*mu(ic+1,jc,kc+1) )+
     *     wgy*( (1-wgx)*mu(ic,jc+1,kc+1) + wgx*mu(ic+1,jc+1,kc+1) ) )
               lambdagrid(i,j,k) = lambdagrid(i,j,k) + 
     *      (1-wgz)*(
     *   (1-wgy)*( (1-wgx)*lambda(ic,jc,kc) + wgx*lambda(ic+1,jc,kc) )+
     *   wgy*( (1-wgx)*lambda(ic,jc+1,kc) + wgx*lambda(ic+1,jc+1,kc) ) ) 
     *      + wgz*(
     * (1-wgy)*( (1-wgx)*lambda(ic,jc,kc+1) + wgx*lambda(ic+1,jc,kc+1))+
     *  wgy*((1-wgx)*lambda(ic,jc+1,kc+1) + wgx*lambda(ic+1,jc+1,kc+1)))
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine GRADIENTS( nx, ny, nz, xmin, ymin, zmin, hx, hy,
     *     hz, gradrho, gradmu, gradlambda, ib, ie, jb, je, kb, ke,
     *     ibact, ieact, jbact, jeact, kbact, keact, gradrhogrid, 
     *     gradmugrid, gradlambdagrid, hf, zmingrid )

*** Chain rule of trilinear interpolation
      implicit none
      integer nx, ny, nz, ib, ie, jb, je, kb, ke, ibact, ieact
      integer jbact, jeact, kbact, keact, i, j, k, ic, jc, kc
      real*8  xmin, ymin, zmin, hx, hy, hz, hf, zmingrid
      real*8  wgx, wgy, wgz, x, y, z
      real*8  gradrho(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradmu(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradlambda(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradrhogrid(ib:ie,jb:je,kb:ke)
      real*8  gradmugrid(ib:ie,jb:je,kb:ke)
      real*8  gradlambdagrid(ib:ie,jb:je,kb:ke)

      do k=kbact,keact
         z = (k-1)*hf + zmingrid
         kc = (z-zmin)/hz
         if( kc.lt.0 )then
            kc  = 0
            wgz = 0
         elseif( kc.gt.nz )then
            kc  = nz
            wgz = 1
         else
            wgz = (z-zmin)/hz-kc
         endif
         do j=jbact, jeact
            y = (j-1)*hf
            jc = (y-ymin)/hy
            if( jc.lt.0 )then
               jc  = 0
               wgy = 0
            elseif( jc.gt.ny )then
               jc  = ny
               wgy = 1
            else
               wgy = (y-ymin)/hy-jc
            endif
            do i=ibact, ieact
               x = (i-1)*hf
               ic = (x-xmin)/hx
               if( ic.lt.0 )then
                  ic  = 0
                  wgx = 0
               elseif( ic.gt.nx )then
                  ic  = nx
                  wgx = 1
               else
                  wgx = (x-xmin)/hx-ic
               endif
               gradrho(ic,jc,kc) = gradrho(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc,kc) = gradrho(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc+1,kc) = gradrho(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc+1,kc) = gradrho(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc,kc+1) = gradrho(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc,kc+1) = gradrho(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc+1,kc+1) = gradrho(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc+1,kc+1) = gradrho(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradrhogrid(i,j,k)
               gradmu(ic,jc,kc) = gradmu(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc,kc) = gradmu(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc+1,kc) = gradmu(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc+1,kc) = gradmu(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc,kc+1) = gradmu(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc,kc+1) = gradmu(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc+1,kc+1) = gradmu(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc+1,kc+1) = gradmu(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradmugrid(i,j,k)
               gradlambda(ic,jc,kc) = gradlambda(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc,kc) = gradlambda(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc+1,kc) = gradlambda(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc+1,kc) = gradlambda(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc,kc+1) = gradlambda(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc,kc+1) = gradlambda(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc+1,kc+1) = gradlambda(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc+1,kc+1) = gradlambda(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradlambdagrid(i,j,k)
            enddo
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine GRADIENTSC( nx, ny, nz, xmin, ymin, zmin, hx, hy,
     *     hz, gradrho, gradmu, gradlambda, ib, ie, jb, je, kb, ke,
     *     ibact, ieact, jbact, jeact, kbact, keact, gradrhogrid, 
     *     gradmugrid, gradlambdagrid, hf, zgrid )

*** Chain rule of trilinear interpolation
      implicit none
      integer nx, ny, nz, ib, ie, jb, je, kb, ke, ibact, ieact
      integer jbact, jeact, kbact, keact, i, j, k, ic, jc, kc
      real*8  xmin, ymin, zmin, hx, hy, hz, hf, zmingrid
      real*8  wgx, wgy, wgz, x, y, z
      real*8  gradrho(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradmu(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradlambda(0:nx+1,0:ny+1,0:nz+1)
      real*8  gradrhogrid(ib:ie,jb:je,kb:ke)
      real*8  gradmugrid(ib:ie,jb:je,kb:ke)
      real*8  gradlambdagrid(ib:ie,jb:je,kb:ke)
      real*8  zgrid(ib:ie,jb:je,kb:ke)

      do k=kbact,keact
         do j=jbact, jeact
            y = (j-1)*hf
            jc = (y-ymin)/hy
            if( jc.lt.0 )then
               jc  = 0
               wgy = 0
            elseif( jc.gt.ny )then
               jc  = ny
               wgy = 1
            else
               wgy = (y-ymin)/hy-jc
            endif
            do i=ibact, ieact
               kc = (zgrid(i,j,k)-zmin)/hz
               if( kc.lt.0 )then
                  kc  = 0
                  wgz = 0
               elseif( kc.gt.nz )then
                  kc  = nz
                  wgz = 1
               else
                  wgz = (z-zmin)/hz-kc
               endif
               x = (i-1)*hf
               ic = (x-xmin)/hx
               if( ic.lt.0 )then
                  ic  = 0
                  wgx = 0
               elseif( ic.gt.nx )then
                  ic  = nx
                  wgx = 1
               else
                  wgx = (x-xmin)/hx-ic
               endif
               gradrho(ic,jc,kc) = gradrho(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc,kc) = gradrho(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc+1,kc) = gradrho(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc+1,kc) = gradrho(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc,kc+1) = gradrho(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc,kc+1) = gradrho(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradrhogrid(i,j,k)
               gradrho(ic,jc+1,kc+1) = gradrho(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradrhogrid(i,j,k)
               gradrho(ic+1,jc+1,kc+1) = gradrho(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradrhogrid(i,j,k)
               gradmu(ic,jc,kc) = gradmu(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc,kc) = gradmu(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc+1,kc) = gradmu(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc+1,kc) = gradmu(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc,kc+1) = gradmu(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc,kc+1) = gradmu(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradmugrid(i,j,k)
               gradmu(ic,jc+1,kc+1) = gradmu(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradmugrid(i,j,k)
               gradmu(ic+1,jc+1,kc+1) = gradmu(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradmugrid(i,j,k)
               gradlambda(ic,jc,kc) = gradlambda(ic,jc,kc) +
     *                 (1-wgz)*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc,kc) = gradlambda(ic+1,jc,kc) +
     *                 (1-wgz)*(1-wgy)*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc+1,kc) = gradlambda(ic,jc+1,kc) +
     *                 (1-wgz)*wgy*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc+1,kc) = gradlambda(ic+1,jc+1,kc) +
     *                 (1-wgz)*wgy*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc,kc+1) = gradlambda(ic,jc,kc+1) +
     *                 wgz*(1-wgy)*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc,kc+1) = gradlambda(ic+1,jc,kc+1) +
     *                 wgz*(1-wgy)*wgx*gradlambdagrid(i,j,k)
               gradlambda(ic,jc+1,kc+1) = gradlambda(ic,jc+1,kc+1) +
     *                 wgz*wgy*(1-wgx)*gradlambdagrid(i,j,k)
               gradlambda(ic+1,jc+1,kc+1) = gradlambda(ic+1,jc+1,kc+1) +
     *                 wgz*wgy*wgx*gradlambdagrid(i,j,k)
            enddo
         enddo
      enddo
      end
