      subroutine ENERGY4( is, ie, js, je, ks, ke, i1, i2, j1, j2, k1,
     *                    k2, onesided, um, u, up, rho, h, energy )
      implicit none
      integer is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2
      integer i, j, k, onesided(6)
      real*8 um(3,is:ie,js:je,ks:ke), u(3,is:ie,js:je,ks:ke)
      real*8 up(3,is:ie,js:je,ks:ke), rho(is:ie,js:je,ks:ke), energy
      real*8 normwgh(4), term, normfact, h
      normwgh(1) = 17d0/48
      normwgh(2) = 59d0/48
      normwgh(3) = 43d0/48
      normwgh(4) = 49d0/48
      energy = 0
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               term =(
     *            (up(1,i,j,k)-u(1,i,j,k))*(up(1,i,j,k)-u(1,i,j,k)) + 
     *            (up(2,i,j,k)-u(2,i,j,k))*(up(2,i,j,k)-u(2,i,j,k)) +
     *            (up(3,i,j,k)-u(3,i,j,k))*(up(3,i,j,k)-u(3,i,j,k)) -
     *              up(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k)) -
     *              up(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k)) -
     *              up(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))
     *              )*rho(i,j,k)
               normfact = 1               
               if( k.le.4 .and. onesided(5).eq.1 )then
                  normfact = normwgh(k)
               endif
               if( k.ge.k2-3 .and. onesided(6).eq.1 )then
                  normfact = normwgh(k2-k+1)
               endif
               energy = energy + normfact*h*h*h*term
            enddo
         enddo
      enddo
      end
c-----------------------------------------------------------------------
      subroutine ENERGY4C( is, ie, js, je, ks, ke, i1, i2, j1, j2, k1,
     *                     k2, onesided, um, u, up, rho, jac, energy )
      implicit none
      integer is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2
      integer i, j, k, onesided(6)
      real*8 um(3,is:ie,js:je,ks:ke), u(3,is:ie,js:je,ks:ke)
      real*8 up(3,is:ie,js:je,ks:ke), rho(is:ie,js:je,ks:ke), energy
      real*8 jac(is:ie,js:je,ks:ke)
      real*8 normwgh(4), term, normfact
      normwgh(1) = 17d0/48
      normwgh(2) = 59d0/48
      normwgh(3) = 43d0/48
      normwgh(4) = 49d0/48
      energy = 0
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               term =(
     *            (up(1,i,j,k)-u(1,i,j,k))*(up(1,i,j,k)-u(1,i,j,k)) + 
     *            (up(2,i,j,k)-u(2,i,j,k))*(up(2,i,j,k)-u(2,i,j,k)) +
     *            (up(3,i,j,k)-u(3,i,j,k))*(up(3,i,j,k)-u(3,i,j,k)) -
     *              up(1,i,j,k)*(up(1,i,j,k)-2*u(1,i,j,k)+um(1,i,j,k)) -
     *              up(2,i,j,k)*(up(2,i,j,k)-2*u(2,i,j,k)+um(2,i,j,k)) -
     *              up(3,i,j,k)*(up(3,i,j,k)-2*u(3,i,j,k)+um(3,i,j,k))
     *              )*rho(i,j,k)*jac(i,j,k)
               normfact = 1               
               if( k.le.4 .and. onesided(5).eq.1 )then
                  normfact = normwgh(k)
               endif
               if( k.ge.k2-3 .and. onesided(6).eq.1 )then
                  normfact = normwgh(k2-k+1)
               endif
               energy = energy + normfact*term
            enddo
         enddo
      enddo
      end
