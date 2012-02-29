      subroutine VELSUM( is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2,
     *                   mu, lambda, rho, cp, cs, npts )
      implicit none
      integer is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2
      integer i, j, k
      real*8 mu(is:ie,js:je,ks:ke), lambda(is:ie,js:je,ks:ke)
      real*8 rho(is:ie,js:je,ks:ke), cp, cs, npts
      cp = 0
      cs = 0
      npts = DBLE(i2-i1+1)*DBLE(j2-j1+1)*DBLE(k2-k1+1)
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               cp = cp + SQRT( (2*mu(i,j,k)+lambda(i,j,k))/rho(i,j,k) )
               cs = cs + SQRT( mu(i,j,k)/rho(i,j,k) )
            enddo
         enddo
      enddo
      end
