      subroutine LINSOLVELU( n, a, b, x )

      integer i, j, k
      real*8  a(n,n), b(n), x(n)

      do k=1,n-1
*** Find pivot element
         piv = a(k,k)
         ipiv = k
         do i=k+1,n
            if( ABS(a(i,k))>abs(piv) )then
               piv = a(i,k)
               ipiv = i
            endif
         enddo
*** Exchange rows
         if( ipiv .ne. k )then
            do j=k,n
               slask = a(k,j)
               a(k,j) = a(ipiv,j)
               a(ipiv,j) = slask
            enddo
            slask = b(k)
            b(k) = b(ipiv)
            b(ipiv) = slask
         endif
*** Eliminate
         pivi = 1/a(k,k)
         do j=k+1,n
            do i=k+1,n
               a(i,j) = a(i,j) - a(i,k)*pivi*a(k,j)
            enddo
            b(j) = b(j) - a(j,k)*pivi*b(k)
         enddo
      enddo
*** Back substitution
      do k=n,1,-1
         x(k) = b(k)
         do j=k+1,n
            x(k) = x(k)-a(k,j)*x(j)
         enddo
         x(k) = x(k)/a(k,k)
      enddo

      return
      end
