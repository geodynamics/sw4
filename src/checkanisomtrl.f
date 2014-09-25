      subroutine CHECKANISOMTRL( ifirst, ilast, jfirst, jlast, kfirst,
     *                           klast, rho, c, 
     *                           rhomin, rhomax, eigmin, eigmax )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer m, info
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rhomin, rhomax, eigmin, eigmax
      real*8 a(21), eig(6), work(18), z

      rhomin =  1d38
      rhomax = -1d38
      eigmin =  1d38
      eigmax = -1d38
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               if( rho(i,j,k) .lt. rhomin )then
                  rhomin = rho(i,j,k)
               endif
               if( rho(i,j,k) .gt. rhomax )then
                  rhomax = rho(i,j,k)
               endif
               do m=1,21
                  a(m) = c(m,i,j,k)
               enddo
               call DSPEV('N', 'L', 6, a, eig, z, 1, work, info )
               if( info .ne. 0 )then
                  write(*,*) 'ERROR in check_anisotropic_material:',
     *        ' could not compute eigenvalues. info from DSPEV = ',
     *         info
               endif
               if( eig(1) .lt. eigmin )then
                  eigmin = eig(1)
               endif
               if( eig(6) .gt. eigmax )then
                  eigmax = eig(6)
               endif
            enddo
         enddo
      enddo
      end
