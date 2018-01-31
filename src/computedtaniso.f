      subroutine COMPUTEDTANISO( ifirst, ilast, jfirst, jlast, kfirst,
     *                           klast, rho, c, cfl, dx, dtloc ) 
***********************************************************************
***
*** Estimate the spectral radius of the discretized operator by
***  (a_1+a_2+....+a_6)/3, where a_i are the eigenvalues
*** of the 6x6 matrix C that occurs in the energy, 
***     energy = (u_t,rho*u_t) +  (e,Ce)
*** (u,v) is the scalar product over the 3D domain,
*** and e are the strains, rho density, u displacements.
*** Only reason for this estimate, is that it comes fairly close
*** to the exact spectral radius in the case of isotropic material.
***
***********************************************************************

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer info
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 dtloc, dtgp, cfl, dx
      real*8 a(21), eig(6), work(18), z, eigestimate

      dtloc = 1d38
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               a(1)  = c(1,i,j,k)
               a(2)  = 2*c(2,i,j,k)
               a(3)  = 2*c(3,i,j,k)
               a(4)  = c(4,i,j,k)
               a(5)  = 2*c(5,i,j,k)
               a(6)  = c(6,i,j,k)
               a(7)  = 4*c(7,i,j,k)
               a(8)  = 4*c(8,i,j,k)
               a(9)  = 2*c(9,i,j,k)
               a(10)  = 4*c(10,i,j,k)
               a(11) = 2*c(11,i,j,k)
               a(12) = 4*c(12,i,j,k)
               a(13) = 2*c(13,i,j,k)
               a(14) = 4*c(14,i,j,k)
               a(15) = 2*c(15,i,j,k)
               a(16) = c(16,i,j,k)
               a(17) = 2*c(17,i,j,k)
               a(18) = c(18,i,j,k)
               a(19) = 4*c(19,i,j,k)
               a(20) = 2*c(20,i,j,k)
               a(21) = c(21,i,j,k)
               call DSPEV('N', 'L', 6, a, eig, z, 1, work, info )
               if( info .ne. 0 )then
                  write(*,*) 'ERROR in computedtaniso:',
     *        ' could not compute eigenvalues. info from DSPEV = ',
     *             info
               endif
               eigestimate=(eig(1)+eig(2)+eig(3)+eig(4)+eig(5)+eig(6))/3
               dtgp = cfl*SQRT(rho(i,j,k)/eigestimate)*dx
               if( dtgp .lt. dtloc )then
                  dtloc = dtgp
               endif
            enddo
         enddo
      enddo
      end
c-----------------------------------------------------------------------
      subroutine COMPUTEDTANISO2( ifirst, ilast, jfirst, jlast, kfirst,
     *                           klast, rho, c, cfl, dx, dtloc ) bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
c      integer nphi, nth
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 dtloc, dtgp, cfl, dx, eigestimate
c      real*8, allocatable, dimension(:) :: ws
c      nphi = 40
c      nth  = 10
c      allocate(ws(2*nth))
      dtloc = 1d38
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
c               call MAXWAVE2( c(1,i,j,k), rho(i,j,k), eigestimate, 
c     *                       nphi, nth, ws )
               call MAXWAVE( c(1,i,j,k), rho(i,j,k), eigestimate ) 
               dtgp = cfl*dx/SQRT(eigestimate)
               if( dtgp .lt. dtloc )then
                  dtloc = dtgp
               endif
            enddo
         enddo
      enddo
c      deallocate(ws)
      end
c-----------------------------------------------------------------------
      subroutine MAXWAVE( c, rho, wavespeed2 )
*** Compute the maximum of the sum of eigenvalues of the matrix 
***   \sum_{i,j} ki* kj* Aij
***   where \sum_{i} ki*ki = 1, scaled by the density.
*** This is equal to (4*mu+lambda)/rho in the isotropic case 
      implicit none
      integer info
      real*8 c(21), rho, wavespeed2
      real*8 eg(3), a(6), work(9), z
      a(1) = c(1)+c(7)+c(12)
      a(2) = c(2)+c(9)+c(14)
      a(3) = c(3)+c(10)+c(15)
      a(4) = c(7)+c(16)+c(19)
      a(5) = c(8)+c(17)+c(20)
      a(6) = c(12)+c(19)+c(21)
      call DSPEV('N', 'L', 3, a, eg, z, 1, work, info )
      if( info .ne. 0 )then
         write(*,*) 'ERROR in maxwave:',
     *        ' could not compute eigenvalues. info from DSPEV = ',
     *        info
      endif
      wavespeed2 = MAX(eg(1),eg(2),eg(3))/rho
      end
c-----------------------------------------------------------------------
      subroutine MAXWAVE2( c, rho, wavespeed2, nphi, nth, ws )
      implicit none
      real*8 pi
      parameter(pi=3.14159265358979d0)
      integer i, j, nth, nphi
      real*8 c(21), rho, wavespeed2, eg(3), wslocal, ws(2,nth)
      real*8 ith, iph, cphi, sphi
      ith = pi/(2*(nth-1))
      iph = 2*pi/(nphi-1)
      do i=1,nth
         ws(1,i) = COS( (i-1d0)*ith )
         ws(2,i) = SIN( (i-1d0)*ith )
      enddo
      wavespeed2 = -1d38
      do i=1,nth
         do j=1,nphi-1
            cphi = COS( (j-1d0)*iph )
            sphi = SIN( (j-1d0)*iph )
            call EIGDIRECTION( ws(2,i)*cphi, ws(2,i)*sphi, ws(1,i), 
     *                         c, eg )
            wslocal = (eg(1)+eg(2)+eg(3))/rho
            if( wslocal.gt.wavespeed2 )then
               wavespeed2 = wslocal
            endif
         enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine EIGDIRECTION( k1, k2, k3, c, eigs )
      implicit none
      integer info
      real*8 k1, k2, k3, c(21)
      real*8 eigs(3), a(6), work(9), z
c      real*8 mat(3,3)
c      mat(1,1) = k1*k1*c(1) +k1*k2*2*c(2) +k1*k3*2*c(3) + 
c     *                   k2*k2*c(7) + k2*k3*2*c(8) + k3*k3*c(12)
c      mat(1,2) = k1*k1*c(2) +k1*k2*(c(4)+c(7)) + k1*k3*(c(5)+c(8)) + 
c     *                   k2*k2*c(9) + k2*k3*(c(10)+c(13)) + k3*k3*c(14)
c      mat(2,1) = mat(1,2)
c      mat(1,3) = k1*k1*c(3) +k1*k2*(c(5)+c(8)) + k1*k3*(c(6)+c(12)) + 
c     *                   k2*k2*c(10) + k2*k3*(c(11)+c(14)) + k3*k3*c(15)
c      mat(3,1) = mat(1,3)
c      mat(2,2) = k1*k1*c(7) +k1*k2*2*c(9) +k1*k3*2*c(10) + 
c     *                   k2*k2*c(16) + k2*k3*2*c(17) + k3*k3*c(19)
c      mat(2,3) = k1*k1*c(8) +k1*k2*(c(10)+c(13)) + k1*k3*(c(11)+c(14)) + 
c     *                   k2*k2*c(17) + k2*k3*(c(18)+c(19)) + k3*k3*c(20)
c      mat(3,2) = mat(2,3)
c      mat(3,3) = k1*k1*c(12) +k1*k2*2*c(14) +k1*k3*2*c(15) + 
c     *                   k2*k2*c(19) + k2*k3*2*c(20) + k3*k3*c(21)
      a(1) = k1*k1*c(1) +k1*k2*2*c(2) +k1*k3*2*c(3) + 
     *                k2*k2*c(7) + k2*k3*2*c(8) + k3*k3*c(12)
      a(2) = k1*k1*c(2) +k1*k2*(c(4)+c(7)) + k1*k3*(c(5)+c(8)) + 
     *                k2*k2*c(9) + k2*k3*(c(10)+c(13)) + k3*k3*c(14)
      a(3) = k1*k1*c(3) +k1*k2*(c(5)+c(8)) + k1*k3*(c(6)+c(12)) + 
     *                k2*k2*c(10) + k2*k3*(c(11)+c(14)) + k3*k3*c(15)
      a(4) = k1*k1*c(7) +k1*k2*2*c(9) +k1*k3*2*c(10) + 
     *                k2*k2*c(16) + k2*k3*2*c(17) + k3*k3*c(19)
      a(5) = k1*k1*c(8) +k1*k2*(c(10)+c(13)) + k1*k3*(c(11)+c(14)) + 
     *                k2*k2*c(17) + k2*k3*(c(18)+c(19)) + k3*k3*c(20)
      a(6) = k1*k1*c(12) +k1*k2*2*c(14) +k1*k3*2*c(15) + 
     *                k2*k2*c(19) + k2*k3*2*c(20) + k3*k3*c(21)
      call DSPEV('N', 'L', 3, a, eigs, z, 1, work, info )
      if( info .ne. 0 )then
         write(*,*) 'ERROR in eigdirection:',
     *        ' could not compute eigenvalues. info from DSPEV = ',
     *        info
      endif
      end
      
c-----------------------------------------------------------------------
      subroutine COMPUTEDTANISO2CURV( ifirst, ilast, jfirst, jlast, 
     *                  kfirst, klast, rho, c, jac, cfl, dtloc )bind(c)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
c      integer nphi, nth
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(45,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 dtloc, dtgp, cfl, eigestimate
c      real*8, allocatable, dimension(:) :: ws
c      nphi = 40
c      nth  = 10
c      allocate(ws(2*nth))
      dtloc = 1d38
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
c               call MAXWAVE2( c(1,i,j,k), rho(i,j,k), eigestimate, 
c     *                       nphi, nth, ws )
               call MAXWAVECURV( c(1,i,j,k), rho(i,j,k), jac(i,j,k), 
     *                           eigestimate )
               dtgp = cfl/SQRT( eigestimate )
               if( dtgp .lt. dtloc )then
                  dtloc = dtgp
               endif
            enddo
         enddo
      enddo
c      deallocate(ws)
      end

c-----------------------------------------------------------------------
      subroutine MAXWAVECURV( c, rho, jac, wavespeed2 )
*** Compute the maximum of the sum of eigenvalues of the matrix 
***   \sum_{i,j} ki* kj* Aij
***   where \sum_{i} ki*ki = 1, scaled by the density.
*** This is equal to (4*mu+lambda)/rho in the isotropic case 
      implicit none
      integer info
      real*8 c(45), rho, wavespeed2, jac
      real*8 eg(3), a(6), work(9), z
*** Traces of matrices, in order xx,xy,xz,yy,yz,zz
      a(1) = c(1)+c(4)+c(6)
      a(2) = c(19)+c(23)+c(27)
      a(3) = c(28)+c(32)+c(36)
      a(4) = c(7)+c(10)+c(12)
      a(5) = c(37)+c(41)+c(45)
      a(6) = c(13)+c(16)+c(18)
      call DSPEV('N', 'L', 3, a, eg, z, 1, work, info )
      if( info .ne. 0 )then
         write(*,*) 'ERROR in maxwave:',
     *        ' could not compute eigenvalues. info from DSPEV = ',
     *        info
      endif
      wavespeed2 = MAX(eg(1),eg(2),eg(3))/(jac*rho)
      end
