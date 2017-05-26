      subroutine bcfreesurfcurvani( ifirst, ilast, jfirst, jlast,kfirst, 
     +     klast, nz, u, c, side, sbop, bforce5, bforce6, strx, stry )
     + bind(c)
      implicit none
      real*8 a1, a2
      parameter( a1=2d0/3, a2=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nz
      integer side, i, j, k, qq, ipiv(3), info
      real*8 sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(45,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 du, dv, dw, rhs1, rhs2, rhs3, x(3), a(9), s0i


      s0i = 1/sbop(0)
!$OMP PARALLEL PRIVATE(i,j,k,qq,du,dv,dw,rhs1,rhs2,rhs3,x,a,ipiv,info)
      if( side.eq.5 )then
         k = 1
!$OMP DO               
         do j=jfirst+2,jlast-2
            do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
               qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
***  x-derivatives
               du = strx(i)*(a2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *                       a1*(u(1,i+1,j,k)-u(1,i-1,j,k)))
               dv = strx(i)*(a2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *                       a1*(u(2,i+1,j,k)-u(2,i-1,j,k)))
               dw = strx(i)*(a2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *                       a1*(u(3,i+1,j,k)-u(3,i-1,j,k)))

               rhs1 =c(28,i,j,k)*du + c(31,i,j,k)*dv + c(34,i,j,k)*dw
               rhs2 =c(29,i,j,k)*du + c(32,i,j,k)*dv + c(35,i,j,k)*dw
               rhs3 =c(30,i,j,k)*du + c(33,i,j,k)*dv + c(36,i,j,k)*dw
***  y-derivatives
               du = stry(j)*(a2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *                       a1*(u(1,i,j+1,k)-u(1,i,j-1,k)))
               dv = stry(j)*(a2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *                       a1*(u(2,i,j+1,k)-u(2,i,j-1,k)))
               dw = stry(j)*(a2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *                       a1*(u(3,i,j+1,k)-u(3,i,j-1,k)))
               rhs1 = rhs1 + 
     *                 c(37,i,j,k)*du + c(40,i,j,k)*dv + c(43,i,j,k)*dw
               rhs2 = rhs2 + 
     *                 c(38,i,j,k)*du + c(41,i,j,k)*dv + c(44,i,j,k)*dw
               rhs3 = rhs3 + 
     *                 c(39,i,j,k)*du + c(42,i,j,k)*dv + c(45,i,j,k)*dw

*** z-derivatives, excepting ghost points.
               du = sbop(1)*u(1,i,j,1)+sbop(2)*u(1,i,j,2)+
     *              sbop(3)*u(1,i,j,3)+sbop(4)*u(1,i,j,4)
               dv = sbop(1)*u(2,i,j,1)+sbop(2)*u(2,i,j,2)+
     *              sbop(3)*u(2,i,j,3)+sbop(4)*u(2,i,j,4)
               dw = sbop(1)*u(3,i,j,1)+sbop(2)*u(3,i,j,2)+
     *              sbop(3)*u(3,i,j,3)+sbop(4)*u(3,i,j,4)
               rhs1 = rhs1 + c(13,i,j,k)*du  + c(14,i,j,k)*dv + 
     *                       c(15,i,j,k)*dw - bforce5(1,qq)
               rhs2 = rhs2 + c(14,i,j,k)*du  + c(16,i,j,k)*dv + 
     *                       c(17,i,j,k)*dw - bforce5(2,qq)
               rhs3 = rhs3 + c(15,i,j,k)*du  + c(17,i,j,k)*dv + 
     *                       c(18,i,j,k)*dw - bforce5(3,qq)
*** Solve symmetric system for ghost point values
               x(1) = rhs1
               x(2) = rhs2
               x(3) = rhs3
               a(1) = c(13,i,j,k)
               a(2) = c(14,i,j,k)
               a(3) = c(15,i,j,k)
               a(4) = c(14,i,j,k)
               a(5) = c(16,i,j,k)
               a(6) = c(17,i,j,k)
               a(7) = c(15,i,j,k)
               a(8) = c(17,i,j,k)
               a(9) = c(18,i,j,k)
               call DGESV( 3, 1, a, 3, ipiv, x, 3, info )
               if( info .ne. 0 )then
                  write(*,*)'ERROR in enforceBCanisotropic: info = ',
     *                      info, ' from routine DGESV'
               endif
               u(1,i,j,0) = -s0i*x(1)
               u(2,i,j,0) = -s0i*x(2)
               u(3,i,j,0) = -s0i*x(3)
            enddo
         enddo
!$OMP ENDDO               
      else
c s=6
         k = nz
!$OMP DO               
         do j=jfirst+2,jlast-2
            do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
               qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
***  x-derivatives
               du = strx(i)*(a2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *                       a1*(u(1,i+1,j,k)-u(1,i-1,j,k)))
               dv = strx(i)*(a2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *                       a1*(u(2,i+1,j,k)-u(2,i-1,j,k)))
               dw = strx(i)*(a2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *                       a1*(u(3,i+1,j,k)-u(3,i-1,j,k)))

               rhs1 =c(28,i,j,k)*du + c(31,i,j,k)*dv + c(34,i,j,k)*dw
               rhs2 =c(29,i,j,k)*du + c(32,i,j,k)*dv + c(35,i,j,k)*dw
               rhs3 =c(30,i,j,k)*du + c(33,i,j,k)*dv + c(36,i,j,k)*dw
***  y-derivatives
               du = stry(j)*(a2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *                       a1*(u(1,i,j+1,k)-u(1,i,j-1,k)))
               dv = stry(j)*(a2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *                       a1*(u(2,i,j+1,k)-u(2,i,j-1,k)))
               dw = stry(j)*(a2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *                       a1*(u(3,i,j+1,k)-u(3,i,j-1,k)))
               rhs1 = rhs1 + 
     *                 c(37,i,j,k)*du + c(40,i,j,k)*dv + c(43,i,j,k)*dw
               rhs2 = rhs2 + 
     *                 c(38,i,j,k)*du + c(41,i,j,k)*dv + c(44,i,j,k)*dw
               rhs3 = rhs3 + 
     *                 c(39,i,j,k)*du + c(42,i,j,k)*dv + c(45,i,j,k)*dw

*** z-derivatives, excepting ghost points.
               du = -(sbop(1)*u(1,i,j,k)  +sbop(2)*u(1,i,j,k-1)+
     *                sbop(3)*u(1,i,j,k-2)+sbop(4)*u(1,i,j,k-3))
               dv = -(sbop(1)*u(2,i,j,k)  +sbop(2)*u(2,i,j,k-1)+
     *                sbop(3)*u(2,i,j,k-2)+sbop(4)*u(2,i,j,k-3))
               dw = -(sbop(1)*u(3,i,j,k)  +sbop(2)*u(3,i,j,k-1)+
     *                sbop(3)*u(3,i,j,k-2)+sbop(4)*u(3,i,j,k-3))
               rhs1 = rhs1 + c(13,i,j,k)*du  + c(14,i,j,k)*dv + 
     *                       c(15,i,j,k)*dw - bforce6(1,qq)
               rhs2 = rhs2 + c(14,i,j,k)*du  + c(16,i,j,k)*dv + 
     *                       c(17,i,j,k)*dw - bforce6(2,qq)
               rhs3 = rhs3 + c(15,i,j,k)*du  + c(17,i,j,k)*dv + 
     *                       c(18,i,j,k)*dw - bforce6(3,qq)

*** Solve symmetric system for ghost point values
               x(1) = rhs1
               x(2) = rhs2
               x(3) = rhs3
               a(1) = c(13,i,j,k)
               a(2) = c(14,i,j,k)
               a(3) = c(15,i,j,k)
               a(4) = c(14,i,j,k)
               a(5) = c(16,i,j,k)
               a(6) = c(17,i,j,k)
               a(7) = c(15,i,j,k)
               a(8) = c(17,i,j,k)
               a(9) = c(18,i,j,k)
               call DGESV( 3, 1, a, 3, ipiv, x, 3, info )
               if( info .ne. 0 )then
                  write(*,*)'ERROR in enforceBCanisotropic: info = ',
     *                      info, ' from routine DGESV'
               endif
               u(1,i,j,k+1) = s0i*x(1)
               u(2,i,j,k+1) = s0i*x(2)
               u(3,i,j,k+1) = s0i*x(3)
            enddo
         enddo
!$OMP ENDDO               
      endif
!$OMP END PARALLEL
      end

