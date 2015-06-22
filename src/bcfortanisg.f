c-----------------------------------------------------------------------
      subroutine bcfortanisg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, wind, nx, ny, nz,
     +     u, h, bccnd, sbop, c,
     *     bforce1, bforce2, bforce3, bforce4, bforce5, bforce6, 
     +     strx, stry )
      implicit none
      real*8 a1, a2
      parameter( a1=2d0/3, a2=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6,6), i, j, k, bccnd(6), qq, ipiv(3), info
      real*8 h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
c note that the numbering of bforce adds one from C (side goes from 1 in Fortran)
      real*8 bforce1(3,*),  bforce2(3,*)
      real*8 bforce3(3,*),  bforce4(3,*)
      real*8 bforce5(3,*),  bforce6(3,*)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 du, dv, dw, rhs1, rhs2, rhs3, x(3), a(9), s0i

c the boundary window 'wind' is now an input argument

c loop over all sides of the 3-D domain
      do s=1,6
*** dirichlet condition, bccnd=1
*** supergrid condition, bccnd=2
c now assigning the forcing arrays outside of this routine!
        if( bccnd(s).eq.1 .or. bccnd(s).eq.2)then
            qq=1
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce1(1,qq)
                    u(2,i,j,k) = bforce1(2,qq)
                    u(3,i,j,k) = bforce1(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce2(1,qq)
                    u(2,i,j,k) = bforce2(2,qq)
                    u(3,i,j,k) = bforce2(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce3(1,qq)
                    u(2,i,j,k) = bforce3(2,qq)
                    u(3,i,j,k) = bforce3(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce4(1,qq)
                    u(2,i,j,k) = bforce4(2,qq)
                    u(3,i,j,k) = bforce4(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce5(1,qq)
                    u(2,i,j,k) = bforce5(2,qq)
                    u(3,i,j,k) = bforce5(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            else if (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = bforce6(1,qq)
                    u(2,i,j,k) = bforce6(2,qq)
                    u(3,i,j,k) = bforce6(3,qq)
                    qq = qq+1
                  enddo
                enddo
              enddo
            endif

          elseif( bccnd(s).eq.3 )then
*** Periodic condition, bccnd=3
            if (s.eq.1) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i+nx,j,k)
                    u(2,i,j,k) = u(2,i+nx,j,k)
                    u(3,i,j,k) = u(3,i+nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.2) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i-nx,j,k)
                    u(2,i,j,k) = u(2,i-nx,j,k)
                    u(3,i,j,k) = u(3,i-nx,j,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.3) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j+ny,k)
                    u(2,i,j,k) = u(2,i,j+ny,k)
                    u(3,i,j,k) = u(3,i,j+ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.4) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j-ny,k)
                    u(2,i,j,k) = u(2,i,j-ny,k)
                    u(3,i,j,k) = u(3,i,j-ny,k)
                  enddo
                enddo
              enddo
            elseif (s.eq.5) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k+nz)
                    u(2,i,j,k) = u(2,i,j,k+nz)
                    u(3,i,j,k) = u(3,i,j,k+nz)
                  enddo
                enddo
              enddo
            elseif (s.eq.6) then
              do k=wind(5,s),wind(6,s)
                do j=wind(3,s),wind(4,s)
                  do i=wind(1,s),wind(2,s)
                    u(1,i,j,k) = u(1,i,j,k-nz)
                    u(2,i,j,k) = u(2,i,j,k-nz)
                    u(3,i,j,k) = u(3,i,j,k-nz)
                  enddo
                enddo
              enddo
            endif
          elseif( bccnd(s).eq.0 )then
*** Free surface condition, bccnd=0
            if( s.ne.5 .and. s.ne.6 )then
              write(*,*) 
     *             'ERROR: Free surface condition ',
     +             'not implemented for side ', s
              stop
            endif

c moved the assignment of bforce5/6 into its own routine

**** Do the free surface condition (
            s0i = 1/sbop(0)
            if( s.eq.5 )then
               k = 1
               do j=jfirst+2,jlast-2
                  do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                  qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
***  x-derivatives
                  du = strx(i)*(a2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *                          a1*(u(1,i+1,j,k)-u(1,i-1,j,k)))
                  dv = strx(i)*(a2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *                          a1*(u(2,i+1,j,k)-u(2,i-1,j,k)))
                  dw = strx(i)*(a2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *                          a1*(u(3,i+1,j,k)-u(3,i-1,j,k)))

                  rhs1 = c(3,i,j,k)*du +  c(8,i,j,k)*dv + c(12,i,j,k)*dw
                  rhs2 = c(5,i,j,k)*du + c(10,i,j,k)*dv + c(14,i,j,k)*dw
                  rhs3 = c(6,i,j,k)*du + c(11,i,j,k)*dv + c(15,i,j,k)*dw
***  y-derivatives
                  du = stry(j)*(a2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *                          a1*(u(1,i,j+1,k)-u(1,i,j-1,k)))
                  dv = stry(j)*(a2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *                          a1*(u(2,i,j+1,k)-u(2,i,j-1,k)))
                  dw = stry(j)*(a2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *                          a1*(u(3,i,j+1,k)-u(3,i,j-1,k)))
                  rhs1 = rhs1 + c(8,i,j,k)*du  + c(13,i,j,k)*dv + 
     *                                                  c(14,i,j,k)*dw
                  rhs2 = rhs2 + c(10,i,j,k)*du + c(17,i,j,k)*dv + 
     *                                                  c(19,i,j,k)*dw
                  rhs3 = rhs3 + c(11,i,j,k)*du + c(18,i,j,k)*dv + 
     *                                                  c(20,i,j,k)*dw

*** z-derivatives, excepting ghost points.
                  du = sbop(1)*u(1,i,j,1)+sbop(2)*u(1,i,j,2)+
     *                            sbop(3)*u(1,i,j,3)+sbop(4)*u(1,i,j,4)
                  dv = sbop(1)*u(2,i,j,1)+sbop(2)*u(2,i,j,2)+
     *                            sbop(3)*u(2,i,j,3)+sbop(4)*u(2,i,j,4)
                  dw = sbop(1)*u(3,i,j,1)+sbop(2)*u(3,i,j,2)+
     *                            sbop(3)*u(3,i,j,3)+sbop(4)*u(3,i,j,4)
                  rhs1 = rhs1 + c(12,i,j,k)*du  + c(14,i,j,k)*dv + 
     *                             c(15,i,j,k)*dw - h*bforce5(1,qq)
                  rhs2 = rhs2 + c(14,i,j,k)*du  + c(19,i,j,k)*dv + 
     *                             c(20,i,j,k)*dw - h*bforce5(2,qq)
                  rhs3 = rhs3 + c(15,i,j,k)*du  + c(20,i,j,k)*dv + 
     *                             c(21,i,j,k)*dw - h*bforce5(3,qq)


*** Solve symmetric system for ghost point values
                  x(1) = rhs1
                  x(2) = rhs2
                  x(3) = rhs3
                  a(1) = c(12,i,j,k)
                  a(2) = c(14,i,j,k)
                  a(3) = c(15,i,j,k)
                  a(4) = c(14,i,j,k)
                  a(5) = c(19,i,j,k)
                  a(6) = c(20,i,j,k)
                  a(7) = c(15,i,j,k)
                  a(8) = c(20,i,j,k)
                  a(9) = c(21,i,j,k)
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
          else
c s=6
            k = nz
               do j=jfirst+2,jlast-2
                  do i=ifirst+2,ilast-2
** compute 1-d index in forcing array
                  qq = i-ifirst+1 + (j-jfirst)*(ilast-ifirst+1)
***  x-derivatives
                  du = strx(i)*(a2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *                          a1*(u(1,i+1,j,k)-u(1,i-1,j,k)))
                  dv = strx(i)*(a2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *                          a1*(u(2,i+1,j,k)-u(2,i-1,j,k)))
                  dw = strx(i)*(a2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *                          a1*(u(3,i+1,j,k)-u(3,i-1,j,k)))

                  rhs1 = c(3,i,j,k)*du +  c(8,i,j,k)*dv + c(12,i,j,k)*dw
                  rhs2 = c(5,i,j,k)*du + c(10,i,j,k)*dv + c(14,i,j,k)*dw
                  rhs3 = c(6,i,j,k)*du + c(11,i,j,k)*dv + c(15,i,j,k)*dw
***  y-derivatives
                  du = stry(j)*(a2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *                          a1*(u(1,i,j+1,k)-u(1,i,j-1,k)))
                  dv = stry(j)*(a2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *                          a1*(u(2,i,j+1,k)-u(2,i,j-1,k)))
                  dw = stry(j)*(a2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *                          a1*(u(3,i,j+1,k)-u(3,i,j-1,k)))
                  rhs1 = rhs1 + c(8,i,j,k)*du  + c(13,i,j,k)*dv + 
     *                                                  c(14,i,j,k)*dw
                  rhs2 = rhs2 + c(10,i,j,k)*du + c(17,i,j,k)*dv + 
     *                                                  c(19,i,j,k)*dw
                  rhs3 = rhs3 + c(11,i,j,k)*du + c(18,i,j,k)*dv + 
     *                                                  c(20,i,j,k)*dw

*** z-derivatives, excepting ghost points.
                  du = -(sbop(1)*u(1,i,j,k)+sbop(2)*u(1,i,j,k-1)+
     *                        sbop(3)*u(1,i,j,k-2)+sbop(4)*u(1,i,j,k-3))
                  dv = -(sbop(1)*u(2,i,j,k)+sbop(2)*u(2,i,j,k-1)+
     *                        sbop(3)*u(2,i,j,k-2)+sbop(4)*u(2,i,j,k-3))
                  dw = -(sbop(1)*u(3,i,j,k)+sbop(2)*u(3,i,j,k-1)+
     *                        sbop(3)*u(3,i,j,k-2)+sbop(4)*u(3,i,j,k-3))
                  rhs1 = rhs1 + c(12,i,j,k)*du  + c(14,i,j,k)*dv + 
     *                             c(15,i,j,k)*dw - h*bforce6(1,qq)
                  rhs2 = rhs2 + c(14,i,j,k)*du  + c(19,i,j,k)*dv + 
     *                             c(20,i,j,k)*dw - h*bforce6(2,qq)
                  rhs3 = rhs3 + c(15,i,j,k)*du  + c(20,i,j,k)*dv + 
     *                             c(21,i,j,k)*dw - h*bforce6(3,qq)
c AP testing
c$$$                  rhs1 = rhs1 + c(12,i,j,k)*du  + c(14,i,j,k)*dv + 
c$$$     *                             c(15,i,j,k)*dw + h*bforce6(1,qq)
c$$$                  rhs2 = rhs2 + c(14,i,j,k)*du  + c(19,i,j,k)*dv + 
c$$$     *                             c(20,i,j,k)*dw + h*bforce6(2,qq)
c$$$                  rhs3 = rhs3 + c(15,i,j,k)*du  + c(20,i,j,k)*dv + 
c$$$     *                             c(21,i,j,k)*dw + h*bforce6(3,qq)
c$$$
*** Solve symmetric system for ghost point values
                  x(1) = rhs1
                  x(2) = rhs2
                  x(3) = rhs3
                  a(1) = c(12,i,j,k)
                  a(2) = c(14,i,j,k)
                  a(3) = c(15,i,j,k)
                  a(4) = c(14,i,j,k)
                  a(5) = c(19,i,j,k)
                  a(6) = c(20,i,j,k)
                  a(7) = c(15,i,j,k)
                  a(8) = c(20,i,j,k)
                  a(9) = c(21,i,j,k)
                  call DGESV( 3, 1, a, 3, ipiv, x, 3, info )
                  if( info .ne. 0 )then
                     write(*,*)'ERROR in enforceBCanisotropic: info = ' 
     *                ,info, ' from routine DGESV'
                  endif
                  u(1,i,j,k+1) = s0i*x(1)
                  u(2,i,j,k+1) = s0i*x(2)
                  u(3,i,j,k+1) = s0i*x(3)
               enddo
            enddo
          endif
        endif
      enddo
      end

