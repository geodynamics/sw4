c-----------------------------------------------------------------------
      subroutine BOUNDARIES( ni, nj, up, u, um, mu, lambda, rho, dt, h, 
     *                       fo, bop, acof, ghcof )

*** 4th order SBP in x, periodic in y
      implicit none
      real*8 d4a, d4b, os, oe, ot, otf, fs, tq
      parameter( d4a=2d0/3, d4b=1d0/12, os=1d0/6, oe=1d0/8, ot=1d0/2 )
      parameter( fs=5d0/6, tq=3d0/4, otf=1d0/24 )

      integer ni, nj, wb, i, j, k, l, ind
      real*8    u(2,-1:ni+1,-1:nj+1),   mu(-1:ni+1,-1:nj+1)
      real*8 lambda(-1:ni+1,-1:nj+1),  rho(-1:ni+1,-1:nj+1)
      real*8   fo(2,-1:ni+1,-1:nj+1), up(2,-1:ni+1,-1:nj+1)
      real*8   um(2,-1:ni+1,-1:nj+1)
      real*8 h, dt, dt2, la2, cof
      real*8 r1l, r2l, c1, c2, f1l, f2l, f3l, f4l
      real*8, allocatable, dimension(:,:) :: r1, r2, f1, f2, f3, f4
c      real*8 f1(8,nj), f2(8,nj), f3(8,nj), f4(8,nj), r1(6,nj), r2(6,nj)
      real*8 acof(6,8,8), bop(4,6), ghcof

      allocate( r1(6,nj), r2(6,nj), f1(8,nj), f2(6,-1:nj+1),
     *      f3(6,-1:nj+1), f4(8,nj) )
      dt2 = dt*dt
      la2 = dt2/(h*h)
      wb = 6

      do j=1,nj-1
         do i=1,6
*** compute xx,yy and f1=lambda*v_y, f2=mu*v_x, f3=lambda*u_x, f4=mu*u_y
            r1l = 0
            r2l = 0
            do k=1,8
               c1 = 0
               c2 = 0
               do l=1,8
                  c1 = c1 + acof(i,k,l)*(2*mu(l,j)+lambda(l,j))
                  c2 = c2 + acof(i,k,l)*mu(l,j)
               enddo
               r1l = r1l + c1*u(1,k,j)
               r2l = r2l + c2*u(2,k,j)
            enddo
            if( i.eq.1 )then
               r1l = r1l + ghcof*(2*mu(1,j)+lambda(1,j))*u(1,0,j)
               r2l = r2l + ghcof*mu(1,j)*u(2,0,j)
            endif

            r1(i,j) = r1l + 
     *       ( os*mu(i,j-1)-oe*(mu(i,j-2)+mu(i,j))         )*u(1,i,j-2)
     *  +( os*(mu(i,j-2)+mu(i,j+1))+ot*(mu(i,j-1)+mu(i,j)) )*u(1,i,j-1)
     *     + ( -otf*(mu(i,j-2)+mu(i,j+2))
     *         -fs *(mu(i,j-1)+mu(i,j+1)) 
     *         -tq *mu(i,j)                                )*u(1,i,j)  
     *  +( os*(mu(i,j-1)+mu(i,j+2))+ot*(mu(i,j)+mu(i,j+1)) )*u(1,i,j+1)
     *  +   (os*mu(i,j+1)-oe*(mu(i,j)+mu(i,j+2))           )*u(1,i,j+2)

            r2(i,j) = r2l + (os*(2*mu(i,j-1)+lambda(i,j-1))-
     * oe*(2*mu(i,j-2)+2*mu(i,j)+lambda(i,j-2)+lambda(i,j)))*u(2,i,j-2)+
     * (os*(2*mu(i,j-2)+lambda(i,j-2)+2*mu(i,j+1)+lambda(i,j+1))
     * +ot*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j)+lambda(i,j)))*u(2,i,j-1)
     *     + (-otf*(2*mu(i,j-2)+lambda(i,j-2)+2*mu(i,j+2)+lambda(i,j+2))
     *  -fs*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j+1)+lambda(i,j+1)) 
     *   -tq*(2*mu(i,j)+lambda(i,j)) )*u(2,i,j) + 
     *   (os*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j+2)+lambda(i,j+2))
     * +ot*(2*mu(i,j)+lambda(i,j)+2*mu(i,j+1)+lambda(i,j+1)))*u(2,i,j+1)
     * +(os*(2*mu(i,j+1)+lambda(i,j+1))-
     * oe*(2*mu(i,j)+2*mu(i,j+2)+lambda(i,j)+lambda(i,j+2)))*u(2,i,j+2)

            f1(i,j) = lambda(i,j)*(d4b*(u(2,i,j-2)-u(2,i,j+2))+
     *                                   d4a*(u(2,i,j+1)-u(2,i,j-1)))
            f4(i,j) = mu(i,j)*(d4b*(u(1,i,j-2)-u(1,i,j+2))+
     *                                   d4a*(u(1,i,j+1)-u(1,i,j-1)))
            if( i.le.4 )then
               f2l = 0
               f3l = 0
               do k=1,wb
                  f2l = f2l + bop(i,k)*u(2,k,j)
                  f3l = f3l + bop(i,k)*u(1,k,j)
               enddo
               f2(i,j) =f2l*mu(i,j)
               f3(i,j) =f3l*lambda(i,j)
            else
               f2(i,j) = mu(i,j)*(d4b*(u(2,i-2,j)-u(2,i+2,j))+
     *              d4a*(u(2,i+1,j)-u(2,i-1,j)))
               f3(i,j) = lambda(i,j)*(d4b*(u(1,i-2,j)-u(1,i+2,j))+
     *              d4a*(u(1,i+1,j)-u(1,i-1,j)))
            endif
         enddo
         do i=7,8
            f1(i,j) = lambda(i,j)*(d4b*(u(2,i,j-2)-u(2,i,j+2))+
     *                                   d4a*(u(2,i,j+1)-u(2,i,j-1)))
            f4(i,j) = mu(i,j)*(d4b*(u(1,i,j-2)-u(1,i,j+2))+
     *                                   d4a*(u(1,i,j+1)-u(1,i,j-1)))
         enddo
      enddo
      do i=1,6
         f2(i,-1)  = f2(i,nj-2)
         f2(i,0)   = f2(i,nj-1)
         f2(i,nj)  = f2(i,1)
         f2(i,nj+1)= f2(i,2)
         f3(i,-1)  = f3(i,nj-2)
         f3(i,0)   = f3(i,nj-1)
         f3(i,nj)  = f3(i,1)
         f3(i,nj+1)= f3(i,2)
      enddo

      do j=1,nj-1
         do i=1,6
***     compute f1_x, f2_y, f3_y, f4_x
            r1(i,j) = r1(i,j) + (d4b*(f2(i,j-2)-f2(i,j+2))+
     *                                   d4a*(f2(i,j+1)-f2(i,j-1)))
            r2(i,j) = r2(i,j) + (d4b*(f3(i,j-2)-f3(i,j+2))+
     *                                   d4a*(f3(i,j+1)-f3(i,j-1)))
            if( i.le.4 )then
               f1l = 0
               f4l = 0
               do k=1,wb
                  f1l = f1l + bop(i,k)*f1(k,j)
                  f4l = f4l + bop(i,k)*f4(k,j)
               enddo
               r1(i,j) = r1(i,j) + f1l
               r2(i,j) = r2(i,j) + f4l
            else
               r1(i,j) = r1(i,j) + d4b*(f1(i-2,j)-f1(i+2,j))+
     *              d4a*(f1(i+1,j)-f1(i-1,j))
               r2(i,j) = r2(i,j) + d4b*(f4(i-2,j)-f4(i+2,j))+
     *              d4a*(f4(i+1,j)-f4(i-1,j))
            endif
*** Update solution
            cof = la2/rho(i,j)
            up(1,i,j) = 2*u(1,i,j)-um(1,i,j)+cof*r1(i,j)+ dt2*fo(1,i,j)
            up(2,i,j) = 2*u(2,i,j)-um(2,i,j)+cof*r2(i,j)+ dt2*fo(2,i,j)
         enddo
      enddo

      do j=1,nj-1
         do i=ni-5,ni
*** compute xx,yy and f1=lambda*v_y, f2=mu*v_x, f3=lambda*u_x, f4=mu*u_y
            r1l = 0
            r2l = 0
            do k=1,8
               c1 = 0
               c2 = 0
               do l=1,8
                  c1 = c1 + acof(ni-i+1,k,l)*
     *                   (2*mu(ni-l+1,j)+lambda(ni-l+1,j))
                  c2 = c2 + acof(ni-i+1,k,l)*mu(ni-l+1,j)
               enddo
               r1l = r1l + c1*u(1,ni-k+1,j)
               r2l = r2l + c2*u(2,ni-k+1,j)
            enddo
            if( i.eq.ni )then
               r1l = r1l + ghcof*(2*mu(ni,j)+lambda(ni,j))*u(1,ni+1,j)
               r2l = r2l + ghcof*mu(ni,j)*u(2,ni+1,j)
            endif
            r1(ni-i+1,j) = r1l + 
     *       (os*mu(i,j-1)- oe*(mu(i,j-2)+mu(i,j)))*u(1,i,j-2)+
     * (os*(mu(i,j-2)+mu(i,j+1))+ot*(mu(i,j-1)+mu(i,j)))*u(1,i,j-1)
     *     + (-otf*(mu(i,j-2)+mu(i,j+2))
     *  -fs*(mu(i,j-1)+mu(i,j+1)) 
     *   -tq*mu(i,j))*u(1,i,j) + 
     *   (os*(mu(i,j-1)+mu(i,j+2))
     * +ot*(mu(i,j)+mu(i,j+1)))*u(1,i,j+1)
     * +(os*mu(i,j+1)-
     * oe*(mu(i,j)+mu(i,j+2)))*u(1,i,j+2)

            r2(ni-i+1,j) = r2l + (os*(2*mu(i,j-1)+lambda(i,j-1))-
     * oe*(2*mu(i,j-2)+2*mu(i,j)+lambda(i,j-2)+lambda(i,j)))*u(2,i,j-2)+
     * (os*(2*mu(i,j-2)+lambda(i,j-2)+2*mu(i,j+1)+lambda(i,j+1))
     * +ot*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j)+lambda(i,j)))*u(2,i,j-1)
     *     + (-otf*(2*mu(i,j-2)+lambda(i,j-2)+2*mu(i,j+2)+lambda(i,j+2))
     *  -fs*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j+1)+lambda(i,j+1)) 
     *   -tq*(2*mu(i,j)+lambda(i,j)) )*u(2,i,j) + 
     *   (os*(2*mu(i,j-1)+lambda(i,j-1)+2*mu(i,j+2)+lambda(i,j+2))
     * +ot*(2*mu(i,j)+lambda(i,j)+2*mu(i,j+1)+lambda(i,j+1)))*u(2,i,j+1)
     * +(os*(2*mu(i,j+1)+lambda(i,j+1))-
     * oe*(2*mu(i,j)+2*mu(i,j+2)+lambda(i,j)+lambda(i,j+2)))*u(2,i,j+2)

            f1(ni-i+1,j) = lambda(i,j)*(d4b*(u(2,i,j-2)-u(2,i,j+2))+
     *                                   d4a*(u(2,i,j+1)-u(2,i,j-1)))
            f4(ni-i+1,j) = mu(i,j)*(d4b*(u(1,i,j-2)-u(1,i,j+2))+
     *                                   d4a*(u(1,i,j+1)-u(1,i,j-1)))
            if( ni-i+1.le.4 )then
               f2l = 0
               f3l = 0
               do k=1,wb
                  f2l = f2l - bop(ni-i+1,k)*u(2,ni-k+1,j)
                  f3l = f3l - bop(ni-i+1,k)*u(1,ni-k+1,j)
               enddo
               f2(ni-i+1,j) =f2l*mu(i,j)
               f3(ni-i+1,j) =f3l*lambda(i,j)
            else
               f2(ni-i+1,j) = mu(i,j)*(d4b*(u(2,i-2,j)-u(2,i+2,j))+
     *              d4a*(u(2,i+1,j)-u(2,i-1,j)))
               f3(ni-i+1,j) = lambda(i,j)*(d4b*(u(1,i-2,j)-u(1,i+2,j))+
     *              d4a*(u(1,i+1,j)-u(1,i-1,j)))
            endif
         enddo
         do i=ni-7,ni-6
            f1(ni-i+1,j) = lambda(i,j)*(d4b*(u(2,i,j-2)-u(2,i,j+2))+
     *                                  d4a*(u(2,i,j+1)-u(2,i,j-1)))
            f4(ni-i+1,j) = mu(i,j)*(d4b*(u(1,i,j-2)-u(1,i,j+2))+
     *                              d4a*(u(1,i,j+1)-u(1,i,j-1)))
         enddo
      enddo
      do i=1,6
         f2(i,-1)  = f2(i,nj-2)
         f2(i,0)   = f2(i,nj-1)
         f2(i,nj)  = f2(i,1)
         f2(i,nj+1)= f2(i,2)
         f3(i,-1)  = f3(i,nj-2)
         f3(i,0)   = f3(i,nj-1)
         f3(i,nj)  = f3(i,1)
         f3(i,nj+1)= f3(i,2)
      enddo

      do j=1,nj-1
         do i=ni-5,ni
***     compute f1_x, f2_y, f3_y, f4_x
            ind = ni-i+1
            r1(ni-i+1,j) = r1(ni-i+1,j) +(d4b*(f2(ind,j-2)-f2(ind,j+2))+
     *                                   d4a*(f2(ind,j+1)-f2(ind,j-1)))
            r2(ni-i+1,j) = r2(ni-i+1,j) +(d4b*(f3(ind,j-2)-f3(ind,j+2))+
     *                                   d4a*(f3(ind,j+1)-f3(ind,j-1)))
            if( ni-i+1.le.4 )then
               f1l = 0
               f4l = 0
               do k=1,wb
                  f1l = f1l - bop(ni-i+1,k)*f1(k,j)
                  f4l = f4l - bop(ni-i+1,k)*f4(k,j)
               enddo
               r1(ni-i+1,j) = r1(ni-i+1,j) + f1l
               r2(ni-i+1,j) = r2(ni-i+1,j) + f4l
            else
               r1(ni-i+1,j) = r1(ni-i+1,j) - (
     *              d4b*(f1(ind-2,j)-f1(ind+2,j))+
     *              d4a*(f1(ind+1,j)-f1(ind-1,j)))
               r2(ni-i+1,j) = r2(ni-i+1,j) - (
     *              d4b*(f4(ind-2,j)-f4(ind+2,j))+
     *              d4a*(f4(ind+1,j)-f4(ind-1,j)))
            endif
*** Update solution
            cof = la2/rho(i,j)
            up(1,i,j) = 2*u(1,i,j)-um(1,i,j)+cof*r1(ni-i+1,j) + 
     *                     dt2*fo(1,i,j)
            up(2,i,j) = 2*u(2,i,j)-um(2,i,j)+cof*r2(ni-i+1,j) +
     *                     dt2*fo(2,i,j)
         enddo
      enddo

      deallocate(r1,r2,f1,f2,f3,f4)

      end
