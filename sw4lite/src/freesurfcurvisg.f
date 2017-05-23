c-----------------------------------------------------------------------
      subroutine FREESURFCURVISG( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, nz, side, u, mu, la, met, s, forcing, strx, stry )
      implicit none
      real*8 c1, c2
      parameter( c1=2d0/3, c2=-1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k, kl, nz, side
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 forcing(3,ifirst:ilast,jfirst:jlast)
      real*8 strx(ifirst:ilast), stry(jfirst:jlast)
      real*8 s(0:4), rhs1, rhs2, rhs3, s0i, ac, bc, cc, dc
      real*8 istry, istrx, xoysqrt, yoxsqrt, isqrtxy

      if( side.eq.5 )then
         k = 1
         kl= 1
      elseif( side.eq.6 )then
         k = nz
         kl= -1
      endif

      s0i = 1/s(0)
      do j=jfirst+2,jlast-2
         istry = 1/stry(j)
         do i=ifirst+2,ilast-2
            istrx = 1/strx(i)

*** First tangential derivatives
            rhs1 = 
*** pr
     *   (2*mu(i,j,k)+la(i,j,k))*met(2,i,j,k)*met(1,i,j,k)*(
     *          c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *          c1*(u(1,i+1,j,k)-u(1,i-1,j,k))  )*strx(i)*istry 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  ) 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*istry   
*** qr
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   )*istrx*stry(j) 
     *  + la(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )  -
     *                 forcing(1,i,j)

*** (v-eq)
            rhs2 = 
*** pr
     *    la(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   ) 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i+2,j,k)-u(2,i-2,j,k)) +
     *        c1*(u(2,i+1,j,k)-u(2,i-1,j,k))  )*strx(i)*istry 
*** qr
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i,j+2,k)-u(1,i,j-2,k)) +
     *        c1*(u(1,i,j+1,k)-u(1,i,j-1,k))   ) 
     * + (2*mu(i,j,k)+la(i,j,k))*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*stry(j)*istrx 
     *  + mu(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*istrx -
     *                  forcing(2,i,j)

*** (w-eq)
            rhs3 = 
*** pr
     *    la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(1,i+2,j,k)-u(1,i-2,j,k)) +
     *        c1*(u(1,i+1,j,k)-u(1,i-1,j,k))   )*istry 
     *  + mu(i,j,k)*met(2,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i+2,j,k)-u(3,i-2,j,k)) +
     *        c1*(u(3,i+1,j,k)-u(3,i-1,j,k))  )*strx(i)*istry
*** qr 
     *  + mu(i,j,k)*met(3,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(3,i,j+2,k)-u(3,i,j-2,k)) +
     *        c1*(u(3,i,j+1,k)-u(3,i,j-1,k))   )*stry(j)*istrx
     *  + la(i,j,k)*met(4,i,j,k)*met(1,i,j,k)*(
     *        c2*(u(2,i,j+2,k)-u(2,i,j-2,k)) +
     *        c1*(u(2,i,j+1,k)-u(2,i,j-1,k))  )*istrx -
     *                  forcing(3,i,j)

*** Normal derivatives
            ac = strx(i)*istry*met(2,i,j,k)**2+
     *         stry(j)*istrx*met(3,i,j,k)**2+met(4,i,j,k)**2*istry*istrx
            bc = 1/(mu(i,j,k)*ac)
            cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac

            xoysqrt = SQRT(strx(i)*istry)
            yoxsqrt = 1/xoysqrt
            isqrtxy = istrx*xoysqrt
            dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
     *           yoxsqrt*met(3,i,j,k)*rhs2 + isqrtxy*met(4,i,j,k)*rhs3)

            u(1,i,j,k-kl) = -s0i*(  s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+
     *           s(3)*u(1,i,j,k+2*kl)+s(4)*u(1,i,j,k+3*kl) + bc*rhs1 - 
     *                                       dc*met(2,i,j,k)*xoysqrt )
            u(2,i,j,k-kl) = -s0i*(  s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+
     *           s(3)*u(2,i,j,k+2*kl)+s(4)*u(2,i,j,k+3*kl) + bc*rhs2 - 
     *                                       dc*met(3,i,j,k)*yoxsqrt )
            u(3,i,j,k-kl) = -s0i*(  s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+
     *           s(3)*u(3,i,j,k+2*kl)+s(4)*u(3,i,j,k+3*kl) + bc*rhs3 - 
     *                                       dc*met(4,i,j,k)*isqrtxy )
c            ac = strx(i)*istry*met(2,i,j,k)**2+
c     *                stry(j)*istrx*met(3,i,j,k)**2+met(4,i,j,k)**2
c            bc = 1/(mu(i,j,k)*ac)
c            cc = (mu(i,j,k)+la(i,j,k))/(2*mu(i,j,k)+la(i,j,k))*bc/ac
c
c            xoysqrt = SQRT(strx(i)*istry)
c            yoxsqrt = 1/xoysqrt
c
c            dc = cc*( xoysqrt*met(2,i,j,k)*rhs1 + 
c     *                yoxsqrt*met(3,i,j,k)*rhs2 + met(4,i,j,k)*rhs3)
c
c            u(1,i,j,k-kl) = -s0i*(  s(1)*u(1,i,j,k)+s(2)*u(1,i,j,k+kl)+
c     *           s(3)*u(1,i,j,k+2*kl)+s(4)*u(1,i,j,k+3*kl) + bc*rhs1 - 
c     *                                       dc*met(2,i,j,k)*xoysqrt )
c            u(2,i,j,k-kl) = -s0i*(  s(1)*u(2,i,j,k)+s(2)*u(2,i,j,k+kl)+
c     *           s(3)*u(2,i,j,k+2*kl)+s(4)*u(2,i,j,k+3*kl) + bc*rhs2 - 
c     *                                       dc*met(3,i,j,k)*yoxsqrt )
c            u(3,i,j,k-kl) = -s0i*(  s(1)*u(3,i,j,k)+s(2)*u(3,i,j,k+kl)+
c     *           s(3)*u(3,i,j,k+2*kl)+s(4)*u(3,i,j,k+3*kl) + bc*rhs3 - 
c     *                                       dc*met(4,i,j,k) )
         enddo
      enddo
      end
