      subroutine rhs4th3fortsgstr( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, nz, onesided, acof, bope, ghcof,
     +     uacc, u, mu, la, h, strx, stry, strz, op )

*** Routine with supergrid stretchings, strx, stry, and strz.
***
*** in the interior: centered approximation of the spatial operator in the elastic wave equation
*** near physical boundaries: one-sided approximation of the spatial operator in the elastic wave equation

      implicit none

      real*8 tf, i6, i144, i12
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144, i12=1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer nz, onesided(6), q, kb, mb, qb, k1, k2
      real*8 acof(6,8,8), bope(6,8), ghcof(6)
      real*8 uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
c      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 muz1, muz2, muz3, muz4
      real*8 mucof, mu1zz, mu2zz, lau2yz, lacof
      real*8 lap2mu, mu3zz, mu3xz, mu3yz, lau1xz
      real*8 lau3zx, u3zim2, u3zim1, u3zip1, u3zip2
      real*8 lau3zy, u3zjm2, u3zjm1, u3zjp1, u3zjp2
      real*8 mu1zx, u1zim2, u1zim1, u1zip1, u1zip2
      real*8 mu2zy, u2zjm2, u2zjm1, u2zjp1, u2zjp2
      real*8 r1, r2, r3, h, cof, d4a, d4b, a1
      real*8 strx(ifirst:ilast), stry(jfirst:jlast), strz(kfirst:klast)
      character*1 op
      parameter( d4a=2d0/3, d4b=-1d0/12 )

      cof = 1d0/(h*h)
      if( op.eq.'=' )then
         a1 = 0
      elseif( op.eq.'+')then
         a1 = 1
      elseif( op.eq.'-')then
         a1 = 1
         cof = -cof
      endif

      k1 = kfirst+2
      if (onesided(5).eq.1) k1 = 7
      k2 = klast-2
      if (onesided(6).eq.1) k2 = nz-6
c the centered stencil can be evaluated 2 points away from the boundary
      do k=k1,k2
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a, 28x3 = 84 ops
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                     tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                      3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                      3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                     tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
c
c            write(*,*) 'F mux 1-4 ',mux1, mux2, mux3, mux4
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                    tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                     3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                     3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                    tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))
c            write(*,*) 'F muy 1-4 ',muy1, muy2, muy3, muy4
c
            muz1 = mu(i,j,k-1)*strz(k-1)-
     *                    tf*(mu(i,j,k)*strz(k)+mu(i,j,k-2)*strz(k-2))
            muz2 = mu(i,j,k-2)*strz(k-2)+mu(i,j,k+1)*strz(k+1)+
     *                     3*(mu(i,j,k)*strz(k)+mu(i,j,k-1)*strz(k-1))
            muz3 = mu(i,j,k-1)*strz(k-1)+mu(i,j,k+2)*strz(k+2)+
     *                     3*(mu(i,j,k+1)*strz(k+1)+mu(i,j,k)*strz(k))
            muz4 = mu(i,j,k+1)*strz(k+1)-
     *                    tf*(mu(i,j,k)*strz(k)+mu(i,j,k+2)*strz(k+2))
c            write(*,*) 'F muz 1-4 ',muz1, muz2, muz3, muz4
*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda in the same way as we did for mu
*** 75 ops
            r1 = i6*( strx(i)*( (2*mux1+la(i-1,j,k)*strx(i-1)-
     *          tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *           3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *           3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *          tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *                muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) + strz(k)*(
     *                muz1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *                muz2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *                muz3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *                muz4*(u(1,i,j,k+2)-u(1,i,j,k)) ) )

*** 75 ops
            r2 = i6*( strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) ) + stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                 tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *               tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) + strz(k)*(
     *                muz1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *                muz2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *                muz3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *                muz4*(u(2,i,j,k+2)-u(2,i,j,k)) ) )

*** 75 ops
            r3 = i6*( strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) + strz(k)*(
     *             (2*muz1+la(i,j,k-1)*strz(k-1)-
     *                 tf*(la(i,j,k)*strz(k)+la(i,j,k-2)*strz(k-2)))*
     *                     (u(3,i,j,k-2)-u(3,i,j,k))+
     *      (2*muz2+la(i,j,k-2)*strz(k-2)+la(i,j,k+1)*strz(k+1)+
     *                 3*(la(i,j,k)*strz(k)+la(i,j,k-1)*strz(k-1)))*
     *                     (u(3,i,j,k-1)-u(3,i,j,k))+ 
     *      (2*muz3+la(i,j,k-1)*strz(k-1)+la(i,j,k+2)*strz(k+2)+
     *                 3*(la(i,j,k+1)*strz(k+1)+la(i,j,k)*strz(k)))*
     *                     (u(3,i,j,k+1)-u(3,i,j,k))+
     *             (2*muz4+la(i,j,k+1)*strz(k+1)-
     *               tf*(la(i,j,k)*strz(k)+la(i,j,k+2)*strz(k+2)))*
     *                     (u(3,i,j,k+2)-u(3,i,j,k)) ) )


*** Mixed derivatives:
*** 29ops /mixed derivative
*** 116 ops for r1
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (la*w_z)_x
     *          + strx(i)*strz(k)*       
     *            i144*( la(i-2,j,k)*(u(3,i-2,j,k-2)-u(3,i-2,j,k+2)+
     *                        8*(-u(3,i-2,j,k-1)+u(3,i-2,j,k+1))) - 8*(
     *                   la(i-1,j,k)*(u(3,i-1,j,k-2)-u(3,i-1,j,k+2)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i-1,j,k+1))) )+8*(
     *                   la(i+1,j,k)*(u(3,i+1,j,k-2)-u(3,i+1,j,k+2)+
     *                        8*(-u(3,i+1,j,k-1)+u(3,i+1,j,k+1))) ) - (
     *                   la(i+2,j,k)*(u(3,i+2,j,k-2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i+2,j,k-1)+u(3,i+2,j,k+1))) )) 
***   (mu*v_x)_y
     *          + strx(i)*stry(j)*       
     *            i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) 
***   (mu*w_x)_z
     *          + strx(i)*strz(k)*       
     *            i144*( mu(i,j,k-2)*(u(3,i-2,j,k-2)-u(3,i+2,j,k-2)+
     *                        8*(-u(3,i-1,j,k-2)+u(3,i+1,j,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i-2,j,k-1)-u(3,i+2,j,k-1)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i+1,j,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i-2,j,k+1)-u(3,i+2,j,k+1)+
     *                        8*(-u(3,i-1,j,k+1)+u(3,i+1,j,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i-2,j,k+2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i-1,j,k+2)+u(3,i+1,j,k+2))) )) 

*** 116 ops for r2
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *         + strx(i)*stry(j)*
     *            i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) 
*** (la*w_z)_y
     *          + stry(j)*strz(k)*
     *            i144*( la(i,j-2,k)*(u(3,i,j-2,k-2)-u(3,i,j-2,k+2)+
     *                        8*(-u(3,i,j-2,k-1)+u(3,i,j-2,k+1))) - 8*(
     *                   la(i,j-1,k)*(u(3,i,j-1,k-2)-u(3,i,j-1,k+2)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j-1,k+1))) )+8*(
     *                   la(i,j+1,k)*(u(3,i,j+1,k-2)-u(3,i,j+1,k+2)+
     *                        8*(-u(3,i,j+1,k-1)+u(3,i,j+1,k+1))) ) - (
     *                   la(i,j+2,k)*(u(3,i,j+2,k-2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j+2,k-1)+u(3,i,j+2,k+1))) ))
*** (mu*w_y)_z
     *          + stry(j)*strz(k)*
     *            i144*( mu(i,j,k-2)*(u(3,i,j-2,k-2)-u(3,i,j+2,k-2)+
     *                        8*(-u(3,i,j-1,k-2)+u(3,i,j+1,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i,j-2,k-1)-u(3,i,j+2,k-1)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j+1,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i,j-2,k+1)-u(3,i,j+2,k+1)+
     *                        8*(-u(3,i,j-1,k+1)+u(3,i,j+1,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i,j-2,k+2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j-1,k+2)+u(3,i,j+1,k+2))) )) 
*** 116 ops for r3
***  (mu*u_z)_x
            r3 = r3 + strx(i)*strz(k)*
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j,k-2)-u(1,i-2,j,k+2)+
     *                        8*(-u(1,i-2,j,k-1)+u(1,i-2,j,k+1))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j,k-2)-u(1,i-1,j,k+2)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i-1,j,k+1))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j,k-2)-u(1,i+1,j,k+2)+
     *                        8*(-u(1,i+1,j,k-1)+u(1,i+1,j,k+1))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j,k-2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i+2,j,k-1)+u(1,i+2,j,k+1))) )) 
*** (mu*v_z)_y
     *         + stry(j)*strz(k)*
     *            i144*( mu(i,j-2,k)*(u(2,i,j-2,k-2)-u(2,i,j-2,k+2)+
     *                        8*(-u(2,i,j-2,k-1)+u(2,i,j-2,k+1))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i,j-1,k-2)-u(2,i,j-1,k+2)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j-1,k+1))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i,j+1,k-2)-u(2,i,j+1,k+2)+
     *                        8*(-u(2,i,j+1,k-1)+u(2,i,j+1,k+1))) ) - (
     *                   mu(i,j+2,k)*(u(2,i,j+2,k-2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j+2,k-1)+u(2,i,j+2,k+1))) ))
***   (la*u_x)_z
     *         + strx(i)*strz(k)*
     *            i144*( la(i,j,k-2)*(u(1,i-2,j,k-2)-u(1,i+2,j,k-2)+
     *                        8*(-u(1,i-1,j,k-2)+u(1,i+1,j,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(1,i-2,j,k-1)-u(1,i+2,j,k-1)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i+1,j,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(1,i-2,j,k+1)-u(1,i+2,j,k+1)+
     *                        8*(-u(1,i-1,j,k+1)+u(1,i+1,j,k+1))) ) - (
     *                   la(i,j,k+2)*(u(1,i-2,j,k+2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i-1,j,k+2)+u(1,i+1,j,k+2))) )) 
*** (la*v_y)_z
     *         + stry(j)*strz(k)*
     *            i144*( la(i,j,k-2)*(u(2,i,j-2,k-2)-u(2,i,j+2,k-2)+
     *                        8*(-u(2,i,j-1,k-2)+u(2,i,j+1,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(2,i,j-2,k-1)-u(2,i,j+2,k-1)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j+1,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(2,i,j-2,k+1)-u(2,i,j+2,k+1)+
     *                        8*(-u(2,i,j-1,k+1)+u(2,i,j+1,k+1))) ) - (
     *                   la(i,j,k+2)*(u(2,i,j-2,k+2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j-1,k+2)+u(2,i,j+1,k+2))) )) 

*** 9 ops
c            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
c            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
c            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3
            uacc(1,i,j,k) = cof*r1
            uacc(2,i,j,k) = cof*r2
            uacc(3,i,j,k) = cof*r3

            enddo
         enddo
      enddo
*** Total number of ops away from boundary loop, per grid point computed:
*** 84+75*3+3*116+9 = 666 ops
*** Memory access mu,lambda,u = 5^3*(1+1+3) = 625
***                uacc = 3
***             strx,stry,strz = 5*3 = 15
***  Memory access read, 643 vars. = 5144 bytes (double prec.)
***                write, 3 vars. = 24 bytes (double prec.)

c low-k boundary modified stencils
      if (onesided(5).eq.1) then
      do k=1,6
c the centered stencil can be used in the x- and y-directions
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                  3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                  3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
c
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                  3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                  3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*(strx(i)*((2*mux1+la(i-1,j,k)*strx(i-1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *                   3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *                   3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) )
c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient, 
c leave out the z- supergrid stretching strz, since it will
c never be used together with the sbp-boundary operator
c computing the second derivative
            mu1zz = 0
            mu2zz = 0
            mu3zz = 0
            do q=1,8
c              mucof =0
c              lap2mu=0
c              do m=1,8
c                mucof = mucof+acof(k,q,m)*mu(i,j,m)
c                lap2mu= lap2mu+acof(k,q,m)*(2*mu(i,j,m)+la(i,j,m))
               mucof = acof(k,q,1)*mu(i,j,1)+acof(k,q,2)*mu(i,j,2)+
     *                 acof(k,q,3)*mu(i,j,3)+acof(k,q,4)*mu(i,j,4)+
     *                 acof(k,q,5)*mu(i,j,5)+acof(k,q,6)*mu(i,j,6)+
     *                 acof(k,q,7)*mu(i,j,7)+acof(k,q,8)*mu(i,j,8)
               lacof = acof(k,q,1)*la(i,j,1)+acof(k,q,2)*la(i,j,2)+
     *                 acof(k,q,3)*la(i,j,3)+acof(k,q,4)*la(i,j,4)+
     *                 acof(k,q,5)*la(i,j,5)+acof(k,q,6)*la(i,j,6)+
     *                 acof(k,q,7)*la(i,j,7)+acof(k,q,8)*la(i,j,8)
c              enddo
              mu1zz = mu1zz + mucof*u(1,i,j,q)
              mu2zz = mu2zz + mucof*u(2,i,j,q)
c              mu3zz = mu3zz + lap2mu*u(3,i,j,q)
              mu3zz = mu3zz + (lacof+2*mucof)*u(3,i,j,q)
            end do

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(k)*mu(i,j,1)*u(1,i,j,0))

            r2 = i6*(strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) )+ stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                   tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                   3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                   3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *                  tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) )

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(k)*mu(i,j,1)*u(2,i,j,0))

            r3 = i6*(strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) )

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(k)*(la(i,j,1)+2*mu(i,j,1))*
     +           u(3,i,j,0))

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*(
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) )
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do q=1,8
              u3zip2 = u3zip2 + bope(k,q)*u(3,i+2,j,q)
              u3zip1 = u3zip1 + bope(k,q)*u(3,i+1,j,q)
              u3zim1 = u3zim1 + bope(k,q)*u(3,i-1,j,q)
              u3zim2 = u3zim2 + bope(k,q)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + strx(i)*lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do q=1,8
              mu3xz = mu3xz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + strx(i)*mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*(
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) )
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do q=1,8
              u3zjp2 = u3zjp2 + bope(k,q)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 + bope(k,q)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 + bope(k,q)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 + bope(k,q)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + stry(j)*lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do q=1,8
              mu3yz = mu3yz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + stry(j)*mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do q=1,8
              u1zip2 = u1zip2 + bope(k,q)*u(1,i+2,j,q)
              u1zip1 = u1zip1 + bope(k,q)*u(1,i+1,j,q)
              u1zim1 = u1zim1 + bope(k,q)*u(1,i-1,j,q)
              u1zim2 = u1zim2 + bope(k,q)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + strx(i)*mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do q=1,8
              u2zjp2 = u2zjp2 + bope(k,q)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 + bope(k,q)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 + bope(k,q)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 + bope(k,q)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + stry(j)*mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do q=1,8
              lau1xz = lau1xz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + strx(i)*lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do q=1,8
              lau2yz = lau2yz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + stry(j)*lau2yz

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
      endif

c high-k boundary
      if (onesided(6).eq.1) then
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)*strx(i-1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i-2,j,k)*strx(i-2))
            mux2 = mu(i-2,j,k)*strx(i-2)+mu(i+1,j,k)*strx(i+1)+
     *                  3*(mu(i,j,k)*strx(i)+mu(i-1,j,k)*strx(i-1))
            mux3 = mu(i-1,j,k)*strx(i-1)+mu(i+2,j,k)*strx(i+2)+
     *                  3*(mu(i+1,j,k)*strx(i+1)+mu(i,j,k)*strx(i))
            mux4 = mu(i+1,j,k)*strx(i+1)-
     *                 tf*(mu(i,j,k)*strx(i)+mu(i+2,j,k)*strx(i+2))
c
            muy1 = mu(i,j-1,k)*stry(j-1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j-2,k)*stry(j-2))
            muy2 = mu(i,j-2,k)*stry(j-2)+mu(i,j+1,k)*stry(j+1)+
     *                  3*(mu(i,j,k)*stry(j)+mu(i,j-1,k)*stry(j-1))
            muy3 = mu(i,j-1,k)*stry(j-1)+mu(i,j+2,k)*stry(j+2)+
     *                  3*(mu(i,j+1,k)*stry(j+1)+mu(i,j,k)*stry(j))
            muy4 = mu(i,j+1,k)*stry(j+1)-
     *                 tf*(mu(i,j,k)*stry(j)+mu(i,j+2,k)*stry(j+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*(strx(i)*((2*mux1+la(i-1,j,k)*strx(i-1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i-2,j,k)*strx(i-2)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)*strx(i-2)+la(i+1,j,k)*strx(i+1)+
     *                   3*(la(i,j,k)*strx(i)+la(i-1,j,k)*strx(i-1)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)*strx(i-1)+la(i+2,j,k)*strx(i+2)+
     *                   3*(la(i+1,j,k)*strx(i+1)+la(i,j,k)*strx(i)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)*strx(i+1)-
     *                  tf*(la(i,j,k)*strx(i)+la(i+2,j,k)*strx(i+2)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) ) + stry(j)*(
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) ) )

c all indices ending with 'b' are indices relative to the boundary, going into the domain (1,2,3,...)
            kb = nz-k+1
c all coefficient arrays (acof, bope, ghcof) should be indexed with these indices
c all solution and material property arrays should be indexed with (i,j,k)

c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            mu1zz = 0
            mu2zz = 0
            mu3zz = 0
            do qb=1,8
              mucof=0
              lap2mu=0
              do mb=1,8
                mucof = mucof+acof(kb,qb,mb)*mu(i,j,nz-mb+1)
                lap2mu = lap2mu+acof(kb,qb,mb)*
     *                   (2*mu(i,j,nz-mb+1)+la(i,j,nz-mb+1))
              enddo
              mu1zz = mu1zz + mucof*u(1,i,j,nz-qb+1)
              mu2zz = mu2zz + mucof*u(2,i,j,nz-qb+1)
              mu3zz = mu3zz + lap2mu*u(3,i,j,nz-qb+1)
            end do

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(kb)*mu(i,j,nz)*u(1,i,j,nz+1))

            r2 = i6*(strx(i)*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) )+ stry(j)*(
     *             (2*muy1+la(i,j-1,k)*stry(j-1)-
     *                   tf*(la(i,j,k)*stry(j)+la(i,j-2,k)*stry(j-2)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)*stry(j-2)+la(i,j+1,k)*stry(j+1)+
     *                   3*(la(i,j,k)*stry(j)+la(i,j-1,k)*stry(j-1)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)*stry(j-1)+la(i,j+2,k)*stry(j+2)+
     *                   3*(la(i,j+1,k)*stry(j+1)+la(i,j,k)*stry(j)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)*stry(j+1)-
     *                  tf*(la(i,j,k)*stry(j)+la(i,j+2,k)*stry(j+2)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) ) )

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(kb)*mu(i,j,nz)*u(2,i,j,nz+1))

            r3 = i6*(strx(i)*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  ) + stry(j)*(
     *                muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) ) )

c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(kb)*(la(i,j,nz)+2*mu(i,j,nz))*
     +           u(3,i,j,nz+1))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + strx(i)*stry(j)*(
     *            i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) )
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do qb=1,8
              q = nz-qb+1
              u3zip2 = u3zip2 - bope(kb,qb)*u(3,i+2,j,q)
              u3zip1 = u3zip1 - bope(kb,qb)*u(3,i+1,j,q)
              u3zim1 = u3zim1 - bope(kb,qb)*u(3,i-1,j,q)
              u3zim2 = u3zim2 - bope(kb,qb)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + strx(i)*lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do qb=1,8
              q = nz-qb+1
              mu3xz = mu3xz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + strx(i)*mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + strx(i)*stry(j)*(
     *            i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) )
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do qb=1,8
              q = nz-qb+1
              u3zjp2 = u3zjp2 - bope(kb,qb)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 - bope(kb,qb)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 - bope(kb,qb)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 - bope(kb,qb)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + stry(j)*lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do qb=1,8
              q = nz-qb+1
              mu3yz = mu3yz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + stry(j)*mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do qb=1,8
              q = nz-qb+1
              u1zip2 = u1zip2 - bope(kb,qb)*u(1,i+2,j,q)
              u1zip1 = u1zip1 - bope(kb,qb)*u(1,i+1,j,q)
              u1zim1 = u1zim1 - bope(kb,qb)*u(1,i-1,j,q)
              u1zim2 = u1zim2 - bope(kb,qb)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + strx(i)*mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do qb=1,8
              q = nz-qb+1
              u2zjp2 = u2zjp2 - bope(kb,qb)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 - bope(kb,qb)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 - bope(kb,qb)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 - bope(kb,qb)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + stry(j)*mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do qb=1,8
              q = nz-qb+1
              lau1xz = lau1xz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + strx(i)*lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do qb=1,8
              q = nz-qb+1
              lau2yz = lau2yz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + stry(j)*lau2yz

            uacc(1,i,j,k) = a1*uacc(1,i,j,k) + cof*r1
            uacc(2,i,j,k) = a1*uacc(2,i,j,k) + cof*r2
            uacc(3,i,j,k) = a1*uacc(3,i,j,k) + cof*r3

            enddo
         enddo
      enddo
      endif

      end
