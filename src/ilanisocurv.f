      subroutine ILANISOCURV( ifirst, ilast, jfirst, jlast, kfirst, 
     *                        klast, nk, u, c, jac, lu, 
     *                        onesided, acof, bop, ghcof, 
     *                        strx, stry, strz ) bind(c)
      implicit none
      real*8 i6, a1, a2
      parameter( i6=1d0/6, a1=2d0/3, a2=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nk, i, j, k
      integer onesided(6), kb, ke, m, q
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 c(45,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 jac(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 strx(ifirst:ilast),stry(jfirst:jlast),strz(kfirst:klast)
      real*8 cm2, cm1, cp1, cp2, r1, r2, r3, dup2, dup1, dum1, dum2
      real*8 acof(6,8,8), ghcof(6), bop(6,8)
      real*8 ijac, ac1, ac2, ac3, ac4, ac5, ac6, du
      if( onesided(5).eq.1 )then
         kb=7
      else
         kb = kfirst+2
      endif
      if( onesided(6).eq.1 )then
         ke = nk-6
      else
         ke = klast-2
      endif
      if( onesided(5).eq.1 )then
         do k=1,6
            do j=jfirst+2,jlast-2
               do i=ifirst+2,ilast-2
         ijac = 1/jac(i,j,k)
         r1 = 0
         r2 = 0
         r3 = 0
      cm2 = c(1,i-1,j,k)*strx(i-1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i-2,j,k)*strx(i-2))
      cm1 = c(1,i-2,j,k)*strx(i-2)+c(1,i+1,j,k)*strx(i+1)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i-1,j,k)*strx(i-1))
      cp1 = c(1,i-1,j,k)*strx(i-1)+c(1,i+2,j,k)*strx(i+2)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i+1,j,k)*strx(i+1))
      cp2 = c(1,i+1,j,k)*strx(i+1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(2,i-1,j,k)*strx(i-1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i-2,j,k)*strx(i-2))
      cm1 = c(2,i-2,j,k)*strx(i-2)+c(2,i+1,j,k)*strx(i+1)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i-1,j,k)*strx(i-1))
      cp1 = c(2,i-1,j,k)*strx(i-1)+c(2,i+2,j,k)*strx(i+2)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i+1,j,k)*strx(i+1))
      cp2 = c(2,i+1,j,k)*strx(i+1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      r2 = r2+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(3,i-1,j,k)*strx(i-1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i-2,j,k)*strx(i-2))
      cm1 = c(3,i-2,j,k)*strx(i-2)+c(3,i+1,j,k)*strx(i+1)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i-1,j,k)*strx(i-1))
      cp1 = c(3,i-1,j,k)*strx(i-1)+c(3,i+2,j,k)*strx(i+2)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i+1,j,k)*strx(i+1))
      cp2 = c(3,i+1,j,k)*strx(i+1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(4,i-1,j,k)*strx(i-1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i-2,j,k)*strx(i-2))
      cm1 = c(4,i-2,j,k)*strx(i-2)+c(4,i+1,j,k)*strx(i+1)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i-1,j,k)*strx(i-1))
      cp1 = c(4,i-1,j,k)*strx(i-1)+c(4,i+2,j,k)*strx(i+2)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i+1,j,k)*strx(i+1))
      cp2 = c(4,i+1,j,k)*strx(i+1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(5,i-1,j,k)*strx(i-1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i-2,j,k)*strx(i-2))
      cm1 = c(5,i-2,j,k)*strx(i-2)+c(5,i+1,j,k)*strx(i+1)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i-1,j,k)*strx(i-1))
      cp1 = c(5,i-1,j,k)*strx(i-1)+c(5,i+2,j,k)*strx(i+2)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i+1,j,k)*strx(i+1))
      cp2 = c(5,i+1,j,k)*strx(i+1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(6,i-1,j,k)*strx(i-1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i-2,j,k)*strx(i-2))
      cm1 = c(6,i-2,j,k)*strx(i-2)+c(6,i+1,j,k)*strx(i+1)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i-1,j,k)*strx(i-1))
      cp1 = c(6,i-1,j,k)*strx(i-1)+c(6,i+2,j,k)*strx(i+2)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i+1,j,k)*strx(i+1))
      cp2 = c(6,i+1,j,k)*strx(i+1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i+2,j,k)*strx(i+2))

      r3 = r3+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      cm2 = c(7,i,j-1,k)*stry(j-1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j-2,k)*stry(j-2))
      cm1 = c(7,i,j-2,k)*stry(j-2)+c(7,i,j+1,k)*stry(j+1)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j-1,k)*stry(j-1))
      cp1 = c(7,i,j-1,k)*stry(j-1)+c(7,i,j+2,k)*stry(j+2)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j+1,k)*stry(j+1))
      cp2 = c(7,i,j+1,k)*stry(j+1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(8,i,j-1,k)*stry(j-1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j-2,k)*stry(j-2))
      cm1 = c(8,i,j-2,k)*stry(j-2)+c(8,i,j+1,k)*stry(j+1)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j-1,k)*stry(j-1))
      cp1 = c(8,i,j-1,k)*stry(j-1)+c(8,i,j+2,k)*stry(j+2)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j+1,k)*stry(j+1))
      cp2 = c(8,i,j+1,k)*stry(j+1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      r2 = r2+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(9,i,j-1,k)*stry(j-1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j-2,k)*stry(j-2))
      cm1 = c(9,i,j-2,k)*stry(j-2)+c(9,i,j+1,k)*stry(j+1)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j-1,k)*stry(j-1))
      cp1 = c(9,i,j-1,k)*stry(j-1)+c(9,i,j+2,k)*stry(j+2)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j+1,k)*stry(j+1))
      cp2 = c(9,i,j+1,k)*stry(j+1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(10,i,j-1,k)*stry(j-1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j-2,k)*stry(j-2))
      cm1 = c(10,i,j-2,k)*stry(j-2)+c(10,i,j+1,k)*stry(j+1)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j-1,k)*stry(j-1))
      cp1 = c(10,i,j-1,k)*stry(j-1)+c(10,i,j+2,k)*stry(j+2)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j+1,k)*stry(j+1))
      cp2 = c(10,i,j+1,k)*stry(j+1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(11,i,j-1,k)*stry(j-1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j-2,k)*stry(j-2))
      cm1 = c(11,i,j-2,k)*stry(j-2)+c(11,i,j+1,k)*stry(j+1)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j-1,k)*stry(j-1))
      cp1 = c(11,i,j-1,k)*stry(j-1)+c(11,i,j+2,k)*stry(j+2)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j+1,k)*stry(j+1))
      cp2 = c(11,i,j+1,k)*stry(j+1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(12,i,j-1,k)*stry(j-1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j-2,k)*stry(j-2))
      cm1 = c(12,i,j-2,k)*stry(j-2)+c(12,i,j+1,k)*stry(j+1)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j-1,k)*stry(j-1))
      cp1 = c(12,i,j-1,k)*stry(j-1)+c(12,i,j+2,k)*stry(j+2)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j+1,k)*stry(j+1))
      cp2 = c(12,i,j+1,k)*stry(j+1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j+2,k)*stry(j+2))

      r3 = r3+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r1 = r1 + ghcof(k)*c(13,i,j,1)*u(1,i,j,0) + 
     &ghcof(k)*c(14,i,j,1)*u(2,i,j,0) + 
     &ghcof(k)*c(15,i,j,1)*u(3,i,j,0)
      r2 = r2 + ghcof(k)*c(14,i,j,1)*u(1,i,j,0) + 
     &ghcof(k)*c(16,i,j,1)*u(2,i,j,0) + 
     &ghcof(k)*c(17,i,j,1)*u(3,i,j,0)
      r3 = r3 + ghcof(k)*c(15,i,j,1)*u(1,i,j,0) + 
     &ghcof(k)*c(17,i,j,1)*u(2,i,j,0) + 
     &ghcof(k)*c(18,i,j,1)*u(3,i,j,0)
      do q=1,8 
      ac1=0
      ac2=0
      ac3=0
      ac4=0
      ac5=0
      ac6=0
      do m=1,8 
      ac1=  ac1+ acof(k,q,m)*c(13,i,j,m) 
      ac2=  ac2+ acof(k,q,m)*c(14,i,j,m) 
      ac3=  ac3+ acof(k,q,m)*c(15,i,j,m) 
      ac4=  ac4+ acof(k,q,m)*c(16,i,j,m) 
      ac5=  ac5+ acof(k,q,m)*c(17,i,j,m) 
      ac6=  ac6+ acof(k,q,m)*c(18,i,j,m) 
      enddo
      r1 = r1 + ac1*u(1,i,j,q) + ac2*u(2,i,j,q) + ac3*u(3,i,j,q)
      r2 = r2 + ac2*u(1,i,j,q) + ac4*u(2,i,j,q) + ac5*u(3,i,j,q)
      r3 = r3 + ac3*u(1,i,j,q) + ac5*u(2,i,j,q) + ac6*u(3,i,j,q)
      enddo



      dup1 = a2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))
      dum1 = a2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k))+
     &a1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k))+
     &a1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))
      dum2 = a2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(19,i+1,j,k)*dup1-
     &c(19,i-1,j,k)*dum1)+a2*(c(19,i+2,j,k)*dup2-c(19,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(22,i+1,j,k)*dup1-
     &c(22,i-1,j,k)*dum1)+a2*(c(22,i+2,j,k)*dup2-c(22,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(25,i+1,j,k)*dup1-
     &c(25,i-1,j,k)*dum1)+a2*(c(25,i+2,j,k)*dup2-c(25,i-2,j,k)*dum2))


      dup1 = a2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))
      dum1 = a2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k))+
     &a1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k))+
     &a1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))
      dum2 = a2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(20,i+1,j,k)*dup1-
     &c(20,i-1,j,k)*dum1)+a2*(c(20,i+2,j,k)*dup2-c(20,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(23,i+1,j,k)*dup1-
     &c(23,i-1,j,k)*dum1)+a2*(c(23,i+2,j,k)*dup2-c(23,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(26,i+1,j,k)*dup1-
     &c(26,i-1,j,k)*dum1)+a2*(c(26,i+2,j,k)*dup2-c(26,i-2,j,k)*dum2))


      dup1 = a2*(u(3,i+1,j+2,k)-u(3,i+1,j-2,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i+1,j-1,k))
      dum1 = a2*(u(3,i-1,j+2,k)-u(3,i-1,j-2,k))+
     &a1*(u(3,i-1,j+1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i+2,j-2,k))+
     &a1*(u(3,i+2,j+1,k)-u(3,i+2,j-1,k))
      dum2 = a2*(u(3,i-2,j+2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i-2,j+1,k)-u(3,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(21,i+1,j,k)*dup1-
     &c(21,i-1,j,k)*dum1)+a2*(c(21,i+2,j,k)*dup2-c(21,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(24,i+1,j,k)*dup1-
     &c(24,i-1,j,k)*dum1)+a2*(c(24,i+2,j,k)*dup2-c(24,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(27,i+1,j,k)*dup1-
     &c(27,i-1,j,k)*dum1)+a2*(c(27,i+2,j,k)*dup2-c(27,i-2,j,k)*dum2))


      dup1 = a2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))
      dum1 = a2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k))+
     &a1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k))+
     &a1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))
      dum2 = a2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(19,i,j+1,k)*dup1-
     &c(19,i,j-1,k)*dum1)+a2*(c(19,i,j+2,k)*dup2-c(19,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(20,i,j+1,k)*dup1-
     &c(20,i,j-1,k)*dum1)+a2*(c(20,i,j+2,k)*dup2-c(20,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(21,i,j+1,k)*dup1-
     &c(21,i,j-1,k)*dum1)+a2*(c(21,i,j+2,k)*dup2-c(21,i,j-2,k)*dum2))


      dup1 = a2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))
      dum1 = a2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k))+
     &a1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k))+
     &a1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))
      dum2 = a2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(22,i,j+1,k)*dup1-
     &c(22,i,j-1,k)*dum1)+a2*(c(22,i,j+2,k)*dup2-c(22,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(23,i,j+1,k)*dup1-
     &c(23,i,j-1,k)*dum1)+a2*(c(23,i,j+2,k)*dup2-c(23,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(24,i,j+1,k)*dup1-
     &c(24,i,j-1,k)*dum1)+a2*(c(24,i,j+2,k)*dup2-c(24,i,j-2,k)*dum2))


      dup1 = a2*(u(3,i+2,j+1,k)-u(3,i-2,j+1,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i-1,j+1,k))
      dum1 = a2*(u(3,i+2,j-1,k)-u(3,i-2,j-1,k))+
     &a1*(u(3,i+1,j-1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i-2,j+2,k))+
     &a1*(u(3,i+1,j+2,k)-u(3,i-1,j+2,k))
      dum2 = a2*(u(3,i+2,j-2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i+1,j-2,k)-u(3,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(25,i,j+1,k)*dup1-
     &c(25,i,j-1,k)*dum1)+a2*(c(25,i,j+2,k)*dup2-c(25,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(26,i,j+1,k)*dup1-
     &c(26,i,j-1,k)*dum1)+a2*(c(26,i,j+2,k)*dup2-c(26,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(27,i,j+1,k)*dup1-
     &c(27,i,j-1,k)*dum1)+a2*(c(27,i,j+2,k)*dup2-c(27,i,j-2,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(1,i+2,j,q)
      dup1 = dup1 +bop(k,q)*u(1,i+1,j,q)
      dum1 = dum1 +bop(k,q)*u(1,i-1,j,q)
      dum2 = dum2 +bop(k,q)*u(1,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(28,i+1,j,k)*dup1-c(28,i-1,j,k)*dum1)
     &+ a2*(c(28,i+2,j,k)*dup2-c(28,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(31,i+1,j,k)*dup1-c(31,i-1,j,k)*dum1)
     &+ a2*(c(31,i+2,j,k)*dup2-c(31,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(34,i+1,j,k)*dup1-c(34,i-1,j,k)*dum1)
     &+ a2*(c(34,i+2,j,k)*dup2-c(34,i-2,j,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(2,i+2,j,q)
      dup1 = dup1 +bop(k,q)*u(2,i+1,j,q)
      dum1 = dum1 +bop(k,q)*u(2,i-1,j,q)
      dum2 = dum2 +bop(k,q)*u(2,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(29,i+1,j,k)*dup1-c(29,i-1,j,k)*dum1)
     &+ a2*(c(29,i+2,j,k)*dup2-c(29,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(32,i+1,j,k)*dup1-c(32,i-1,j,k)*dum1)
     &+ a2*(c(32,i+2,j,k)*dup2-c(32,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(35,i+1,j,k)*dup1-c(35,i-1,j,k)*dum1)
     &+ a2*(c(35,i+2,j,k)*dup2-c(35,i-2,j,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(3,i+2,j,q)
      dup1 = dup1 +bop(k,q)*u(3,i+1,j,q)
      dum1 = dum1 +bop(k,q)*u(3,i-1,j,q)
      dum2 = dum2 +bop(k,q)*u(3,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(30,i+1,j,k)*dup1-c(30,i-1,j,k)*dum1)
     &+ a2*(c(30,i+2,j,k)*dup2-c(30,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(33,i+1,j,k)*dup1-c(33,i-1,j,k)*dum1)
     &+ a2*(c(33,i+2,j,k)*dup2-c(33,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(36,i+1,j,k)*dup1-c(36,i-1,j,k)*dum1)
     &+ a2*(c(36,i+2,j,k)*dup2-c(36,i-2,j,k)*dum2))


      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(1,i+2,j,q)-u(1,i-2,j,q))+
     &a1*( u(1,i+1,j,q)-u(1,i-1,j,q))
         ac1 = ac1+bop(k,q)*c(28,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(29,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(30,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3



      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(2,i+2,j,q)-u(2,i-2,j,q))+
     &a1*( u(2,i+1,j,q)-u(2,i-1,j,q))
         ac1 = ac1+bop(k,q)*c(31,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(32,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(33,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3



      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(3,i+2,j,q)-u(3,i-2,j,q))+
     &a1*( u(3,i+1,j,q)-u(3,i-1,j,q))
         ac1 = ac1+bop(k,q)*c(34,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(35,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(36,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3



      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(1,i,j+2,q)
      dup1 = dup1 +bop(k,q)*u(1,i,j+1,q)
      dum1 = dum1 +bop(k,q)*u(1,i,j-1,q)
      dum2 = dum2 +bop(k,q)*u(1,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(37,i,j+1,k)*dup1-c(37,i,j-1,k)*dum1)
     &+ a2*(c(37,i,j+2,k)*dup2-c(37,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(40,i,j+1,k)*dup1-c(40,i,j-1,k)*dum1)
     &+ a2*(c(40,i,j+2,k)*dup2-c(40,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(43,i,j+1,k)*dup1-c(43,i,j-1,k)*dum1)
     &+ a2*(c(43,i,j+2,k)*dup2-c(43,i,j-2,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(2,i,j+2,q)
      dup1 = dup1 +bop(k,q)*u(2,i,j+1,q)
      dum1 = dum1 +bop(k,q)*u(2,i,j-1,q)
      dum2 = dum2 +bop(k,q)*u(2,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(38,i,j+1,k)*dup1-c(38,i,j-1,k)*dum1)
     &+ a2*(c(38,i,j+2,k)*dup2-c(38,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(41,i,j+1,k)*dup1-c(41,i,j-1,k)*dum1)
     &+ a2*(c(41,i,j+2,k)*dup2-c(41,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(44,i,j+1,k)*dup1-c(44,i,j-1,k)*dum1)
     &+ a2*(c(44,i,j+2,k)*dup2-c(44,i,j-2,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=1,8 
      dup2 = dup2 +bop(k,q)*u(3,i,j+2,q)
      dup1 = dup1 +bop(k,q)*u(3,i,j+1,q)
      dum1 = dum1 +bop(k,q)*u(3,i,j-1,q)
      dum2 = dum2 +bop(k,q)*u(3,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(39,i,j+1,k)*dup1-c(39,i,j-1,k)*dum1)
     &+ a2*(c(39,i,j+2,k)*dup2-c(39,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(42,i,j+1,k)*dup1-c(42,i,j-1,k)*dum1)
     &+ a2*(c(42,i,j+2,k)*dup2-c(42,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(45,i,j+1,k)*dup1-c(45,i,j-1,k)*dum1)
     &+ a2*(c(45,i,j+2,k)*dup2-c(45,i,j-2,k)*dum2))
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(1,i,j+2,q)-u(1,i,j-2,q))+
     &a1*( u(1,i,j+1,q)-u(1,i,j-1,q))
         ac1 = ac1+bop(k,q)*c(37,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(38,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(39,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(2,i,j+2,q)-u(2,i,j-2,q))+
     &a1*( u(2,i,j+1,q)-u(2,i,j-1,q))
         ac1 = ac1+bop(k,q)*c(40,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(41,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(42,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=1,8 
         du = a2*(u(3,i,j+2,q)-u(3,i,j-2,q))+
     &a1*( u(3,i,j+1,q)-u(3,i,j-1,q))
         ac1 = ac1+bop(k,q)*c(43,i,j,q)*du
         ac2 = ac2+bop(k,q)*c(44,i,j,q)*du
         ac3 = ac3+bop(k,q)*c(45,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3

      lu(1,i,j,k) = r1*ijac
      lu(2,i,j,k) = r2*ijac
      lu(3,i,j,k) = r3*ijac
               enddo
            enddo
         enddo
      endif
      do k=kb,ke
      do j=jfirst+2,jlast-2
      do i=ifirst+2,ilast-2
         ijac = 1/jac(i,j,k)
         r1 = 0
         r2 = 0
         r3 = 0
      cm2 = c(1,i-1,j,k)*strx(i-1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i-2,j,k)*strx(i-2))
      cm1 = c(1,i-2,j,k)*strx(i-2)+c(1,i+1,j,k)*strx(i+1)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i-1,j,k)*strx(i-1))
      cp1 = c(1,i-1,j,k)*strx(i-1)+c(1,i+2,j,k)*strx(i+2)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i+1,j,k)*strx(i+1))
      cp2 = c(1,i+1,j,k)*strx(i+1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(2,i-1,j,k)*strx(i-1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i-2,j,k)*strx(i-2))
      cm1 = c(2,i-2,j,k)*strx(i-2)+c(2,i+1,j,k)*strx(i+1)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i-1,j,k)*strx(i-1))
      cp1 = c(2,i-1,j,k)*strx(i-1)+c(2,i+2,j,k)*strx(i+2)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i+1,j,k)*strx(i+1))
      cp2 = c(2,i+1,j,k)*strx(i+1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      r2 = r2+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(3,i-1,j,k)*strx(i-1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i-2,j,k)*strx(i-2))
      cm1 = c(3,i-2,j,k)*strx(i-2)+c(3,i+1,j,k)*strx(i+1)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i-1,j,k)*strx(i-1))
      cp1 = c(3,i-1,j,k)*strx(i-1)+c(3,i+2,j,k)*strx(i+2)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i+1,j,k)*strx(i+1))
      cp2 = c(3,i+1,j,k)*strx(i+1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(4,i-1,j,k)*strx(i-1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i-2,j,k)*strx(i-2))
      cm1 = c(4,i-2,j,k)*strx(i-2)+c(4,i+1,j,k)*strx(i+1)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i-1,j,k)*strx(i-1))
      cp1 = c(4,i-1,j,k)*strx(i-1)+c(4,i+2,j,k)*strx(i+2)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i+1,j,k)*strx(i+1))
      cp2 = c(4,i+1,j,k)*strx(i+1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(5,i-1,j,k)*strx(i-1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i-2,j,k)*strx(i-2))
      cm1 = c(5,i-2,j,k)*strx(i-2)+c(5,i+1,j,k)*strx(i+1)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i-1,j,k)*strx(i-1))
      cp1 = c(5,i-1,j,k)*strx(i-1)+c(5,i+2,j,k)*strx(i+2)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i+1,j,k)*strx(i+1))
      cp2 = c(5,i+1,j,k)*strx(i+1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(6,i-1,j,k)*strx(i-1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i-2,j,k)*strx(i-2))
      cm1 = c(6,i-2,j,k)*strx(i-2)+c(6,i+1,j,k)*strx(i+1)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i-1,j,k)*strx(i-1))
      cp1 = c(6,i-1,j,k)*strx(i-1)+c(6,i+2,j,k)*strx(i+2)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i+1,j,k)*strx(i+1))
      cp2 = c(6,i+1,j,k)*strx(i+1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i+2,j,k)*strx(i+2))

      r3 = r3+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      cm2 = c(7,i,j-1,k)*stry(j-1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j-2,k)*stry(j-2))
      cm1 = c(7,i,j-2,k)*stry(j-2)+c(7,i,j+1,k)*stry(j+1)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j-1,k)*stry(j-1))
      cp1 = c(7,i,j-1,k)*stry(j-1)+c(7,i,j+2,k)*stry(j+2)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j+1,k)*stry(j+1))
      cp2 = c(7,i,j+1,k)*stry(j+1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(8,i,j-1,k)*stry(j-1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j-2,k)*stry(j-2))
      cm1 = c(8,i,j-2,k)*stry(j-2)+c(8,i,j+1,k)*stry(j+1)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j-1,k)*stry(j-1))
      cp1 = c(8,i,j-1,k)*stry(j-1)+c(8,i,j+2,k)*stry(j+2)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j+1,k)*stry(j+1))
      cp2 = c(8,i,j+1,k)*stry(j+1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      r2 = r2+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(9,i,j-1,k)*stry(j-1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j-2,k)*stry(j-2))
      cm1 = c(9,i,j-2,k)*stry(j-2)+c(9,i,j+1,k)*stry(j+1)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j-1,k)*stry(j-1))
      cp1 = c(9,i,j-1,k)*stry(j-1)+c(9,i,j+2,k)*stry(j+2)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j+1,k)*stry(j+1))
      cp2 = c(9,i,j+1,k)*stry(j+1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(10,i,j-1,k)*stry(j-1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j-2,k)*stry(j-2))
      cm1 = c(10,i,j-2,k)*stry(j-2)+c(10,i,j+1,k)*stry(j+1)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j-1,k)*stry(j-1))
      cp1 = c(10,i,j-1,k)*stry(j-1)+c(10,i,j+2,k)*stry(j+2)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j+1,k)*stry(j+1))
      cp2 = c(10,i,j+1,k)*stry(j+1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(11,i,j-1,k)*stry(j-1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j-2,k)*stry(j-2))
      cm1 = c(11,i,j-2,k)*stry(j-2)+c(11,i,j+1,k)*stry(j+1)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j-1,k)*stry(j-1))
      cp1 = c(11,i,j-1,k)*stry(j-1)+c(11,i,j+2,k)*stry(j+2)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j+1,k)*stry(j+1))
      cp2 = c(11,i,j+1,k)*stry(j+1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(12,i,j-1,k)*stry(j-1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j-2,k)*stry(j-2))
      cm1 = c(12,i,j-2,k)*stry(j-2)+c(12,i,j+1,k)*stry(j+1)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j-1,k)*stry(j-1))
      cp1 = c(12,i,j-1,k)*stry(j-1)+c(12,i,j+2,k)*stry(j+2)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j+1,k)*stry(j+1))
      cp2 = c(12,i,j+1,k)*stry(j+1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j+2,k)*stry(j+2))

      r3 = r3+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      cm2 = c(13,i,j,k-1)*strz(k-1)-0.75d0*(c(13,i,j,k)*strz(k)+
     &c(13,i,j,k-2)*strz(k-2))
      cm1 = c(13,i,j,k-2)*strz(k-2)+c(13,i,j,k+1)*strz(k+1)+3*(
     &c(13,i,j,k)*strz(k)+c(13,i,j,k-1)*strz(k-1))
      cp1 = c(13,i,j,k-1)*strz(k-1)+c(13,i,j,k+2)*strz(k+2)+3*(
     &c(13,i,j,k)*strz(k)+c(13,i,j,k+1)*strz(k+1))
      cp2 = c(13,i,j,k+1)*strz(k+1)-0.75d0*(c(13,i,j,k)*strz(k)+
     &c(13,i,j,k+2)*strz(k+2))

      r1 = r1+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
      cm2 = c(14,i,j,k-1)*strz(k-1)-0.75d0*(c(14,i,j,k)*strz(k)+
     &c(14,i,j,k-2)*strz(k-2))
      cm1 = c(14,i,j,k-2)*strz(k-2)+c(14,i,j,k+1)*strz(k+1)+3*(
     &c(14,i,j,k)*strz(k)+c(14,i,j,k-1)*strz(k-1))
      cp1 = c(14,i,j,k-1)*strz(k-1)+c(14,i,j,k+2)*strz(k+2)+3*(
     &c(14,i,j,k)*strz(k)+c(14,i,j,k+1)*strz(k+1))
      cp2 = c(14,i,j,k+1)*strz(k+1)-0.75d0*(c(14,i,j,k)*strz(k)+
     &c(14,i,j,k+2)*strz(k+2))

      r1 = r1+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
      r2 = r2+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
      cm2 = c(15,i,j,k-1)*strz(k-1)-0.75d0*(c(15,i,j,k)*strz(k)+
     &c(15,i,j,k-2)*strz(k-2))
      cm1 = c(15,i,j,k-2)*strz(k-2)+c(15,i,j,k+1)*strz(k+1)+3*(
     &c(15,i,j,k)*strz(k)+c(15,i,j,k-1)*strz(k-1))
      cp1 = c(15,i,j,k-1)*strz(k-1)+c(15,i,j,k+2)*strz(k+2)+3*(
     &c(15,i,j,k)*strz(k)+c(15,i,j,k+1)*strz(k+1))
      cp2 = c(15,i,j,k+1)*strz(k+1)-0.75d0*(c(15,i,j,k)*strz(k)+
     &c(15,i,j,k+2)*strz(k+2))

      r1 = r1+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )
      r3 = r3+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
      cm2 = c(16,i,j,k-1)*strz(k-1)-0.75d0*(c(16,i,j,k)*strz(k)+
     &c(16,i,j,k-2)*strz(k-2))
      cm1 = c(16,i,j,k-2)*strz(k-2)+c(16,i,j,k+1)*strz(k+1)+3*(
     &c(16,i,j,k)*strz(k)+c(16,i,j,k-1)*strz(k-1))
      cp1 = c(16,i,j,k-1)*strz(k-1)+c(16,i,j,k+2)*strz(k+2)+3*(
     &c(16,i,j,k)*strz(k)+c(16,i,j,k+1)*strz(k+1))
      cp2 = c(16,i,j,k+1)*strz(k+1)-0.75d0*(c(16,i,j,k)*strz(k)+
     &c(16,i,j,k+2)*strz(k+2))

      r2 = r2+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
      cm2 = c(17,i,j,k-1)*strz(k-1)-0.75d0*(c(17,i,j,k)*strz(k)+
     &c(17,i,j,k-2)*strz(k-2))
      cm1 = c(17,i,j,k-2)*strz(k-2)+c(17,i,j,k+1)*strz(k+1)+3*(
     &c(17,i,j,k)*strz(k)+c(17,i,j,k-1)*strz(k-1))
      cp1 = c(17,i,j,k-1)*strz(k-1)+c(17,i,j,k+2)*strz(k+2)+3*(
     &c(17,i,j,k)*strz(k)+c(17,i,j,k+1)*strz(k+1))
      cp2 = c(17,i,j,k+1)*strz(k+1)-0.75d0*(c(17,i,j,k)*strz(k)+
     &c(17,i,j,k+2)*strz(k+2))

      r2 = r2+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )
      r3 = r3+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
      cm2 = c(18,i,j,k-1)*strz(k-1)-0.75d0*(c(18,i,j,k)*strz(k)+
     &c(18,i,j,k-2)*strz(k-2))
      cm1 = c(18,i,j,k-2)*strz(k-2)+c(18,i,j,k+1)*strz(k+1)+3*(
     &c(18,i,j,k)*strz(k)+c(18,i,j,k-1)*strz(k-1))
      cp1 = c(18,i,j,k-1)*strz(k-1)+c(18,i,j,k+2)*strz(k+2)+3*(
     &c(18,i,j,k)*strz(k)+c(18,i,j,k+1)*strz(k+1))
      cp2 = c(18,i,j,k+1)*strz(k+1)-0.75d0*(c(18,i,j,k)*strz(k)+
     &c(18,i,j,k+2)*strz(k+2))

      r3 = r3+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )


      dup1 = a2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))
      dum1 = a2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k))+
     &a1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k))+
     &a1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))
      dum2 = a2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(19,i+1,j,k)*dup1-
     &c(19,i-1,j,k)*dum1)+a2*(c(19,i+2,j,k)*dup2-c(19,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(22,i+1,j,k)*dup1-
     &c(22,i-1,j,k)*dum1)+a2*(c(22,i+2,j,k)*dup2-c(22,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(25,i+1,j,k)*dup1-
     &c(25,i-1,j,k)*dum1)+a2*(c(25,i+2,j,k)*dup2-c(25,i-2,j,k)*dum2))


      dup1 = a2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))
      dum1 = a2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k))+
     &a1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k))+
     &a1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))
      dum2 = a2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(20,i+1,j,k)*dup1-
     &c(20,i-1,j,k)*dum1)+a2*(c(20,i+2,j,k)*dup2-c(20,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(23,i+1,j,k)*dup1-
     &c(23,i-1,j,k)*dum1)+a2*(c(23,i+2,j,k)*dup2-c(23,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(26,i+1,j,k)*dup1-
     &c(26,i-1,j,k)*dum1)+a2*(c(26,i+2,j,k)*dup2-c(26,i-2,j,k)*dum2))


      dup1 = a2*(u(3,i+1,j+2,k)-u(3,i+1,j-2,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i+1,j-1,k))
      dum1 = a2*(u(3,i-1,j+2,k)-u(3,i-1,j-2,k))+
     &a1*(u(3,i-1,j+1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i+2,j-2,k))+
     &a1*(u(3,i+2,j+1,k)-u(3,i+2,j-1,k))
      dum2 = a2*(u(3,i-2,j+2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i-2,j+1,k)-u(3,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(21,i+1,j,k)*dup1-
     &c(21,i-1,j,k)*dum1)+a2*(c(21,i+2,j,k)*dup2-c(21,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(24,i+1,j,k)*dup1-
     &c(24,i-1,j,k)*dum1)+a2*(c(24,i+2,j,k)*dup2-c(24,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(27,i+1,j,k)*dup1-
     &c(27,i-1,j,k)*dum1)+a2*(c(27,i+2,j,k)*dup2-c(27,i-2,j,k)*dum2))


      dup1 = a2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))
      dum1 = a2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k))+
     &a1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k))+
     &a1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))
      dum2 = a2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(19,i,j+1,k)*dup1-
     &c(19,i,j-1,k)*dum1)+a2*(c(19,i,j+2,k)*dup2-c(19,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(20,i,j+1,k)*dup1-
     &c(20,i,j-1,k)*dum1)+a2*(c(20,i,j+2,k)*dup2-c(20,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(21,i,j+1,k)*dup1-
     &c(21,i,j-1,k)*dum1)+a2*(c(21,i,j+2,k)*dup2-c(21,i,j-2,k)*dum2))


      dup1 = a2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))
      dum1 = a2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k))+
     &a1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k))+
     &a1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))
      dum2 = a2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(22,i,j+1,k)*dup1-
     &c(22,i,j-1,k)*dum1)+a2*(c(22,i,j+2,k)*dup2-c(22,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(23,i,j+1,k)*dup1-
     &c(23,i,j-1,k)*dum1)+a2*(c(23,i,j+2,k)*dup2-c(23,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(24,i,j+1,k)*dup1-
     &c(24,i,j-1,k)*dum1)+a2*(c(24,i,j+2,k)*dup2-c(24,i,j-2,k)*dum2))


      dup1 = a2*(u(3,i+2,j+1,k)-u(3,i-2,j+1,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i-1,j+1,k))
      dum1 = a2*(u(3,i+2,j-1,k)-u(3,i-2,j-1,k))+
     &a1*(u(3,i+1,j-1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i-2,j+2,k))+
     &a1*(u(3,i+1,j+2,k)-u(3,i-1,j+2,k))
      dum2 = a2*(u(3,i+2,j-2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i+1,j-2,k)-u(3,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(25,i,j+1,k)*dup1-
     &c(25,i,j-1,k)*dum1)+a2*(c(25,i,j+2,k)*dup2-c(25,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(26,i,j+1,k)*dup1-
     &c(26,i,j-1,k)*dum1)+a2*(c(26,i,j+2,k)*dup2-c(26,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(27,i,j+1,k)*dup1-
     &c(27,i,j-1,k)*dum1)+a2*(c(27,i,j+2,k)*dup2-c(27,i,j-2,k)*dum2))


      dup1 = a2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2))+
     &a1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1))
      dum1 = a2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2))+
     &a1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1))
      dup2 = a2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2))+
     &a1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))
      dum2 = a2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2))+
     &a1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))
      r1 = r1 +strz(k)*strx(i)*(a1*(c(28,i+1,j,k)*dup1-
     &c(28,i-1,j,k)*dum1)+a2*(c(28,i+2,j,k)*dup2-c(28,i-2,j,k)*dum2))
      r2 = r2 +strz(k)*strx(i)*(a1*(c(31,i+1,j,k)*dup1-
     &c(31,i-1,j,k)*dum1)+a2*(c(31,i+2,j,k)*dup2-c(31,i-2,j,k)*dum2))
      r3 = r3 +strz(k)*strx(i)*(a1*(c(34,i+1,j,k)*dup1-
     &c(34,i-1,j,k)*dum1)+a2*(c(34,i+2,j,k)*dup2-c(34,i-2,j,k)*dum2))


      dup1 = a2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2))+
     &a1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1))
      dum1 = a2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2))+
     &a1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1))
      dup2 = a2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2))+
     &a1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))
      dum2 = a2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2))+
     &a1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))
      r1 = r1 +strz(k)*strx(i)*(a1*(c(29,i+1,j,k)*dup1-
     &c(29,i-1,j,k)*dum1)+a2*(c(29,i+2,j,k)*dup2-c(29,i-2,j,k)*dum2))
      r2 = r2 +strz(k)*strx(i)*(a1*(c(32,i+1,j,k)*dup1-
     &c(32,i-1,j,k)*dum1)+a2*(c(32,i+2,j,k)*dup2-c(32,i-2,j,k)*dum2))
      r3 = r3 +strz(k)*strx(i)*(a1*(c(35,i+1,j,k)*dup1-
     &c(35,i-1,j,k)*dum1)+a2*(c(35,i+2,j,k)*dup2-c(35,i-2,j,k)*dum2))


      dup1 = a2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2))+
     &a1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1))
      dum1 = a2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2))+
     &a1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1))
      dup2 = a2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2))+
     &a1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1))
      dum2 = a2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2))+
     &a1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1))
      r1 = r1 +strz(k)*strx(i)*(a1*(c(30,i+1,j,k)*dup1-
     &c(30,i-1,j,k)*dum1)+a2*(c(30,i+2,j,k)*dup2-c(30,i-2,j,k)*dum2))
      r2 = r2 +strz(k)*strx(i)*(a1*(c(33,i+1,j,k)*dup1-
     &c(33,i-1,j,k)*dum1)+a2*(c(33,i+2,j,k)*dup2-c(33,i-2,j,k)*dum2))
      r3 = r3 +strz(k)*strx(i)*(a1*(c(36,i+1,j,k)*dup1-
     &c(36,i-1,j,k)*dum1)+a2*(c(36,i+2,j,k)*dup2-c(36,i-2,j,k)*dum2))


      dup1 = a2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1))+
     &a1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1))
      dum1 = a2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1))+
     &a1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1))
      dup2 = a2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2))+
     &a1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))
      dum2 = a2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2))+
     &a1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))
      r1 = r1 +strx(i)*strz(k)*(a1*(c(28,i,j,k+1)*dup1-
     &c(28,i,j,k-1)*dum1)+a2*(c(28,i,j,k+2)*dup2-c(28,i,j,k-2)*dum2))
      r2 = r2 +strx(i)*strz(k)*(a1*(c(29,i,j,k+1)*dup1-
     &c(29,i,j,k-1)*dum1)+a2*(c(29,i,j,k+2)*dup2-c(29,i,j,k-2)*dum2))
      r3 = r3 +strx(i)*strz(k)*(a1*(c(30,i,j,k+1)*dup1-
     &c(30,i,j,k-1)*dum1)+a2*(c(30,i,j,k+2)*dup2-c(30,i,j,k-2)*dum2))


      dup1 = a2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1))+
     &a1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1))
      dum1 = a2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1))+
     &a1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1))
      dup2 = a2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2))+
     &a1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))
      dum2 = a2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2))+
     &a1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2))
      r1 = r1 +strx(i)*strz(k)*(a1*(c(31,i,j,k+1)*dup1-
     &c(31,i,j,k-1)*dum1)+a2*(c(31,i,j,k+2)*dup2-c(31,i,j,k-2)*dum2))
      r2 = r2 +strx(i)*strz(k)*(a1*(c(32,i,j,k+1)*dup1-
     &c(32,i,j,k-1)*dum1)+a2*(c(32,i,j,k+2)*dup2-c(32,i,j,k-2)*dum2))
      r3 = r3 +strx(i)*strz(k)*(a1*(c(33,i,j,k+1)*dup1-
     &c(33,i,j,k-1)*dum1)+a2*(c(33,i,j,k+2)*dup2-c(33,i,j,k-2)*dum2))


      dup1 = a2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1))+
     &a1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1))
      dum1 = a2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1))+
     &a1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1))
      dup2 = a2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2))+
     &a1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))
      dum2 = a2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2))+
     &a1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2))
      r1 = r1 +strx(i)*strz(k)*(a1*(c(34,i,j,k+1)*dup1-
     &c(34,i,j,k-1)*dum1)+a2*(c(34,i,j,k+2)*dup2-c(34,i,j,k-2)*dum2))
      r2 = r2 +strx(i)*strz(k)*(a1*(c(35,i,j,k+1)*dup1-
     &c(35,i,j,k-1)*dum1)+a2*(c(35,i,j,k+2)*dup2-c(35,i,j,k-2)*dum2))
      r3 = r3 +strx(i)*strz(k)*(a1*(c(36,i,j,k+1)*dup1-
     &c(36,i,j,k-1)*dum1)+a2*(c(36,i,j,k+2)*dup2-c(36,i,j,k-2)*dum2))


      dup1 = a2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2))+
     &a1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1))
      dum1 = a2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2))+
     &a1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1))
      dup2 = a2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2))+
     &a1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))
      dum2 = a2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2))+
     &a1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))
      r1 = r1 +strz(k)*stry(j)*(a1*(c(37,i,j+1,k)*dup1-
     &c(37,i,j-1,k)*dum1)+a2*(c(37,i,j+2,k)*dup2-c(37,i,j-2,k)*dum2))
      r2 = r2 +strz(k)*stry(j)*(a1*(c(40,i,j+1,k)*dup1-
     &c(40,i,j-1,k)*dum1)+a2*(c(40,i,j+2,k)*dup2-c(40,i,j-2,k)*dum2))
      r3 = r3 +strz(k)*stry(j)*(a1*(c(43,i,j+1,k)*dup1-
     &c(43,i,j-1,k)*dum1)+a2*(c(43,i,j+2,k)*dup2-c(43,i,j-2,k)*dum2))


      dup1 = a2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2))+
     &a1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1))
      dum1 = a2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2))+
     &a1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1))
      dup2 = a2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2))+
     &a1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))
      dum2 = a2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2))+
     &a1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))
      r1 = r1 +strz(k)*stry(j)*(a1*(c(38,i,j+1,k)*dup1-
     &c(38,i,j-1,k)*dum1)+a2*(c(38,i,j+2,k)*dup2-c(38,i,j-2,k)*dum2))
      r2 = r2 +strz(k)*stry(j)*(a1*(c(41,i,j+1,k)*dup1-
     &c(41,i,j-1,k)*dum1)+a2*(c(41,i,j+2,k)*dup2-c(41,i,j-2,k)*dum2))
      r3 = r3 +strz(k)*stry(j)*(a1*(c(44,i,j+1,k)*dup1-
     &c(44,i,j-1,k)*dum1)+a2*(c(44,i,j+2,k)*dup2-c(44,i,j-2,k)*dum2))


      dup1 = a2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2))+
     &a1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1))
      dum1 = a2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2))+
     &a1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1))
      dup2 = a2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2))+
     &a1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))
      dum2 = a2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2))+
     &a1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))
      r1 = r1 +strz(k)*stry(j)*(a1*(c(39,i,j+1,k)*dup1-
     &c(39,i,j-1,k)*dum1)+a2*(c(39,i,j+2,k)*dup2-c(39,i,j-2,k)*dum2))
      r2 = r2 +strz(k)*stry(j)*(a1*(c(42,i,j+1,k)*dup1-
     &c(42,i,j-1,k)*dum1)+a2*(c(42,i,j+2,k)*dup2-c(42,i,j-2,k)*dum2))
      r3 = r3 +strz(k)*stry(j)*(a1*(c(45,i,j+1,k)*dup1-
     &c(45,i,j-1,k)*dum1)+a2*(c(45,i,j+2,k)*dup2-c(45,i,j-2,k)*dum2))


      dup1 = a2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1))+
     &a1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1))
      dum1 = a2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1))+
     &a1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1))
      dup2 = a2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2))+
     &a1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))
      dum2 = a2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2))+
     &a1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))
      r1 = r1 +stry(j)*strz(k)*(a1*(c(37,i,j,k+1)*dup1-
     &c(37,i,j,k-1)*dum1)+a2*(c(37,i,j,k+2)*dup2-c(37,i,j,k-2)*dum2))
      r2 = r2 +stry(j)*strz(k)*(a1*(c(38,i,j,k+1)*dup1-
     &c(38,i,j,k-1)*dum1)+a2*(c(38,i,j,k+2)*dup2-c(38,i,j,k-2)*dum2))
      r3 = r3 +stry(j)*strz(k)*(a1*(c(39,i,j,k+1)*dup1-
     &c(39,i,j,k-1)*dum1)+a2*(c(39,i,j,k+2)*dup2-c(39,i,j,k-2)*dum2))


      dup1 = a2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1))+
     &a1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1))
      dum1 = a2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1))+
     &a1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1))
      dup2 = a2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2))+
     &a1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))
      dum2 = a2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2))+
     &a1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))
      r1 = r1 +stry(j)*strz(k)*(a1*(c(40,i,j,k+1)*dup1-
     &c(40,i,j,k-1)*dum1)+a2*(c(40,i,j,k+2)*dup2-c(40,i,j,k-2)*dum2))
      r2 = r2 +stry(j)*strz(k)*(a1*(c(41,i,j,k+1)*dup1-
     &c(41,i,j,k-1)*dum1)+a2*(c(41,i,j,k+2)*dup2-c(41,i,j,k-2)*dum2))
      r3 = r3 +stry(j)*strz(k)*(a1*(c(42,i,j,k+1)*dup1-
     &c(42,i,j,k-1)*dum1)+a2*(c(42,i,j,k+2)*dup2-c(42,i,j,k-2)*dum2))


      dup1 = a2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1))+
     &a1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1))
      dum1 = a2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1))+
     &a1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1))
      dup2 = a2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2))+
     &a1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))
      dum2 = a2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2))+
     &a1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))
      r1 = r1 +stry(j)*strz(k)*(a1*(c(43,i,j,k+1)*dup1-
     &c(43,i,j,k-1)*dum1)+a2*(c(43,i,j,k+2)*dup2-c(43,i,j,k-2)*dum2))
      r2 = r2 +stry(j)*strz(k)*(a1*(c(44,i,j,k+1)*dup1-
     &c(44,i,j,k-1)*dum1)+a2*(c(44,i,j,k+2)*dup2-c(44,i,j,k-2)*dum2))
      r3 = r3 +stry(j)*strz(k)*(a1*(c(45,i,j,k+1)*dup1-
     &c(45,i,j,k-1)*dum1)+a2*(c(45,i,j,k+2)*dup2-c(45,i,j,k-2)*dum2))


c$$$
c$$$      cm2 = c(1,i-1,j,k)*strx(i-1)-0.75d0*(c(1,i,j,k)*strx(i)+
c$$$     &c(1,i-2,j,k)*strx(i-2))
c$$$      cm1 = c(1,i-2,j,k)*strx(i-2)+c(1,i+1,j,k)*strx(i+1)+3*(
c$$$     &c(1,i,j,k)*strx(i)+c(1,i-1,j,k)*strx(i-1))
c$$$      cp1 = c(1,i-1,j,k)*strx(i-1)+c(1,i+2,j,k)*strx(i+2)+3*(
c$$$     &c(1,i,j,k)*strx(i)+c(1,i+1,j,k)*strx(i+1))
c$$$      cp2 = c(1,i+1,j,k)*strx(i+1)-0.75d0*(c(1,i,j,k)*strx(i)+
c$$$     &c(1,i+2,j,k)*strx(i+2))
c$$$
c$$$      r1 = r1+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
c$$$      cm2 = c(2,i-1,j,k)*strx(i-1)-0.75d0*(c(2,i,j,k)*strx(i)+
c$$$     &c(2,i-2,j,k)*strx(i-2))
c$$$      cm1 = c(2,i-2,j,k)*strx(i-2)+c(2,i+1,j,k)*strx(i+1)+3*(
c$$$     &c(2,i,j,k)*strx(i)+c(2,i-1,j,k)*strx(i-1))
c$$$      cp1 = c(2,i-1,j,k)*strx(i-1)+c(2,i+2,j,k)*strx(i+2)+3*(
c$$$     &c(2,i,j,k)*strx(i)+c(2,i+1,j,k)*strx(i+1))
c$$$      cp2 = c(2,i+1,j,k)*strx(i+1)-0.75d0*(c(2,i,j,k)*strx(i)+
c$$$     &c(2,i+2,j,k)*strx(i+2))
c$$$
c$$$      r1 = r1+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
c$$$      r2 = r2+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
c$$$      cm2 = c(3,i-1,j,k)*strx(i-1)-0.75d0*(c(3,i,j,k)*strx(i)+
c$$$     &c(3,i-2,j,k)*strx(i-2))
c$$$      cm1 = c(3,i-2,j,k)*strx(i-2)+c(3,i+1,j,k)*strx(i+1)+3*(
c$$$     &c(3,i,j,k)*strx(i)+c(3,i-1,j,k)*strx(i-1))
c$$$      cp1 = c(3,i-1,j,k)*strx(i-1)+c(3,i+2,j,k)*strx(i+2)+3*(
c$$$     &c(3,i,j,k)*strx(i)+c(3,i+1,j,k)*strx(i+1))
c$$$      cp2 = c(3,i+1,j,k)*strx(i+1)-0.75d0*(c(3,i,j,k)*strx(i)+
c$$$     &c(3,i+2,j,k)*strx(i+2))
c$$$
c$$$      r1 = r1+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
c$$$      cm2 = c(4,i-1,j,k)*strx(i-1)-0.75d0*(c(4,i,j,k)*strx(i)+
c$$$     &c(4,i-2,j,k)*strx(i-2))
c$$$      cm1 = c(4,i-2,j,k)*strx(i-2)+c(4,i+1,j,k)*strx(i+1)+3*(
c$$$     &c(4,i,j,k)*strx(i)+c(4,i-1,j,k)*strx(i-1))
c$$$      cp1 = c(4,i-1,j,k)*strx(i-1)+c(4,i+2,j,k)*strx(i+2)+3*(
c$$$     &c(4,i,j,k)*strx(i)+c(4,i+1,j,k)*strx(i+1))
c$$$      cp2 = c(4,i+1,j,k)*strx(i+1)-0.75d0*(c(4,i,j,k)*strx(i)+
c$$$     &c(4,i+2,j,k)*strx(i+2))
c$$$
c$$$      r2 = r2+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
c$$$      cm2 = c(5,i-1,j,k)*strx(i-1)-0.75d0*(c(5,i,j,k)*strx(i)+
c$$$     &c(5,i-2,j,k)*strx(i-2))
c$$$      cm1 = c(5,i-2,j,k)*strx(i-2)+c(5,i+1,j,k)*strx(i+1)+3*(
c$$$     &c(5,i,j,k)*strx(i)+c(5,i-1,j,k)*strx(i-1))
c$$$      cp1 = c(5,i-1,j,k)*strx(i-1)+c(5,i+2,j,k)*strx(i+2)+3*(
c$$$     &c(5,i,j,k)*strx(i)+c(5,i+1,j,k)*strx(i+1))
c$$$      cp2 = c(5,i+1,j,k)*strx(i+1)-0.75d0*(c(5,i,j,k)*strx(i)+
c$$$     &c(5,i+2,j,k)*strx(i+2))
c$$$
c$$$      r2 = r2+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
c$$$      cm2 = c(6,i,j,k-1)*strz(k-1)-0.75d0*(c(6,i,j,k)*strz(k)+
c$$$     &c(6,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(6,i,j,k-2)*strz(k-2)+c(6,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(6,i,j,k)*strz(k)+c(6,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(6,i,j,k-1)*strz(k-1)+c(6,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(6,i,j,k)*strz(k)+c(6,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(6,i,j,k+1)*strz(k+1)-0.75d0*(c(6,i,j,k)*strz(k)+
c$$$     &c(6,i,j,k+2)*strz(k+2))
c$$$
c$$$      r3 = r3+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
c$$$      cm2 = c(7,i,j-1,k)*stry(j-1)-0.75d0*(c(7,i,j,k)*stry(j)+
c$$$     &c(7,i,j-2,k)*stry(j-2))
c$$$      cm1 = c(7,i,j-2,k)*stry(j-2)+c(7,i,j+1,k)*stry(j+1)+3*(
c$$$     &c(7,i,j,k)*stry(j)+c(7,i,j-1,k)*stry(j-1))
c$$$      cp1 = c(7,i,j-1,k)*stry(j-1)+c(7,i,j+2,k)*stry(j+2)+3*(
c$$$     &c(7,i,j,k)*stry(j)+c(7,i,j+1,k)*stry(j+1))
c$$$      cp2 = c(7,i,j+1,k)*stry(j+1)-0.75d0*(c(7,i,j,k)*stry(j)+
c$$$     &c(7,i,j+2,k)*stry(j+2))
c$$$
c$$$      r1 = r1+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
c$$$      cm2 = c(8,i,j-1,k)*stry(j-1)-0.75d0*(c(8,i,j,k)*stry(j)+
c$$$     &c(8,i,j-2,k)*stry(j-2))
c$$$      cm1 = c(8,i,j-2,k)*stry(j-2)+c(8,i,j+1,k)*stry(j+1)+3*(
c$$$     &c(8,i,j,k)*stry(j)+c(8,i,j-1,k)*stry(j-1))
c$$$      cp1 = c(8,i,j-1,k)*stry(j-1)+c(8,i,j+2,k)*stry(j+2)+3*(
c$$$     &c(8,i,j,k)*stry(j)+c(8,i,j+1,k)*stry(j+1))
c$$$      cp2 = c(8,i,j+1,k)*stry(j+1)-0.75d0*(c(8,i,j,k)*stry(j)+
c$$$     &c(8,i,j+2,k)*stry(j+2))
c$$$
c$$$      r1 = r1+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
c$$$      r2 = r2+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
c$$$      cm2 = c(9,i,j-1,k)*stry(j-1)-0.75d0*(c(9,i,j,k)*stry(j)+
c$$$     &c(9,i,j-2,k)*stry(j-2))
c$$$      cm1 = c(9,i,j-2,k)*stry(j-2)+c(9,i,j+1,k)*stry(j+1)+3*(
c$$$     &c(9,i,j,k)*stry(j)+c(9,i,j-1,k)*stry(j-1))
c$$$      cp1 = c(9,i,j-1,k)*stry(j-1)+c(9,i,j+2,k)*stry(j+2)+3*(
c$$$     &c(9,i,j,k)*stry(j)+c(9,i,j+1,k)*stry(j+1))
c$$$      cp2 = c(9,i,j+1,k)*stry(j+1)-0.75d0*(c(9,i,j,k)*stry(j)+
c$$$     &c(9,i,j+2,k)*stry(j+2))
c$$$
c$$$      r1 = r1+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
c$$$      cm2 = c(10,i,j-1,k)*stry(j-1)-0.75d0*(c(10,i,j,k)*stry(j)+
c$$$     &c(10,i,j-2,k)*stry(j-2))
c$$$      cm1 = c(10,i,j-2,k)*stry(j-2)+c(10,i,j+1,k)*stry(j+1)+3*(
c$$$     &c(10,i,j,k)*stry(j)+c(10,i,j-1,k)*stry(j-1))
c$$$      cp1 = c(10,i,j-1,k)*stry(j-1)+c(10,i,j+2,k)*stry(j+2)+3*(
c$$$     &c(10,i,j,k)*stry(j)+c(10,i,j+1,k)*stry(j+1))
c$$$      cp2 = c(10,i,j+1,k)*stry(j+1)-0.75d0*(c(10,i,j,k)*stry(j)+
c$$$     &c(10,i,j+2,k)*stry(j+2))
c$$$
c$$$      r2 = r2+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
c$$$      cm2 = c(11,i,j-1,k)*stry(j-1)-0.75d0*(c(11,i,j,k)*stry(j)+
c$$$     &c(11,i,j-2,k)*stry(j-2))
c$$$      cm1 = c(11,i,j-2,k)*stry(j-2)+c(11,i,j+1,k)*stry(j+1)+3*(
c$$$     &c(11,i,j,k)*stry(j)+c(11,i,j-1,k)*stry(j-1))
c$$$      cp1 = c(11,i,j-1,k)*stry(j-1)+c(11,i,j+2,k)*stry(j+2)+3*(
c$$$     &c(11,i,j,k)*stry(j)+c(11,i,j+1,k)*stry(j+1))
c$$$      cp2 = c(11,i,j+1,k)*stry(j+1)-0.75d0*(c(11,i,j,k)*stry(j)+
c$$$     &c(11,i,j+2,k)*stry(j+2))
c$$$
c$$$      r2 = r2+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
c$$$      cm2 = c(12,i,j,k-1)*strz(k-1)-0.75d0*(c(12,i,j,k)*strz(k)+
c$$$     &c(12,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(12,i,j,k-2)*strz(k-2)+c(12,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(12,i,j,k)*strz(k)+c(12,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(12,i,j,k-1)*strz(k-1)+c(12,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(12,i,j,k)*strz(k)+c(12,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(12,i,j,k+1)*strz(k+1)-0.75d0*(c(12,i,j,k)*strz(k)+
c$$$     &c(12,i,j,k+2)*strz(k+2))
c$$$
c$$$      r3 = r3+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
c$$$      cm2 = c(13,i,j,k-1)*strz(k-1)-0.75d0*(c(13,i,j,k)*strz(k)+
c$$$     &c(13,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(13,i,j,k-2)*strz(k-2)+c(13,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(13,i,j,k)*strz(k)+c(13,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(13,i,j,k-1)*strz(k-1)+c(13,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(13,i,j,k)*strz(k)+c(13,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(13,i,j,k+1)*strz(k+1)-0.75d0*(c(13,i,j,k)*strz(k)+
c$$$     &c(13,i,j,k+2)*strz(k+2))
c$$$
c$$$      r1 = r1+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
c$$$      cm2 = c(14,i,j,k-1)*strz(k-1)-0.75d0*(c(14,i,j,k)*strz(k)+
c$$$     &c(14,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(14,i,j,k-2)*strz(k-2)+c(14,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(14,i,j,k)*strz(k)+c(14,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(14,i,j,k-1)*strz(k-1)+c(14,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(14,i,j,k)*strz(k)+c(14,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(14,i,j,k+1)*strz(k+1)-0.75d0*(c(14,i,j,k)*strz(k)+
c$$$     &c(14,i,j,k+2)*strz(k+2))
c$$$
c$$$      r1 = r1+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
c$$$      r2 = r2+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
c$$$      cm2 = c(15,i,j,k-1)*strz(k-1)-0.75d0*(c(15,i,j,k)*strz(k)+
c$$$     &c(15,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(15,i,j,k-2)*strz(k-2)+c(15,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(15,i,j,k)*strz(k)+c(15,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(15,i,j,k-1)*strz(k-1)+c(15,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(15,i,j,k)*strz(k)+c(15,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(15,i,j,k+1)*strz(k+1)-0.75d0*(c(15,i,j,k)*strz(k)+
c$$$     &c(15,i,j,k+2)*strz(k+2))
c$$$
c$$$      r1 = r1+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*strz(k)*(cm2*(u(1,i,j,k-2)-u(1,i,j,k)) + 
c$$$     &cm1*(u(1,i,j,k-1)-u(1,i,j,k)) + 
c$$$     &cp1*(u(1,i,j,k+1)-u(1,i,j,k)) + 
c$$$     &cp2*(u(1,i,j,k+2)-u(1,i,j,k)) )
c$$$      cm2 = c(16,i,j,k-1)*strz(k-1)-0.75d0*(c(16,i,j,k)*strz(k)+
c$$$     &c(16,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(16,i,j,k-2)*strz(k-2)+c(16,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(16,i,j,k)*strz(k)+c(16,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(16,i,j,k-1)*strz(k-1)+c(16,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(16,i,j,k)*strz(k)+c(16,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(16,i,j,k+1)*strz(k+1)-0.75d0*(c(16,i,j,k)*strz(k)+
c$$$     &c(16,i,j,k+2)*strz(k+2))
c$$$
c$$$      r2 = r2+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
c$$$      cm2 = c(17,i,j,k-1)*strz(k-1)-0.75d0*(c(17,i,j,k)*strz(k)+
c$$$     &c(17,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(17,i,j,k-2)*strz(k-2)+c(17,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(17,i,j,k)*strz(k)+c(17,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(17,i,j,k-1)*strz(k-1)+c(17,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(17,i,j,k)*strz(k)+c(17,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(17,i,j,k+1)*strz(k+1)-0.75d0*(c(17,i,j,k)*strz(k)+
c$$$     &c(17,i,j,k+2)*strz(k+2))
c$$$
c$$$      r2 = r2+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )
c$$$      r3 = r3+i6*strz(k)*(cm2*(u(2,i,j,k-2)-u(2,i,j,k)) + 
c$$$     &cm1*(u(2,i,j,k-1)-u(2,i,j,k)) + 
c$$$     &cp1*(u(2,i,j,k+1)-u(2,i,j,k)) + 
c$$$     &cp2*(u(2,i,j,k+2)-u(2,i,j,k)) )
c$$$      cm2 = c(18,i,j,k-1)*strz(k-1)-0.75d0*(c(18,i,j,k)*strz(k)+
c$$$     &c(18,i,j,k-2)*strz(k-2))
c$$$      cm1 = c(18,i,j,k-2)*strz(k-2)+c(18,i,j,k+1)*strz(k+1)+3*(
c$$$     &c(18,i,j,k)*strz(k)+c(18,i,j,k-1)*strz(k-1))
c$$$      cp1 = c(18,i,j,k-1)*strz(k-1)+c(18,i,j,k+2)*strz(k+2)+3*(
c$$$     &c(18,i,j,k)*strz(k)+c(18,i,j,k+1)*strz(k+1))
c$$$      cp2 = c(18,i,j,k+1)*strz(k+1)-0.75d0*(c(18,i,j,k)*strz(k)+
c$$$     &c(18,i,j,k+2)*strz(k+2))
c$$$
c$$$      r3 = r3+i6*strz(k)*(cm2*(u(3,i,j,k-2)-u(3,i,j,k)) + 
c$$$     &cm1*(u(3,i,j,k-1)-u(3,i,j,k)) + 
c$$$     &cp1*(u(3,i,j,k+1)-u(3,i,j,k)) + 
c$$$     &cp2*(u(3,i,j,k+2)-u(3,i,j,k)) )
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k))+
c$$$     &a1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))
c$$$      dum1 = a2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k))+
c$$$     &a1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))
c$$$      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k))+
c$$$     &a1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))
c$$$      dum2 = a2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
c$$$     &a1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))
c$$$      r1 =r1+stry(j)*strx(i)*(a1*(c(19,i+1,j,k)*dup1-c(19,i-1,j,k)*dum1)
c$$$     &+a2*(c(19,i+2,j,k)*dup2-c(19,i-2,j,k)*dum2))
c$$$      r2 =r2+stry(j)*strx(i)*(a1*(c(22,i+1,j,k)*dup1-c(22,i-1,j,k)*dum1)
c$$$     &+a2*(c(22,i+2,j,k)*dup2-c(22,i-2,j,k)*dum2))
c$$$      r3 =r3+stry(j)*strx(i)*(a1*(c(25,i+1,j,k)*dup1-c(25,i-1,j,k)*dum1)
c$$$     &+a2*(c(25,i+2,j,k)*dup2-c(25,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k))+
c$$$     &a1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))
c$$$      dum1 = a2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k))+
c$$$     &a1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))
c$$$      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k))+
c$$$     &a1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))
c$$$      dum2 = a2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
c$$$     &a1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))
c$$$      r1 =r1+stry(j)*strx(i)*(a1*(c(20,i+1,j,k)*dup1-c(20,i-1,j,k)*dum1)
c$$$     &+a2*(c(20,i+2,j,k)*dup2-c(20,i-2,j,k)*dum2))
c$$$      r2 =r2+stry(j)*strx(i)*(a1*(c(23,i+1,j,k)*dup1-c(23,i-1,j,k)*dum1)
c$$$     &+a2*(c(23,i+2,j,k)*dup2-c(23,i-2,j,k)*dum2))
c$$$      r3 =r3+stry(j)*strx(i)*(a1*(c(26,i+1,j,k)*dup1-c(26,i-1,j,k)*dum1)
c$$$     &+a2*(c(26,i+2,j,k)*dup2-c(26,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i+1,j+2,k)-u(3,i+1,j-2,k))+
c$$$     &a1*(u(3,i+1,j+1,k)-u(3,i+1,j-1,k))
c$$$      dum1 = a2*(u(3,i-1,j+2,k)-u(3,i-1,j-2,k))+
c$$$     &a1*(u(3,i-1,j+1,k)-u(3,i-1,j-1,k))
c$$$      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i+2,j-2,k))+
c$$$     &a1*(u(3,i+2,j+1,k)-u(3,i+2,j-1,k))
c$$$      dum2 = a2*(u(3,i-2,j+2,k)-u(3,i-2,j-2,k))+
c$$$     &a1*(u(3,i-2,j+1,k)-u(3,i-2,j-1,k))
c$$$      r1 =r1+stry(j)*strx(i)*(a1*(c(21,i+1,j,k)*dup1-c(21,i-1,j,k)*dum1)
c$$$     &+a2*(c(21,i+2,j,k)*dup2-c(21,i-2,j,k)*dum2))
c$$$      r2 =r2+stry(j)*strx(i)*(a1*(c(24,i+1,j,k)*dup1-c(24,i-1,j,k)*dum1)
c$$$     &+a2*(c(24,i+2,j,k)*dup2-c(24,i-2,j,k)*dum2))
c$$$      r3 =r3+stry(j)*strx(i)*(a1*(c(27,i+1,j,k)*dup1-c(27,i-1,j,k)*dum1)
c$$$     &+a2*(c(27,i+2,j,k)*dup2-c(27,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k))+
c$$$     &a1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))
c$$$      dum1 = a2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k))+
c$$$     &a1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))
c$$$      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k))+
c$$$     &a1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))
c$$$      dum2 = a2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
c$$$     &a1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))
c$$$      r1 =r1+strx(i)*stry(j)*(a1*(c(19,i,j+1,k)*dup1-c(19,i,j-1,k)*dum1)
c$$$     &+a2*(c(19,i,j+2,k)*dup2-c(19,i,j-2,k)*dum2))
c$$$      r2 =r2+strx(i)*stry(j)*(a1*(c(20,i,j+1,k)*dup1-c(20,i,j-1,k)*dum1)
c$$$     &+a2*(c(20,i,j+2,k)*dup2-c(20,i,j-2,k)*dum2))
c$$$      r3 =r3+strx(i)*stry(j)*(a1*(c(21,i,j+1,k)*dup1-c(21,i,j-1,k)*dum1)
c$$$     &+a2*(c(21,i,j+2,k)*dup2-c(21,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k))+
c$$$     &a1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))
c$$$      dum1 = a2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k))+
c$$$     &a1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))
c$$$      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k))+
c$$$     &a1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))
c$$$      dum2 = a2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
c$$$     &a1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))
c$$$      r1 =r1+strx(i)*stry(j)*(a1*(c(22,i,j+1,k)*dup1-c(22,i,j-1,k)*dum1)
c$$$     &+a2*(c(22,i,j+2,k)*dup2-c(22,i,j-2,k)*dum2))
c$$$      r2 =r2+strx(i)*stry(j)*(a1*(c(23,i,j+1,k)*dup1-c(23,i,j-1,k)*dum1)
c$$$     &+a2*(c(23,i,j+2,k)*dup2-c(23,i,j-2,k)*dum2))
c$$$      r3 =r3+strx(i)*stry(j)*(a1*(c(24,i,j+1,k)*dup1-c(24,i,j-1,k)*dum1)
c$$$     &+a2*(c(24,i,j+2,k)*dup2-c(24,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i+2,j+1,k)-u(3,i-2,j+1,k))+
c$$$     &a1*(u(3,i+1,j+1,k)-u(3,i-1,j+1,k))
c$$$      dum1 = a2*(u(3,i+2,j-1,k)-u(3,i-2,j-1,k))+
c$$$     &a1*(u(3,i+1,j-1,k)-u(3,i-1,j-1,k))
c$$$      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i-2,j+2,k))+
c$$$     &a1*(u(3,i+1,j+2,k)-u(3,i-1,j+2,k))
c$$$      dum2 = a2*(u(3,i+2,j-2,k)-u(3,i-2,j-2,k))+
c$$$     &a1*(u(3,i+1,j-2,k)-u(3,i-1,j-2,k))
c$$$      r1 =r1+strx(i)*stry(j)*(a1*(c(25,i,j+1,k)*dup1-c(25,i,j-1,k)*dum1)
c$$$     &+a2*(c(25,i,j+2,k)*dup2-c(25,i,j-2,k)*dum2))
c$$$      r2 =r2+strx(i)*stry(j)*(a1*(c(26,i,j+1,k)*dup1-c(26,i,j-1,k)*dum1)
c$$$     &+a2*(c(26,i,j+2,k)*dup2-c(26,i,j-2,k)*dum2))
c$$$      r3 =r3+strx(i)*stry(j)*(a1*(c(27,i,j+1,k)*dup1-c(27,i,j-1,k)*dum1)
c$$$     &+a2*(c(27,i,j+2,k)*dup2-c(27,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i+1,j,k+2)-u(1,i+1,j,k-2))+
c$$$     &a1*(u(1,i+1,j,k+1)-u(1,i+1,j,k-1))
c$$$      dum1 = a2*(u(1,i-1,j,k+2)-u(1,i-1,j,k-2))+
c$$$     &a1*(u(1,i-1,j,k+1)-u(1,i-1,j,k-1))
c$$$      dup2 = a2*(u(1,i+2,j,k+2)-u(1,i+2,j,k-2))+
c$$$     &a1*(u(1,i+2,j,k+1)-u(1,i+2,j,k-1))
c$$$      dum2 = a2*(u(1,i-2,j,k+2)-u(1,i-2,j,k-2))+
c$$$     &a1*(u(1,i-2,j,k+1)-u(1,i-2,j,k-1))
c$$$      r1 =r1+strz(k)*strx(i)*(a1*(c(28,i+1,j,k)*dup1-c(28,i-1,j,k)*dum1)
c$$$     &+a2*(c(28,i+2,j,k)*dup2-c(28,i-2,j,k)*dum2))
c$$$      r2 =r2+strz(k)*strx(i)*(a1*(c(31,i+1,j,k)*dup1-c(31,i-1,j,k)*dum1)
c$$$     &+a2*(c(31,i+2,j,k)*dup2-c(31,i-2,j,k)*dum2))
c$$$      r3 =r3+strz(k)*strx(i)*(a1*(c(34,i+1,j,k)*dup1-c(34,i-1,j,k)*dum1)
c$$$     &+a2*(c(34,i+2,j,k)*dup2-c(34,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i+1,j,k+2)-u(2,i+1,j,k-2))+
c$$$     &a1*(u(2,i+1,j,k+1)-u(2,i+1,j,k-1))
c$$$      dum1 = a2*(u(2,i-1,j,k+2)-u(2,i-1,j,k-2))+
c$$$     &a1*(u(2,i-1,j,k+1)-u(2,i-1,j,k-1))
c$$$      dup2 = a2*(u(2,i+2,j,k+2)-u(2,i+2,j,k-2))+
c$$$     &a1*(u(2,i+2,j,k+1)-u(2,i+2,j,k-1))
c$$$      dum2 = a2*(u(2,i-2,j,k+2)-u(2,i-2,j,k-2))+
c$$$     &a1*(u(2,i-2,j,k+1)-u(2,i-2,j,k-1))
c$$$      r1 =r1+strz(k)*strx(i)*(a1*(c(29,i+1,j,k)*dup1-c(29,i-1,j,k)*dum1)
c$$$     &+a2*(c(29,i+2,j,k)*dup2-c(29,i-2,j,k)*dum2))
c$$$      r2 =r2+strz(k)*strx(i)*(a1*(c(32,i+1,j,k)*dup1-c(32,i-1,j,k)*dum1)
c$$$     &+a2*(c(32,i+2,j,k)*dup2-c(32,i-2,j,k)*dum2))
c$$$      r3 =r3+strz(k)*strx(i)*(a1*(c(35,i+1,j,k)*dup1-c(35,i-1,j,k)*dum1)
c$$$     &+a2*(c(35,i+2,j,k)*dup2-c(35,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i+1,j,k+2)-u(3,i+1,j,k-2))+
c$$$     &a1*(u(3,i+1,j,k+1)-u(3,i+1,j,k-1))
c$$$      dum1 = a2*(u(3,i-1,j,k+2)-u(3,i-1,j,k-2))+
c$$$     &a1*(u(3,i-1,j,k+1)-u(3,i-1,j,k-1))
c$$$      dup2 = a2*(u(3,i+2,j,k+2)-u(3,i+2,j,k-2))+
c$$$     &a1*(u(3,i+2,j,k+1)-u(3,i+2,j,k-1))
c$$$      dum2 = a2*(u(3,i-2,j,k+2)-u(3,i-2,j,k-2))+
c$$$     &a1*(u(3,i-2,j,k+1)-u(3,i-2,j,k-1))
c$$$      r1 =r1+strz(k)*strx(i)*(a1*(c(30,i+1,j,k)*dup1-c(30,i-1,j,k)*dum1)
c$$$     &+a2*(c(30,i+2,j,k)*dup2-c(30,i-2,j,k)*dum2))
c$$$      r2 =r2+strz(k)*strx(i)*(a1*(c(33,i+1,j,k)*dup1-c(33,i-1,j,k)*dum1)
c$$$     &+a2*(c(33,i+2,j,k)*dup2-c(33,i-2,j,k)*dum2))
c$$$      r3 =r3+strz(k)*strx(i)*(a1*(c(36,i+1,j,k)*dup1-c(36,i-1,j,k)*dum1)
c$$$     &+a2*(c(36,i+2,j,k)*dup2-c(36,i-2,j,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i+2,j,k+1)-u(1,i-2,j,k+1))+
c$$$     &a1*(u(1,i+1,j,k+1)-u(1,i-1,j,k+1))
c$$$      dum1 = a2*(u(1,i+2,j,k-1)-u(1,i-2,j,k-1))+
c$$$     &a1*(u(1,i+1,j,k-1)-u(1,i-1,j,k-1))
c$$$      dup2 = a2*(u(1,i+2,j,k+2)-u(1,i-2,j,k+2))+
c$$$     &a1*(u(1,i+1,j,k+2)-u(1,i-1,j,k+2))
c$$$      dum2 = a2*(u(1,i+2,j,k-2)-u(1,i-2,j,k-2))+
c$$$     &a1*(u(1,i+1,j,k-2)-u(1,i-1,j,k-2))
c$$$      r1 =r1+strx(i)*strz(k)*(a1*(c(28,i,j,k+1)*dup1-c(28,i,j,k-1)*dum1)
c$$$     &+a2*(c(28,i,j,k+2)*dup2-c(28,i,j,k-2)*dum2))
c$$$      r2 =r2+strx(i)*strz(k)*(a1*(c(29,i,j,k+1)*dup1-c(29,i,j,k-1)*dum1)
c$$$     &+a2*(c(29,i,j,k+2)*dup2-c(29,i,j,k-2)*dum2))
c$$$      r3 =r3+strx(i)*strz(k)*(a1*(c(30,i,j,k+1)*dup1-c(30,i,j,k-1)*dum1)
c$$$     &+a2*(c(30,i,j,k+2)*dup2-c(30,i,j,k-2)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i+2,j,k+1)-u(2,i-2,j,k+1))+
c$$$     &a1*(u(2,i+1,j,k+1)-u(2,i-1,j,k+1))
c$$$      dum1 = a2*(u(2,i+2,j,k-1)-u(2,i-2,j,k-1))+
c$$$     &a1*(u(2,i+1,j,k-1)-u(2,i-1,j,k-1))
c$$$      dup2 = a2*(u(2,i+2,j,k+2)-u(2,i-2,j,k+2))+
c$$$     &a1*(u(2,i+1,j,k+2)-u(2,i-1,j,k+2))
c$$$      dum2 = a2*(u(2,i+2,j,k-2)-u(2,i-2,j,k-2))+
c$$$     &a1*(u(2,i+1,j,k-2)-u(2,i-1,j,k-2))
c$$$      r1 =r1+strx(i)*strz(k)*(a1*(c(31,i,j,k+1)*dup1-c(31,i,j,k-1)*dum1)
c$$$     &+a2*(c(31,i,j,k+2)*dup2-c(31,i,j,k-2)*dum2))
c$$$      r2 =r2+strx(i)*strz(k)*(a1*(c(32,i,j,k+1)*dup1-c(32,i,j,k-1)*dum1)
c$$$     &+a2*(c(32,i,j,k+2)*dup2-c(32,i,j,k-2)*dum2))
c$$$      r3 =r3+strx(i)*strz(k)*(a1*(c(33,i,j,k+1)*dup1-c(33,i,j,k-1)*dum1)
c$$$     &+a2*(c(33,i,j,k+2)*dup2-c(33,i,j,k-2)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i+2,j,k+1)-u(3,i-2,j,k+1))+
c$$$     &a1*(u(3,i+1,j,k+1)-u(3,i-1,j,k+1))
c$$$      dum1 = a2*(u(3,i+2,j,k-1)-u(3,i-2,j,k-1))+
c$$$     &a1*(u(3,i+1,j,k-1)-u(3,i-1,j,k-1))
c$$$      dup2 = a2*(u(3,i+2,j,k+2)-u(3,i-2,j,k+2))+
c$$$     &a1*(u(3,i+1,j,k+2)-u(3,i-1,j,k+2))
c$$$      dum2 = a2*(u(3,i+2,j,k-2)-u(3,i-2,j,k-2))+
c$$$     &a1*(u(3,i+1,j,k-2)-u(3,i-1,j,k-2))
c$$$      r1 =r1+strx(i)*strz(k)*(a1*(c(34,i,j,k+1)*dup1-c(34,i,j,k-1)*dum1)
c$$$     &+a2*(c(34,i,j,k+2)*dup2-c(34,i,j,k-2)*dum2))
c$$$      r2 =r2+strx(i)*strz(k)*(a1*(c(35,i,j,k+1)*dup1-c(35,i,j,k-1)*dum1)
c$$$     &+a2*(c(35,i,j,k+2)*dup2-c(35,i,j,k-2)*dum2))
c$$$      r3 =r3+strx(i)*strz(k)*(a1*(c(36,i,j,k+1)*dup1-c(36,i,j,k-1)*dum1)
c$$$     &+a2*(c(36,i,j,k+2)*dup2-c(36,i,j,k-2)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i,j+1,k+2)-u(1,i,j+1,k-2))+
c$$$     &a1*(u(1,i,j+1,k+1)-u(1,i,j+1,k-1))
c$$$      dum1 = a2*(u(1,i,j-1,k+2)-u(1,i,j-1,k-2))+
c$$$     &a1*(u(1,i,j-1,k+1)-u(1,i,j-1,k-1))
c$$$      dup2 = a2*(u(1,i,j+2,k+2)-u(1,i,j+2,k-2))+
c$$$     &a1*(u(1,i,j+2,k+1)-u(1,i,j+2,k-1))
c$$$      dum2 = a2*(u(1,i,j-2,k+2)-u(1,i,j-2,k-2))+
c$$$     &a1*(u(1,i,j-2,k+1)-u(1,i,j-2,k-1))
c$$$      r1 =r1+strz(k)*stry(j)*(a1*(c(37,i,j+1,k)*dup1-c(37,i,j-1,k)*dum1)
c$$$     &+a2*(c(37,i,j+2,k)*dup2-c(37,i,j-2,k)*dum2))
c$$$      r2 =r2+strz(k)*stry(j)*(a1*(c(40,i,j+1,k)*dup1-c(40,i,j-1,k)*dum1)
c$$$     &+a2*(c(40,i,j+2,k)*dup2-c(40,i,j-2,k)*dum2))
c$$$      r3 =r3+strz(k)*stry(j)*(a1*(c(43,i,j+1,k)*dup1-c(43,i,j-1,k)*dum1)
c$$$     &+a2*(c(43,i,j+2,k)*dup2-c(43,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i,j+1,k+2)-u(2,i,j+1,k-2))+
c$$$     &a1*(u(2,i,j+1,k+1)-u(2,i,j+1,k-1))
c$$$      dum1 = a2*(u(2,i,j-1,k+2)-u(2,i,j-1,k-2))+
c$$$     &a1*(u(2,i,j-1,k+1)-u(2,i,j-1,k-1))
c$$$      dup2 = a2*(u(2,i,j+2,k+2)-u(2,i,j+2,k-2))+
c$$$     &a1*(u(2,i,j+2,k+1)-u(2,i,j+2,k-1))
c$$$      dum2 = a2*(u(2,i,j-2,k+2)-u(2,i,j-2,k-2))+
c$$$     &a1*(u(2,i,j-2,k+1)-u(2,i,j-2,k-1))
c$$$      r1 =r1+strz(k)*stry(j)*(a1*(c(38,i,j+1,k)*dup1-c(38,i,j-1,k)*dum1)
c$$$     &+a2*(c(38,i,j+2,k)*dup2-c(38,i,j-2,k)*dum2))
c$$$      r2 =r2+strz(k)*stry(j)*(a1*(c(41,i,j+1,k)*dup1-c(41,i,j-1,k)*dum1)
c$$$     &+a2*(c(41,i,j+2,k)*dup2-c(41,i,j-2,k)*dum2))
c$$$      r3 =r3+strz(k)*stry(j)*(a1*(c(44,i,j+1,k)*dup1-c(44,i,j-1,k)*dum1)
c$$$     &+a2*(c(44,i,j+2,k)*dup2-c(44,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i,j+1,k+2)-u(3,i,j+1,k-2))+
c$$$     &a1*(u(3,i,j+1,k+1)-u(3,i,j+1,k-1))
c$$$      dum1 = a2*(u(3,i,j-1,k+2)-u(3,i,j-1,k-2))+
c$$$     &a1*(u(3,i,j-1,k+1)-u(3,i,j-1,k-1))
c$$$      dup2 = a2*(u(3,i,j+2,k+2)-u(3,i,j+2,k-2))+
c$$$     &a1*(u(3,i,j+2,k+1)-u(3,i,j+2,k-1))
c$$$      dum2 = a2*(u(3,i,j-2,k+2)-u(3,i,j-2,k-2))+
c$$$     &a1*(u(3,i,j-2,k+1)-u(3,i,j-2,k-1))
c$$$      r1 =r1+strz(k)*stry(j)*(a1*(c(39,i,j+1,k)*dup1-c(39,i,j-1,k)*dum1)
c$$$     &+a2*(c(39,i,j+2,k)*dup2-c(39,i,j-2,k)*dum2))
c$$$      r2 =r2+strz(k)*stry(j)*(a1*(c(42,i,j+1,k)*dup1-c(42,i,j-1,k)*dum1)
c$$$     &+a2*(c(42,i,j+2,k)*dup2-c(42,i,j-2,k)*dum2))
c$$$      r3 =r3+strz(k)*stry(j)*(a1*(c(45,i,j+1,k)*dup1-c(45,i,j-1,k)*dum1)
c$$$     &+a2*(c(45,i,j+2,k)*dup2-c(45,i,j-2,k)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(1,i,j+2,k+1)-u(1,i,j-2,k+1))+
c$$$     &a1*(u(1,i,j+1,k+1)-u(1,i,j-1,k+1))
c$$$      dum1 = a2*(u(1,i,j+2,k-1)-u(1,i,j-2,k-1))+
c$$$     &a1*(u(1,i,j+1,k-1)-u(1,i,j-1,k-1))
c$$$      dup2 = a2*(u(1,i,j+2,k+2)-u(1,i,j-2,k+2))+
c$$$     &a1*(u(1,i,j+1,k+2)-u(1,i,j-1,k+2))
c$$$      dum2 = a2*(u(1,i,j+2,k-2)-u(1,i,j-2,k-2))+
c$$$     &a1*(u(1,i,j+1,k-2)-u(1,i,j-1,k-2))
c$$$      r1 =r1+stry(j)*strz(k)*(a1*(c(37,i,j,k+1)*dup1-c(37,i,j,k-1)*dum1)
c$$$     &+a2*(c(37,i,j,k+2)*dup2-c(37,i,j,k-2)*dum2))
c$$$      r2 =r2+stry(j)*strz(k)*(a1*(c(38,i,j,k+1)*dup1-c(38,i,j,k-1)*dum1)
c$$$     &+a2*(c(38,i,j,k+2)*dup2-c(38,i,j,k-2)*dum2))
c$$$      r3 =r3+stry(j)*strz(k)*(a1*(c(39,i,j,k+1)*dup1-c(39,i,j,k-1)*dum1)
c$$$     &+a2*(c(39,i,j,k+2)*dup2-c(39,i,j,k-2)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(2,i,j+2,k+1)-u(2,i,j-2,k+1))+
c$$$     &a1*(u(2,i,j+1,k+1)-u(2,i,j-1,k+1))
c$$$      dum1 = a2*(u(2,i,j+2,k-1)-u(2,i,j-2,k-1))+
c$$$     &a1*(u(2,i,j+1,k-1)-u(2,i,j-1,k-1))
c$$$      dup2 = a2*(u(2,i,j+2,k+2)-u(2,i,j-2,k+2))+
c$$$     &a1*(u(2,i,j+1,k+2)-u(2,i,j-1,k+2))
c$$$      dum2 = a2*(u(2,i,j+2,k-2)-u(2,i,j-2,k-2))+
c$$$     &a1*(u(2,i,j+1,k-2)-u(2,i,j-1,k-2))
c$$$      r1 =r1+stry(j)*strz(k)*(a1*(c(40,i,j,k+1)*dup1-c(40,i,j,k-1)*dum1)
c$$$     &+a2*(c(40,i,j,k+2)*dup2-c(40,i,j,k-2)*dum2))
c$$$      r2 =r2+stry(j)*strz(k)*(a1*(c(41,i,j,k+1)*dup1-c(41,i,j,k-1)*dum1)
c$$$     &+a2*(c(41,i,j,k+2)*dup2-c(41,i,j,k-2)*dum2))
c$$$      r3 =r3+stry(j)*strz(k)*(a1*(c(42,i,j,k+1)*dup1-c(42,i,j,k-1)*dum1)
c$$$     &+a2*(c(42,i,j,k+2)*dup2-c(42,i,j,k-2)*dum2))
c$$$
c$$$
c$$$      dup1 = a2*(u(3,i,j+2,k+1)-u(3,i,j-2,k+1))+
c$$$     &a1*(u(3,i,j+1,k+1)-u(3,i,j-1,k+1))
c$$$      dum1 = a2*(u(3,i,j+2,k-1)-u(3,i,j-2,k-1))+
c$$$     &a1*(u(3,i,j+1,k-1)-u(3,i,j-1,k-1))
c$$$      dup2 = a2*(u(3,i,j+2,k+2)-u(3,i,j-2,k+2))+
c$$$     &a1*(u(3,i,j+1,k+2)-u(3,i,j-1,k+2))
c$$$      dum2 = a2*(u(3,i,j+2,k-2)-u(3,i,j-2,k-2))+
c$$$     &a1*(u(3,i,j+1,k-2)-u(3,i,j-1,k-2))
c$$$      r1 =r1+stry(j)*strz(k)*(a1*(c(43,i,j,k+1)*dup1-c(43,i,j,k-1)*dum1)
c$$$     &+a2*(c(43,i,j,k+2)*dup2-c(43,i,j,k-2)*dum2))
c$$$      r2 =r2+stry(j)*strz(k)*(a1*(c(44,i,j,k+1)*dup1-c(44,i,j,k-1)*dum1)
c$$$     &+a2*(c(44,i,j,k+2)*dup2-c(44,i,j,k-2)*dum2))
c$$$      r3 =r3+stry(j)*strz(k)*(a1*(c(45,i,j,k+1)*dup1-c(45,i,j,k-1)*dum1)
c$$$     &+a2*(c(45,i,j,k+2)*dup2-c(45,i,j,k-2)*dum2))

      lu(1,i,j,k) = r1*ijac
      lu(2,i,j,k) = r2*ijac
      lu(3,i,j,k) = r3*ijac
      enddo
      enddo
      enddo
      if( onesided(6).eq.1 )then
         do k=nk-5,nk
            do j=jfirst+2,jlast-2
               do i=ifirst+2,ilast-2
                  ijac = 1/jac(i,j,k)
                  r1 = 0
                  r2 = 0
                  r3 = 0
      cm2 = c(1,i-1,j,k)*strx(i-1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i-2,j,k)*strx(i-2))
      cm1 = c(1,i-2,j,k)*strx(i-2)+c(1,i+1,j,k)*strx(i+1)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i-1,j,k)*strx(i-1))
      cp1 = c(1,i-1,j,k)*strx(i-1)+c(1,i+2,j,k)*strx(i+2)+3*(
     &c(1,i,j,k)*strx(i)+c(1,i+1,j,k)*strx(i+1))
      cp2 = c(1,i+1,j,k)*strx(i+1)-0.75d0*(c(1,i,j,k)*strx(i)+
     &c(1,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(2,i-1,j,k)*strx(i-1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i-2,j,k)*strx(i-2))
      cm1 = c(2,i-2,j,k)*strx(i-2)+c(2,i+1,j,k)*strx(i+1)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i-1,j,k)*strx(i-1))
      cp1 = c(2,i-1,j,k)*strx(i-1)+c(2,i+2,j,k)*strx(i+2)+3*(
     &c(2,i,j,k)*strx(i)+c(2,i+1,j,k)*strx(i+1))
      cp2 = c(2,i+1,j,k)*strx(i+1)-0.75d0*(c(2,i,j,k)*strx(i)+
     &c(2,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      r2 = r2+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(3,i-1,j,k)*strx(i-1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i-2,j,k)*strx(i-2))
      cm1 = c(3,i-2,j,k)*strx(i-2)+c(3,i+1,j,k)*strx(i+1)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i-1,j,k)*strx(i-1))
      cp1 = c(3,i-1,j,k)*strx(i-1)+c(3,i+2,j,k)*strx(i+2)+3*(
     &c(3,i,j,k)*strx(i)+c(3,i+1,j,k)*strx(i+1))
      cp2 = c(3,i+1,j,k)*strx(i+1)-0.75d0*(c(3,i,j,k)*strx(i)+
     &c(3,i+2,j,k)*strx(i+2))

      r1 = r1+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(1,i-2,j,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i-1,j,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i+1,j,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i+2,j,k)-u(1,i,j,k)) )
      cm2 = c(4,i-1,j,k)*strx(i-1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i-2,j,k)*strx(i-2))
      cm1 = c(4,i-2,j,k)*strx(i-2)+c(4,i+1,j,k)*strx(i+1)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i-1,j,k)*strx(i-1))
      cp1 = c(4,i-1,j,k)*strx(i-1)+c(4,i+2,j,k)*strx(i+2)+3*(
     &c(4,i,j,k)*strx(i)+c(4,i+1,j,k)*strx(i+1))
      cp2 = c(4,i+1,j,k)*strx(i+1)-0.75d0*(c(4,i,j,k)*strx(i)+
     &c(4,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(5,i-1,j,k)*strx(i-1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i-2,j,k)*strx(i-2))
      cm1 = c(5,i-2,j,k)*strx(i-2)+c(5,i+1,j,k)*strx(i+1)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i-1,j,k)*strx(i-1))
      cp1 = c(5,i-1,j,k)*strx(i-1)+c(5,i+2,j,k)*strx(i+2)+3*(
     &c(5,i,j,k)*strx(i)+c(5,i+1,j,k)*strx(i+1))
      cp2 = c(5,i+1,j,k)*strx(i+1)-0.75d0*(c(5,i,j,k)*strx(i)+
     &c(5,i+2,j,k)*strx(i+2))

      r2 = r2+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      r3 = r3+i6*strx(i)*(cm2*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i+1,j,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i+2,j,k)-u(2,i,j,k)) )
      cm2 = c(6,i-1,j,k)*strx(i-1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i-2,j,k)*strx(i-2))
      cm1 = c(6,i-2,j,k)*strx(i-2)+c(6,i+1,j,k)*strx(i+1)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i-1,j,k)*strx(i-1))
      cp1 = c(6,i-1,j,k)*strx(i-1)+c(6,i+2,j,k)*strx(i+2)+3*(
     &c(6,i,j,k)*strx(i)+c(6,i+1,j,k)*strx(i+1))
      cp2 = c(6,i+1,j,k)*strx(i+1)-0.75d0*(c(6,i,j,k)*strx(i)+
     &c(6,i+2,j,k)*strx(i+2))

      r3 = r3+i6*strx(i)*(cm2*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i+1,j,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i+2,j,k)-u(3,i,j,k)) )
      cm2 = c(7,i,j-1,k)*stry(j-1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j-2,k)*stry(j-2))
      cm1 = c(7,i,j-2,k)*stry(j-2)+c(7,i,j+1,k)*stry(j+1)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j-1,k)*stry(j-1))
      cp1 = c(7,i,j-1,k)*stry(j-1)+c(7,i,j+2,k)*stry(j+2)+3*(
     &c(7,i,j,k)*stry(j)+c(7,i,j+1,k)*stry(j+1))
      cp2 = c(7,i,j+1,k)*stry(j+1)-0.75d0*(c(7,i,j,k)*stry(j)+
     &c(7,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(8,i,j-1,k)*stry(j-1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j-2,k)*stry(j-2))
      cm1 = c(8,i,j-2,k)*stry(j-2)+c(8,i,j+1,k)*stry(j+1)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j-1,k)*stry(j-1))
      cp1 = c(8,i,j-1,k)*stry(j-1)+c(8,i,j+2,k)*stry(j+2)+3*(
     &c(8,i,j,k)*stry(j)+c(8,i,j+1,k)*stry(j+1))
      cp2 = c(8,i,j+1,k)*stry(j+1)-0.75d0*(c(8,i,j,k)*stry(j)+
     &c(8,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      r2 = r2+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(9,i,j-1,k)*stry(j-1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j-2,k)*stry(j-2))
      cm1 = c(9,i,j-2,k)*stry(j-2)+c(9,i,j+1,k)*stry(j+1)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j-1,k)*stry(j-1))
      cp1 = c(9,i,j-1,k)*stry(j-1)+c(9,i,j+2,k)*stry(j+2)+3*(
     &c(9,i,j,k)*stry(j)+c(9,i,j+1,k)*stry(j+1))
      cp2 = c(9,i,j+1,k)*stry(j+1)-0.75d0*(c(9,i,j,k)*stry(j)+
     &c(9,i,j+2,k)*stry(j+2))

      r1 = r1+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     &cm1*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     &cp1*(u(1,i,j+1,k)-u(1,i,j,k)) + 
     &cp2*(u(1,i,j+2,k)-u(1,i,j,k)) )
      cm2 = c(10,i,j-1,k)*stry(j-1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j-2,k)*stry(j-2))
      cm1 = c(10,i,j-2,k)*stry(j-2)+c(10,i,j+1,k)*stry(j+1)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j-1,k)*stry(j-1))
      cp1 = c(10,i,j-1,k)*stry(j-1)+c(10,i,j+2,k)*stry(j+2)+3*(
     &c(10,i,j,k)*stry(j)+c(10,i,j+1,k)*stry(j+1))
      cp2 = c(10,i,j+1,k)*stry(j+1)-0.75d0*(c(10,i,j,k)*stry(j)+
     &c(10,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(11,i,j-1,k)*stry(j-1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j-2,k)*stry(j-2))
      cm1 = c(11,i,j-2,k)*stry(j-2)+c(11,i,j+1,k)*stry(j+1)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j-1,k)*stry(j-1))
      cp1 = c(11,i,j-1,k)*stry(j-1)+c(11,i,j+2,k)*stry(j+2)+3*(
     &c(11,i,j,k)*stry(j)+c(11,i,j+1,k)*stry(j+1))
      cp2 = c(11,i,j+1,k)*stry(j+1)-0.75d0*(c(11,i,j,k)*stry(j)+
     &c(11,i,j+2,k)*stry(j+2))

      r2 = r2+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r3 = r3+i6*stry(j)*(cm2*(u(2,i,j-2,k)-u(2,i,j,k)) + 
     &cm1*(u(2,i,j-1,k)-u(2,i,j,k)) + 
     &cp1*(u(2,i,j+1,k)-u(2,i,j,k)) + 
     &cp2*(u(2,i,j+2,k)-u(2,i,j,k)) )
      cm2 = c(12,i,j-1,k)*stry(j-1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j-2,k)*stry(j-2))
      cm1 = c(12,i,j-2,k)*stry(j-2)+c(12,i,j+1,k)*stry(j+1)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j-1,k)*stry(j-1))
      cp1 = c(12,i,j-1,k)*stry(j-1)+c(12,i,j+2,k)*stry(j+2)+3*(
     &c(12,i,j,k)*stry(j)+c(12,i,j+1,k)*stry(j+1))
      cp2 = c(12,i,j+1,k)*stry(j+1)-0.75d0*(c(12,i,j,k)*stry(j)+
     &c(12,i,j+2,k)*stry(j+2))

      r3 = r3+i6*stry(j)*(cm2*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     &cm1*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     &cp1*(u(3,i,j+1,k)-u(3,i,j,k)) + 
     &cp2*(u(3,i,j+2,k)-u(3,i,j,k)) )
      r1 = r1 + ghcof(nk-k+1)*c(13,i,j,nk)*u(1,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(14,i,j,nk)*u(2,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(15,i,j,nk)*u(3,i,j,nk+1)
      r2 = r2 + ghcof(nk-k+1)*c(14,i,j,nk)*u(1,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(16,i,j,nk)*u(2,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(17,i,j,nk)*u(3,i,j,nk+1)
      r3 = r3 + ghcof(nk-k+1)*c(15,i,j,nk)*u(1,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(17,i,j,nk)*u(2,i,j,nk+1) + 
     &ghcof(nk-k+1)*c(18,i,j,nk)*u(3,i,j,nk+1)
      do q=nk-7,nk
      ac1=0
      ac2=0
      ac3=0
      ac4=0
      ac5=0
      ac6=0
      do m=nk-7,nk
      ac1=  ac1+ acof(nk-k+1,nk-q+1,nk-m+1)*c(13,i,j,m) 
      ac2=  ac2+ acof(nk-k+1,nk-q+1,nk-m+1)*c(14,i,j,m) 
      ac3=  ac3+ acof(nk-k+1,nk-q+1,nk-m+1)*c(15,i,j,m) 
      ac4=  ac4+ acof(nk-k+1,nk-q+1,nk-m+1)*c(16,i,j,m) 
      ac5=  ac5+ acof(nk-k+1,nk-q+1,nk-m+1)*c(17,i,j,m) 
      ac6=  ac6+ acof(nk-k+1,nk-q+1,nk-m+1)*c(18,i,j,m) 
      enddo
      r1 = r1 + ac1*u(1,i,j,q) + ac2*u(2,i,j,q) + ac3*u(3,i,j,q)
      r2 = r2 + ac2*u(1,i,j,q) + ac4*u(2,i,j,q) + ac5*u(3,i,j,q)
      r3 = r3 + ac3*u(1,i,j,q) + ac5*u(2,i,j,q) + ac6*u(3,i,j,q)
      enddo



      dup1 = a2*(u(1,i+1,j+2,k)-u(1,i+1,j-2,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i+1,j-1,k))
      dum1 = a2*(u(1,i-1,j+2,k)-u(1,i-1,j-2,k))+
     &a1*(u(1,i-1,j+1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i+2,j-2,k))+
     &a1*(u(1,i+2,j+1,k)-u(1,i+2,j-1,k))
      dum2 = a2*(u(1,i-2,j+2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i-2,j+1,k)-u(1,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(19,i+1,j,k)*dup1-
     &c(19,i-1,j,k)*dum1)+a2*(c(19,i+2,j,k)*dup2-c(19,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(22,i+1,j,k)*dup1-
     &c(22,i-1,j,k)*dum1)+a2*(c(22,i+2,j,k)*dup2-c(22,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(25,i+1,j,k)*dup1-
     &c(25,i-1,j,k)*dum1)+a2*(c(25,i+2,j,k)*dup2-c(25,i-2,j,k)*dum2))


      dup1 = a2*(u(2,i+1,j+2,k)-u(2,i+1,j-2,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i+1,j-1,k))
      dum1 = a2*(u(2,i-1,j+2,k)-u(2,i-1,j-2,k))+
     &a1*(u(2,i-1,j+1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i+2,j-2,k))+
     &a1*(u(2,i+2,j+1,k)-u(2,i+2,j-1,k))
      dum2 = a2*(u(2,i-2,j+2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i-2,j+1,k)-u(2,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(20,i+1,j,k)*dup1-
     &c(20,i-1,j,k)*dum1)+a2*(c(20,i+2,j,k)*dup2-c(20,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(23,i+1,j,k)*dup1-
     &c(23,i-1,j,k)*dum1)+a2*(c(23,i+2,j,k)*dup2-c(23,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(26,i+1,j,k)*dup1-
     &c(26,i-1,j,k)*dum1)+a2*(c(26,i+2,j,k)*dup2-c(26,i-2,j,k)*dum2))


      dup1 = a2*(u(3,i+1,j+2,k)-u(3,i+1,j-2,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i+1,j-1,k))
      dum1 = a2*(u(3,i-1,j+2,k)-u(3,i-1,j-2,k))+
     &a1*(u(3,i-1,j+1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i+2,j-2,k))+
     &a1*(u(3,i+2,j+1,k)-u(3,i+2,j-1,k))
      dum2 = a2*(u(3,i-2,j+2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i-2,j+1,k)-u(3,i-2,j-1,k))
      r1 = r1 +stry(j)*strx(i)*(a1*(c(21,i+1,j,k)*dup1-
     &c(21,i-1,j,k)*dum1)+a2*(c(21,i+2,j,k)*dup2-c(21,i-2,j,k)*dum2))
      r2 = r2 +stry(j)*strx(i)*(a1*(c(24,i+1,j,k)*dup1-
     &c(24,i-1,j,k)*dum1)+a2*(c(24,i+2,j,k)*dup2-c(24,i-2,j,k)*dum2))
      r3 = r3 +stry(j)*strx(i)*(a1*(c(27,i+1,j,k)*dup1-
     &c(27,i-1,j,k)*dum1)+a2*(c(27,i+2,j,k)*dup2-c(27,i-2,j,k)*dum2))


      dup1 = a2*(u(1,i+2,j+1,k)-u(1,i-2,j+1,k))+
     &a1*(u(1,i+1,j+1,k)-u(1,i-1,j+1,k))
      dum1 = a2*(u(1,i+2,j-1,k)-u(1,i-2,j-1,k))+
     &a1*(u(1,i+1,j-1,k)-u(1,i-1,j-1,k))
      dup2 = a2*(u(1,i+2,j+2,k)-u(1,i-2,j+2,k))+
     &a1*(u(1,i+1,j+2,k)-u(1,i-1,j+2,k))
      dum2 = a2*(u(1,i+2,j-2,k)-u(1,i-2,j-2,k))+
     &a1*(u(1,i+1,j-2,k)-u(1,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(19,i,j+1,k)*dup1-
     &c(19,i,j-1,k)*dum1)+a2*(c(19,i,j+2,k)*dup2-c(19,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(20,i,j+1,k)*dup1-
     &c(20,i,j-1,k)*dum1)+a2*(c(20,i,j+2,k)*dup2-c(20,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(21,i,j+1,k)*dup1-
     &c(21,i,j-1,k)*dum1)+a2*(c(21,i,j+2,k)*dup2-c(21,i,j-2,k)*dum2))


      dup1 = a2*(u(2,i+2,j+1,k)-u(2,i-2,j+1,k))+
     &a1*(u(2,i+1,j+1,k)-u(2,i-1,j+1,k))
      dum1 = a2*(u(2,i+2,j-1,k)-u(2,i-2,j-1,k))+
     &a1*(u(2,i+1,j-1,k)-u(2,i-1,j-1,k))
      dup2 = a2*(u(2,i+2,j+2,k)-u(2,i-2,j+2,k))+
     &a1*(u(2,i+1,j+2,k)-u(2,i-1,j+2,k))
      dum2 = a2*(u(2,i+2,j-2,k)-u(2,i-2,j-2,k))+
     &a1*(u(2,i+1,j-2,k)-u(2,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(22,i,j+1,k)*dup1-
     &c(22,i,j-1,k)*dum1)+a2*(c(22,i,j+2,k)*dup2-c(22,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(23,i,j+1,k)*dup1-
     &c(23,i,j-1,k)*dum1)+a2*(c(23,i,j+2,k)*dup2-c(23,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(24,i,j+1,k)*dup1-
     &c(24,i,j-1,k)*dum1)+a2*(c(24,i,j+2,k)*dup2-c(24,i,j-2,k)*dum2))


      dup1 = a2*(u(3,i+2,j+1,k)-u(3,i-2,j+1,k))+
     &a1*(u(3,i+1,j+1,k)-u(3,i-1,j+1,k))
      dum1 = a2*(u(3,i+2,j-1,k)-u(3,i-2,j-1,k))+
     &a1*(u(3,i+1,j-1,k)-u(3,i-1,j-1,k))
      dup2 = a2*(u(3,i+2,j+2,k)-u(3,i-2,j+2,k))+
     &a1*(u(3,i+1,j+2,k)-u(3,i-1,j+2,k))
      dum2 = a2*(u(3,i+2,j-2,k)-u(3,i-2,j-2,k))+
     &a1*(u(3,i+1,j-2,k)-u(3,i-1,j-2,k))
      r1 = r1 +strx(i)*stry(j)*(a1*(c(25,i,j+1,k)*dup1-
     &c(25,i,j-1,k)*dum1)+a2*(c(25,i,j+2,k)*dup2-c(25,i,j-2,k)*dum2))
      r2 = r2 +strx(i)*stry(j)*(a1*(c(26,i,j+1,k)*dup1-
     &c(26,i,j-1,k)*dum1)+a2*(c(26,i,j+2,k)*dup2-c(26,i,j-2,k)*dum2))
      r3 = r3 +strx(i)*stry(j)*(a1*(c(27,i,j+1,k)*dup1-
     &c(27,i,j-1,k)*dum1)+a2*(c(27,i,j+2,k)*dup2-c(27,i,j-2,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(1,i+2,j,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(1,i+1,j,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(1,i-1,j,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(1,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(28,i+1,j,k)*dup1-c(28,i-1,j,k)*dum1)
     &+ a2*(c(28,i+2,j,k)*dup2-c(28,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(31,i+1,j,k)*dup1-c(31,i-1,j,k)*dum1)
     &+ a2*(c(31,i+2,j,k)*dup2-c(31,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(34,i+1,j,k)*dup1-c(34,i-1,j,k)*dum1)
     &+ a2*(c(34,i+2,j,k)*dup2-c(34,i-2,j,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(2,i+2,j,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(2,i+1,j,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(2,i-1,j,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(2,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(29,i+1,j,k)*dup1-c(29,i-1,j,k)*dum1)
     &+ a2*(c(29,i+2,j,k)*dup2-c(29,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(32,i+1,j,k)*dup1-c(32,i-1,j,k)*dum1)
     &+ a2*(c(32,i+2,j,k)*dup2-c(32,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(35,i+1,j,k)*dup1-c(35,i-1,j,k)*dum1)
     &+ a2*(c(35,i+2,j,k)*dup2-c(35,i-2,j,k)*dum2))


      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(3,i+2,j,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(3,i+1,j,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(3,i-1,j,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(3,i-2,j,q)
      enddo
      r1 = r1 +strx(i)*( a1*(c(30,i+1,j,k)*dup1-c(30,i-1,j,k)*dum1)
     &+ a2*(c(30,i+2,j,k)*dup2-c(30,i-2,j,k)*dum2))
      r2 = r2 +strx(i)*( a1*(c(33,i+1,j,k)*dup1-c(33,i-1,j,k)*dum1)
     &+ a2*(c(33,i+2,j,k)*dup2-c(33,i-2,j,k)*dum2))
      r3 = r3 +strx(i)*( a1*(c(36,i+1,j,k)*dup1-c(36,i-1,j,k)*dum1)
     &+ a2*(c(36,i+2,j,k)*dup2-c(36,i-2,j,k)*dum2))


      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(1,i+2,j,q)-u(1,i-2,j,q))+
     &a1*( u(1,i+1,j,q)-u(1,i-1,j,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(28,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(29,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(30,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3



      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(2,i+2,j,q)-u(2,i-2,j,q))+
     &a1*( u(2,i+1,j,q)-u(2,i-1,j,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(31,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(32,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(33,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3



      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(3,i+2,j,q)-u(3,i-2,j,q))+
     &a1*( u(3,i+1,j,q)-u(3,i-1,j,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(34,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(35,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(36,i,j,q)*du
      enddo
      r1 = r1 + strx(i)*ac1
      r2 = r2 + strx(i)*ac2
      r3 = r3 + strx(i)*ac3

      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(1,i,j+2,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(1,i,j+1,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(1,i,j-1,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(1,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(37,i,j+1,k)*dup1-c(37,i,j-1,k)*dum1)
     &+ a2*(c(37,i,j+2,k)*dup2-c(37,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(40,i,j+1,k)*dup1-c(40,i,j-1,k)*dum1)
     &+ a2*(c(40,i,j+2,k)*dup2-c(40,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(43,i,j+1,k)*dup1-c(43,i,j-1,k)*dum1)
     &+ a2*(c(43,i,j+2,k)*dup2-c(43,i,j-2,k)*dum2))

      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(2,i,j+2,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(2,i,j+1,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(2,i,j-1,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(2,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(38,i,j+1,k)*dup1-c(38,i,j-1,k)*dum1)
     &+ a2*(c(38,i,j+2,k)*dup2-c(38,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(41,i,j+1,k)*dup1-c(41,i,j-1,k)*dum1)
     &+ a2*(c(41,i,j+2,k)*dup2-c(41,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(44,i,j+1,k)*dup1-c(44,i,j-1,k)*dum1)
     &+ a2*(c(44,i,j+2,k)*dup2-c(44,i,j-2,k)*dum2))
      dum2=0
      dum1=0
      dup1=0
      dup2=0
      do q=nk-7,nk
      dup2 = dup2 -bop(nk-k+1,nk-q+1)*u(3,i,j+2,q)
      dup1 = dup1 -bop(nk-k+1,nk-q+1)*u(3,i,j+1,q)
      dum1 = dum1 -bop(nk-k+1,nk-q+1)*u(3,i,j-1,q)
      dum2 = dum2 -bop(nk-k+1,nk-q+1)*u(3,i,j-2,q)
      enddo
      r1 = r1 +stry(j)*( a1*(c(39,i,j+1,k)*dup1-c(39,i,j-1,k)*dum1)
     &+ a2*(c(39,i,j+2,k)*dup2-c(39,i,j-2,k)*dum2))
      r2 = r2 +stry(j)*( a1*(c(42,i,j+1,k)*dup1-c(42,i,j-1,k)*dum1)
     &+ a2*(c(42,i,j+2,k)*dup2-c(42,i,j-2,k)*dum2))
      r3 = r3 +stry(j)*( a1*(c(45,i,j+1,k)*dup1-c(45,i,j-1,k)*dum1)
     &+ a2*(c(45,i,j+2,k)*dup2-c(45,i,j-2,k)*dum2))
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(1,i,j+2,q)-u(1,i,j-2,q))+
     &a1*( u(1,i,j+1,q)-u(1,i,j-1,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(37,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(38,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(39,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(2,i,j+2,q)-u(2,i,j-2,q))+
     &a1*( u(2,i,j+1,q)-u(2,i,j-1,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(40,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(41,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(42,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3
      ac1 = 0
      ac2 = 0
      ac3 = 0
      do q=nk-7,nk
         du = a2*(u(3,i,j+2,q)-u(3,i,j-2,q))+
     &a1*( u(3,i,j+1,q)-u(3,i,j-1,q))
         ac1 = ac1-bop(nk-k+1,nk-q+1)*c(43,i,j,q)*du
         ac2 = ac2-bop(nk-k+1,nk-q+1)*c(44,i,j,q)*du
         ac3 = ac3-bop(nk-k+1,nk-q+1)*c(45,i,j,q)*du
      enddo
      r1 = r1 + stry(j)*ac1
      r2 = r2 + stry(j)*ac2
      r3 = r3 + stry(j)*ac3

      lu(1,i,j,k) = r1*ijac
      lu(2,i,j,k) = r2*ijac
      lu(3,i,j,k) = r3*ijac
               enddo
            enddo
         enddo
      endif   
      end
