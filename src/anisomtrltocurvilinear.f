      subroutine ANISOMTRLTOCURVILINEAR( ifirst, ilast, jfirst, jlast, 
     *        kfirst, klast, met, c, cnew )

***********************************************************************
***
*** Set up matrices, C_{i,j} for the anisotropic elastic wave equation
*** in curvilinear coordinates:
***
***  div(sigma) = 1/J*(
***               (C_{1,1}u_p)_p + (C_{1,2}u_q)_p + (C_{1,3}u_r)_p + 
***               (C_{2,1}u_p)_q + (C_{2,2}u_q)_q + (C_{2,3}u_r)_q + 
***               (C_{3,1}u_p)_r + (C_{3,2}u_q)_r + (C_{3,3}u_r)_r  )
***
***   where (p,q,r) are the curvilinear coordinates.
***
***   Input: ifirst, ilast, jfirst, jlast, kfirst, klast - Declared
***              dimensions of arrays.
***          met - Metric derivatives, as stored in SW4
***          c   - Anisotropic tensor in Cartesian formulation
***   Output: cnew - Anisotropic tensor in Curvilinear formulation (C_{i,j} above).
***       The matrices C_{i,i} are symmetric and 
***            C_{i,j} = C_{j,i}^T when i is not equal to j.
***       Hence, only need to store 6+6+6+9+9+9 = 45 coefficients at each grid point.
***
*** The input Cartesian stress tensor has 21 elements, these are assumed enumerated
*** in the non-standard way used internally in SW4.
***
*** The output stress tensor has 45 elements. These are ordered according to
***  A11=[c1 c2 c3]  A22=[c7 c8 c9  ]  A33=[c13 c14 c15]  A21=[c19 c20 c21]
***      [c2 c4 c5]      [c8 c10 c11]      [c14 c16 c17]      [c22 c23 c24]
***      [c3 c5 c6]      [c9 c11 c12]      [c15 c17 c18]      [c25 c26 c27]
***
***  A31=[c28 c29 c30]   A32=[c37 c38 c39]
***      [c31 c32 c33]       [c40 c41 c42]
***      [c34 c35 c36]       [c43 c44 c45]
*** 
*** where the divergence of stress operator is of the form
***
***  L(u) = 1/J*( Dp(A11*Dp(u)) + Dq(A22*Dq(u)) + Dr(A33*Dr(u)) + Dp(A21*Dq(u))+
***               Dp(A31*Dr(u)) + Dq(A32*Dr(u)) + Dq(A12*Dp(u)) + Dr(A13*Dp(u)) +
***               Dr(A23*Dq(u)) )
***
***  Matrices A12, A13, and A23 are obtained by 
***  transposing A21, A31, and A32 respectively.
***********************************************************************

      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer is, js, ks, i, j, k, l, n, p, cind(3,3,3,3), i45
      real*8 c(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 met(4,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 cnew(45,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mder(3,3)

*** xx-matrix
      cind(1,1,1,1) = 1
      cind(1,2,1,1) = 2
      cind(1,3,1,1) = 3
      cind(2,1,1,1) = 2
      cind(2,2,1,1) = 7
      cind(2,3,1,1) = 8
      cind(3,1,1,1) = 3
      cind(3,2,1,1) = 8
      cind(3,3,1,1) = 12

*** yy-matrix
      cind(1,1,2,2) = 7
      cind(1,2,2,2) = 9 
      cind(1,3,2,2) = 10
      cind(2,1,2,2) = 9
      cind(2,2,2,2) = 16
      cind(2,3,2,2) = 17
      cind(3,1,2,2) = 10
      cind(3,2,2,2) = 17
      cind(3,3,2,2) = 19

*** zz-matrix
      cind(1,1,3,3) = 12
      cind(1,2,3,3) = 14
      cind(1,3,3,3) = 15
      cind(2,1,3,3) = 14
      cind(2,2,3,3) = 19
      cind(2,3,3,3) = 20
      cind(3,1,3,3) = 15
      cind(3,2,3,3) = 20
      cind(3,3,3,3) = 21

*** yx-matrix
      cind(1,1,1,2) = 2
      cind(1,2,1,2) = 4
      cind(1,3,1,2) = 5
      cind(2,1,1,2) = 7
      cind(2,2,1,2) = 9
      cind(2,3,1,2) = 10
      cind(3,1,1,2) = 8
      cind(3,2,1,2) = 13
      cind(3,3,1,2) = 14
      do j=1,3
         do i=1,3
            cind(i,j,2,1)=cind(j,i,1,2)
         enddo
      enddo
         
*** zx-matrix
      cind(1,1,1,3) = 3
      cind(1,2,1,3) = 5
      cind(1,3,1,3) = 6
      cind(2,1,1,3) = 8
      cind(2,2,1,3) = 10
      cind(2,3,1,3) = 11
      cind(3,1,1,3) = 12
      cind(3,2,1,3) = 14
      cind(3,3,1,3) = 15
      do j=1,3
         do i=1,3
            cind(i,j,3,1)=cind(j,i,1,3)
         enddo
      enddo

*** zy-matrix
      cind(1,1,2,3) = 8
      cind(1,2,2,3) = 10
      cind(1,3,2,3) = 11
      cind(2,1,2,3) = 13
      cind(2,2,2,3) = 17
      cind(2,3,2,3) = 18
      cind(3,1,2,3) = 14
      cind(3,2,2,3) = 19
      cind(3,3,2,3) = 20
      do j=1,3
         do i=1,3
            cind(i,j,3,2)=cind(j,i,2,3)
         enddo
      enddo

      do ks=kfirst,klast
         do js=jfirst,jlast
            do is=ifirst,ilast
*** General expression for metric. Input a full 3x3 mder-matrix instead of met 
*** will work for general curvilinear grids.
               mder(1,1) = met(1,is,js,ks)
               mder(1,2) = 0
               mder(1,3) = met(2,is,js,ks)
               mder(2,1) = 0
               mder(2,2) = met(1,is,js,ks)
               mder(2,3) = met(3,is,js,ks)
               mder(3,1) = 0
               mder(3,2) = 0
               mder(3,3) = met(4,is,js,ks)
               i45 = 1
               k = 1
               l = 1
*** matrix A in (A u_x)_x 
               do n=1,3
                  do p=n,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
               k = 2
               l = 2
*** matrix A in (A u_y)_y 
               do n=1,3
                  do p=n,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                    c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
               k = 3
               l = 3
***   matrix A in (A u_z)_z 
               do n=1,3
                  do p=n,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
               k = 1
               l = 2
*** matrix A in (A u_y)_x 
               do n=1,3
                  do p=1,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
               k = 1
               l = 3
*** matrix A in (A u_z)_x 
               do n=1,3
                  do p=1,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
               k = 2
               l = 3
*** matrix A in (A u_z)_y
               do n=1,3
                  do p=1,3
                     cnew(i45,is,js,ks) = 0
                     do j=1,3
                        do i=1,3
                           cnew(i45,is,js,ks) = cnew(i45,is,js,ks) + 
     *                c(cind(n,p,i,j),is,js,ks)*mder(i,k)*mder(j,l)
                        enddo
                     enddo
                     i45 = i45 + 1
                  enddo
               enddo
            enddo
         enddo
      enddo

      end
