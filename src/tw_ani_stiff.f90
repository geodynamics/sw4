subroutine tw_ani_stiff( ifirst, ilast, jfirst, jlast, kfirst, klast, h, zmin, omm, phm, amprho, rho, phc, cm) bind(c)
  use iso_c_binding
  implicit none
  integer, parameter:: dp=c_double
!
  integer, intent(in), value:: ifirst, ilast, jfirst, jlast, kfirst, klast
  real(dp), intent(in), value:: h, zmin, omm, phm, amprho 
  real(dp), intent(out):: rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
  real(dp), intent(out):: cm(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real(dp), intent(in):: phc(21)

! local variables
  integer i, j, k
  real(dp) x, y, z

!$OMP PARALLEL PRIVATE(i,j,k,x,y,z)
!$OMP DO
  do k=kfirst,klast
     do j=jfirst,jlast
        do i=ifirst,ilast
           x=(i-1)*h
           y=(j-1)*h
           z=zmin+(k-1)*h
! note that the constant in the diagonal elements is 10 instead of 2           
           rho(i,j,k) = amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
           cm(1,i,j,k)= 10+sin(omm*x+phc(1))*cos(omm*y+phc(1))*sin(omm*z+phc(1))
           cm(2,i,j,k)= 2+sin(omm*x+phc(2))*cos(omm*y+phc(1))*sin(omm*z+phc(2))
           cm(3,i,j,k)= 2+sin(omm*x+phc(3))*cos(omm*y+phc(2))*sin(omm*z+phc(3))
           cm(4,i,j,k)= 2+sin(omm*x+phc(4))*cos(omm*y+phc(2))*sin(omm*z+phc(4))
           cm(5,i,j,k)= 2+sin(omm*x+phc(5))*cos(omm*y+phc(3))*sin(omm*z+phc(5))
           cm(6,i,j,k)= 2+sin(omm*x+phc(6))*cos(omm*y+phc(3))*sin(omm*z+phc(6))
           cm(7,i,j,k)= 10+sin(omm*x+phc(7))*cos(omm*y+phc(4))*sin(omm*z+phc(7))
           cm(8,i,j,k)= 2+sin(omm*x+phc(8))*cos(omm*y+phc(4))*sin(omm*z+phc(8))
           cm(9,i,j,k)= 2+sin(omm*x+phc(9))*cos(omm*y+phc(5))*sin(omm*z+phc(9))
           cm(10,i,j,k)=2+sin(omm*x+phc(10))*cos(omm*y+phc(5))*sin(omm*z+phc(1))
           cm(11,i,j,k)=2+sin(omm*x+phc(11))*cos(omm*y+phc(6))*sin(omm*z+phc(2))
           cm(12,i,j,k)=10+sin(omm*x+phc(12))*cos(omm*y+phc(6))*sin(omm*z+phc(3))
           cm(13,i,j,k)=2+sin(omm*x+phc(13))*cos(omm*y+phc(7))*sin(omm*z+phc(4))
           cm(14,i,j,k)=2+sin(omm*x+phc(14))*cos(omm*y+phc(7))*sin(omm*z+phc(5))
           cm(15,i,j,k)=2+sin(omm*x+phc(15))*cos(omm*y+phc(8))*sin(omm*z+phc(6))
           cm(16,i,j,k)=10+sin(omm*x+phc(16))*cos(omm*y+phc(8))*sin(omm*z+phc(7))
           cm(17,i,j,k)=2+sin(omm*x+phc(17))*cos(omm*y+phc(9))*sin(omm*z+phc(8))
           cm(18,i,j,k)=2+sin(omm*x+phc(18))*cos(omm*y+phc(9))*sin(omm*z+phc(9))
           cm(19,i,j,k)=10+sin(omm*x+phc(19))*cos(omm*y+phc(10))*sin(omm*z+phc(10))
           cm(20,i,j,k)=2+sin(omm*x+phc(20))*cos(omm*y+phc(10))*sin(omm*z+phc(11))
           cm(21,i,j,k)=10+sin(omm*x+phc(21))*cos(omm*y+phc(10))*sin(omm*z+phc(12))
        enddo
     enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL
end subroutine tw_ani_stiff


subroutine tw_ani_curvi_stiff( ifirst, ilast, jfirst, jlast, kfirst, klast, xx, yy, zz, omm, phm, amprho, rho, phc, cm) bind(c)
  use iso_c_binding
  implicit none
  integer, parameter:: dp=c_double
!
  integer, intent(in), value:: ifirst, ilast, jfirst, jlast, kfirst, klast
  real(dp), intent(in), value:: omm, phm, amprho 
  real(dp), intent(out):: rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
  real(dp), intent(in), dimension(ifirst:ilast,jfirst:jlast,kfirst:klast):: xx, yy, zz
  real(dp), intent(out):: cm(21,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real(dp), intent(in):: phc(21)

! local variables
  integer i, j, k
  real(dp) x, y, z

!$OMP PARALLEL PRIVATE(i,j,k,x,y,z)
!$OMP DO
  do k=kfirst,klast
     do j=jfirst,jlast
        do i=ifirst,ilast
           x=xx(i,j,k)
           y=yy(i,j,k)
           z=zz(i,j,k)
! note that the constant in the diagonal elements is 10 instead of 2           
           rho(i,j,k) = amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
           cm(1,i,j,k)= 10+sin(omm*x+phc(1))*cos(omm*y+phc(1))*sin(omm*z+phc(1))
           cm(2,i,j,k)= 2+sin(omm*x+phc(2))*cos(omm*y+phc(1))*sin(omm*z+phc(2))
           cm(3,i,j,k)= 2+sin(omm*x+phc(3))*cos(omm*y+phc(2))*sin(omm*z+phc(3))
           cm(4,i,j,k)= 2+sin(omm*x+phc(4))*cos(omm*y+phc(2))*sin(omm*z+phc(4))
           cm(5,i,j,k)= 2+sin(omm*x+phc(5))*cos(omm*y+phc(3))*sin(omm*z+phc(5))
           cm(6,i,j,k)= 2+sin(omm*x+phc(6))*cos(omm*y+phc(3))*sin(omm*z+phc(6))
           cm(7,i,j,k)= 10+sin(omm*x+phc(7))*cos(omm*y+phc(4))*sin(omm*z+phc(7))
           cm(8,i,j,k)= 2+sin(omm*x+phc(8))*cos(omm*y+phc(4))*sin(omm*z+phc(8))
           cm(9,i,j,k)= 2+sin(omm*x+phc(9))*cos(omm*y+phc(5))*sin(omm*z+phc(9))
           cm(10,i,j,k)=2+sin(omm*x+phc(10))*cos(omm*y+phc(5))*sin(omm*z+phc(1))
           cm(11,i,j,k)=2+sin(omm*x+phc(11))*cos(omm*y+phc(6))*sin(omm*z+phc(2))
           cm(12,i,j,k)=10+sin(omm*x+phc(12))*cos(omm*y+phc(6))*sin(omm*z+phc(3))
           cm(13,i,j,k)=2+sin(omm*x+phc(13))*cos(omm*y+phc(7))*sin(omm*z+phc(4))
           cm(14,i,j,k)=2+sin(omm*x+phc(14))*cos(omm*y+phc(7))*sin(omm*z+phc(5))
           cm(15,i,j,k)=2+sin(omm*x+phc(15))*cos(omm*y+phc(8))*sin(omm*z+phc(6))
           cm(16,i,j,k)=10+sin(omm*x+phc(16))*cos(omm*y+phc(8))*sin(omm*z+phc(7))
           cm(17,i,j,k)=2+sin(omm*x+phc(17))*cos(omm*y+phc(9))*sin(omm*z+phc(8))
           cm(18,i,j,k)=2+sin(omm*x+phc(18))*cos(omm*y+phc(9))*sin(omm*z+phc(9))
           cm(19,i,j,k)=10+sin(omm*x+phc(19))*cos(omm*y+phc(10))*sin(omm*z+phc(10))
           cm(20,i,j,k)=2+sin(omm*x+phc(20))*cos(omm*y+phc(10))*sin(omm*z+phc(11))
           cm(21,i,j,k)=10+sin(omm*x+phc(21))*cos(omm*y+phc(10))*sin(omm*z+phc(12))
        enddo
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine tw_ani_curvi_stiff
