subroutine dgetrf_wrap(i1,i2,A,i3,piv,info) bind(c)
  
  use iso_fortran_env
  use iso_c_binding
  !use problemsetup_new_3d
  implicit none
  integer, parameter :: dp=c_double
  integer :: i1,i2,i3,info
  real(dp):: A(1:3,1:3)
  integer:: piv(1:3)
  call dgetrf(i1,i2,A,i3,piv,info)
  if (info.ne.0) then
     write(*,*)"LU FAILED",A(INFO,INFO)
  endif
end subroutine dgetrf_wrap

subroutine dgetrs_wrap(i1,i2,A,i3,piv,rhs,i4,info) bind(c)
  
  use iso_fortran_env
  use iso_c_binding
  !use problemsetup_new_3d
  implicit none
  integer, parameter :: dp=c_double
  integer,value,intent(in) :: i1,i2,i3,i4
  integer :: info
  real(dp):: A(1:3,1:3)
  integer:: piv(1:3)
  integer,dimension(3) :: rhs
  call dgetrs('N',i1,i2,A,i3,piv,rhs,i4,info)
  if (info.ne.0) then
     write(*,*)"DGETRS FAILED",A(INFO,INFO)
  endif
end subroutine dgetrs_wrap

subroutine print_uf(u_f) bind(c)
!!$    use iso_fortran_env
!!$    use iso_c_binding
!!$    !use problemsetup_new_3d
!!$    implicit none
!!$    
!!$    real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
!!$    open (unit = 7, file = "compare.dat")
!!$    write(7,*)u_f
!!$    close(7)
  end subroutine print_uf
  
subroutine print_f(f) bind(c)
!!$    use iso_fortran_env
!!$    use iso_c_binding
!!$    use problemsetup_new_3d
!!$    implicit none
!!$    real(dp) ::  f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim)
!!$    open (unit = 7, file = "compare.dat")
!!$    write(7,*)f
!!$    close(7)
  end subroutine print_f
