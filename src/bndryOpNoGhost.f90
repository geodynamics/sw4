subroutine bndryOpNoGhost( acof_no_gp, ghcof_no_gp, sbop_no_gp) bind(c, name="bndryOpNoGhost")
  use iso_fortran_env
  implicit none
  integer, parameter:: dp=real64;
  
!** acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j)
!** ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.

  real(dp):: ghcof_no_gp(6), sbop_no_gp(0:5);
  real(dp):: acof_no_gp(6,8,8);
! local variables
  real(dp):: d5(0:8), w0;
  real(dp):: acof(6,8,8), ghcof(6);
  integer:: i, j, k;
  interface
      subroutine VARCOEFFS4( acof, ghcof ) bind(c)
      implicit none
      real*8 acof(6,8,8), ghcof(6)
      end subroutine VARCOEFFS4
  end interface

  call varcoeffs4( acof, ghcof );
!  print *, "Called varcoeffs4"
  ! modified coefficients for d(a(x) du): NOT using the ghost point
  acof_no_gp = acof
  d5 = 0;
  ! Use 5th divided difference to cancel the ghost point contribution
  d5(0) = -1.0;
  d5(1) = 5.0;
  d5(2) = -10.0;
  d5(3) = 10.0;
  d5(4) = -5.0;
  d5(5) = 1.0;
  w0 = 17.0/48.0_dp
  i = 1;
  k=1; ! only depends on the coefficient a(1)
  do j=1,8
     acof_no_gp(i,j,k) = acof(i,j,k) + d5(j)/(4*w0)
  end do

  ! the coeff for all ghost points are zero (don't divided by them!)
  ghcof_no_gp = 0;

  ! boundary normal derivative, not using ghost points
  ! sb = (-25*f(1)/12 + 4*f(2) - 3*f(3) + 4*f(4)/3 - f(5)/4)/h(q);
  sbop_no_gp(0) = 0;
  sbop_no_gp(1) = -25.0/12.0_dp;
  sbop_no_gp(2) = 4.0_dp;
  sbop_no_gp(3) = -3.0_dp;
  sbop_no_gp(4) = 4.0/3.0_dp;
  sbop_no_gp(5) = -1.0/4.0_dp;
  
end subroutine BNDRYOPNOGHOST
