subroutine interpolation_restriction(P,Rop,RPop) bind(c)
  use iso_fortran_env
  use iso_c_binding
  implicit none
  integer, parameter :: dp=c_double
  
  
  real(dp) :: RPop(-3:3), P(-1:2), Rop(-4:2)
  P(-1) = -1.d0/16.d0
  P(0) = 9.d0/16.d0
  P(1) = 9.d0/16.d0
  P(2) = -1.d0/16.d0
  
  Rop(-4) = -1.d0/32.d0
  Rop(-3) = 0.d0
  Rop(-2) = 9.d0/32.d0
  Rop(-1) = 1.d0/2.d0
  Rop(0) = 9.d0/32.d0
  Rop(1) = 0.d0
  Rop(2) = -1.d0/32.d0
  
  RPop(-3) = 1.d0/512.d0
  RPop(-2) = -9.d0/256.d0
  RPop(-1) = 63.d0/512.d0
  RPop(0)  = 105.d0/128.d0
  RPop(1)  = 63.d0/512.d0
  RPop(2)  = -9.d0/256.d0
  RPop(3)  = 1.d0/512.d0
end subroutine interpolation_restriction

subroutine central_difference(ux_cof,uxx_cof) bind(c)
  use iso_fortran_env
  use iso_c_binding
  implicit none
  integer, parameter :: dp=c_double
  
  real(dp) :: ux_cof(-2:2), uxx_cof(-2:2)
  ux_cof(-2) = 1.d0/12.d0
  ux_cof(-1) = -2.d0/3.d0
  ux_cof(0) = 0.d0
  ux_cof(1) = 2.d0/3.d0
  ux_cof(2) = -1.d0/12.d0
  
  uxx_cof(-2) = -1.d0/12.d0
  uxx_cof(-1) = 4.d0/3.d0
  uxx_cof(0) = -5.d0/2.d0
  uxx_cof(1) = 4.d0/3.d0
  uxx_cof(2) = -1.d0/12.d0
end subroutine central_difference

subroutine VARCOEFFS4(acof, ghcof, Sb) bind(c)
  use iso_fortran_env
  use iso_c_binding
  implicit none
  integer, parameter :: dp=c_double
    real(dp):: acof(6,8,8), ghcof(6), Sb(0:4)
    ! acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j)
    ! ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.
    ghcof(1) = 12.d0/17.d0
    ghcof(2) = 0
    ghcof(3) = 0
    ghcof(4) = 0
    ghcof(5) = 0
    ghcof(6) = 0
    acof(1,1,1) = 104.D0/289.D0
    acof(1,1,2) = -2476335.D0/2435692.D0
    acof(1,1,3) = -16189.D0/84966.D0
    acof(1,1,4) = -9.D0/3332.D0
    acof(1,1,5) = 0
    acof(1,1,6) = 0
    acof(1,1,7) = 0
    acof(1,1,8) = 0
    acof(1,2,1) = -516.D0/289.D0
    acof(1,2,2) = 544521.D0/1217846.D0
    acof(1,2,3) = 2509879.D0/3653538.D0
    acof(1,2,4) = 0
    acof(1,2,5) = 0
    acof(1,2,6) = 0
    acof(1,2,7) = 0
    acof(1,2,8) = 0
    acof(1,3,1) = 312.D0/289.D0
    acof(1,3,2) = 1024279.D0/2435692.D0
    acof(1,3,3) = -687797.D0/1217846.D0
    acof(1,3,4) = 177.D0/3332.D0
    acof(1,3,5) = 0
    acof(1,3,6) = 0
    acof(1,3,7) = 0
    acof(1,3,8) = 0
    acof(1,4,1) = -104.D0/289.D0
    acof(1,4,2) = 181507.D0/1217846.D0
    acof(1,4,3) = 241309.D0/3653538.D0
    acof(1,4,4) = 0
    acof(1,4,5) = 0
    acof(1,4,6) = 0
    acof(1,4,7) = 0
    acof(1,4,8) = 0
    acof(1,5,1) = 0
    acof(1,5,2) = 0
    acof(1,5,3) = 5.D0/2193.D0
    acof(1,5,4) = -48.D0/833.D0
    acof(1,5,5) = 0
    acof(1,5,6) = 0
    acof(1,5,7) = 0
    acof(1,5,8) = 0
    acof(1,6,1) = 0
    acof(1,6,2) = 0
    acof(1,6,3) = 0
    acof(1,6,4) = 6.D0/833.D0
    acof(1,6,5) = 0
    acof(1,6,6) = 0
    acof(1,6,7) = 0
    acof(1,6,8) = 0
    acof(1,7,1) = 0
    acof(1,7,2) = 0
    acof(1,7,3) = 0
    acof(1,7,4) = 0
    acof(1,7,5) = 0
    acof(1,7,6) = 0
    acof(1,7,7) = 0
    acof(1,7,8) = 0
    acof(1,8,1) = 0
    acof(1,8,2) = 0
    acof(1,8,3) = 0
    acof(1,8,4) = 0
    acof(1,8,5) = 0
    acof(1,8,6) = 0
    acof(1,8,7) = 0
    acof(1,8,8) = 0
    acof(2,1,1) = 12.D0/17.D0
    acof(2,1,2) = 544521.D0/4226642.D0
    acof(2,1,3) = 2509879.D0/12679926.D0
    acof(2,1,4) = 0
    acof(2,1,5) = 0
    acof(2,1,6) = 0
    acof(2,1,7) = 0
    acof(2,1,8) = 0
    acof(2,2,1) = -59.D0/68.D0
    acof(2,2,2) = -1633563.D0/4226642.D0
    acof(2,2,3) = -21510077.D0/25359852.D0
    acof(2,2,4) = -12655.D0/372939.D0
    acof(2,2,5) = 0
    acof(2,2,6) = 0
    acof(2,2,7) = 0
    acof(2,2,8) = 0
    acof(2,3,1) = 2.D0/17.D0
    acof(2,3,2) = 1633563.D0/4226642.D0
    acof(2,3,3) = 2565299.D0/4226642.D0
    acof(2,3,4) = 40072.D0/372939.D0
    acof(2,3,5) = 0
    acof(2,3,6) = 0
    acof(2,3,7) = 0
    acof(2,3,8) = 0
    acof(2,4,1) = 3.D0/68.D0
    acof(2,4,2) = -544521.D0/4226642.D0
    acof(2,4,3) = 987685.D0/25359852.D0
    acof(2,4,4) = -14762.D0/124313.D0
    acof(2,4,5) = 0
    acof(2,4,6) = 0
    acof(2,4,7) = 0
    acof(2,4,8) = 0
    acof(2,5,1) = 0
    acof(2,5,2) = 0
    acof(2,5,3) = 1630.D0/372939.D0
    acof(2,5,4) = 18976.D0/372939.D0
    acof(2,5,5) = 0
    acof(2,5,6) = 0
    acof(2,5,7) = 0
    acof(2,5,8) = 0
    acof(2,6,1) = 0
    acof(2,6,2) = 0
    acof(2,6,3) = 0
    acof(2,6,4) = -1.D0/177.D0
    acof(2,6,5) = 0
    acof(2,6,6) = 0
    acof(2,6,7) = 0
    acof(2,6,8) = 0
    acof(2,7,1) = 0
    acof(2,7,2) = 0
    acof(2,7,3) = 0
    acof(2,7,4) = 0
    acof(2,7,5) = 0
    acof(2,7,6) = 0
    acof(2,7,7) = 0
    acof(2,7,8) = 0
    acof(2,8,1) = 0
    acof(2,8,2) = 0
    acof(2,8,3) = 0
    acof(2,8,4) = 0
    acof(2,8,5) = 0
    acof(2,8,6) = 0
    acof(2,8,7) = 0
    acof(2,8,8) = 0
    acof(3,1,1) = -96.D0/731.D0
    acof(3,1,2) = 1024279.D0/6160868.D0
    acof(3,1,3) = -687797.D0/3080434.D0
    acof(3,1,4) = 177.D0/8428.D0
    acof(3,1,5) = 0
    acof(3,1,6) = 0
    acof(3,1,7) = 0
    acof(3,1,8) = 0
    acof(3,2,1) = 118.D0/731.D0
    acof(3,2,2) = 1633563.D0/3080434.D0
    acof(3,2,3) = 2565299.D0/3080434.D0
    acof(3,2,4) = 40072.D0/271803.D0
    acof(3,2,5) = 0
    acof(3,2,6) = 0
    acof(3,2,7) = 0
    acof(3,2,8) = 0
    acof(3,3,1) = -16.D0/731.D0
    acof(3,3,2) = -5380447.D0/6160868.D0
    acof(3,3,3) = -3569115.D0/3080434.D0
    acof(3,3,4) = -331815.D0/362404.D0
    acof(3,3,5) = -283.D0/6321.D0
    acof(3,3,6) = 0
    acof(3,3,7) = 0
    acof(3,3,8) = 0
    acof(3,4,1) = -6.D0/731.D0
    acof(3,4,2) = 544521.D0/3080434.D0
    acof(3,4,3) = 2193521.D0/3080434.D0
    acof(3,4,4) = 8065.D0/12943.D0
    acof(3,4,5) = 381.D0/2107.D0
    acof(3,4,6) = 0
    acof(3,4,7) = 0
    acof(3,4,8) = 0
    acof(3,5,1) = 0
    acof(3,5,2) = 0
    acof(3,5,3) = -14762.D0/90601.D0
    acof(3,5,4) = 32555.D0/271803.D0
    acof(3,5,5) = -283.D0/2107.D0
    acof(3,5,6) = 0
    acof(3,5,7) = 0
    acof(3,5,8) = 0
    acof(3,6,1) = 0
    acof(3,6,2) = 0
    acof(3,6,3) = 0
    acof(3,6,4) = 9.D0/2107.D0
    acof(3,6,5) = -11.D0/6321.D0
    acof(3,6,6) = 0
    acof(3,6,7) = 0
    acof(3,6,8) = 0
    acof(3,7,1) = 0
    acof(3,7,2) = 0
    acof(3,7,3) = 0
    acof(3,7,4) = 0
    acof(3,7,5) = 0
    acof(3,7,6) = 0
    acof(3,7,7) = 0
    acof(3,7,8) = 0
    acof(3,8,1) = 0
    acof(3,8,2) = 0
    acof(3,8,3) = 0
    acof(3,8,4) = 0
    acof(3,8,5) = 0
    acof(3,8,6) = 0
    acof(3,8,7) = 0
    acof(3,8,8) = 0
    acof(4,1,1) = -36.D0/833.D0
    acof(4,1,2) = 181507.D0/3510262.D0
    acof(4,1,3) = 241309.D0/10530786.D0
    acof(4,1,4) = 0
    acof(4,1,5) = 0
    acof(4,1,6) = 0
    acof(4,1,7) = 0
    acof(4,1,8) = 0
    acof(4,2,1) = 177.D0/3332.D0
    acof(4,2,2) = -544521.D0/3510262.D0
    acof(4,2,3) = 987685.D0/21061572.D0
    acof(4,2,4) = -14762.D0/103243.D0
    acof(4,2,5) = 0
    acof(4,2,6) = 0
    acof(4,2,7) = 0
    acof(4,2,8) = 0
    acof(4,3,1) = -6.D0/833.D0
    acof(4,3,2) = 544521.D0/3510262.D0
    acof(4,3,3) = 2193521.D0/3510262.D0
    acof(4,3,4) = 8065.D0/14749.D0
    acof(4,3,5) = 381.D0/2401.D0
    acof(4,3,6) = 0
    acof(4,3,7) = 0
    acof(4,3,8) = 0
    acof(4,4,1) = -9.D0/3332.D0
    acof(4,4,2) = -181507.D0/3510262.D0
    acof(4,4,3) = -2647979.D0/3008796.D0
    acof(4,4,4) = -80793.D0/103243.D0
    acof(4,4,5) = -1927.D0/2401.D0
    acof(4,4,6) = -2.D0/49.D0
    acof(4,4,7) = 0
    acof(4,4,8) = 0
    acof(4,5,1) = 0
    acof(4,5,2) = 0
    acof(4,5,3) = 57418.D0/309729.D0
    acof(4,5,4) = 51269.D0/103243.D0
    acof(4,5,5) = 1143.D0/2401.D0
    acof(4,5,6) = 8.D0/49.D0
    acof(4,5,7) = 0
    acof(4,5,8) = 0
    acof(4,6,1) = 0
    acof(4,6,2) = 0
    acof(4,6,3) = 0
    acof(4,6,4) = -283.D0/2401.D0
    acof(4,6,5) = 403.D0/2401.D0
    acof(4,6,6) = -6.D0/49.D0
    acof(4,6,7) = 0
    acof(4,6,8) = 0
    acof(4,7,1) = 0
    acof(4,7,2) = 0
    acof(4,7,3) = 0
    acof(4,7,4) = 0
    acof(4,7,5) = 0
    acof(4,7,6) = 0
    acof(4,7,7) = 0
    acof(4,7,8) = 0
    acof(4,8,1) = 0
    acof(4,8,2) = 0
    acof(4,8,3) = 0
    acof(4,8,4) = 0
    acof(4,8,5) = 0
    acof(4,8,6) = 0
    acof(4,8,7) = 0
    acof(4,8,8) = 0
    acof(5,1,1) = 0
    acof(5,1,2) = 0
    acof(5,1,3) = 5.D0/6192.D0
    acof(5,1,4) = -1.D0/49.D0
    acof(5,1,5) = 0
    acof(5,1,6) = 0
    acof(5,1,7) = 0
    acof(5,1,8) = 0
    acof(5,2,1) = 0
    acof(5,2,2) = 0
    acof(5,2,3) = 815.D0/151704.D0
    acof(5,2,4) = 1186.D0/18963.D0
    acof(5,2,5) = 0
    acof(5,2,6) = 0
    acof(5,2,7) = 0
    acof(5,2,8) = 0
    acof(5,3,1) = 0
    acof(5,3,2) = 0
    acof(5,3,3) = -7381.D0/50568.D0
    acof(5,3,4) = 32555.D0/303408.D0
    acof(5,3,5) = -283.D0/2352.D0
    acof(5,3,6) = 0
    acof(5,3,7) = 0
    acof(5,3,8) = 0
    acof(5,4,1) = 0
    acof(5,4,2) = 0
    acof(5,4,3) = 28709.D0/151704.D0
    acof(5,4,4) = 51269.D0/101136.D0
    acof(5,4,5) = 381.D0/784.D0
    acof(5,4,6) = 1.D0/6.D0
    acof(5,4,7) = 0
    acof(5,4,8) = 0
    acof(5,5,1) = 0
    acof(5,5,2) = 0
    acof(5,5,3) = -349.D0/7056.D0
    acof(5,5,4) = -247951.D0/303408.D0
    acof(5,5,5) = -577.D0/784.D0
    acof(5,5,6) = -5.D0/6.D0
    acof(5,5,7) = -1.D0/24.D0
    acof(5,5,8) = 0
    acof(5,6,1) = 0
    acof(5,6,2) = 0
    acof(5,6,3) = 0
    acof(5,6,4) = 1135.D0/7056.D0
    acof(5,6,5) = 1165.D0/2352.D0
    acof(5,6,6) = 1.D0/2.D0
    acof(5,6,7) = 1.D0/6.D0
    acof(5,6,8) = 0
    acof(5,7,1) = 0
    acof(5,7,2) = 0
    acof(5,7,3) = 0
    acof(5,7,4) = 0
    acof(5,7,5) = -1.D0/8.D0
    acof(5,7,6) = 1.D0/6.D0
    acof(5,7,7) = -1.D0/8.D0
    acof(5,7,8) = 0
    acof(5,8,1) = 0
    acof(5,8,2) = 0
    acof(5,8,3) = 0
    acof(5,8,4) = 0
    acof(5,8,5) = 0
    acof(5,8,6) = 0
    acof(5,8,7) = 0
    acof(5,8,8) = 0
    acof(6,1,1) = 0
    acof(6,1,2) = 0
    acof(6,1,3) = 0
    acof(6,1,4) = 1.D0/392.D0
    acof(6,1,5) = 0
    acof(6,1,6) = 0
    acof(6,1,7) = 0
    acof(6,1,8) = 0
    acof(6,2,1) = 0
    acof(6,2,2) = 0
    acof(6,2,3) = 0
    acof(6,2,4) = -1.D0/144.D0
    acof(6,2,5) = 0
    acof(6,2,6) = 0
    acof(6,2,7) = 0
    acof(6,2,8) = 0
    acof(6,3,1) = 0
    acof(6,3,2) = 0
    acof(6,3,3) = 0
    acof(6,3,4) = 3.D0/784.D0
    acof(6,3,5) = -11.D0/7056.D0
    acof(6,3,6) = 0
    acof(6,3,7) = 0
    acof(6,3,8) = 0
    acof(6,4,1) = 0
    acof(6,4,2) = 0
    acof(6,4,3) = 0
    acof(6,4,4) = -283.D0/2352.D0
    acof(6,4,5) = 403.D0/2352.D0
    acof(6,4,6) = -1.D0/8.D0
    acof(6,4,7) = 0
    acof(6,4,8) = 0
    acof(6,5,1) = 0
    acof(6,5,2) = 0
    acof(6,5,3) = 0
    acof(6,5,4) = 1135.D0/7056.D0
    acof(6,5,5) = 1165.D0/2352.D0
    acof(6,5,6) = 1.D0/2.D0
    acof(6,5,7) = 1.D0/6.D0
    acof(6,5,8) = 0
    acof(6,6,1) = 0
    acof(6,6,2) = 0
    acof(6,6,3) = 0
    acof(6,6,4) = -47.D0/1176.D0
    acof(6,6,5) = -5869.D0/7056.D0
    acof(6,6,6) = -3.D0/4.D0
    acof(6,6,7) = -5.D0/6.D0
    acof(6,6,8) = -1.D0/24.D0
    acof(6,7,1) = 0
    acof(6,7,2) = 0
    acof(6,7,3) = 0
    acof(6,7,4) = 0
    acof(6,7,5) = 1.D0/6.D0
    acof(6,7,6) = 1.D0/2.D0
    acof(6,7,7) = 1.D0/2.D0
    acof(6,7,8) = 1.D0/6.D0
    acof(6,8,1) = 0
    acof(6,8,2) = 0
    acof(6,8,3) = 0
    acof(6,8,4) = 0
    acof(6,8,5) = 0
    acof(6,8,6) = -1.D0/8.D0
    acof(6,8,7) = 1.D0/6.D0
    acof(6,8,8) = -1.D0/8.D0
    ! 129 non-zero out of 384.

  Sb(0) = -3.d0/12.d0
  Sb(1) = -10.d0/12.d0
  Sb(2) = 18.d0/12.d0
  Sb(3) = -6.d0/12.d0
  Sb(4) = 1.d0/12.d0
  end subroutine VARCOEFFS4


subroutine varcoeff_NoGhost( acof_no_gp, ghcof_no_gp, sbop_no_gp ) bind(c)
  use iso_fortran_env
  use iso_c_binding
  implicit none
  integer, parameter :: dp=c_double
  
  real(dp):: ghcof_no_gp(6), sbop_no_gp(0:5)
  real(dp):: acof_no_gp(6,8,8)
  ! local variables
  real(dp):: d5(0:8), w0;
  real(dp):: acof(6,8,8), ghcof(6)
  integer:: i, j, k;
  call VARCOEFFS4( acof, ghcof,sbop_no_gp(0:4) )
  acof_no_gp = acof
  d5 = 0.d0;
  ! Use 5th divided difference to cancel the ghost point contribution
  d5(0) = -1.0d0;
  d5(1) = 5.0d0;
  d5(2) = -10.0d0;
  d5(3) = 10.0d0;
  d5(4) = -5.0d0;
  d5(5) = 1.0d0;
  w0 = 17.0/48.0d0
  i = 1;
  k=1; ! only depends on the coefficient a(1)
  do j=1,8
     acof_no_gp(i,j,k) = acof(i,j,k) + d5(j)/(4.d0*w0)
  end do
  
  ! the coeff for all ghost points are zero (don't divided by them!)
  ghcof_no_gp = 0.d0;
  
  ! boundary normal derivative, not using ghost points
  ! sb = (-25*f(1)/12 + 4*f(2) - 3*f(3) + 4*f(4)/3 - f(5)/4)/h(q);
  sbop_no_gp(0) = 0.d0;
  sbop_no_gp(1) = -25.d0/12.d0;
  sbop_no_gp(2) = 4.d0;
  sbop_no_gp(3) = -3.d0;
  sbop_no_gp(4) = 4.0/3.d0;
  sbop_no_gp(5) = -1.0/4.d0;
end subroutine varcoeff_NoGhost

subroutine dx_46(bop) bind(c)
  use iso_fortran_env
  use iso_c_binding
  implicit none
  integer, parameter :: dp=c_double
    real(dp):: bop(4,6)
    bop(1,1) = -24.D0/17.D0
    bop(1,2) = 59.D0/34.D0
    bop(1,3) = -4.D0/17.D0
    bop(1,4) = -3.D0/34.D0
    bop(1,5) = 0
    bop(1,6) = 0
    bop(2,1) = -1.D0/2.D0
    bop(2,2) = 0
    bop(2,3) = 1.D0/2.D0
    bop(2,4) = 0
    bop(2,5) = 0
    bop(2,6) = 0
    bop(3,1) = 4.D0/43.D0
    bop(3,2) = -59.D0/86.D0
    bop(3,3) = 0
    bop(3,4) = 59.D0/86.D0
    bop(3,5) = -4.D0/43.D0
    bop(3,6) = 0
    bop(4,1) = 3.D0/98.D0
    bop(4,2) = 0
    bop(4,3) = -59.D0/98.D0
    bop(4,4) = 0
    bop(4,5) = 32.D0/49.D0
    bop(4,6) = -4.D0/49.D0
  end subroutine dx_46




  
  subroutine equation_cof(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3, &
       XI13_c,XI23_c,XI33_c,Jacobian_c,rho_c,rho_f,XI13_f,XI23_f,XI33_f,Jacobian_f,&
       mu_c,mu_f,lambda_c,lambda_f) bind(c)
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
!!$    integer, parameter :: dp=c_double
!!$    integer :: nrg,n1_c,n2_c,n3_c,n1_f,n2_f,n3_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: rho_c,mu_c,lambda_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: rho_f,mu_f,lambda_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: Xgrid_c_3
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: Xgrid_f_3
    real(dp), dimension (1-nrg:n1_c+nrg) :: Xgrid_c_1
    real(dp), dimension (1-nrg:n1_f+nrg) :: Xgrid_f_1
    real(dp), dimension (1-nrg:n2_c+nrg) :: Xgrid_c_2
    real(dp), dimension (1-nrg:n2_f+nrg) :: Xgrid_f_2
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI13_c,XI23_c,XI33_c,Jacobian_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: XI13_f,XI23_f,XI33_f,Jacobian_f
    integer:: i,j,k

    call generate_grid(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3)


    ! compute metric derivatives, coarse
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             call metric_derivative(dble(k-1)/(n1_c-1),dble(j-1)/(n2_c-1),dble(i-1)/(n3_c-1),&
                 XI13_c(k,j,i),XI23_c(k,j,i),XI33_c(k,j,i),Jacobian_c(k,j,i),0)
          end do
       end do
    end do
    ! compute metric derivatives, fine
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
              call metric_derivative(dble(k-1)/(n1_f-1),dble(j-1)/(n2_f-1),dble(i-1)/(n3_f-1),&
                  XI13_f(k,j,i),XI23_f(k,j,i),XI33_f(k,j,i),Jacobian_f(k,j,i),1)
          end do
       end do
    end do
    ! variable coefficients
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             mu_c(k,j,i) = 3.d0 + sin(3.d0*Xgrid_c_1(k)+0.1d0)*sin(3.d0*Xgrid_c_2(j)+0.1d0)*sin(Xgrid_c_3(k,j,i))
             lambda_c(k,j,i) = 21.d0+ cos(Xgrid_c_1(k)+0.1d0)*cos(Xgrid_c_2(j)+0.1d0)*sin(3.d0*Xgrid_c_3(k,j,i))**2
             rho_c(k,j,i) = 2.d0 + sin(Xgrid_c_1(k)+0.3d0)*sin(Xgrid_c_2(j)+0.3d0)*sin(Xgrid_c_3(k,j,i)-0.2d0)
          end do
       end do
    end do
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
             mu_f(k,j,i) = 3.d0 + sin(3.d0*Xgrid_f_1(k)+0.1d0)*sin(3.d0*Xgrid_f_2(j)+0.1d0)*sin(Xgrid_f_3(k,j,i))
             lambda_f(k,j,i) = 21.d0+ cos(Xgrid_f_1(k)+0.1d0)*cos(Xgrid_f_2(j)+0.1d0)*sin(3.d0*Xgrid_f_3(k,j,i))**2
             rho_f(k,j,i) = 2.d0 + sin(Xgrid_f_1(k)+0.3d0)*sin(Xgrid_f_2(j)+0.3d0)*sin(Xgrid_f_3(k,j,i)-0.2d0)
          end do
       end do
    end do

    rho_c = rho_c*Jacobian_c
    rho_f = rho_f*Jacobian_f

  end subroutine  equation_cof
 
  subroutine kappa(mu_f,lambda_f,XI13_f,XI23_f,XI33_f,Jacobian_f,rho_f,kappa2) bind(c)

    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none

    
    real(dp), dimension (3) :: kappa1
    real(dp), dimension (3,3) :: mat_t
    real(dp), dimension (8) :: work
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: rho_f,mu_f,lambda_f
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: XI13_f,XI23_f,XI33_f,Jacobian_f
    integer :: info
    real(dp) :: kappa2
    integer :: i,j,k
    
    kappa2 = 0.d0
    do k = 1, n3_f
       do j = 1, n2_f
          do i = 1,n1_f
            mat_t(1,1) = (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))/l1**2
            mat_t(1,2) = 0.d0
            mat_t(1,3) = (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))/l1*XI13_f(i,j,k)
            mat_t(2,1) = 0.d0
            mat_t(2,2) = (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))/l2**2
            mat_t(2,3) = (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))/l2*XI23_f(i,j,k)
            mat_t(3,1) = mat_t(1,3)
            mat_t(3,2) = mat_t(2,3)
            mat_t(3,3) = (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))*XI13_f(i,j,k)*XI13_f(i,j,k) &
                       + (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))*XI23_f(i,j,k)*XI23_f(i,j,k) &
                       + (4.d0*mu_f(i,j,k)+lambda_f(i,j,k))*XI33_f(i,j,k)*XI33_f(i,j,k)
            mat_t = Jacobian_f(i,j,k)*mat_t/rho_f(i,j,k)
            !print *,"mat_t before ",mat_t
            !print *," inputs ",Jacobian_f(i,j,k),rho_f(i,j,k),l1
            call dsyev('N','U',3,mat_t,3,kappa1,work,8,info)
            !print *,"kappa",kappa1
            !print *,"mat_t after ",mat_t
             if (kappa1(3) > kappa2) then
                kappa2 = kappa1(3)
             end if
          end do
       end do
    end do
  end subroutine kappa
  !
  !exact solution
  ! We want the time-term to be the same in u1,u2 and u3, so that we don't need to call this function
  ! in every time step. We just need to scale the initial condition.
  subroutine exact_solution(x1,x2,x3,t,u1,u2,u3,flag) bind(c)
    use iso_fortran_env
    use iso_c_binding
 
    implicit none
    integer, parameter :: dp=c_double
    real(dp) :: x1, x2, x3, t, u1, u2, u3
    integer flag
    u1 = cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    u2 = sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    u3 = sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
  end subroutine exact_solution


  
  subroutine Interface_block(Rop,ghcof,Sb,rho_c,lambda_c,rho_f,Jacobian_c,mu_c,P,XI13_c,XI23_c,XI33_c,Mass_block) bind(c)
    
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
    
    
    !
    real(dp), dimension (1:3,1:3,1:n1_c,1:n2_c) :: Mass_block
    real(dp) :: int_cof
    real(dp), dimension (-4:2):: Rop
    real(dp), dimension (-1:2) :: P
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: rho_f
    real(dp),dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: rho_c,lambda_c,mu_c,Jacobian_c,XI13_c,XI23_c,XI33_c
    real(dp):: ghcof(6)
    real(dp), dimension (0:4):: Sb
    integer :: i,j,k,l
    !
    int_cof = 17.d0/48.d0*h3_f*ghcof(1)/h3_c**2
    !
    Mass_block = 0.d0
    do l = 1,n2_c
       do k = 1,n1_c
          !
          do j = -4,2,2
             do i = -4,2,2
                ! first set equation w.r.t the first component
                Mass_block(1,1,k,l) = Mass_block(1,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
                  *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
                ! first set equation w.r.t the second component
                Mass_block(1,2,k,l) = Mass_block(1,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! first set equation w.r.t the third component
                Mass_block(1,3,k,l) = Mass_block(1,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! second set equation w.r.t the first component
                Mass_block(2,1,k,l) = Mass_block(2,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! second set equation w.r.t the second component
                Mass_block(2,2,k,l) = Mass_block(2,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
                  *(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
                ! second set equation w.r.t the third component
                Mass_block(2,3,k,l) = Mass_block(2,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                  *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! third set equation w.r.t the first component
                Mass_block(3,1,k,l) = Mass_block(3,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! third set equation w.r.t the second component
                Mass_block(3,2,k,l) = Mass_block(3,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
                ! third set equation w.r.t the third component
                Mass_block(3,3,k,l) = Mass_block(3,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
                   *(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
             end do
          end do
          !
          do j = -4,2,2
             ! first set equation w.r.t the first component
             Mass_block(1,1,k,l) = Mass_block(1,1,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
             ! first set equation w.r.t the second component
             Mass_block(1,2,k,l) = Mass_block(1,2,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! first set equation w.r.t the third component
             Mass_block(1,3,k,l) = Mass_block(1,3,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the first component
             Mass_block(2,1,k,l) = Mass_block(2,1,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the second component
             Mass_block(2,2,k,l) = Mass_block(2,2,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the third component
             Mass_block(2,3,k,l) = Mass_block(2,3,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the first component
             Mass_block(3,1,k,l) = Mass_block(3,1,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the second component
             Mass_block(3,2,k,l) = Mass_block(3,2,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the third component
             Mass_block(3,3,k,l) = Mass_block(3,3,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(-j/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
               +XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
          end do
          !
          do i = -4,2,2
             ! first set equation w.r.t the first component
             Mass_block(1,1,k,l) = Mass_block(1,1,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
             ! first set equation w.r.t the second component
             Mass_block(1,2,k,l) = Mass_block(1,2,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! first set equation w.r.t the third component
             Mass_block(1,3,k,l) = Mass_block(1,3,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the first component
             Mass_block(2,1,k,l) = Mass_block(2,1,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the second component
             Mass_block(2,2,k,l) = Mass_block(2,2,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
               +XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
             ! second set equation w.r.t the third component
             Mass_block(2,3,k,l) = Mass_block(2,3,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the first component
             Mass_block(3,1,k,l) = Mass_block(3,1,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the second component
             Mass_block(3,2,k,l) = Mass_block(3,2,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
             ! third set equation w.r.t the third component
             Mass_block(3,3,k,l) = Mass_block(3,3,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
          end do
          !
          ! first set equation w.r.t the first component
          Mass_block(1,1,k,l) = Mass_block(1,1,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
            *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
          ! first set equation w.r.t the second component
          Mass_block(1,2,k,l) = Mass_block(1,2,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! first set equation w.r.t the third component
          Mass_block(1,3,k,l) = Mass_block(1,3,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! second set equation w.r.t the first component
          Mass_block(2,1,k,l) = Mass_block(2,1,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! second set equation w.r.t the second component
          Mass_block(2,2,k,l) = Mass_block(2,2,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
            *(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
          ! second set equation w.r.t the third component
          Mass_block(2,3,k,l) = Mass_block(2,3,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! third set equation w.r.t the first component
          Mass_block(3,1,k,l) = Mass_block(3,1,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! third set equation w.r.t the second component
          Mass_block(3,2,k,l) = Mass_block(3,2,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof
          ! third set equation w.r.t the third component
          Mass_block(3,3,k,l) = Mass_block(3,3,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
            *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
            *(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof
          ! from the norm derivative
          ! first set equation w.r.t the first component
          Mass_block(1,1,k,l) = Mass_block(1,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
            *XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c
          ! first set equation w.r.t the second component
          Mass_block(1,2,k,l) = Mass_block(1,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c
          ! first set equation w.r.t the third component
          Mass_block(1,3,k,l) = Mass_block(1,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c
          ! second set equation w.r.t the first component
          Mass_block(2,1,k,l) = Mass_block(2,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c
          ! second set equation w.r.t the second component
          Mass_block(2,2,k,l) = Mass_block(2,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
            *XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c
          ! second set equation w.r.t the third component
          Mass_block(2,3,k,l) = Mass_block(2,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c
          ! third set equation w.r.t the first component
          Mass_block(3,1,k,l) = Mass_block(3,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c
          ! third set equation w.r.t the second component
          Mass_block(3,2,k,l) = Mass_block(3,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
            *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c
          ! third set equation w.r.t the third component
          Mass_block(3,3,k,l) = Mass_block(3,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
            *XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/h3_c
       end do
    end do
!!$    print *,"M_CALC",Mass_block(:,:,1,1)
    !open (unit = 7, file = "compare.dat")
    !write(7,*)u_f
    !write(7,*)mu_f
    !close(7)
  !
end subroutine Interface_block

subroutine dgetrf_wrap(i1,i2,A,i3,piv,info) bind(c)
  
  use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
  !integer, parameter :: dp=c_double
  integer :: i1,i2,i3,info
  real(dp):: A(1:3,1:3)
  integer:: piv(1:3)
  call dgetrf(i1,i2,A,i3,piv,info)
  if (info.ne.0) then
     write(*,*)"LU FAILED",A(INFO,INFO)
  endif
end subroutine dgetrf_wrap

subroutine Injection(u_f,u_c,P,index) bind(c)
   use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
  real(dp) :: u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
  real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
  real(dp), dimension (-1:2) :: P
  integer,value,intent(in)::  index
  integer :: i,j
  ! Injection at the interface
  !write(*,*)"Index in injection is ",index
  do j = 1, n2_c
     do i = 1, n1_c
        u_f(2*i-1,2*j-1,1,:,index) = u_c(i,j,n3_c,:,index)
        u_f(2*i,2*j-1,1,:,index) = P(-1)*u_c(i-1,j,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
             +P(1)*u_c(i+1,j,n3_c,:,index)+P(2)*u_c(i+2,j,n3_c,:,index)
        u_f(2*i-1,2*j,1,:,index) = P(-1)*u_c(i,j-1,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
             +P(1)*u_c(i,j+1,n3_c,:,index)+P(2)*u_c(i,j+2,n3_c,:,index)
        u_f(2*i,2*j,1,:,index) = P(-1)*(P(-1)*u_c(i-1,j-1,n3_c,:,index)+P(0)*u_c(i,j-1,n3_c,:,index)&
             +P(1)*u_c(i+1,j-1,n3_c,:,index)+P(2)*u_c(i+2,j-1,n3_c,:,index))&
             +P(0)*(P(-1)*u_c(i-1,j,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
             +P(1)*u_c(i+1,j,n3_c,:,index)+P(2)*u_c(i+2,j,n3_c,:,index))&
             +P(1)*(P(-1)*u_c(i-1,j+1,n3_c,:,index)+P(0)*u_c(i,j+1,n3_c,:,index)&
             +P(1)*u_c(i+1,j+1,n3_c,:,index)+P(2)*u_c(i+2,j+1,n3_c,:,index))&
             +P(2)*(P(-1)*u_c(i-1,j+2,n3_c,:,index)+P(0)*u_c(i,j+2,n3_c,:,index)&
             +P(1)*u_c(i+1,j+2,n3_c,:,index)+P(2)*u_c(i+2,j+2,n3_c,:,index))
     end do
  end do
!!$  open (unit = 7, file = "compare.dat")
!!$  write(7,*)u_f
!!$  !write(7,*)mu_f
!!$  close(7)
end subroutine Injection


subroutine Interface_RHS(Vass,lh_c,lh_f,Jacobian_c,Jacobian_f, mu_c,mu_f,lambda_c,lambda_f,rho_c,rho_f,&
  XI13_c,XI23_c,XI33_c,XI13_f,XI23_f,XI33_f,P,Sb,Rop,sbop_no_gp,acof_no_gp,u_c,u_f,Mass_f1,ux_cof,ghcof,acof,bof,index) bind(c)
  use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
  !integer, parameter :: dp=c_double
  integer,value,intent(in) ::  index
  real(dp) :: Vass(1:n1_c*n2_c*3)
  real(dp) :: lh_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim)
  real(dp) :: lh_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim)
  real(dp) :: Jacobian_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp) :: mu_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp) :: lambda_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp),dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI13_c,XI23_c,XI33_c,rho_c
  real(dp), dimension (0:4):: Sb
  real(dp) :: u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
  real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
  real(dp) :: Mass_f1(-2:n1_f+3,-2:n2_f+3),ux_cof(-2:2),ghcof(6)
  real(dp) :: acof(6,8,8), bof(4,6),sbop_no_gp(0:5),acof_no_gp(6,8,8)
  !real(dp) :: rho_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)
  real(dp) ,dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: lambda_f,mu_f,rho_f,Jacobian_f,XI13_f,XI23_f,XI33_f
  real(dp), dimension (-1:2) :: P
  real(dp), dimension (-4:2):: Rop
   
  integer :: i,j,k,k1,l,m
    !
    Vass = 0.d0
    lh_c = 0.d0
    lh_f = 0.d0
    !
    ! term 1
    do k = 1,n2_c
       do i = 1,n1_c
          do j = 1,4
             ! 33
             ! first set equation
             Vass((k-1)*3*n1_c+(i-1)*3+1) = Vass((k-1)*3*n1_c+(i-1)*3+1)+Jacobian_c(i,k,n3_c)*((2.d0*mu_c(i,k,n3_c)&
               +lambda_c(i,k,n3_c))*XI13_c(i,k,n3_c)**2+mu_c(i,k,n3_c)*(XI23_c(i,k,n3_c)**2+XI33_c(i,k,n3_c)**2))*Sb(j)&
               *u_c(i,k,n3_c+1-j,1,index)/h3_c+Jacobian_c(i,k,n3_c)*(lambda_c(i,k,n3_c)+mu_c(i,k,n3_c))*XI13_c(i,k,n3_c)&
               *XI23_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,2,index)/h3_c+Jacobian_c(i,k,n3_c)*(lambda_c(i,k,n3_c)+mu_c(i,k,n3_c))&
               *XI13_c(i,k,n3_c)*XI33_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,3,index)/h3_c
             ! second set equation
             Vass((k-1)*3*n1_c+(i-1)*3+2) = Vass((k-1)*3*n1_c+(i-1)*3+2)+Jacobian_c(i,k,n3_c)*(lambda_c(i,k,n3_c)&
               +mu_c(i,k,n3_c))*XI13_c(i,k,n3_c)*XI23_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,1,index)/h3_c+Jacobian_c(i,k,n3_c)&
               *((2.d0*mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))*XI23_c(i,k,n3_c)**2+mu_c(i,k,n3_c)*(XI13_c(i,k,n3_c)**2&
               +XI33_c(i,k,n3_c)**2))*Sb(j)*u_c(i,k,n3_c+1-j,2,index)/h3_c+Jacobian_c(i,k,n3_c)*(lambda_c(i,k,n3_c)&
               +mu_c(i,k,n3_c))*XI23_c(i,k,n3_c)*XI33_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,3,index)/h3_c
             ! third set equation
             Vass((k-1)*3*n1_c+(i-1)*3+3) = Vass((k-1)*3*n1_c+(i-1)*3+3)+Jacobian_c(i,k,n3_c)*(lambda_c(i,k,n3_c)&
               +mu_c(i,k,n3_c))*XI13_c(i,k,n3_c)*XI33_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,1,index)/h3_c+Jacobian_c(i,k,n3_c)&
               *(lambda_c(i,k,n3_c)+mu_c(i,k,n3_c))*XI23_c(i,k,n3_c)*XI33_c(i,k,n3_c)*Sb(j)*u_c(i,k,n3_c+1-j,2,index)/h3_c&
               +Jacobian_c(i,k,n3_c)*((2.d0*mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))*XI33_c(i,k,n3_c)**2&
               +mu_c(i,k,n3_c)*(XI13_c(i,k,n3_c)**2+XI23_c(i,k,n3_c)**2))*Sb(j)*u_c(i,k,n3_c+1-j,3,index)/h3_c
          end do
       end do
    end do
    !
    do k = 1,n2_c
       do i = 1,n1_c
          do j = -2,2
             ! 31
             ! first set equation
             Vass((k-1)*3*n1_c+(i-1)*3+1) = Vass((k-1)*3*n1_c+(i-1)*3+1)-Jacobian_c(i,k,n3_c)*(2.d0*mu_c(i,k,n3_c)&
               +lambda_c(i,k,n3_c))/l1*XI13_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,1,index)/h1_c-Jacobian_c(i,k,n3_c)&
               *mu_c(i,k,n3_c)/l1*XI23_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,2,index)/h1_c-Jacobian_c(i,k,n3_c)&
               *mu_c(i,k,n3_c)/l1*XI33_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,3,index)/h1_c
             ! second set equation
             Vass((k-1)*3*n1_c+(i-1)*3+2) = Vass((k-1)*3*n1_c+(i-1)*3+2)-Jacobian_c(i,k,n3_c)*lambda_c(i,k,n3_c)/l1&
               *XI23_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,1,index)/h1_c-Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l1&
               *XI13_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,2,index)/h1_c
             ! third set equation
             Vass((k-1)*3*n1_c+(i-1)*3+3) = Vass((k-1)*3*n1_c+(i-1)*3+3)-Jacobian_c(i,k,n3_c)*lambda_c(i,k,n3_c)/l1&
               *XI33_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,1,index)/h1_c-Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l1&
               *XI13_c(i,k,n3_c)*ux_cof(j)*u_c(i+j,k,n3_c,3,index)/h1_c
             ! 32
             ! first set equation
             Vass((k-1)*3*n1_c+(i-1)*3+1) = Vass((k-1)*3*n1_c+(i-1)*3+1)-Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l2&
               *XI23_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,1,index)/h2_c&
               -Jacobian_c(i,k,n3_c)*lambda_c(i,k,n3_c)/l2*XI13_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,2,index)/h2_c
             ! second set equation
             Vass((k-1)*3*n1_c+(i-1)*3+2) = Vass((k-1)*3*n1_c+(i-1)*3+2)-Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l2&
               *XI13_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,1,index)/h2_c-Jacobian_c(i,k,n3_c)*(2.d0*mu_c(i,k,n3_c)&
               +lambda_c(i,k,n3_c))/l2*XI23_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,2,index)/h2_c&
               -Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l2*XI33_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,3,index)/h2_c
             ! third set equation
             Vass((k-1)*3*n1_c+(i-1)*3+3) = Vass((k-1)*3*n1_c+(i-1)*3+3)-Jacobian_c(i,k,n3_c)*lambda_c(i,k,n3_c)/l2&
               *XI33_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,2,index)/h2_c&
               -Jacobian_c(i,k,n3_c)*mu_c(i,k,n3_c)/l2*XI23_c(i,k,n3_c)*ux_cof(j)*u_c(i,k+j,n3_c,3,index)/h2_c
          end do
       end do
    end do
    ! term 2
    ! interior
    do j = -2,n2_c+3
       do i = -2,n1_c+3
          ! second derivative 11 & 22 & 12 & 21
          ! first set
          lh_c(i,j,n3_c,1) = lh_c(i,j,n3_c,1)+((-Jacobian_c(i-2,j,n3_c)*(2.d0*mu_c(i-2,j,n3_c)+lambda_c(i-2,j,n3_c))/8.d0 &
              +Jacobian_c(i-1,j,n3_c)*(2.d0*mu_c(i-1,j,n3_c)+lambda_c(i-1,j,n3_c))/6.d0-Jacobian_c(i,j,n3_c)&
              *(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/8.d0)*u_c(i-2,j,n3_c,1,index)+(Jacobian_c(i-2,j,n3_c)&
              *(2.d0*mu_c(i-2,j,n3_c)+lambda_c(i-2,j,n3_c))/6.d0+Jacobian_c(i-1,j,n3_c)*(2.d0*mu_c(i-1,j,n3_c)&
              +lambda_c(i-1,j,n3_c))/2.d0+Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/2.d0 &
              +Jacobian_c(i+1,j,n3_c)*(2.d0*mu_c(i+1,j,n3_c)+lambda_c(i+1,j,n3_c))/6.d0)*u_c(i-1,j,n3_c,1,index)&
              +(-Jacobian_c(i-2,j,n3_c)*(2.d0*mu_c(i-2,j,n3_c)+lambda_c(i-2,j,n3_c))/24.d0-Jacobian_c(i-1,j,n3_c)&
              *(2.d0*mu_c(i-1,j,n3_c)+lambda_c(i-1,j,n3_c))*5.d0/6.d0-Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)&
              +lambda_c(i,j,n3_c))*3.d0/4.d0-Jacobian_c(i+1,j,n3_c)*(2.d0*mu_c(i+1,j,n3_c)+lambda_c(i+1,j,n3_c))*5.d0/6.d0 &
              -Jacobian_c(i+2,j,n3_c)*(2.d0*mu_c(i+2,j,n3_c)+lambda_c(i+2,j,n3_c))/24.d0)*u_c(i-0,j,n3_c,1,index)&
              +(Jacobian_c(i-1,j,n3_c)*(2.d0*mu_c(i-1,j,n3_c)+lambda_c(i-1,j,n3_c))/6.d0+Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)&
              +lambda_c(i,j,n3_c))/2.d0+Jacobian_c(i+1,j,n3_c)*(2.d0*mu_c(i+1,j,n3_c)+lambda_c(i+1,j,n3_c))/2.d0 &
              + Jacobian_c(i+2,j,n3_c)*(2.d0*mu_c(i+2,j,n3_c)+lambda_c(i+2,j,n3_c))/6.d0)*u_c(i+1,j,n3_c,1,index) &
              +(-Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/8.d0+Jacobian_c(i+1,j,n3_c)*(2.d0*mu_c(i+1,j,n3_c)&
              +lambda_c(i+1,j,n3_c))/6.d0-Jacobian_c(i+2,j,n3_c)*(2.d0*mu_c(i+2,j,n3_c)+lambda_c(i+2,j,n3_c))/8.d0)&
              *u_c(i+2,j,n3_c,1,index))/h1_c**2/l1**2+((-Jacobian_c(i,j-2,n3_c)*mu_c(i,j-2,n3_c)/8.d0+Jacobian_c(i,j-1,n3_c)&
              *mu_c(i,j-1,n3_c)/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0)*u_c(i,j-2,n3_c,1,index)+(Jacobian_c(i,j-2,n3_c)&
              *mu_c(i,j-2,n3_c)/6.d0 + Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)/2.d0+Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
              +Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/6.d0)*u_c(i,j-1,n3_c,1,index)+(-Jacobian_c(i,j-2,n3_c)*mu_c(i,j-2,n3_c)&
              /24.d0-Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)*5.d0/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)*3.d0/4.d0&
              -Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)*5.d0/6.d0 -Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)/24.d0)&
              *u_c(i,j-0,n3_c,1,index)+(Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)/6.d0 + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
              +Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/2.d0+Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)/6.d0)*u_c(i,j+1,n3_c,1,index)&
              +(-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0 + Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/6.d0-Jacobian_c(i,j+2,n3_c)&
              *mu_c(i,j+2,n3_c)/8.d0)*u_c(i,j+2,n3_c,1,index))/h2_c**2/l2**2+(Jacobian_c(i-2,j,n3_c)*lambda_c(i-2,j,n3_c)&
              *(u_c(i-2,j-2,n3_c,2,index)/12.d0-u_c(i-2,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i-2,j+1,n3_c,2,index)*2.d0/3.d0&
              -u_c(i-2,j+2,n3_c,2,index)/12.d0)/12.d0-Jacobian_c(i-1,j,n3_c)*lambda_c(i-1,j,n3_c)*(u_c(i-1,j-2,n3_c,2,index)/12.d0&
              -u_c(i-1,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i-1,j+1,n3_c,2,index)*2.d0/3.d0 - u_c(i-1,j+2,n3_c,2,index)/12.d0)*2.d0/3.d0&
              +Jacobian_c(i+1,j,n3_c)*lambda_c(i+1,j,n3_c)*(u_c(i+1,j-2,n3_c,2,index)/12.d0-u_c(i+1,j-1,n3_c,2,index)*2.d0/3.d0 &
              +u_c(i+1,j+1,n3_c,2,index)*2.d0/3.d0 - u_c(i+1,j+2,n3_c,2,index)/12.d0)*2.d0/3.d0-Jacobian_c(i+2,j,n3_c)&
              *lambda_c(i+2,j,n3_c)*(u_c(i+2,j-2,n3_c,2,index)/12.d0-u_c(i+2,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i+2,j+1,n3_c,2,index)&
              *2.d0/3.d0 - u_c(i+2,j+2,n3_c,2,index)/12.d0)/12.d0)/l1/l2/h1_c/h2_c+(Jacobian_c(i,j-2,n3_c)*mu_c(i,j-2,n3_c)&
              *(u_c(i-2,j-2,n3_c,2,index)/12.d0-u_c(i-1,j-2,n3_c,2,index)*2.d0/3.d0+u_c(i+1,j-2,n3_c,2,index)*2.d0/3.d0&
              -u_c(i+2,j-2,n3_c,2,index)/12.d0)/12.d0-Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)*(u_c(i-2,j-1,n3_c,2,index)/12.d0&
              -u_c(i-1,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i+1,j-1,n3_c,2,index)*2.d0/3.d0 - u_c(i+2,j-1,n3_c,2,index)/12.d0)*2.d0/3.d0&
              +Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)*(u_c(i-2,j+1,n3_c,2,index)/12.d0-u_c(i-1,j+1,n3_c,2,index)*2.d0/3.d0 &
              +u_c(i+1,j+1,n3_c,2,index)*2.d0/3.d0 - u_c(i+2,j+1,n3_c,2,index)/12.d0)*2.d0/3.d0 &
              -Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)*(u_c(i-2,j+2,n3_c,2,index)/12.d0-u_c(i-1,j+2,n3_c,2,index)*2.d0/3.d0 &
              +u_c(i+1,j+2,n3_c,2,index)*2.d0/3.d0 - u_c(i+2,j+2,n3_c,2,index)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
         !
         ! second set
         lh_c(i,j,n3_c,2) = lh_c(i,j,n3_c,2)+((-Jacobian_c(i-2,j,n3_c)*mu_c(i-2,j,n3_c)/8.d0+Jacobian_c(i-1,j,n3_c)&
             *mu_c(i-1,j,n3_c)/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0)*u_c(i-2,j,n3_c,2,index)+(Jacobian_c(i-2,j,n3_c)&
             *mu_c(i-2,j,n3_c)/6.d0 + Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)/2.d0+Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
             +Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/6.d0)*u_c(i-1,j,n3_c,2,index)+(-Jacobian_c(i-2,j,n3_c)*mu_c(i-2,j,n3_c)/24.d0&
             -Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)*5.d0/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)*3.d0/4.d0&
             -Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)*5.d0/6.d0 -Jacobian_c(i+2,j,n3_c)*mu_c(i+2,j,n3_c)/24.d0)&
             *u_c(i-0,j,n3_c,2,index)+(Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)/6.d0 + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0&
             +Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/2.d0+Jacobian_c(i+2,j,n3_c)*mu_c(i+2,j,n3_c)/6.d0)*u_c(i+1,j,n3_c,2,index)&
             +(-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0 + Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/6.d0-Jacobian_c(i+2,j,n3_c)&
             *mu_c(i+2,j,n3_c)/8.d0)*u_c(i+2,j,n3_c,2,index))/h1_c**2/l1**2+((-Jacobian_c(i,j-2,n3_c)*(2.d0*mu_c(i,j-2,n3_c)&
             +lambda_c(i,j-2,n3_c))/8.d0+Jacobian_c(i,j-1,n3_c)*(2.d0*mu_c(i,j-1,n3_c)+lambda_c(i,j-1,n3_c))/6.d0&
             -Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/8.d0)*u_c(i,j-2,n3_c,2,index)+(Jacobian_c(i,j-2,n3_c)&
             *(2.d0*mu_c(i,j-2,n3_c)+lambda_c(i,j-2,n3_c))/6.d0+Jacobian_c(i,j-1,n3_c)*(2.d0*mu_c(i,j-1,n3_c)+lambda_c(i,j-1,n3_c))&
             /2.d0+Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/2.d0+Jacobian_c(i,j+1,n3_c)*(2.d0*mu_c(i,j+1,n3_c)&
             +lambda_c(i,j+1,n3_c))/6.d0)*u_c(i,j-1,n3_c,2,index)+(-Jacobian_c(i,j-2,n3_c)*(2.d0*mu_c(i,j-2,n3_c)&
             +lambda_c(i,j-2,n3_c))/24.d0-Jacobian_c(i,j-1,n3_c)*(2.d0*mu_c(i,j-1,n3_c)+lambda_c(i,j-1,n3_c))*5.d0/6.d0&
             -Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*3.d0/4.d0-Jacobian_c(i,j+1,n3_c)*(2.d0*mu_c(i,j+1,n3_c)&
             +lambda_c(i,j+1,n3_c))*5.d0/6.d0-Jacobian_c(i,j+2,n3_c)*(2.d0*mu_c(i,j+2,n3_c)+lambda_c(i,j+2,n3_c))/24.d0)&
             *u_c(i,j-0,n3_c,2,index)+(Jacobian_c(i,j-1,n3_c)*(2.d0*mu_c(i,j-1,n3_c)+lambda_c(i,j-1,n3_c))/6.d0&
             +Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/2.d0+Jacobian_c(i,j+1,n3_c)*(2.d0*mu_c(i,j+1,n3_c)&
             +lambda_c(i,j+1,n3_c))/2.d0+Jacobian_c(i,j+2,n3_c)*(2.d0*mu_c(i,j+2,n3_c)+lambda_c(i,j+2,n3_c))/6.d0)&
             *u_c(i,j+1,n3_c,2,index)+(-Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/8.d0+Jacobian_c(i,j+1,n3_c)&
             *(2.d0*mu_c(i,j+1,n3_c)+lambda_c(i,j+1,n3_c))/6.d0-Jacobian_c(i,j+2,n3_c)*(2.d0*mu_c(i,j+2,n3_c)+lambda_c(i,j+2,n3_c))&
             /8.d0)*u_c(i,j+2,n3_c,2,index))/h2_c**2/l2**2+(Jacobian_c(i-2,j,n3_c)*mu_c(i-2,j,n3_c)*(u_c(i-2,j-2,n3_c,1,index)&
             /12.d0-u_c(i-2,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i-2,j+1,n3_c,1,index)*2.d0/3.d0-u_c(i-2,j+2,n3_c,1,index)/12.d0)/12.d0&
             -Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)*(u_c(i-1,j-2,n3_c,1,index)/12.d0-u_c(i-1,j-1,n3_c,1,index)*2.d0/3.d0&
             +u_c(i-1,j+1,n3_c,1,index)*2.d0/3.d0-u_c(i-1,j+2,n3_c,1,index)/12.d0)*2.d0/3.d0+Jacobian_c(i+1,j,n3_c)&
             *mu_c(i+1,j,n3_c)*(u_c(i+1,j-2,n3_c,1,index)/12.d0-u_c(i+1,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i+1,j+1,n3_c,1,index)&
             *2.d0/3.d0-u_c(i+1,j+2,n3_c,1,index)/12.d0)*2.d0/3.d0-Jacobian_c(i+2,j,n3_c)*mu_c(i+2,j,n3_c)&
             *(u_c(i+2,j-2,n3_c,1,index)/12.d0-u_c(i+2,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i+2,j+1,n3_c,1,index)*2.d0/3.d0&
             -u_c(i+2,j+2,n3_c,1,index)/12.d0)/12.d0)/l1/l2/h1_c/h2_c+(Jacobian_c(i,j-2,n3_c)*lambda_c(i,j-2,n3_c)&
             *(u_c(i-2,j-2,n3_c,1,index)/12.d0-u_c(i-1,j-2,n3_c,1,index)*2.d0/3.d0+u_c(i+1,j-2,n3_c,1,index)*2.d0/3.d0 &
             -u_c(i+2,j-2,n3_c,1,index)/12.d0)/12.d0-Jacobian_c(i,j-1,n3_c)*lambda_c(i,j-1,n3_c)*(u_c(i-2,j-1,n3_c,1,index)/12.d0&
             -u_c(i-1,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i+1,j-1,n3_c,1,index)*2.d0/3.d0-u_c(i+2,j-1,n3_c,1,index)/12.d0)*2.d0/3.d0&
             +Jacobian_c(i,j+1,n3_c)*lambda_c(i,j+1,n3_c)*(u_c(i-2,j+1,n3_c,1,index)/12.d0-u_c(i-1,j+1,n3_c,1,index)*2.d0/3.d0&
             +u_c(i+1,j+1,n3_c,1,index)*2.d0/3.d0 - u_c(i+2,j+1,n3_c,1,index)/12.d0)*2.d0/3.d0-Jacobian_c(i,j+2,n3_c)&
             *lambda_c(i,j+2,n3_c)*(u_c(i-2,j+2,n3_c,1,index)/12.d0-u_c(i-1,j+2,n3_c,1,index)*2.d0/3.d0 &
             +u_c(i+1,j+2,n3_c,1,index)*2.d0/3.d0 - u_c(i+2,j+2,n3_c,1,index)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
         ! third set
         lh_c(i,j,n3_c,3) = lh_c(i,j,n3_c,3)+((-Jacobian_c(i-2,j,n3_c)*mu_c(i-2,j,n3_c)/8.d0+Jacobian_c(i-1,j,n3_c)&
             *mu_c(i-1,j,n3_c)/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0)*u_c(i-2,j,n3_c,3,index)+(Jacobian_c(i-2,j,n3_c)&
             *mu_c(i-2,j,n3_c)/6.d0 + Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)/2.d0+Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
             + Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/6.d0)*u_c(i-1,j,n3_c,3,index)+(-Jacobian_c(i-2,j,n3_c)*mu_c(i-2,j,n3_c)&
             /24.d0-Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)*5.d0/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)*3.d0/4.d0 &
             -Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)*5.d0/6.d0 -Jacobian_c(i+2,j,n3_c)*mu_c(i+2,j,n3_c)/24.d0)&
             *u_c(i-0,j,n3_c,3,index)+(Jacobian_c(i-1,j,n3_c)*mu_c(i-1,j,n3_c)/6.d0 + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
             +Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/2.d0+Jacobian_c(i+2,j,n3_c)*mu_c(i+2,j,n3_c)/6.d0)*u_c(i+1,j,n3_c,3,index)&
             +(-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0 + Jacobian_c(i+1,j,n3_c)*mu_c(i+1,j,n3_c)/6.d0-Jacobian_c(i+2,j,n3_c)&
             *mu_c(i+2,j,n3_c)/8.d0)*u_c(i+2,j,n3_c,3,index))/h1_c**2/l1**2+((-Jacobian_c(i,j-2,n3_c)*mu_c(i,j-2,n3_c)/8.d0 &
             +Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0)*u_c(i,j-2,n3_c,3,index) &
             +(Jacobian_c(i,j-2,n3_c)*mu_c(i,j-2,n3_c)/6.d0 + Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)/2.d0+Jacobian_c(i,j,n3_c)&
             *mu_c(i,j,n3_c)/2.d0+Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/6.d0)*u_c(i,j-1,n3_c,3,index)+(-Jacobian_c(i,j-2,n3_c)&
             *mu_c(i,j-2,n3_c)/24.d0 - Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)*5.d0/6.d0-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)&
             *3.d0/4.d0-Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)*5.d0/6.d0 -Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)/24.d0)&
             *u_c(i,j-0,n3_c,3,index)+(Jacobian_c(i,j-1,n3_c)*mu_c(i,j-1,n3_c)/6.d0 + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/2.d0 &
             +Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/2.d0+Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)/6.d0)*u_c(i,j+1,n3_c,3,index)&
             +(-Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/8.d0 + Jacobian_c(i,j+1,n3_c)*mu_c(i,j+1,n3_c)/6.d0 &
             -Jacobian_c(i,j+2,n3_c)*mu_c(i,j+2,n3_c)/8.d0)*u_c(i,j+2,n3_c,3,index))/h2_c**2/l2**2
       end do
    end do
    !
    do j = -2,n2_c+3
       do k = -2,n1_c+3
          do k1 = 1,8
             do m = 1,8
                ! second derivative 33
                ! first set equation
                lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) +(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI23_c(k,j,n3_c+1-m)**2&
                     +XI33_c(k,j,n3_c+1-m)**2))*u_c(k,j,n3_c+1-k1,1,index))/h3_c**2+(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)&
                     *(mu_c(k,j,n3_c+1-m)+lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)&
                     *u_c(k,j,n3_c+1-k1,2,index))/h3_c**2+(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c(k,j,n3_c+1-k1,3,index))/h3_c**2
                ! second set equation
                lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) +(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)*u_c(k,j,n3_c+1-k1,1,index))/h3_c**2&
                     +(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)+lambda_c(k,j,n3_c+1-m))&
                     *XI23_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI13_c(k,j,n3_c+1-m)**2+XI33_c(k,j,n3_c+1-m)**2))&
                     *u_c(k,j,n3_c+1-k1,2,index))/h3_c**2+(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c(k,j,n3_c+1-k1,3,index))/h3_c**2
                ! third set equation
                lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) +(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c(k,j,n3_c+1-k1,1,index))/h3_c**2&
                     +(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)+lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)&
                     *XI33_c(k,j,n3_c+1-m)*u_c(k,j,n3_c+1-k1,2,index))/h3_c**2+(acof(1,k1,m)*Jacobian_c(k,j,n3_c+1-m)&
                     *((2.d0*mu_c(k,j,n3_c+1-m)+lambda_c(k,j,n3_c+1-m))*XI33_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)&
                     *(XI13_c(k,j,n3_c+1-m)**2+XI23_c(k,j,n3_c+1-m)**2))*u_c(k,j,n3_c+1-k1,3,index))/h3_c**2
             end do
          end do
       end do
    end do
    ! ghost points
    do j = -2,n2_c+3
       do k = -2,0
          ! first set equation
          lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               +(u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)+lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)&
               *XI23_c(k,j,n3_c))/h3_c**2+(u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! second set equation
          lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2+(u_c(k,j,n3_c+1,2,index)*ghcof(1)&
               *Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)+lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)&
               *(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2+(u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)&
               *(mu_c(k,j,n3_c)+lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! third set equation
          lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2+(u_c(k,j,n3_c+1,2,index)*ghcof(1)&
               *Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)+lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               +(u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
    !
    do j = -2,n2_c+3
       do k = n1_c+1,n1_c+3
          ! first set equation
          lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! second set equation
          lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! third set equation
          lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
    !
    do j = -2,0
       do k = 1,n1_c
          ! first set equation
          lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2 &
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! second set equation
          lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! third set equation
          lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
    !
    do j = n2_c+1,n2_c+3
      do k = 1,n1_c
         ! first set equation
         lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
              + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2 &
              + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
         ! second set equation
         lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
              + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
              + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
         ! third set equation
         lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c(k,j,n3_c+1,1,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
              + (u_c(k,j,n3_c+1,2,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
              + (u_c(k,j,n3_c+1,3,index)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
              +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
    !
    do i = -2,n2_c+3
       do j = -2,n1_c+3
          do k1 = 1,6
             ! mixed derivative 13 & 23 & 31 & 32
             ! first set equation
             lh_c(j,i,n3_c,1) = lh_c(j,i,n3_c,1)+(-Jacobian_c(j-2,i,n3_c)*(2.d0*mu_c(j-2,i,n3_c)+lambda_c(j-2,i,n3_c))&
                 *XI13_c(j-2,i,n3_c)*bof(1,k1)*u_c(j-2,i,n3_c+1-k1,1,index)/12.d0+Jacobian_c(j-1,i,n3_c)*(2.d0*mu_c(j-1,i,n3_c)&
                 +lambda_c(j-1,i,n3_c))*XI13_c(j-1,i,n3_c)*bof(1,k1)*u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)&
                 *(2.d0*mu_c(j+1,i,n3_c)+lambda_c(j+1,i,n3_c))*XI13_c(j+1,i,n3_c)*bof(1,k1)*u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0&
                 +Jacobian_c(j+2,i,n3_c)*(2.d0*mu_c(j+2,i,n3_c)+lambda_c(j+2,i,n3_c))*XI13_c(j+2,i,n3_c)*bof(1,k1)&
                 *u_c(j+2,i,n3_c+1-k1,1,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j-2,i,n3_c)*lambda_c(j-2,i,n3_c)*XI23_c(j-2,i,n3_c)&
                 *bof(1,k1)*u_c(j-2,i,n3_c+1-k1,2,index)/12.d0+Jacobian_c(j-1,i,n3_c)*lambda_c(j-1,i,n3_c)*XI23_c(j-1,i,n3_c)&
                 *bof(1,k1)*u_c(j-1,i,n3_c+1-k1,2,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*lambda_c(j+1,i,n3_c)*XI23_c(j+1,i,n3_c)&
                 *bof(1,k1)*u_c(j+1,i,n3_c+1-k1,2,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*lambda_c(j+2,i,n3_c)*XI23_c(j+2,i,n3_c)&
                 *bof(1,k1)*u_c(j+2,i,n3_c+1-k1,2,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j-2,i,n3_c)*lambda_c(j-2,i,n3_c)&
                 *XI33_c(j-2,i,n3_c)*bof(1,k1)*u_c(j-2,i,n3_c+1-k1,3,index)/12.d0+Jacobian_c(j-1,i,n3_c)*lambda_c(j-1,i,n3_c)&
                 *XI33_c(j-1,i,n3_c)*bof(1,k1)*u_c(j-1,i,n3_c+1-k1,3,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*lambda_c(j+1,i,n3_c)&
                 *XI33_c(j+1,i,n3_c)*bof(1,k1)*u_c(j+1,i,n3_c+1-k1,3,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*lambda_c(j+2,i,n3_c)&
                 *XI33_c(j+2,i,n3_c)*bof(1,k1)*u_c(j+2,i,n3_c+1-k1,3,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j,i-2,n3_c)&
                 *mu_c(j,i-2,n3_c)*XI23_c(j,i-2,n3_c)*bof(1,k1)*u_c(j,i-2,n3_c+1-k1,1,index)/12.d0+Jacobian_c(j,i-1,n3_c)&
                 *mu_c(j,i-1,n3_c)*XI23_c(j,i-1,n3_c)*bof(1,k1)*u_c(j,i-1,n3_c+1-k1,1,index)*2.d0/3.d0-Jacobian_c(j,i+1,n3_c)&
                 *mu_c(j,i+1,n3_c)*XI23_c(j,i+1,n3_c)*bof(1,k1)*u_c(j,i+1,n3_c+1-k1,1,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)&
                 *mu_c(j,i+2,n3_c)*XI23_c(j,i+2,n3_c)*bof(1,k1)*u_c(j,i+2,n3_c+1-k1,1,index)/12.d0)/l2/h2_c/h3_c&
                 +(-Jacobian_c(j,i-2,n3_c)*mu_c(j,i-2,n3_c)*XI13_c(j,i-2,n3_c)*bof(1,k1)*u_c(j,i-2,n3_c+1-k1,2,index)/12.d0 &
                 +Jacobian_c(j,i-1,n3_c)*mu_c(j,i-1,n3_c)*XI13_c(j,i-1,n3_c)*bof(1,k1)*u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0 &
                 -Jacobian_c(j,i+1,n3_c)*mu_c(j,i+1,n3_c)*XI13_c(j,i+1,n3_c)*bof(1,k1)*u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0 &
                 +Jacobian_c(j,i+2,n3_c)*mu_c(j,i+2,n3_c)*XI13_c(j,i+2,n3_c)*bof(1,k1)*u_c(j,i+2,n3_c+1-k1,2,index)/12.d0)&
                 /l2/h2_c/h3_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*(2.d0*mu_c(j,i,n3_c+1-k1)+lambda_c(j,i,n3_c+1-k1))&
                 *XI13_c(j,i,n3_c+1-k1)*(u_c(j-2,i,n3_c+1-k1,1,index)/12.d0-u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0&
                 +u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0-u_c(j+2,i,n3_c+1-k1,1,index)/12.d0))/l1/h3_c/h1_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)*XI23_c(j,i,n3_c+1-k1)*(u_c(j-2,i,n3_c+1-k1,2,index)/12.d0 &
                 -u_c(j-1,i,n3_c+1-k1,2,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,2,index)*2.d0/3.d0-u_c(j+2,i,n3_c+1-k1,2,index)&
                 /12.d0))/l1/h3_c/h1_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)*XI33_c(j,i,n3_c+1-k1)&
                 *(u_c(j-2,i,n3_c+1-k1,3,index)/12.d0-u_c(j-1,i,n3_c+1-k1,3,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,3,index)*2.d0/3.d0&
                 -u_c(j+2,i,n3_c+1-k1,3,index)/12.d0))/l1/h1_c/h3_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)&
                 *XI23_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,1,index)/12.d0-u_c(j,i-1,n3_c+1-k1,1,index)*2.d0/3.d0&
                 +u_c(j,i+1,n3_c+1-k1,1,index)*2.d0/3.d0-u_c(j,i+2,n3_c+1-k1,1,index)/12.d0))/l2/h3_c/h2_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*lambda_c(j,i,n3_c+1-k1)*XI13_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,2,index)/12.d0 &
                 -u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0+u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0&
                 -u_c(j,i+2,n3_c+1-k1,2,index)/12.d0))/l2/h3_c/h2_c
             ! second set equation
             lh_c(j,i,n3_c,2) = lh_c(j,i,n3_c,2)+(-Jacobian_c(j-2,i,n3_c)*mu_c(j-2,i,n3_c)*XI23_c(j-2,i,n3_c)*bof(1,k1)&
                 *u_c(j-2,i,n3_c+1-k1,1,index)/12.d0+Jacobian_c(j-1,i,n3_c)*mu_c(j-1,i,n3_c)*XI23_c(j-1,i,n3_c)*bof(1,k1)&
                 *u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*mu_c(j+1,i,n3_c)*XI23_c(j+1,i,n3_c)*bof(1,k1)&
                 *u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*mu_c(j+2,i,n3_c)*XI23_c(j+2,i,n3_c)*bof(1,k1)&
                 *u_c(j+2,i,n3_c+1-k1,1,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j-2,i,n3_c)*mu_c(j-2,i,n3_c)*XI13_c(j-2,i,n3_c)&
                 *bof(1,k1)*u_c(j-2,i,n3_c+1-k1,2,index)/12.d0+Jacobian_c(j-1,i,n3_c)*mu_c(j-1,i,n3_c)*XI13_c(j-1,i,n3_c)*bof(1,k1)&
                 *u_c(j-1,i,n3_c+1-k1,2,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*mu_c(j+1,i,n3_c)*XI13_c(j+1,i,n3_c)*bof(1,k1)&
                 *u_c(j+1,i,n3_c+1-k1,2,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*mu_c(j+2,i,n3_c)*XI13_c(j+2,i,n3_c)*bof(1,k1)&
                 *u_c(j+2,i,n3_c+1-k1,2,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j,i-2,n3_c)*lambda_c(j,i-2,n3_c)*XI13_c(j,i-2,n3_c)&
                 *bof(1,k1)*u_c(j,i-2,n3_c+1-k1,1,index)/12.d0+Jacobian_c(j,i-1,n3_c)*lambda_c(j,i-1,n3_c)*XI13_c(j,i-1,n3_c)&
                 *bof(1,k1)*u_c(j,i-1,n3_c+1-k1,1,index)*2.d0/3.d0-Jacobian_c(j,i+1,n3_c)*lambda_c(j,i+1,n3_c)*XI13_c(j,i+1,n3_c)&
                 *bof(1,k1)*u_c(j,i+1,n3_c+1-k1,1,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)*lambda_c(j,i+2,n3_c)*XI13_c(j,i+2,n3_c)&
                 *bof(1,k1)*u_c(j,i+2,n3_c+1-k1,1,index)/12.d0)/l2/h2_c/h3_c+(-Jacobian_c(j,i-2,n3_c)*(2.d0*mu_c(j,i-2,n3_c)&
                 +lambda_c(j,i-2,n3_c))*XI23_c(j,i-2,n3_c)*bof(1,k1)*u_c(j,i-2,n3_c+1-k1,2,index)/12.d0+Jacobian_c(j,i-1,n3_c)&
                 *(2.d0*mu_c(j,i-1,n3_c)+lambda_c(j,i-1,n3_c))*XI23_c(j,i-1,n3_c)*bof(1,k1)*u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0&
                 -Jacobian_c(j,i+1,n3_c)*(2.d0*mu_c(j,i+1,n3_c)+lambda_c(j,i+1,n3_c))*XI23_c(j,i+1,n3_c)*bof(1,k1)&
                 *u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)*(2.d0*mu_c(j,i+2,n3_c)+lambda_c(j,i+2,n3_c))&
                 *XI23_c(j,i+2,n3_c)*bof(1,k1)*u_c(j,i+2,n3_c+1-k1,2,index)/12.d0)/l2/h2_c/h3_c+(-Jacobian_c(j,i-2,n3_c)&
                 *lambda_c(j,i-2,n3_c)*XI33_c(j,i-2,n3_c)*bof(1,k1)*u_c(j,i-2,n3_c+1-k1,3,index)/12.d0+Jacobian_c(j,i-1,n3_c)&
                 *lambda_c(j,i-1,n3_c)*XI33_c(j,i-1,n3_c)*bof(1,k1)*u_c(j,i-1,n3_c+1-k1,3,index)*2.d0/3.d0-Jacobian_c(j,i+1,n3_c)&
                 *lambda_c(j,i+1,n3_c)*XI33_c(j,i+1,n3_c)*bof(1,k1)*u_c(j,i+1,n3_c+1-k1,3,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)&
                 *lambda_c(j,i+2,n3_c)*XI33_c(j,i+2,n3_c)*bof(1,k1)*u_c(j,i+2,n3_c+1-k1,3,index)/12.d0)/l2/h2_c/h3_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*lambda_c(j,i,n3_c+1-k1)*XI23_c(j,i,n3_c+1-k1)*(u_c(j-2,i,n3_c+1-k1,1,index)/12.d0 &
                 -u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0-u_c(j+2,i,n3_c+1-k1,1,index)&
                 /12.d0))/l1/h1_c/h3_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)*XI13_c(j,i,n3_c+1-k1)&
                 *(u_c(j-2,i,n3_c+1-k1,2,index)/12.d0-u_c(j-1,i,n3_c+1-k1,2,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,2,index)*2.d0/3.d0&
                 -u_c(j+2,i,n3_c+1-k1,2,index)/12.d0))/l1/h1_c/h3_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)&
                 *XI13_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,1,index)/12.d0-u_c(j,i-1,n3_c+1-k1,1,index)*2.d0/3.d0&
                 +u_c(j,i+1,n3_c+1-k1,1,index)*2.d0/3.d0-u_c(j,i+2,n3_c+1-k1,1,index)/12.d0))/l2/h3_c/h2_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*(2.d0*mu_c(j,i,n3_c+1-k1)+lambda_c(j,i,n3_c+1-k1))*XI23_c(j,i,n3_c+1-k1)&
                 *(u_c(j,i-2,n3_c+1-k1,2,index)/12.d0-u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0+u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0&
                 -u_c(j,i+2,n3_c+1-k1,2,index)/12.d0))/l2/h3_c/h2_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)&
                 *XI33_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,3,index)/12.d0-u_c(j,i-1,n3_c+1-k1,3,index)*2.d0/3.d0&
                 +u_c(j,i+1,n3_c+1-k1,3,index)*2.d0/3.d0-u_c(j,i+2,n3_c+1-k1,3,index)/12.d0))/l2/h3_c/h2_c
             ! third set equation
             lh_c(j,i,n3_c,3) = lh_c(j,i,n3_c,3)+(-Jacobian_c(j-2,i,n3_c)*mu_c(j-2,i,n3_c)*XI33_c(j-2,i,n3_c)*bof(1,k1)&
                 *u_c(j-2,i,n3_c+1-k1,1,index)/12.d0+Jacobian_c(j-1,i,n3_c)*mu_c(j-1,i,n3_c)*XI33_c(j-1,i,n3_c)*bof(1,k1)&
                 *u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*mu_c(j+1,i,n3_c)*XI33_c(j+1,i,n3_c)*bof(1,k1)&
                 *u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*mu_c(j+2,i,n3_c)*XI33_c(j+2,i,n3_c)*bof(1,k1)&
                 *u_c(j+2,i,n3_c+1-k1,1,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j-2,i,n3_c)*mu_c(j-2,i,n3_c)*XI13_c(j-2,i,n3_c)&
                 *bof(1,k1)*u_c(j-2,i,n3_c+1-k1,3,index)/12.d0+Jacobian_c(j-1,i,n3_c)*mu_c(j-1,i,n3_c)*XI13_c(j-1,i,n3_c)*bof(1,k1)&
                 *u_c(j-1,i,n3_c+1-k1,3,index)*2.d0/3.d0-Jacobian_c(j+1,i,n3_c)*mu_c(j+1,i,n3_c)*XI13_c(j+1,i,n3_c)*bof(1,k1)&
                 *u_c(j+1,i,n3_c+1-k1,3,index)*2.d0/3.d0+Jacobian_c(j+2,i,n3_c)*mu_c(j+2,i,n3_c)*XI13_c(j+2,i,n3_c)*bof(1,k1)&
                 *u_c(j+2,i,n3_c+1-k1,3,index)/12.d0)/l1/h1_c/h3_c+(-Jacobian_c(j,i-2,n3_c)*mu_c(j,i-2,n3_c)*XI33_c(j,i-2,n3_c)&
                 *bof(1,k1)*u_c(j,i-2,n3_c+1-k1,2,index)/12.d0+Jacobian_c(j,i-1,n3_c)*mu_c(j,i-1,n3_c)*XI33_c(j,i-1,n3_c)&
                 *bof(1,k1)*u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0-Jacobian_c(j,i+1,n3_c)*mu_c(j,i+1,n3_c)*XI33_c(j,i+1,n3_c)&
                 *bof(1,k1)*u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)*mu_c(j,i+2,n3_c)*XI33_c(j,i+2,n3_c)&
                 *bof(1,k1)*u_c(j,i+2,n3_c+1-k1,2,index)/12.d0)/l2/h2_c/h3_c+(-Jacobian_c(j,i-2,n3_c)*mu_c(j,i-2,n3_c)&
                 *XI23_c(j,i-2,n3_c)*bof(1,k1)*u_c(j,i-2,n3_c+1-k1,3,index)/12.d0+Jacobian_c(j,i-1,n3_c)*mu_c(j,i-1,n3_c)&
                 *XI23_c(j,i-1,n3_c)*bof(1,k1)*u_c(j,i-1,n3_c+1-k1,3,index)*2.d0/3.d0-Jacobian_c(j,i+1,n3_c)*mu_c(j,i+1,n3_c)&
                 *XI23_c(j,i+1,n3_c)*bof(1,k1)*u_c(j,i+1,n3_c+1-k1,3,index)*2.d0/3.d0+Jacobian_c(j,i+2,n3_c)*mu_c(j,i+2,n3_c)&
                 *XI23_c(j,i+2,n3_c)*bof(1,k1)*u_c(j,i+2,n3_c+1-k1,3,index)/12.d0)/l2/h2_c/h3_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*lambda_c(j,i,n3_c+1-k1)*XI33_c(j,i,n3_c+1-k1)*(u_c(j-2,i,n3_c+1-k1,1,index)/12.d0 &
                 -u_c(j-1,i,n3_c+1-k1,1,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,1,index)*2.d0/3.d0-u_c(j+2,i,n3_c+1-k1,1,index)&
                 /12.d0))/l1/h3_c/h1_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)*XI13_c(j,i,n3_c+1-k1)&
                 *(u_c(j-2,i,n3_c+1-k1,3,index)/12.d0-u_c(j-1,i,n3_c+1-k1,3,index)*2.d0/3.d0+u_c(j+1,i,n3_c+1-k1,3,index)*2.d0/3.d0&
                 -u_c(j+2,i,n3_c+1-k1,3,index)/12.d0))/l1/h3_c/h1_c+(-bof(1,k1)*Jacobian_c(j,i,n3_c+1-k1)*lambda_c(j,i,n3_c+1-k1)&
                 *XI33_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,2,index)/12.d0-u_c(j,i-1,n3_c+1-k1,2,index)*2.d0/3.d0&
                 +u_c(j,i+1,n3_c+1-k1,2,index)*2.d0/3.d0-u_c(j,i+2,n3_c+1-k1,2,index)/12.d0))/l2/h2_c/h3_c+(-bof(1,k1)&
                 *Jacobian_c(j,i,n3_c+1-k1)*mu_c(j,i,n3_c+1-k1)*XI23_c(j,i,n3_c+1-k1)*(u_c(j,i-2,n3_c+1-k1,3,index)/12.d0 &
                 -u_c(j,i-1,n3_c+1-k1,3,index)*2.d0/3.d0+u_c(j,i+1,n3_c+1-k1,3,index)*2.d0/3.d0&
                 -u_c(j,i+2,n3_c+1-k1,3,index)/12.d0))/l2/h2_c/h3_c
          end do
       end do
    end do
    ! project
    ! first set
    Mass_f1 = 0.d0
    do k = -1,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             do l = -1,2
                Mass_f1(2*i,2*k)=Mass_f1(2*i,2*k)+P(j)*(P(l)*lh_c(i+l,k+j,n3_c,1)/rho_c(i+l,k+j,n3_c))
             end do
          end do
      end do
    end do
    !
    do k = -1,n2_c+1
       do i = 0,n1_c+1
          do j = -1,2
             Mass_f1(2*i-1,2*k)=Mass_f1(2*i-1,2*k)+P(j)*lh_c(i,k+j,n3_c,1)/rho_c(i,j+k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             Mass_f1(2*i,2*k-1)=Mass_f1(2*i,2*k-1)+P(j)*lh_c(i+j,k,n3_c,1)/rho_c(i+j,k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = 0,n1_c+1
          Mass_f1(2*i-1,2*k-1)=lh_c(i,k,n3_c,1)/rho_c(i,k,n3_c)
       end do
    end do
    ! restrict
    ! first set
    do k = 1,n2_c
       do i = 1,n1_c
          do j = -4,2
             do l = -4,2
                Vass((k-1)*3*n1_c+(i-1)*3+1)=Vass((k-1)*3*n1_c+(i-1)*3+1)&
                     -17.d0/48.d0*h3_f*Rop(j)*(Rop(l)*rho_f(2*i+l,2*k+j,1)*Mass_f1(2*i+l,2*k+j)*1.d0)
             end do
          end do
       end do
    end do
    ! second set
    Mass_f1 = 0.d0
    do k = -1,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             do l = -1,2
                Mass_f1(2*i,2*k)=Mass_f1(2*i,2*k)+P(j)*(P(l)*lh_c(i+l,k+j,n3_c,2)/rho_c(i+l,k+j,n3_c))
             end do
          end do
      end do
    end do
    !
    do k = -1,n2_c+1
       do i = 0,n1_c+1
          do j = -1,2
             Mass_f1(2*i-1,2*k)=Mass_f1(2*i-1,2*k)+P(j)*lh_c(i,k+j,n3_c,2)/rho_c(i,j+k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             Mass_f1(2*i,2*k-1)=Mass_f1(2*i,2*k-1)+P(j)*lh_c(i+j,k,n3_c,2)/rho_c(i+j,k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = 0,n1_c+1
          Mass_f1(2*i-1,2*k-1)=lh_c(i,k,n3_c,2)/rho_c(i,k,n3_c)
       end do
    end do
    ! restriction
    ! second set
    do k = 1,n2_c
       do i = 1,n1_c
          do j = -4,2
             do l = -4,2
                Vass((k-1)*3*n1_c+(i-1)*3+2)=Vass((k-1)*3*n1_c+(i-1)*3+2) &
                     -17.d0/48.d0*h3_f*Rop(j)*(Rop(l)*rho_f(2*i+l,2*k+j,1)*Mass_f1(2*i+l,2*k+j)*1.d0)
             end do
          end do
       end do
    end do
    ! third set
    Mass_f1 = 0.d0
    do k = -1,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             do l = -1,2
                Mass_f1(2*i,2*k)=Mass_f1(2*i,2*k)+P(j)*(P(l)*lh_c(i+l,k+j,n3_c,3)/rho_c(i+l,k+j,n3_c))
             end do
          end do
      end do
    end do
    !
    do k = -1,n2_c+1
       do i = 0,n1_c+1
          do j = -1,2
             Mass_f1(2*i-1,2*k)=Mass_f1(2*i-1,2*k)+P(j)*lh_c(i,k+j,n3_c,3)/rho_c(i,j+k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = -1,n1_c+1
          do j = -1,2
             Mass_f1(2*i,2*k-1)=Mass_f1(2*i,2*k-1)+P(j)*lh_c(i+j,k,n3_c,3)/rho_c(i+j,k,n3_c)
          end do
       end do
    end do
    !
    do k = 0,n2_c+1
       do i = 0,n1_c+1
          Mass_f1(2*i-1,2*k-1)=lh_c(i,k,n3_c,3)/rho_c(i,k,n3_c)
       end do
    end do
    ! restrict
    ! third set
    do k = 1,n2_c
       do i = 1,n1_c
          do j = -4,2
             do l = -4,2
                Vass((k-1)*3*n1_c+(i-1)*3+3)=Vass((k-1)*3*n1_c+(i-1)*3+3) &
                     -17.d0/48.d0*h3_f*Rop(j)*(Rop(l)*rho_f(2*i+l,2*k+j,1)*Mass_f1(2*i+l,2*k+j)*1.d0)
             end do
          end do
       end do
    end do
    ! term 3
    do j = -2,n2_f+3
       do i = -2,n1_f+3
          ! second derivative 11 & 22 & 12 & 21
          ! first set
          lh_f(i,j,1,1) = lh_f(i,j,1,1)+((-Jacobian_f(i-2,j,1)*(2.d0*mu_f(i-2,j,1)+lambda_f(i-2,j,1))/8.d0+Jacobian_f(i-1,j,1)&
            *(2.d0*mu_f(i-1,j,1)+lambda_f(i-1,j,1))/6.d0-Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))/8.d0)&
            *u_f(i-2,j,1,1,index)+(Jacobian_f(i-2,j,1)*(2.d0*mu_f(i-2,j,1)+lambda_f(i-2,j,1))/6.d0+Jacobian_f(i-1,j,1)&
            *(2.d0*mu_f(i-1,j,1)+lambda_f(i-1,j,1))/2.d0+Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))/2.d0&
            +Jacobian_f(i+1,j,1)*(2.d0*mu_f(i+1,j,1)+lambda_f(i+1,j,1))/6.d0)*u_f(i-1,j,1,1,index)+(-Jacobian_f(i-2,j,1)&
            *(2.d0*mu_f(i-2,j,1)+lambda_f(i-2,j,1))/24.d0-Jacobian_f(i-1,j,1)*(2.d0*mu_f(i-1,j,1)+lambda_f(i-1,j,1))*5.d0/6.d0 &
            -Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))*3.d0/4.d0-Jacobian_f(i+1,j,1)*(2.d0*mu_f(i+1,j,1)&
            +lambda_f(i+1,j,1))*5.d0/6.d0-Jacobian_f(i+2,j,1)*(2.d0*mu_f(i+2,j,1)+lambda_f(i+2,j,1))/24.d0)*u_f(i-0,j,1,1,index)&
            +(Jacobian_f(i-1,j,1)*(2.d0*mu_f(i-1,j,1)+lambda_f(i-1,j,1))/6.d0+Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)&
            +lambda_f(i,j,1))/2.d0+Jacobian_f(i+1,j,1)*(2.d0*mu_f(i+1,j,1)+lambda_f(i+1,j,1))/2.d0+Jacobian_f(i+2,j,1)&
            *(2.d0*mu_f(i+2,j,1)+lambda_f(i+2,j,1))/6.d0)*u_f(i+1,j,1,1,index)+(-Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)&
            +lambda_f(i,j,1))/8.d0+Jacobian_f(i+1,j,1)*(2.d0*mu_f(i+1,j,1)+lambda_f(i+1,j,1))/6.d0-Jacobian_f(i+2,j,1)&
            *(2.d0*mu_f(i+2,j,1)+lambda_f(i+2,j,1))/8.d0)*u_f(i+2,j,1,1,index))/h1_f**2/l1**2+((-Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)&
            /8.d0+ Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)/6.d0-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0)*u_f(i,j-2,1,1,index) &
            +(Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)/6.d0 + Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)/2.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0 &
            +Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)/6.d0)*u_f(i,j-1,1,1,index)+(-Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)/24.d0&
            -Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)*5.d0/6.d0-Jacobian_f(i,j,1)*mu_f(i,j,1)*3.d0/4.d0-Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)&
            *5.d0/6.d0 -Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/24.d0)*u_f(i,j-0,1,1,index)+(Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)/6.d0&
            +Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0+Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)/2.d0+Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/6.d0)&
            *u_f(i,j+1,1,1,index)+(-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0 + Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)/6.d0 &
            -Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/8.d0)*u_f(i,j+2,1,1,index))/h2_f**2/l2**2+(Jacobian_f(i-2,j,1)*lambda_f(i-2,j,1)&
            *(u_f(i-2,j-2,1,2,index)/12.d0-u_f(i-2,j-1,1,2,index)*2.d0/3.d0+u_f(i-2,j+1,1,2,index)*2.d0/3.d0-u_f(i-2,j+2,1,2,index)&
            /12.d0)/12.d0-Jacobian_f(i-1,j,1)*lambda_f(i-1,j,1)*(u_f(i-1,j-2,1,2,index)/12.d0-u_f(i-1,j-1,1,2,index)*2.d0/3.d0 &
            +u_f(i-1,j+1,1,2,index)*2.d0/3.d0 - u_f(i-1,j+2,1,2,index)/12.d0)*2.d0/3.d0+Jacobian_f(i+1,j,1)*lambda_f(i+1,j,1)&
            *(u_f(i+1,j-2,1,2,index)/12.d0-u_f(i+1,j-1,1,2,index)*2.d0/3.d0+u_f(i+1,j+1,1,2,index)*2.d0/3.d0-u_f(i+1,j+2,1,2,index)&
            /12.d0)*2.d0/3.d0-Jacobian_f(i+2,j,1)*lambda_f(i+2,j,1)*(u_f(i+2,j-2,1,2,index)/12.d0-u_f(i+2,j-1,1,2,index)*2.d0/3.d0 &
            + u_f(i+2,j+1,1,2,index)*2.d0/3.d0 - u_f(i+2,j+2,1,2,index)/12.d0)/12.d0)/l1/l2/h1_f/h2_f+(Jacobian_f(i,j-2,1)&
            *mu_f(i,j-2,1)*(u_f(i-2,j-2,1,2,index)/12.d0-u_f(i-1,j-2,1,2,index)*2.d0/3.d0+u_f(i+1,j-2,1,2,index)*2.d0/3.d0&
            -u_f(i+2,j-2,1,2,index)/12.d0)/12.d0-Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)*(u_f(i-2,j-1,1,2,index)/12.d0&
            -u_f(i-1,j-1,1,2,index)*2.d0/3.d0+u_f(i+1,j-1,1,2,index)*2.d0/3.d0 - u_f(i+2,j-1,1,2,index)/12.d0)*2.d0/3.d0 &
            +Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)*(u_f(i-2,j+1,1,2,index)/12.d0-u_f(i-1,j+1,1,2,index)*2.d0/3.d0&
            +u_f(i+1,j+1,1,2,index)*2.d0/3.d0 - u_f(i+2,j+1,1,2,index)/12.d0)*2.d0/3.d0-Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)&
            *(u_f(i-2,j+2,1,2,index)/12.d0-u_f(i-1,j+2,1,2,index)*2.d0/3.d0+u_f(i+1,j+2,1,2,index)*2.d0/3.d0 &
            -u_f(i+2,j+2,1,2,index)/12.d0)/12.d0)/l1/l2/h1_f/h2_f
          ! second set
          lh_f(i,j,1,2) = lh_f(i,j,1,2)+((-Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/8.d0 + Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/6.d0 &
            -Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0)*u_f(i-2,j,1,2,index)+(Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/6.d0&
            +Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/2.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0+Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/6.d0)&
            *u_f(i-1,j,1,2,index)+(-Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/24.d0 - Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)*5.d0/6.d0 &
            -Jacobian_f(i,j,1)*mu_f(i,j,1)*3.d0/4.d0-Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)*5.d0/6.d0 -Jacobian_f(i+2,j,1)&
            *mu_f(i+2,j,1)/24.d0)*u_f(i-0,j,1,2,index)+(Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/6.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0&
            +Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/2.d0+Jacobian_f(i+2,j,1)*mu_f(i+2,j,1)/6.d0)*u_f(i+1,j,1,2,index) &
            +(-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0 + Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/6.d0-Jacobian_f(i+2,j,1)*mu_f(i+2,j,1)/8.d0)&
            *u_f(i+2,j,1,2,index))/h1_f**2/l1**2+((-Jacobian_f(i,j-2,1)*(2.d0*mu_f(i,j-2,1)+lambda_f(i,j-2,1))/8.d0&
            +Jacobian_f(i,j-1,1)*(2.d0*mu_f(i,j-1,1)+lambda_f(i,j-1,1))/6.d0-Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))&
            /8.d0)*u_f(i,j-2,1,2,index)+(Jacobian_f(i,j-2,1)*(2.d0*mu_f(i,j-2,1)+lambda_f(i,j-2,1))/6.d0+Jacobian_f(i,j-1,1)&
            *(2.d0*mu_f(i,j-1,1)+lambda_f(i,j-1,1))/2.d0+Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))/2.d0&
            +Jacobian_f(i,j+1,1)*(2.d0*mu_f(i,j+1,1)+lambda_f(i,j+1,1))/6.d0)*u_f(i,j-1,1,2,index)+(-Jacobian_f(i,j-2,1)&
            *(2.d0*mu_f(i,j-2,1)+lambda_f(i,j-2,1))/24.d0-Jacobian_f(i,j-1,1)*(2.d0*mu_f(i,j-1,1)+lambda_f(i,j-1,1))*5.d0/6.d0 &
            -Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)+lambda_f(i,j,1))*3.d0/4.d0-Jacobian_f(i,j+1,1)*(2.d0*mu_f(i,j+1,1)&
            +lambda_f(i,j+1,1))*5.d0/6.d0-Jacobian_f(i,j+2,1)*(2.d0*mu_f(i,j+2,1)+lambda_f(i,j+2,1))/24.d0)*u_f(i,j-0,1,2,index) &
            +(Jacobian_f(i,j-1,1)*(2.d0*mu_f(i,j-1,1)+lambda_f(i,j-1,1))/6.d0+Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)&
            +lambda_f(i,j,1))/2.d0+Jacobian_f(i,j+1,1)*(2.d0*mu_f(i,j+1,1)+lambda_f(i,j+1,1))/2.d0+Jacobian_f(i,j+2,1)&
            *(2.d0*mu_f(i,j+2,1)+lambda_f(i,j+2,1))/6.d0)*u_f(i,j+1,1,2,index)+(-Jacobian_f(i,j,1)*(2.d0*mu_f(i,j,1)&
            +lambda_f(i,j,1))/8.d0+Jacobian_f(i,j+1,1)*(2.d0*mu_f(i,j+1,1)+lambda_f(i,j+1,1))/6.d0-Jacobian_f(i,j+2,1)&
            *(2.d0*mu_f(i,j+2,1)+lambda_f(i,j+2,1))/8.d0)*u_f(i,j+2,1,2,index))/h2_f**2/l2**2+(Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)&
            *(u_f(i-2,j-2,1,1,index)/12.d0-u_f(i-2,j-1,1,1,index)*2.d0/3.d0+u_f(i-2,j+1,1,1,index)*2.d0/3.d0&
            -u_f(i-2,j+2,1,1,index)/12.d0)/12.d0-Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)*(u_f(i-1,j-2,1,1,index)/12.d0&
            -u_f(i-1,j-1,1,1,index)*2.d0/3.d0+u_f(i-1,j+1,1,1,index)*2.d0/3.d0 - u_f(i-1,j+2,1,1,index)/12.d0)*2.d0/3.d0&
            +Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)*(u_f(i+1,j-2,1,1,index)/12.d0-u_f(i+1,j-1,1,1,index)*2.d0/3.d0 &
            +u_f(i+1,j+1,1,1,index)*2.d0/3.d0 - u_f(i+1,j+2,1,1,index)/12.d0)*2.d0/3.d0-Jacobian_f(i+2,j,1)*mu_f(i+2,j,1)&
            *(u_f(i+2,j-2,1,1,index)/12.d0-u_f(i+2,j-1,1,1,index)*2.d0/3.d0+u_f(i+2,j+1,1,1,index)*2.d0/3.d0&
            -u_f(i+2,j+2,1,1,index)/12.d0)/12.d0)/l1/l2/h1_f/h2_f+(Jacobian_f(i,j-2,1)*lambda_f(i,j-2,1)&
            *(u_f(i-2,j-2,1,1,index)/12.d0-u_f(i-1,j-2,1,1,index)*2.d0/3.d0+u_f(i+1,j-2,1,1,index)*2.d0/3.d0&
            -u_f(i+2,j-2,1,1,index)/12.d0)/12.d0-Jacobian_f(i,j-1,1)*lambda_f(i,j-1,1)*(u_f(i-2,j-1,1,1,index)/12.d0&
            -u_f(i-1,j-1,1,1,index)*2.d0/3.d0+u_f(i+1,j-1,1,1,index)*2.d0/3.d0 - u_f(i+2,j-1,1,1,index)/12.d0)*2.d0/3.d0 &
            +Jacobian_f(i,j+1,1)*lambda_f(i,j+1,1)*(u_f(i-2,j+1,1,1,index)/12.d0-u_f(i-1,j+1,1,1,index)*2.d0/3.d0 &
            +u_f(i+1,j+1,1,1,index)*2.d0/3.d0 - u_f(i+2,j+1,1,1,index)/12.d0)*2.d0/3.d0-Jacobian_f(i,j+2,1)*lambda_f(i,j+2,1)&
            *(u_f(i-2,j+2,1,1,index)/12.d0-u_f(i-1,j+2,1,1,index)*2.d0/3.d0+u_f(i+1,j+2,1,1,index)*2.d0/3.d0&
            -u_f(i+2,j+2,1,1,index)/12.d0)/12.d0)/l1/l2/h1_f/h2_f
          ! third set
          lh_f(i,j,1,3) = lh_f(i,j,1,3)+((-Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/8.d0 + Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/6.d0 &
            -Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0)*u_f(i-2,j,1,3,index)+(Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/6.d0&
            +Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/2.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0+Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/6.d0)&
            *u_f(i-1,j,1,3,index)+(-Jacobian_f(i-2,j,1)*mu_f(i-2,j,1)/24.d0-Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)*5.d0/6.d0 &
            -Jacobian_f(i,j,1)*mu_f(i,j,1)*3.d0/4.d0-Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)*5.d0/6.d0-Jacobian_f(i+2,j,1)&
            *mu_f(i+2,j,1)/24.d0)*u_f(i-0,j,1,3,index)+(Jacobian_f(i-1,j,1)*mu_f(i-1,j,1)/6.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)&
            /2.d0+Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/2.d0+Jacobian_f(i+2,j,1)*mu_f(i+2,j,1)/6.d0)*u_f(i+1,j,1,3,index) &
            +(-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0 + Jacobian_f(i+1,j,1)*mu_f(i+1,j,1)/6.d0-Jacobian_f(i+2,j,1)*mu_f(i+2,j,1)&
            /8.d0)*u_f(i+2,j,1,3,index))/h1_f**2/l1**2+((-Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)/8.d0+Jacobian_f(i,j-1,1)&
            *mu_f(i,j-1,1)/6.d0-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0)*u_f(i,j-2,1,3,index)+(Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)&
            /6.d0 + Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)/2.d0+Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0+Jacobian_f(i,j+1,1)&
            *mu_f(i,j+1,1)/6.d0)*u_f(i,j-1,1,3,index)+(-Jacobian_f(i,j-2,1)*mu_f(i,j-2,1)/24.d0-Jacobian_f(i,j-1,1)&
            *mu_f(i,j-1,1)*5.d0/6.d0-Jacobian_f(i,j,1)*mu_f(i,j,1)*3.d0/4.d0-Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)*5.d0/6.d0&
            -Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/24.d0)*u_f(i,j-0,1,3,index)+(Jacobian_f(i,j-1,1)*mu_f(i,j-1,1)/6.d0&
            +Jacobian_f(i,j,1)*mu_f(i,j,1)/2.d0+Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)/2.d0+Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/6.d0)&
            *u_f(i,j+1,1,3,index)+(-Jacobian_f(i,j,1)*mu_f(i,j,1)/8.d0 + Jacobian_f(i,j+1,1)*mu_f(i,j+1,1)/6.d0 &
            -Jacobian_f(i,j+2,1)*mu_f(i,j+2,1)/8.d0)*u_f(i,j+2,1,3,index))/h2_f**2/l2**2
       end do
    end do
    !
    do j = -2,n2_f+3
       do k = -2,n1_f+3
          do k1 = 1,8
             do m = 1,8
                ! second derivative 33
                ! first set equation
                lh_f(k,j,1,1) = lh_f(k,j,1,1) +(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*((2.d0*mu_f(k,j,m)+lambda_f(k,j,m))&
                  *XI13_f(k,j,m)**2+mu_f(k,j,m)*(XI23_f(k,j,m)**2+XI33_f(k,j,m)**2))*u_f(k,j,k1,1,index))/h3_f**2&
                  +(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)*XI23_f(k,j,m)&
                  *u_f(k,j,k1,2,index))/h3_f**2+(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                  *XI33_f(k,j,m)*u_f(k,j,k1,3,index))/h3_f**2
                ! second set equation
                lh_f(k,j,1,2) = lh_f(k,j,1,2) +(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                  *XI23_f(k,j,m)*u_f(k,j,k1,1,index))/h3_f**2+(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*((2.d0*mu_f(k,j,m)&
                  +lambda_f(k,j,m))*XI23_f(k,j,m)**2+mu_f(k,j,m)*(XI13_f(k,j,m)**2+XI33_f(k,j,m)**2))*u_f(k,j,k1,2,index))/h3_f**2&
                  +(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI23_f(k,j,m)*XI33_f(k,j,m)&
                  *u_f(k,j,k1,3,index))/h3_f**2
                ! third set equation
                lh_f(k,j,1,3) = lh_f(k,j,1,3) +(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                  *XI33_f(k,j,m)*u_f(k,j,k1,1,index))/h3_f**2+(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)&
                  +lambda_f(k,j,m))*XI23_f(k,j,m)*XI33_f(k,j,m)*u_f(k,j,k1,2,index))/h3_f**2+(acof_no_gp(1,k1,m)*Jacobian_f(k,j,m)&
                  *((2.d0*mu_f(k,j,m)+lambda_f(k,j,m))*XI33_f(k,j,m)**2+mu_f(k,j,m)*(XI13_f(k,j,m)**2+XI23_f(k,j,m)**2))&
                  *u_f(k,j,k1,3,index))/h3_f**2
            end do
         end do
     end do
  end do
  !
  do i = -2,n2_f+3
     do j = -2,n1_f+3
        do k1 = 1,6
           ! mixed derivative 13 & 23 & 31 & 32
           ! first set equation
           lh_f(j,i,1,1) = lh_f(j,i,1,1)+(Jacobian_f(j-2,i,1)*(2.d0*mu_f(j-2,i,1)+lambda_f(j-2,i,1))*XI13_f(j-2,i,1)*bof(1,k1)&
             *u_f(j-2,i,k1,1,index)/12.d0-Jacobian_f(j-1,i,1)*(2.d0*mu_f(j-1,i,1)+lambda_f(j-1,i,1))*XI13_f(j-1,i,1)*bof(1,k1)&
             *u_f(j-1,i,k1,1,index)*2.d0/3.d0+Jacobian_f(j+1,i,1)*(2.d0*mu_f(j+1,i,1)+lambda_f(j+1,i,1))*XI13_f(j+1,i,1)*bof(1,k1)&
             *u_f(j+1,i,k1,1,index)*2.d0/3.d0-Jacobian_f(j+2,i,1)*(2.d0*mu_f(j+2,i,1)+lambda_f(j+2,i,1))*XI13_f(j+2,i,1)*bof(1,k1)&
             *u_f(j+2,i,k1,1,index)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j-2,i,1)*lambda_f(j-2,i,1)*XI23_f(j-2,i,1)*bof(1,k1)&
             *u_f(j-2,i,k1,2,index)/12.d0-Jacobian_f(j-1,i,1)*lambda_f(j-1,i,1)*XI23_f(j-1,i,1)*bof(1,k1)*u_f(j-1,i,k1,2,index)&
             *2.d0/3.d0+Jacobian_f(j+1,i,1)*lambda_f(j+1,i,1)*XI23_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,2,index)*2.d0/3.d0 &
             -Jacobian_f(j+2,i,1)*lambda_f(j+2,i,1)*XI23_f(j+2,i,1)*bof(1,k1)*u_f(j+2,i,k1,2,index)/12.d0)/l1/h1_f/h3_f&
             +(Jacobian_f(j-2,i,1)*lambda_f(j-2,i,1)*XI33_f(j-2,i,1)*bof(1,k1)*u_f(j-2,i,k1,3,index)/12.d0 -Jacobian_f(j-1,i,1)&
             *lambda_f(j-1,i,1)*XI33_f(j-1,i,1)*bof(1,k1)*u_f(j-1,i,k1,3,index)*2.d0/3.d0+Jacobian_f(j+1,i,1)*lambda_f(j+1,i,1)&
             *XI33_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,3,index)*2.d0/3.d0-Jacobian_f(j+2,i,1)*lambda_f(j+2,i,1)*XI33_f(j+2,i,1)&
             *bof(1,k1)*u_f(j+2,i,k1,3,index)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j,i-2,1)*mu_f(j,i-2,1)*XI23_f(j,i-2,1)*bof(1,k1)&
             *u_f(j,i-2,k1,1,index)/12.d0-Jacobian_f(j,i-1,1)*mu_f(j,i-1,1)*XI23_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,1,index)&
             *2.d0/3.d0+Jacobian_f(j,i+1,1)*mu_f(j,i+1,1)*XI23_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,1,index)*2.d0/3.d0 &
             -Jacobian_f(j,i+2,1)*mu_f(j,i+2,1)*XI23_f(j,i+2,1)*bof(1,k1)*u_f(j,i+2,k1,1,index)/12.d0)/l2/h2_f/h3_f&
             +(Jacobian_f(j,i-2,1)*mu_f(j,i-2,1)*XI13_f(j,i-2,1)*bof(1,k1)*u_f(j,i-2,k1,2,index)/12.d0-Jacobian_f(j,i-1,1)&
             *mu_f(j,i-1,1)*XI13_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,2,index)*2.d0/3.d0+Jacobian_f(j,i+1,1)*mu_f(j,i+1,1)&
             *XI13_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,2,index)*2.d0/3.d0-Jacobian_f(j,i+2,1)*mu_f(j,i+2,1)*XI13_f(j,i+2,1)&
             *bof(1,k1)*u_f(j,i+2,k1,2,index)/12.d0)/l2/h2_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)*(2.d0*mu_f(j,i,k1)&
             +lambda_f(j,i,k1))*XI13_f(j,i,k1)*(u_f(j-2,i,k1,1,index)/12.d0-u_f(j-1,i,k1,1,index)*2.d0/3.d0+u_f(j+1,i,k1,1,index)&
             *2.d0/3.d0-u_f(j+2,i,k1,1,index)/12.d0))/l1/h3_f/h1_f+(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI23_f(j,i,k1)&
             *(u_f(j-2,i,k1,2,index)/12.d0-u_f(j-1,i,k1,2,index)*2.d0/3.d0+u_f(j+1,i,k1,2,index)*2.d0/3.d0-u_f(j+2,i,k1,2,index)&
             /12.d0))/l1/h3_f/h1_f+(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI33_f(j,i,k1)*(u_f(j-2,i,k1,3,index)/12.d0&
             -u_f(j-1,i,k1,3,index)*2.d0/3.d0+u_f(j+1,i,k1,3,index)*2.d0/3.d0-u_f(j+2,i,k1,3,index)/12.d0))/l1/h1_f/h3_f&
             +(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI23_f(j,i,k1)*(u_f(j,i-2,k1,1,index)/12.d0-u_f(j,i-1,k1,1,index)&
             *2.d0/3.d0+u_f(j,i+1,k1,1,index)*2.d0/3.d0-u_f(j,i+2,k1,1,index)/12.d0))/l2/h3_f/h2_f+(bof(1,k1)*Jacobian_f(j,i,k1)&
             *lambda_f(j,i,k1)*XI13_f(j,i,k1)*(u_f(j,i-2,k1,2,index)/12.d0-u_f(j,i-1,k1,2,index)*2.d0/3.d0+u_f(j,i+1,k1,2,index)&
             *2.d0/3.d0-u_f(j,i+2,k1,2,index)/12.d0))/l2/h3_f/h2_f
           ! second set equation
           lh_f(j,i,1,2) = lh_f(j,i,1,2) &
             +(Jacobian_f(j-2,i,1)*mu_f(j-2,i,1)*XI23_f(j-2,i,1)*bof(1,k1)*u_f(j-2,i,k1,1,index)/12.d0-Jacobian_f(j-1,i,1)&
             *mu_f(j-1,i,1)*XI23_f(j-1,i,1)*bof(1,k1)*u_f(j-1,i,k1,1,index)*2.d0/3.d0+Jacobian_f(j+1,i,1)*mu_f(j+1,i,1)&
             *XI23_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,1,index)*2.d0/3.d0-Jacobian_f(j+2,i,1)*mu_f(j+2,i,1)*XI23_f(j+2,i,1)&
             *bof(1,k1)*u_f(j+2,i,k1,1,index)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j-2,i,1)*mu_f(j-2,i,1)*XI13_f(j-2,i,1)*bof(1,k1)&
             *u_f(j-2,i,k1,2,index)/12.d0-Jacobian_f(j-1,i,1)*mu_f(j-1,i,1)*XI13_f(j-1,i,1)*bof(1,k1)*u_f(j-1,i,k1,2,index)&
             *2.d0/3.d0+Jacobian_f(j+1,i,1)*mu_f(j+1,i,1)*XI13_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,2,index)*2.d0/3.d0 &
             -Jacobian_f(j+2,i,1)*mu_f(j+2,i,1)*XI13_f(j+2,i,1)*bof(1,k1)*u_f(j+2,i,k1,2,index)/12.d0)/l1/h1_f/h3_f&
             +(Jacobian_f(j,i-2,1)*lambda_f(j,i-2,1)*XI13_f(j,i-2,1)*bof(1,k1)*u_f(j,i-2,k1,1,index)/12.d0-Jacobian_f(j,i-1,1)&
             *lambda_f(j,i-1,1)*XI13_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,1,index)*2.d0/3.d0+Jacobian_f(j,i+1,1)*lambda_f(j,i+1,1)&
             *XI13_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,1,index)*2.d0/3.d0-Jacobian_f(j,i+2,1)*lambda_f(j,i+2,1)*XI13_f(j,i+2,1)&
             *bof(1,k1)*u_f(j,i+2,k1,1,index)/12.d0)/l2/h2_f/h3_f+(Jacobian_f(j,i-2,1)*(2.d0*mu_f(j,i-2,1)+lambda_f(j,i-2,1))&
             *XI23_f(j,i-2,1)*bof(1,k1)*u_f(j,i-2,k1,2,index)/12.d0-Jacobian_f(j,i-1,1)*(2.d0*mu_f(j,i-1,1)+lambda_f(j,i-1,1))&
             *XI23_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,2,index)*2.d0/3.d0+Jacobian_f(j,i+1,1)*(2.d0*mu_f(j,i+1,1)&
             +lambda_f(j,i+1,1))*XI23_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,2,index)*2.d0/3.d0-Jacobian_f(j,i+2,1)&
             *(2.d0*mu_f(j,i+2,1)+lambda_f(j,i+2,1))*XI23_f(j,i+2,1)*bof(1,k1)*u_f(j,i+2,k1,2,index)/12.d0)/l2/h2_f/h3_f&
             +(Jacobian_f(j,i-2,1)*lambda_f(j,i-2,1)*XI33_f(j,i-2,1)*bof(1,k1)*u_f(j,i-2,k1,3,index)/12.d0-Jacobian_f(j,i-1,1)&
             *lambda_f(j,i-1,1)*XI33_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,3,index)*2.d0/3.d0+Jacobian_f(j,i+1,1)*lambda_f(j,i+1,1)&
             *XI33_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,3,index)*2.d0/3.d0-Jacobian_f(j,i+2,1)*lambda_f(j,i+2,1)*XI33_f(j,i+2,1)&
             *bof(1,k1)*u_f(j,i+2,k1,3,index)/12.d0)/l2/h2_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)*lambda_f(j,i,k1)*XI23_f(j,i,k1)&
             *(u_f(j-2,i,k1,1,index)/12.d0-u_f(j-1,i,k1,1,index)*2.d0/3.d0+u_f(j+1,i,k1,1,index)*2.d0/3.d0&
             -u_f(j+2,i,k1,1,index)/12.d0))/l1/h1_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI13_f(j,i,k1)&
             *(u_f(j-2,i,k1,2,index)/12.d0-u_f(j-1,i,k1,2,index)*2.d0/3.d0+u_f(j+1,i,k1,2,index)*2.d0/3.d0-u_f(j+2,i,k1,2,index)&
             /12.d0))/l1/h1_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI13_f(j,i,k1)*(u_f(j,i-2,k1,1,index)/12.d0&
             -u_f(j,i-1,k1,1,index)*2.d0/3.d0+u_f(j,i+1,k1,1,index)*2.d0/3.d0-u_f(j,i+2,k1,1,index)/12.d0))/l2/h3_f/h2_f&
             +(bof(1,k1)*Jacobian_f(j,i,k1)*(2.d0*mu_f(j,i,k1)+lambda_f(j,i,k1))*XI23_f(j,i,k1)*(u_f(j,i-2,k1,2,index)/12.d0&
             -u_f(j,i-1,k1,2,index)*2.d0/3.d0+u_f(j,i+1,k1,2,index)*2.d0/3.d0-u_f(j,i+2,k1,2,index)/12.d0))/l2/h3_f/h2_f&
             +(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI33_f(j,i,k1)*(u_f(j,i-2,k1,3,index)/12.d0-u_f(j,i-1,k1,3,index)&
             *2.d0/3.d0+u_f(j,i+1,k1,3,index)*2.d0/3.d0-u_f(j,i+2,k1,3,index)/12.d0))/l2/h3_f/h2_f
           ! third set equation
           lh_f(j,i,1,3) = lh_f(j,i,1,3)+(Jacobian_f(j-2,i,1)*mu_f(j-2,i,1)*XI33_f(j-2,i,1)*bof(1,k1)*u_f(j-2,i,k1,1,index)/12.d0&
             -Jacobian_f(j-1,i,1)*mu_f(j-1,i,1)*XI33_f(j-1,i,1)*bof(1,k1)*u_f(j-1,i,k1,1,index)*2.d0/3.d0+Jacobian_f(j+1,i,1)&
             *mu_f(j+1,i,1)*XI33_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,1,index)*2.d0/3.d0-Jacobian_f(j+2,i,1)*mu_f(j+2,i,1)&
             *XI33_f(j+2,i,1)*bof(1,k1)*u_f(j+2,i,k1,1,index)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j-2,i,1)*mu_f(j-2,i,1)&
             *XI13_f(j-2,i,1)*bof(1,k1)*u_f(j-2,i,k1,3,index)/12.d0-Jacobian_f(j-1,i,1)*mu_f(j-1,i,1)*XI13_f(j-1,i,1)*bof(1,k1)&
             *u_f(j-1,i,k1,3,index)*2.d0/3.d0+Jacobian_f(j+1,i,1)*mu_f(j+1,i,1)*XI13_f(j+1,i,1)*bof(1,k1)*u_f(j+1,i,k1,3,index)&
             *2.d0/3.d0-Jacobian_f(j+2,i,1)*mu_f(j+2,i,1)*XI13_f(j+2,i,1)*bof(1,k1)*u_f(j+2,i,k1,3,index)/12.d0)/l1/h1_f/h3_f&
             +(Jacobian_f(j,i-2,1)*mu_f(j,i-2,1)*XI33_f(j,i-2,1)*bof(1,k1)*u_f(j,i-2,k1,2,index)/12.d0-Jacobian_f(j,i-1,1)&
             *mu_f(j,i-1,1)*XI33_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,2,index)*2.d0/3.d0+Jacobian_f(j,i+1,1)*mu_f(j,i+1,1)&
             *XI33_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,2,index)*2.d0/3.d0-Jacobian_f(j,i+2,1)*mu_f(j,i+2,1)*XI33_f(j,i+2,1)&
             *bof(1,k1)*u_f(j,i+2,k1,2,index)/12.d0)/l2/h2_f/h3_f+(Jacobian_f(j,i-2,1)*mu_f(j,i-2,1)*XI23_f(j,i-2,1)*bof(1,k1)&
             *u_f(j,i-2,k1,3,index)/12.d0-Jacobian_f(j,i-1,1)*mu_f(j,i-1,1)*XI23_f(j,i-1,1)*bof(1,k1)*u_f(j,i-1,k1,3,index)&
             *2.d0/3.d0+Jacobian_f(j,i+1,1)*mu_f(j,i+1,1)*XI23_f(j,i+1,1)*bof(1,k1)*u_f(j,i+1,k1,3,index)*2.d0/3.d0 &
             -Jacobian_f(j,i+2,1)*mu_f(j,i+2,1)*XI23_f(j,i+2,1)*bof(1,k1)*u_f(j,i+2,k1,3,index)/12.d0)/l2/h2_f/h3_f+(bof(1,k1)&
             *Jacobian_f(j,i,k1)*lambda_f(j,i,k1)*XI33_f(j,i,k1)*(u_f(j-2,i,k1,1,index)/12.d0-u_f(j-1,i,k1,1,index)*2.d0/3.d0 &
             +u_f(j+1,i,k1,1,index)*2.d0/3.d0-u_f(j+2,i,k1,1,index)/12.d0))/l1/h1_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)&
             *mu_f(j,i,k1)*XI13_f(j,i,k1)*(u_f(j-2,i,k1,3,index)/12.d0-u_f(j-1,i,k1,3,index)*2.d0/3.d0+u_f(j+1,i,k1,3,index)&
             *2.d0/3.d0-u_f(j+2,i,k1,3,index)/12.d0))/l1/h3_f/h1_f+(bof(1,k1)*Jacobian_f(j,i,k1)*lambda_f(j,i,k1)*XI33_f(j,i,k1)&
             *(u_f(j,i-2,k1,2,index)/12.d0-u_f(j,i-1,k1,2,index)*2.d0/3.d0+u_f(j,i+1,k1,2,index)*2.d0/3.d0-u_f(j,i+2,k1,2,index)&
             /12.d0))/l2/h2_f/h3_f+(bof(1,k1)*Jacobian_f(j,i,k1)*mu_f(j,i,k1)*XI23_f(j,i,k1)*(u_f(j,i-2,k1,3,index)/12.d0&
             -u_f(j,i-1,k1,3,index)*2.d0/3.d0+u_f(j,i+1,k1,3,index)*2.d0/3.d0-u_f(j,i+2,k1,3,index)/12.d0))/l2/h2_f/h3_f
        end do
     end do
  end do

  ! scale L_f
  do j = -2,n2_f+3
     do i = -2,n1_f+3
        lh_f(i,j,1,1) = lh_f(i,j,1,1)*17.d0/48.d0*h3_f
        lh_f(i,j,1,2) = lh_f(i,j,1,2)*17.d0/48.d0*h3_f
        lh_f(i,j,1,3) = lh_f(i,j,1,3)*17.d0/48.d0*h3_f
     end do
  end do
  !
  do j = -2,n2_f+3
     do i = -2,n1_f+3
        do k = 1,5
           ! first set equation
           lh_f(i,j,1,1) = lh_f(i,j,1,1)+Jacobian_f(i,j,1)*((2.d0*mu_f(i,j,1)+lambda_f(i,j,1))*XI13_f(i,j,1)**2+mu_f(i,j,1)&
             *(XI23_f(i,j,1)**2+XI33_f(i,j,1)**2))*sbop_no_gp(k)*u_f(i,j,k,1,index)/h3_f+Jacobian_f(i,j,1)*(lambda_f(i,j,1)&
             +mu_f(i,j,1))*XI13_f(i,j,1)*XI23_f(i,j,1)*sbop_no_gp(k)*u_f(i,j,k,2,index)/h3_f+Jacobian_f(i,j,1)*(lambda_f(i,j,1)&
             +mu_f(i,j,1))*XI13_f(i,j,1)*XI33_f(i,j,1)*sbop_no_gp(k)*u_f(i,j,k,3,index)/h3_f
           ! second set equation
           lh_f(i,j,1,2) = lh_f(i,j,1,2)+Jacobian_f(i,j,1)*(lambda_f(i,j,1)+mu_f(i,j,1))*XI13_f(i,j,1)*XI23_f(i,j,1)&
             *sbop_no_gp(k)*u_f(i,j,k,1,index)/h3_f+Jacobian_f(i,j,1)*((2.d0*mu_f(i,j,1)+lambda_f(i,j,1))*XI23_f(i,j,1)**2&
             +mu_f(i,j,1)*(XI13_f(i,j,1)**2+XI33_f(i,j,1)**2))*sbop_no_gp(k)*u_f(i,j,k,2,index)/h3_f+Jacobian_f(i,j,1)&
             *(lambda_f(i,j,1)+mu_f(i,j,1))*XI23_f(i,j,1)*XI33_f(i,j,1)*sbop_no_gp(k)*u_f(i,j,k,3,index)/h3_f
           ! third set equation
           lh_f(i,j,1,3) = lh_f(i,j,1,3)+Jacobian_f(i,j,1)*(lambda_f(i,j,1)+mu_f(i,j,1))*XI13_f(i,j,1)*XI33_f(i,j,1)*sbop_no_gp(k)&
             *u_f(i,j,k,1,index)/h3_f+Jacobian_f(i,j,1)*(lambda_f(i,j,1)+mu_f(i,j,1))*XI23_f(i,j,1)*XI33_f(i,j,1)*sbop_no_gp(k)&
             *u_f(i,j,k,2,index)/h3_f+Jacobian_f(i,j,1)*((2.d0*mu_f(i,j,1)+lambda_f(i,j,1))*XI33_f(i,j,1)**2+mu_f(i,j,1)&
             *(XI13_f(i,j,1)**2+XI23_f(i,j,1)**2))*sbop_no_gp(k)*u_f(i,j,k,3,index)/h3_f
         end do
      end do
   end do
   !
   do k = -2,n2_f+3
      do i = -2,n1_f+3
         do j = -2,2
            ! 31 & 32
            ! first set equation
            lh_f(i,k,1,1) = lh_f(i,k,1,1)+Jacobian_f(i,k,1)*(2.d0*mu_f(i,k,1)+lambda_f(i,k,1))/l1*XI13_f(i,k,1)*ux_cof(j)&
              *u_f(i+j,k,1,1,index)/h1_f+Jacobian_f(i,k,1)*mu_f(i,k,1)/l1*XI23_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,2,index)/h1_f&
              +Jacobian_f(i,k,1)*mu_f(i,k,1)/l1*XI33_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,3,index)/h1_f+Jacobian_f(i,k,1)*mu_f(i,k,1)/l2&
              *XI23_f(i,k,1)*ux_cof(j)*u_f(i,k+j,1,1,index)/h2_f+Jacobian_f(i,k,1)*lambda_f(i,k,1)/l2*XI13_f(i,k,1)*ux_cof(j)&
              *u_f(i,k+j,1,2,index)/h2_f
            ! second set equation
            lh_f(i,k,1,2) = lh_f(i,k,1,2)+Jacobian_f(i,k,1)*lambda_f(i,k,1)/l1*XI23_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,1,index)/h1_f&
              +Jacobian_f(i,k,1)*mu_f(i,k,1)/l1*XI13_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,2,index)/h1_f+Jacobian_f(i,k,1)*mu_f(i,k,1)/l2&
              *XI13_f(i,k,1)*ux_cof(j)*u_f(i,k+j,1,1,index)/h2_f+Jacobian_f(i,k,1)*(2.d0*mu_f(i,k,1)+lambda_f(i,k,1))/l2&
              *XI23_f(i,k,1)*ux_cof(j)*u_f(i,k+j,1,2,index)/h2_f+Jacobian_f(i,k,1)*mu_f(i,k,1)/l2*XI33_f(i,k,1)*ux_cof(j)&
              *u_f(i,k+j,1,3,index)/h2_f
            ! third set equation
            lh_f(i,k,1,3) = lh_f(i,k,1,3)+Jacobian_f(i,k,1)*lambda_f(i,k,1)/l1*XI33_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,1,index)/h1_f&
              +Jacobian_f(i,k,1)*mu_f(i,k,1)/l1*XI13_f(i,k,1)*ux_cof(j)*u_f(i+j,k,1,3,index)/h1_f+Jacobian_f(i,k,1)&
              *lambda_f(i,k,1)/l2*XI33_f(i,k,1)*ux_cof(j)*u_f(i,k+j,1,2,index)/h2_f+Jacobian_f(i,k,1)*mu_f(i,k,1)/l2*XI23_f(i,k,1)&
              *ux_cof(j)*u_f(i,k+j,1,3,index)/h2_f
         end do
      end do
   end do
   ! now restrict it to the coarse grid
   do k = 1,n2_c
      do i = 1,n1_c
         do j = -4,2
            do l = -4,2
               ! first set
               Vass((k-1)*3*n1_c+(i-1)*3+1) = Vass((k-1)*3*n1_c+(i-1)*3+1)+Rop(j)*(Rop(l)*lh_f(2*i+l,2*k+j,1,1))
               ! second set
               Vass((k-1)*3*n1_c+(i-1)*3+2) = Vass((k-1)*3*n1_c+(i-1)*3+2)+Rop(j)*(Rop(l)*lh_f(2*i+l,2*k+j,1,2))
               ! third set
               Vass((k-1)*3*n1_c+(i-1)*3+3) = Vass((k-1)*3*n1_c+(i-1)*3+3)+Rop(j)*(Rop(l)*lh_f(2*i+l,2*k+j,1,3))
            end do
        end do
      end do
   end do
   !
!!$   open (unit = 7, file = "Vass.dat")
!!$   write(7,*)Vass
!!$   close(7)
  end subroutine Interface_RHS


    subroutine Interface_LHS(LHS,lh_c,lh_f,Jacobian_c,Jacobian_f, mu_c,mu_f,lambda_c,lambda_f,rho_c,rho_f,&
  XI13_c,XI23_c,XI33_c,XI13_f,XI23_f,XI33_f,P,Sb,Rop,sbop_no_gp,acof_no_gp,u_c,u_f,Mass_f1,ux_cof,ghcof,acof,bof,index) bind(c)
  use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
    integer,value,intent(in) :: index
    real(dp) :: int_cof
    real(dp), dimension (:,:,:), allocatable :: u_temp
    real(dp) :: LHS(1:n1_c*n2_c*3)
  real(dp) :: lh_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim)
  real(dp) :: lh_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim)
  real(dp) :: Jacobian_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp) :: mu_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp) :: lambda_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
  real(dp),dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI13_c,XI23_c,XI33_c,rho_c
  real(dp), dimension (0:4):: Sb
  real(dp) :: u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
  real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
  real(dp) :: Mass_f1(-2:n1_f+3,-2:n2_f+3),ux_cof(-2:2),ghcof(6)
  real(dp) :: acof(6,8,8), bof(4,6),sbop_no_gp(0:5),acof_no_gp(6,8,8)
  !real(dp) :: rho_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)
  real(dp) ,dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: lambda_f,mu_f,rho_f,Jacobian_f,XI13_f,XI23_f,XI33_f
  real(dp), dimension (-1:2) :: P
  real(dp), dimension (-4:2):: Rop
  integer :: i,j,k,i1,j1,k1,l
    ! allocate memory for temporary arrays
    allocate(u_temp(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1:dim))
    !
    u_temp = 0.d0
    do j = 1,n2_c
       do i = 1,n1_c
          u_temp(i,j,1) = u_c(i,j,n3_c+1,1,index)
          u_temp(i,j,2) = u_c(i,j,n3_c+1,2,index)
          u_temp(i,j,3) = u_c(i,j,n3_c+1,3,index)
       end do
    end do
    !
    int_cof = 17.d0/48.d0*h3_f*ghcof(1)/h3_c**2
    !
    LHS = 0.d0
    do l = 1,n2_c
       do k = 1,n1_c
          !
          do j = -4,2,2
             do i = -4,2,2
                do j1 = -1,2
                   do i1 = -1,2
                      ! first set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof &
                       +Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI13_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof+Rop(j)&
                       *Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI13_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof
                      ! second set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+ Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)*XI23_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof &
                       + Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI23_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof
                      ! third set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)*XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)*XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof
                   end do
                end do
             end do
          end do
          !
          do j = -4,2,2
             do j1 = -1,2
                ! first set equation
                LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI13_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI23_c(k,l+j/2+j1,n3_c)**2+XI33_c(k,l+j/2+j1,n3_c)**2))&
                 /rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI23_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,3)*int_cof
                ! second set equation
                LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI23_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI23_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI13_c(k,l+j/2+j1,n3_c)**2+XI33_c(k,l+j/2+j1,n3_c)**2))&
                 /rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI23_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,3)*int_cof
                ! third set equation
                LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 * P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI23_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI33_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI13_c(k,l+j/2+j1,n3_c)**2&
                 +XI23_c(k,l+j/2+j1,n3_c)**2))/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,3)*int_cof
             end do
          end do
          !
          do i = -4,2,2
             do i1 = -1,2
                ! first set equation
                LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)**2&
                 +mu_c(k+i/2+i1,l,n3_c)*(XI23_c(k+i/2+i1,l,n3_c)**2+XI33_c(k+i/2+i1,l,n3_c)**2))/rho_c(k+i/2+i1,l,n3_c)&
                 *u_temp(k+i/2+i1,l,1)*int_cof+ Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)&
                 *(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)*XI23_c(k+i/2+i1,l,n3_c)&
                 /rho_c(k+i/2+i1,l,n3_c)* u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)* u_temp(k+i/2+i1,l,3)*int_cof
                ! second set equation
                LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI23_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,1)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))&
                 *XI23_c(k+i/2+i1,l,n3_c)**2+mu_c(k+i/2+i1,l,n3_c)*(XI13_c(k+i/2+i1,l,n3_c)**2+XI33_c(k+i/2+i1,l,n3_c)**2))&
                 /rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI23_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,3)*int_cof
                ! third set equation
                LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,1)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI23_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))&
                 *XI33_c(k+i/2+i1,l,n3_c)**2+mu_c(k+i/2+i1,l,n3_c)*(XI13_c(k+i/2+i1,l,n3_c)**2+XI23_c(k+i/2+i1,l,n3_c)**2))&
                 /rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,3)*int_cof
             end do
          end do
          !
          ! first set equation
          LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
           +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)&
           *XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof
          ! second set equation
          LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof &
           +Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2&
           +mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof
          ! third set equation
          LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof &
           +Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)&
           *XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof
          !
          ! first set equation
          LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,1) &
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c*u_temp(k,l,2) &
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,3)
          ! second set equation
          LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c*u_temp(k,l,1)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,2)&
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,3)
          ! third set equation
          LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,1)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,2)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,3)
       end do
    end do
    !
!!$    open (unit = 7, file = "LHS.dat")
!!$    write(7,*)u_f
!!$    close(7)
  end subroutine Interface_LHS


  subroutine dgetrs_wrap(i1,i2,A,i3,piv,rhs,i4,info) bind(c)
  
  use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
  !integer, parameter :: dp=c_double
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

subroutine Update_interior(u_c_t,u_f_t,bof,ghcof,acof,acof_no_gp,&
     lh_f,Jacobian_f,mu_f,lambda_f,XI13_f,XI23_f,XI33_f,&
     lh_c,Jacobian_c,mu_c,lambda_c,XI13_c,XI23_c,XI33_c) bind(c)
  use iso_fortran_env
  use iso_c_binding
  use problemsetup_new_3d
  implicit none
  real(dp) ghcof(6), bof(4,6),acof(6,8,8),acof_no_gp(6,8,8)
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t,lh_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: u_f_t,lh_f
  real(dp),dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: Jacobian_f,mu_f,lambda_f,XI13_f,XI23_f,XI33_f
  real(dp),dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) ::  Jacobian_c,mu_c,lambda_c,XI13_c,XI23_c,XI33_c
  
  integer :: k1,k,j,i,m

    ! Difference operators in the interior of the domains
    ! fine mesh
    lh_f = 0.d0
    do k = 1,n3_f
       do j = 1,n2_f
          do i = 1,n1_f
             ! second derivative 11 & 22 & 12 & 21
             ! first set
             lh_f(i,j,k,1) = lh_f(i,j,k,1)+((-Jacobian_f(i-2,j,k)*(2.d0*mu_f(i-2,j,k)+lambda_f(i-2,j,k))/8.d0+Jacobian_f(i-1,j,k)&
               *(2.d0*mu_f(i-1,j,k)+lambda_f(i-1,j,k))/6.d0-Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/8.d0)&
               *u_f_t(i-2,j,k,1)+(Jacobian_f(i-2,j,k)*(2.d0*mu_f(i-2,j,k)+lambda_f(i-2,j,k))/6.d0+Jacobian_f(i-1,j,k)&
               *(2.d0*mu_f(i-1,j,k)+lambda_f(i-1,j,k))/2.d0+Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/2.d0 &
               + Jacobian_f(i+1,j,k)*(2.d0*mu_f(i+1,j,k)+lambda_f(i+1,j,k))/6.d0)*u_f_t(i-1,j,k,1)+(-Jacobian_f(i-2,j,k)&
               *(2.d0*mu_f(i-2,j,k)+lambda_f(i-2,j,k))/24.d0-Jacobian_f(i-1,j,k)*(2.d0*mu_f(i-1,j,k)+lambda_f(i-1,j,k))*5.d0/6.d0&
               -Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))*3.d0/4.d0-Jacobian_f(i+1,j,k)*(2.d0*mu_f(i+1,j,k)&
               +lambda_f(i+1,j,k))*5.d0/6.d0-Jacobian_f(i+2,j,k)*(2.d0*mu_f(i+2,j,k)+lambda_f(i+2,j,k))/24.d0)*u_f_t(i-0,j,k,1) &
               +(Jacobian_f(i-1,j,k)*(2.d0*mu_f(i-1,j,k)+lambda_f(i-1,j,k))/6.d0+Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)&
               +lambda_f(i,j,k))/2.d0+Jacobian_f(i+1,j,k)*(2.d0*mu_f(i+1,j,k)+lambda_f(i+1,j,k))/2.d0+Jacobian_f(i+2,j,k)&
               *(2.d0*mu_f(i+2,j,k)+lambda_f(i+2,j,k))/6.d0)*u_f_t(i+1,j,k,1)+(-Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)&
               +lambda_f(i,j,k))/8.d0+Jacobian_f(i+1,j,k)*(2.d0*mu_f(i+1,j,k)+lambda_f(i+1,j,k))/6.d0-Jacobian_f(i+2,j,k)&
               *(2.d0*mu_f(i+2,j,k)+lambda_f(i+2,j,k))/8.d0)*u_f_t(i+2,j,k,1))/h1_f**2/l1**2+((-Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)&
               /8.d0+Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)/6.d0-Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0)*u_f_t(i,j-2,k,1)&
               +(Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)/6.d0 + Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)/2.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0&
               +Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)/6.d0)*u_f_t(i,j-1,k,1)+(-Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)/24.d0&
               -Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)*5.d0/6.d0-Jacobian_f(i,j,k)*mu_f(i,j,k)*3.d0/4.d0-Jacobian_f(i,j+1,k)&
               *mu_f(i,j+1,k)*5.d0/6.d0-Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)/24.d0)*u_f_t(i,j-0,k,1)+(Jacobian_f(i,j-1,k)&
               *mu_f(i,j-1,k)/6.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)/2.d0+Jacobian_f(i,j+2,k)&
               *mu_f(i,j+2,k)/6.d0)*u_f_t(i,j+1,k,1)+(-Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0 + Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)/6.d0&
               -Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)/8.d0)*u_f_t(i,j+2,k,1))/h2_f**2/l2**2+(Jacobian_f(i-2,j,k)*lambda_f(i-2,j,k)&
               *(u_f_t(i-2,j-2,k,2)/12.d0-u_f_t(i-2,j-1,k,2)*2.d0/3.d0+u_f_t(i-2,j+1,k,2)*2.d0/3.d0 - u_f_t(i-2,j+2,k,2)/12.d0)&
               /12.d0-Jacobian_f(i-1,j,k)*lambda_f(i-1,j,k)*(u_f_t(i-1,j-2,k,2)/12.d0-u_f_t(i-1,j-1,k,2)*2.d0/3.d0 &
               +u_f_t(i-1,j+1,k,2)*2.d0/3.d0 - u_f_t(i-1,j+2,k,2)/12.d0)*2.d0/3.d0+Jacobian_f(i+1,j,k)*lambda_f(i+1,j,k)&
               *(u_f_t(i+1,j-2,k,2)/12.d0-u_f_t(i+1,j-1,k,2)*2.d0/3.d0+u_f_t(i+1,j+1,k,2)*2.d0/3.d0 - u_f_t(i+1,j+2,k,2)/12.d0)&
               *2.d0/3.d0-Jacobian_f(i+2,j,k)*lambda_f(i+2,j,k)*(u_f_t(i+2,j-2,k,2)/12.d0-u_f_t(i+2,j-1,k,2)*2.d0/3.d0 &
               +u_f_t(i+2,j+1,k,2)*2.d0/3.d0 - u_f_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_f/h2_f+(Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)&
               *(u_f_t(i-2,j-2,k,2)/12.d0-u_f_t(i-1,j-2,k,2)*2.d0/3.d0+u_f_t(i+1,j-2,k,2)*2.d0/3.d0 - u_f_t(i+2,j-2,k,2)/12.d0)&
               /12.d0-Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)*(u_f_t(i-2,j-1,k,2)/12.d0-u_f_t(i-1,j-1,k,2)*2.d0/3.d0+u_f_t(i+1,j-1,k,2)&
               *2.d0/3.d0-u_f_t(i+2,j-1,k,2)/12.d0)*2.d0/3.d0+Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)*(u_f_t(i-2,j+1,k,2)/12.d0&
               -u_f_t(i-1,j+1,k,2)*2.d0/3.d0+u_f_t(i+1,j+1,k,2)*2.d0/3.d0 - u_f_t(i+2,j+1,k,2)/12.d0)*2.d0/3.d0 &
               -Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)*(u_f_t(i-2,j+2,k,2)/12.d0-u_f_t(i-1,j+2,k,2)*2.d0/3.d0+u_f_t(i+1,j+2,k,2)&
               *2.d0/3.d0-u_f_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_f/h2_f
             ! second set
             lh_f(i,j,k,2) = lh_f(i,j,k,2)+((-Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/8.d0 + Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)/6.d0 &
               -Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0)*u_f_t(i-2,j,k,2)+(Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/6.d0+Jacobian_f(i-1,j,k)&
               *mu_f(i-1,j,k)/2.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)/6.d0)*u_f_t(i-1,j,k,2) &
               +(-Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/24.d0-Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)*5.d0/6.d0-Jacobian_f(i,j,k)&
               *mu_f(i,j,k)*3.d0/4.d0-Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)*5.d0/6.d0 -Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/24.d0)&
               *u_f_t(i-0,j,k,2)+(Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)/6.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i+1,j,k)&
               *mu_f(i+1,j,k)/2.d0+Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/6.d0)*u_f_t(i+1,j,k,2)+(-Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0&
               +Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)/6.d0-Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/8.d0)*u_f_t(i+2,j,k,2))/h1_f**2/l1**2&
               +((-Jacobian_f(i,j-2,k)*(2.d0*mu_f(i,j-2,k)+lambda_f(i,j-2,k))/8.d0+Jacobian_f(i,j-1,k)*(2.d0*mu_f(i,j-1,k)&
               +lambda_f(i,j-1,k))/6.d0-Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/8.d0)*u_f_t(i,j-2,k,2)&
               +(Jacobian_f(i,j-2,k)*(2.d0*mu_f(i,j-2,k)+lambda_f(i,j-2,k))/6.d0+Jacobian_f(i,j-1,k)*(2.d0*mu_f(i,j-1,k)&
               +lambda_f(i,j-1,k))/2.d0+Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/2.d0+Jacobian_f(i,j+1,k)&
               *(2.d0*mu_f(i,j+1,k)+lambda_f(i,j+1,k))/6.d0)*u_f_t(i,j-1,k,2)+(-Jacobian_f(i,j-2,k)*(2.d0*mu_f(i,j-2,k)&
               +lambda_f(i,j-2,k))/24.d0-Jacobian_f(i,j-1,k)*(2.d0*mu_f(i,j-1,k)+lambda_f(i,j-1,k))*5.d0/6.d0-Jacobian_f(i,j,k)&
               *(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))*3.d0/4.d0-Jacobian_f(i,j+1,k)*(2.d0*mu_f(i,j+1,k)+lambda_f(i,j+1,k))*5.d0/6.d0&
               -Jacobian_f(i,j+2,k)*(2.d0*mu_f(i,j+2,k)+lambda_f(i,j+2,k))/24.d0)*u_f_t(i,j-0,k,2)+(Jacobian_f(i,j-1,k)&
               *(2.d0*mu_f(i,j-1,k)+lambda_f(i,j-1,k))/6.d0+Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/2.d0 &
               +Jacobian_f(i,j+1,k)*(2.d0*mu_f(i,j+1,k)+lambda_f(i,j+1,k))/2.d0+Jacobian_f(i,j+2,k)*(2.d0*mu_f(i,j+2,k)&
               +lambda_f(i,j+2,k))/6.d0)*u_f_t(i,j+1,k,2)+(-Jacobian_f(i,j,k)*(2.d0*mu_f(i,j,k)+lambda_f(i,j,k))/8.d0 &
               +Jacobian_f(i,j+1,k)*(2.d0*mu_f(i,j+1,k)+lambda_f(i,j+1,k))/6.d0-Jacobian_f(i,j+2,k)*(2.d0*mu_f(i,j+2,k)&
               +lambda_f(i,j+2,k))/8.d0)*u_f_t(i,j+2,k,2))/h2_f**2/l2**2+(Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)*(u_f_t(i-2,j-2,k,1)&
               /12.d0-u_f_t(i-2,j-1,k,1)*2.d0/3.d0+u_f_t(i-2,j+1,k,1)*2.d0/3.d0 - u_f_t(i-2,j+2,k,1)/12.d0)/12.d0 &
               -Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)*(u_f_t(i-1,j-2,k,1)/12.d0-u_f_t(i-1,j-1,k,1)*2.d0/3.d0+u_f_t(i-1,j+1,k,1)&
               *2.d0/3.d0-u_f_t(i-1,j+2,k,1)/12.d0)*2.d0/3.d0+Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)*(u_f_t(i+1,j-2,k,1)/12.d0&
               -u_f_t(i+1,j-1,k,1)*2.d0/3.d0+u_f_t(i+1,j+1,k,1)*2.d0/3.d0-u_f_t(i+1,j+2,k,1)/12.d0)*2.d0/3.d0-Jacobian_f(i+2,j,k)&
               *mu_f(i+2,j,k)*(u_f_t(i+2,j-2,k,1)/12.d0-u_f_t(i+2,j-1,k,1)*2.d0/3.d0+u_f_t(i+2,j+1,k,1)*2.d0/3.d0&
               -u_f_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_f/h2_f+(Jacobian_f(i,j-2,k)*lambda_f(i,j-2,k)*(u_f_t(i-2,j-2,k,1)/12.d0&
               -u_f_t(i-1,j-2,k,1)*2.d0/3.d0+u_f_t(i+1,j-2,k,1)*2.d0/3.d0-u_f_t(i+2,j-2,k,1)/12.d0)/12.d0-Jacobian_f(i,j-1,k)&
               *lambda_f(i,j-1,k)*(u_f_t(i-2,j-1,k,1)/12.d0-u_f_t(i-1,j-1,k,1)*2.d0/3.d0+u_f_t(i+1,j-1,k,1)*2.d0/3.d0&
               -u_f_t(i+2,j-1,k,1)/12.d0)*2.d0/3.d0+Jacobian_f(i,j+1,k)*lambda_f(i,j+1,k)*(u_f_t(i-2,j+1,k,1)/12.d0&
               -u_f_t(i-1,j+1,k,1)*2.d0/3.d0+u_f_t(i+1,j+1,k,1)*2.d0/3.d0-u_f_t(i+2,j+1,k,1)/12.d0)*2.d0/3.d0-Jacobian_f(i,j+2,k)&
               *lambda_f(i,j+2,k)*(u_f_t(i-2,j+2,k,1)/12.d0-u_f_t(i-1,j+2,k,1)*2.d0/3.d0+u_f_t(i+1,j+2,k,1)*2.d0/3.d0&
               -u_f_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_f/h2_f
             ! third set
             lh_f(i,j,k,3) = lh_f(i,j,k,3)+((-Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/8.d0 + Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)/6.d0 &
               -Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0)*u_f_t(i-2,j,k,3)+(Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/6.d0+Jacobian_f(i-1,j,k)&
               *mu_f(i-1,j,k)/2.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)/6.d0)*u_f_t(i-1,j,k,3) &
               +(-Jacobian_f(i-2,j,k)*mu_f(i-2,j,k)/24.d0-Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)*5.d0/6.d0-Jacobian_f(i,j,k)&
               *mu_f(i,j,k)*3.d0/4.d0-Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)*5.d0/6.d0-Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/24.d0)&
               *u_f_t(i-0,j,k,3)+(Jacobian_f(i-1,j,k)*mu_f(i-1,j,k)/6.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i+1,j,k)&
               *mu_f(i+1,j,k)/2.d0+Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/6.d0)*u_f_t(i+1,j,k,3)+(-Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0&
               +Jacobian_f(i+1,j,k)*mu_f(i+1,j,k)/6.d0-Jacobian_f(i+2,j,k)*mu_f(i+2,j,k)/8.d0)*u_f_t(i+2,j,k,3))/h1_f**2/l1**2&
               +((-Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)/8.d0+Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)/6.d0-Jacobian_f(i,j,k)*mu_f(i,j,k)&
               /8.d0)*u_f_t(i,j-2,k,3)+(Jacobian_f(i,j-2,k)*mu_f(i,j-2,k)/6.d0 + Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)/2.d0&
               +Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)/6.d0)*u_f_t(i,j-1,k,3)+(-Jacobian_f(i,j-2,k)&
               *mu_f(i,j-2,k)/24.d0 - Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)*5.d0/6.d0-Jacobian_f(i,j,k)*mu_f(i,j,k)*3.d0/4.d0 &
               -Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)*5.d0/6.d0 -Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)/24.d0)*u_f_t(i,j-0,k,3) &
               +(Jacobian_f(i,j-1,k)*mu_f(i,j-1,k)/6.d0+Jacobian_f(i,j,k)*mu_f(i,j,k)/2.d0+Jacobian_f(i,j+1,k)*mu_f(i,j+1,k)/2.d0&
               +Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)/6.d0)*u_f_t(i,j+1,k,3)+(-Jacobian_f(i,j,k)*mu_f(i,j,k)/8.d0+Jacobian_f(i,j+1,k)&
               *mu_f(i,j+1,k)/6.d0-Jacobian_f(i,j+2,k)*mu_f(i,j+2,k)/8.d0)*u_f_t(i,j+2,k,3))/h2_f**2/l2**2
          end do
       end do
    end do
    !
    do i = 7,n3_f-6
       do j = 1,n2_f
          do k = 1,n1_f
             ! second derivative 33
             ! first set equation
             lh_f(k,j,i,1) = lh_f(k,j,i,1)+((-Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)**2&
               +mu_f(k,j,i-2)*(XI23_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))/8.d0 + Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI13_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI23_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))/6.d0&
               -Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)**2+mu_f(k,j,i)*(XI23_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/8.d0)*u_f_t(k,j,i-2,1)+(Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))&
               *XI13_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI23_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))/6.d0+Jacobian_f(k,j,i-1)&
               *((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI23_f(k,j,i-1)**2&
               +XI33_f(k,j,i-1)**2))/2.d0+Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)**2+mu_f(k,j,i)&
               *(XI23_f(k,j,i)**2+XI33_f(k,j,i)**2))/2.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))&
               *XI13_f(k,j,i+1)**2+mu_f(k,j,i+1)*(XI23_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/6.d0)*u_f_t(k,j,i-1,1) &
               +(-Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)**2+mu_f(k,j,i-2)&
               *(XI23_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))/24.d0 - Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))&
               *XI13_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI23_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))*5.d0/6.d0-Jacobian_f(k,j,i)&
               *((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)**2+mu_f(k,j,i)*(XI23_f(k,j,i)**2+XI33_f(k,j,i)**2))*3.d0/4.d0&
               -Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)**2+mu_f(k,j,i+1)*(XI23_f(k,j,i+1)**2&
               +XI33_f(k,j,i+1)**2))*5.d0/6.d0 -Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)**2&
               +mu_f(k,j,i+2)*(XI23_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/24.d0)*u_f_t(k,j,i,1)+(Jacobian_f(k,j,i-1)&
               *((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI23_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))&
               /6.d0+Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)**2+mu_f(k,j,i)*(XI23_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/2.d0 + Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI23_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/2.d0+Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI13_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI23_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/6.d0)*u_f_t(k,j,i+1,1)&
               +(-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)**2+mu_f(k,j,i)*(XI23_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/8.d0 + Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI23_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/6.d0-Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI13_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI23_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/8.d0)&
               *u_f_t(k,j,i+2,1))/h3_f**2+((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI23_f(k,j,i-2)&
               /8.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/6.d0-Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,2)+(Jacobian_f(k,j,i-2)&
               *(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI23_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI13_f(k,j,i)*XI23_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
               *XI23_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,2)+(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)&
               *XI23_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)&
               *5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)*3.d0/4.d0-Jacobian_f(k,j,i+1)&
               *(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI23_f(k,j,i+1)*5.d0/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,2)+(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI13_f(k,j,i)*XI23_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
               *XI23_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/6.d0)&
               *u_f_t(k,j,i+1,2)+(-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)/8.d0 &
               +Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI23_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)&
               *(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/8.d0)*u_f_t(k,j,i+2,2))/h3_f**2&
               +((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI33_f(k,j,i-2)/8.d0+Jacobian_f(k,j,i-1)&
               *(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,3)+(Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)&
               +lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI33_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))&
               *XI13_f(k,j,i-1)*XI33_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)&
               /2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,3) &
               +(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI33_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)&
               *(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)*5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)*3.d0/4.d0-Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))&
               *XI13_f(k,j,i+1)*XI33_f(k,j,i+1)*5.d0/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)&
               *XI33_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,3)+(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)&
               *XI33_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)/2.d0 &
               +Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI33_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)&
               *(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI33_f(k,j,i+2)/6.d0)*u_f_t(k,j,i+1,3)+(-Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)/8.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)&
               +lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))&
               *XI13_f(k,j,i+2)*XI33_f(k,j,i+2)/8.d0)*u_f_t(k,j,i+2,3))/h3_f**2
             ! second set of equation
             lh_f(k,j,i,2) = lh_f(k,j,i,2)+((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)&
               *XI23_f(k,j,i-2)/8.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/6.d0&
               -Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,1)&
               +(Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI23_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)&
               *(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))&
               *XI13_f(k,j,i+1)*XI23_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,1)+(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))&
               *XI13_f(k,j,i-2)*XI23_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)&
               *XI23_f(k,j,i-1)*5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)*3.d0/4.d0&
               -Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI23_f(k,j,i+1)*5.d0/6.d0&
               -Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,1)&
               +(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI23_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI23_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)&
               +lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI23_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))&
               *XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/6.d0)*u_f_t(k,j,i+1,1)+(-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI13_f(k,j,i)*XI23_f(k,j,i)/8.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
               *XI23_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI23_f(k,j,i+2)/8.d0)&
               *u_f_t(k,j,i+2,1))/h3_f**2+((-Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)**2&
               +mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))/8.d0 + Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI23_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))/6.d0&
               -Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/8.d0)*u_f_t(k,j,i-2,2)+(Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))&
               *XI23_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))/6.d0+Jacobian_f(k,j,i-1)&
               *((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))&
               /2.d0 + Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/2.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/6.d0)*u_f_t(k,j,i-1,2)+(-Jacobian_f(k,j,i-2)&
               *((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2+XI33_f(k,j,i-2)**2))&
               /24.d0-Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)**2+mu_f(k,j,i-1)&
               *(XI13_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))*5.d0/6.d0-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI23_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2+XI33_f(k,j,i)**2))*3.d0/4.d0-Jacobian_f(k,j,i+1)&
               *((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)**2+mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))&
               *5.d0/6.d0-Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)**2+mu_f(k,j,i+2)&
               *(XI13_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/24.d0)*u_f_t(k,j,i,2)+(Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI23_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2+XI33_f(k,j,i-1)**2))/6.d0&
               +Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/2.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/2.d0+Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI23_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI13_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/6.d0)*u_f_t(k,j,i+1,2)&
               +(-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI33_f(k,j,i)**2))/8.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI33_f(k,j,i+1)**2))/6.d0-Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI23_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI13_f(k,j,i+2)**2+XI33_f(k,j,i+2)**2))/8.d0)&
               *u_f_t(k,j,i+2,2))/h3_f**2+((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)*XI33_f(k,j,i-2)&
               /8.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0-Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,3)+(Jacobian_f(k,j,i-2)&
               *(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)*XI33_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI23_f(k,j,i)*XI33_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)&
               *XI33_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,3)+(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)&
               *XI33_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)&
               *5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)*3.d0/4.d0-Jacobian_f(k,j,i+1)&
               *(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)*XI33_f(k,j,i+1)*5.d0/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI23_f(k,j,i+2)*XI33_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,3)+(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI23_f(k,j,i)*XI33_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)&
               *XI33_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)*XI33_f(k,j,i+2)/6.d0)&
               *u_f_t(k,j,i+1,3)+(-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)/8.d0 &
               +Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)&
               *(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)*XI33_f(k,j,i+2)/8.d0)*u_f_t(k,j,i+2,3))/h3_f**2
             ! third set equation
             lh_f(k,j,i,3) = lh_f(k,j,i,3)+((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI33_f(k,j,i-2)&
               /8.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0-Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,1)+(Jacobian_f(k,j,i-2)&
               *(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)*XI33_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI13_f(k,j,i)*XI33_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
               *XI33_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,1)+(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI13_f(k,j,i-2)&
               *XI33_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)&
               *5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)*3.d0/4.d0-Jacobian_f(k,j,i+1)&
               *(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI33_f(k,j,i+1)*5.d0/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI33_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,1)+(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI13_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))&
               *XI13_f(k,j,i)*XI33_f(k,j,i)/2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
               *XI33_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI33_f(k,j,i+2)/6.d0)&
               *u_f_t(k,j,i+1,1)+(-Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI13_f(k,j,i)*XI33_f(k,j,i)/8.d0&
               +Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)&
               *(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)*XI33_f(k,j,i+2)/8.d0)*u_f_t(k,j,i+2,1))/h3_f**2&
               +((-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)*XI33_f(k,j,i-2)/8.d0+Jacobian_f(k,j,i-1)&
               *(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)/8.d0)*u_f_t(k,j,i-2,2)+(Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)&
               +lambda_f(k,j,i-2))*XI23_f(k,j,i-2)*XI33_f(k,j,i-2)/6.d0+Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))&
               *XI23_f(k,j,i-1)*XI33_f(k,j,i-1)/2.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)&
               /2.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0)*u_f_t(k,j,i-1,2)&
               +(-Jacobian_f(k,j,i-2)*(mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI23_f(k,j,i-2)*XI33_f(k,j,i-2)/24.d0-Jacobian_f(k,j,i-1)&
               *(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)*XI33_f(k,j,i-1)*5.d0/6.d0-Jacobian_f(k,j,i)*(mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)*3.d0/4.d0-Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))&
               *XI23_f(k,j,i+1)*XI33_f(k,j,i+1)*5.d0/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)&
               *XI33_f(k,j,i+2)/24.d0)*u_f_t(k,j,i,2)+(Jacobian_f(k,j,i-1)*(mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)&
               *XI33_f(k,j,i-1)/6.d0+Jacobian_f(k,j,i)*(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)/2.d0 &
               +Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)*XI33_f(k,j,i+1)/2.d0+Jacobian_f(k,j,i+2)&
               *(mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)*XI33_f(k,j,i+2)/6.d0)*u_f_t(k,j,i+1,2)+(-Jacobian_f(k,j,i)&
               *(mu_f(k,j,i)+lambda_f(k,j,i))*XI23_f(k,j,i)*XI33_f(k,j,i)/8.d0+Jacobian_f(k,j,i+1)*(mu_f(k,j,i+1)&
               +lambda_f(k,j,i+1))*XI23_f(k,j,i+1)*XI33_f(k,j,i+1)/6.d0-Jacobian_f(k,j,i+2)*(mu_f(k,j,i+2)+lambda_f(k,j,i+2))&
               *XI23_f(k,j,i+2)*XI33_f(k,j,i+2)/8.d0)*u_f_t(k,j,i+2,2))/h3_f**2+((-Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)&
               +lambda_f(k,j,i-2))*XI33_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2+XI23_f(k,j,i-2)**2))/8.d0&
               +Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI33_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2&
               +XI23_f(k,j,i-1)**2))/6.d0-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI33_f(k,j,i)**2+mu_f(k,j,i)&
               *(XI13_f(k,j,i)**2+XI23_f(k,j,i)**2))/8.d0)*u_f_t(k,j,i-2,3)+(Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)&
               +lambda_f(k,j,i-2))*XI33_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2+XI23_f(k,j,i-2)**2))/6.d0&
               +Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI33_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2&
               +XI23_f(k,j,i-1)**2))/2.d0+Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI33_f(k,j,i)**2+mu_f(k,j,i)&
               *(XI13_f(k,j,i)**2+XI23_f(k,j,i)**2))/2.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))&
               *XI33_f(k,j,i+1)**2+mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI23_f(k,j,i+1)**2))/6.d0)*u_f_t(k,j,i-1,3)&
               +(-Jacobian_f(k,j,i-2)*((2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))*XI33_f(k,j,i-2)**2+mu_f(k,j,i-2)*(XI13_f(k,j,i-2)**2&
               +XI23_f(k,j,i-2)**2))/24.d0-Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI33_f(k,j,i-1)**2&
               +mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2+XI23_f(k,j,i-1)**2))*5.d0/6.d0-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)&
               +lambda_f(k,j,i))*XI33_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2+XI23_f(k,j,i)**2))*3.d0/4.d0-Jacobian_f(k,j,i+1)&
               *((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI33_f(k,j,i+1)**2+mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI23_f(k,j,i+1)**2))&
               *5.d0/6.d0-Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI33_f(k,j,i+2)**2+mu_f(k,j,i+2)&
               *(XI13_f(k,j,i+2)**2+XI23_f(k,j,i+2)**2))/24.d0)*u_f_t(k,j,i,3)+(Jacobian_f(k,j,i-1)*((2.d0*mu_f(k,j,i-1)&
               +lambda_f(k,j,i-1))*XI33_f(k,j,i-1)**2+mu_f(k,j,i-1)*(XI13_f(k,j,i-1)**2+XI23_f(k,j,i-1)**2))/6.d0&
               +Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI33_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI23_f(k,j,i)**2))/2.d0+Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI33_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI23_f(k,j,i+1)**2))/2.d0+Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI33_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI13_f(k,j,i+2)**2+XI23_f(k,j,i+2)**2))/6.d0)*u_f_t(k,j,i+1,3)&
               +(-Jacobian_f(k,j,i)*((2.d0*mu_f(k,j,i)+lambda_f(k,j,i))*XI33_f(k,j,i)**2+mu_f(k,j,i)*(XI13_f(k,j,i)**2&
               +XI23_f(k,j,i)**2))/8.d0 + Jacobian_f(k,j,i+1)*((2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI33_f(k,j,i+1)**2&
               +mu_f(k,j,i+1)*(XI13_f(k,j,i+1)**2+XI23_f(k,j,i+1)**2))/6.d0-Jacobian_f(k,j,i+2)*((2.d0*mu_f(k,j,i+2)&
               +lambda_f(k,j,i+2))*XI33_f(k,j,i+2)**2+mu_f(k,j,i+2)*(XI13_f(k,j,i+2)**2+XI23_f(k,j,i+2)**2))/8.d0)&
               *u_f_t(k,j,i+2,3))/h3_f**2
          end do
       end do
    end do
    do j = 1,n2_f
       do k = 1,n1_f
          do i = 1,6
             do k1 = 1,8
                do m = 1,8
                   ! second derivative 33
                   ! first set equation
                   lh_f(k,j,i,1) = lh_f(k,j,i,1) +(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*((2.d0*mu_f(k,j,m)+lambda_f(k,j,m))&
                     *XI13_f(k,j,m)**2+mu_f(k,j,m)*(XI23_f(k,j,m)**2+XI33_f(k,j,m)**2))*u_f_t(k,j,k1,1))/h3_f**2&
                     +(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)*XI23_f(k,j,m)&
                     *u_f_t(k,j,k1,2))/h3_f**2+(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                     *XI33_f(k,j,m)*u_f_t(k,j,k1,3))/h3_f**2

                   lh_f(k,j,n3_f+1-i,1) = lh_f(k,j,n3_f+1-i,1) +(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*((2.d0*mu_f(k,j,n3_f+1-m)&
                     +lambda_f(k,j,n3_f+1-m))*XI13_f(k,j,n3_f+1-m)**2+mu_f(k,j,n3_f+1-m)*(XI23_f(k,j,n3_f+1-m)**2&
                     +XI33_f(k,j,n3_f+1-m)**2))*u_f_t(k,j,n3_f+1-k1,1))/h3_f**2+(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)&
                     *(mu_f(k,j,n3_f+1-m)+lambda_f(k,j,n3_f+1-m))*XI13_f(k,j,n3_f+1-m)*XI23_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,2))&
                     /h3_f**2+(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*(mu_f(k,j,n3_f+1-m)+lambda_f(k,j,n3_f+1-m))&
                     *XI13_f(k,j,n3_f+1-m)*XI33_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,3))/h3_f**2
                   ! second set equation
                   lh_f(k,j,i,2) = lh_f(k,j,i,2) +(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                     *XI23_f(k,j,m)*u_f_t(k,j,k1,1))/h3_f**2+(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*((2.d0*mu_f(k,j,m)&
                     +lambda_f(k,j,m))*XI23_f(k,j,m)**2+mu_f(k,j,m)*(XI13_f(k,j,m)**2+XI33_f(k,j,m)**2))*u_f_t(k,j,k1,2))/h3_f**2&
                     +(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI23_f(k,j,m)*XI33_f(k,j,m)&
                     *u_f_t(k,j,k1,3))/h3_f**2

                   lh_f(k,j,n3_f+1-i,2) = lh_f(k,j,n3_f+1-i,2) +(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*(mu_f(k,j,n3_f+1-m)&
                     +lambda_f(k,j,n3_f+1-m))*XI13_f(k,j,n3_f+1-m)*XI23_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,1))/h3_f**2&
                     +(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*((2.d0*mu_f(k,j,n3_f+1-m)+lambda_f(k,j,n3_f+1-m))&
                     *XI23_f(k,j,n3_f+1-m)**2+mu_f(k,j,n3_f+1-m)*(XI13_f(k,j,n3_f+1-m)**2+XI33_f(k,j,n3_f+1-m)**2))&
                     *u_f_t(k,j,n3_f+1-k1,2))/h3_f**2+(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*(mu_f(k,j,n3_f+1-m)&
                     +lambda_f(k,j,n3_f+1-m))*XI23_f(k,j,n3_f+1-m)*XI33_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,3))/h3_f**2
                   ! third set equation
                   lh_f(k,j,i,3) = lh_f(k,j,i,3) +(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))*XI13_f(k,j,m)&
                     *XI33_f(k,j,m)*u_f_t(k,j,k1,1))/h3_f**2+(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*(mu_f(k,j,m)+lambda_f(k,j,m))&
                     *XI23_f(k,j,m)*XI33_f(k,j,m)*u_f_t(k,j,k1,2))/h3_f**2+(acof_no_gp(i,k1,m)*Jacobian_f(k,j,m)*((2.d0*mu_f(k,j,m)&
                     +lambda_f(k,j,m))*XI33_f(k,j,m)**2+mu_f(k,j,m)*(XI13_f(k,j,m)**2+XI23_f(k,j,m)**2))*u_f_t(k,j,k1,3))/h3_f**2

                   lh_f(k,j,n3_f+1-i,3) = lh_f(k,j,n3_f+1-i,3) +(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*(mu_f(k,j,n3_f+1-m)&
                     +lambda_f(k,j,n3_f+1-m))*XI13_f(k,j,n3_f+1-m)*XI33_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,1))/h3_f**2&
                     +(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)*(mu_f(k,j,n3_f+1-m)+lambda_f(k,j,n3_f+1-m))*XI23_f(k,j,n3_f+1-m)&
                     *XI33_f(k,j,n3_f+1-m)*u_f_t(k,j,n3_f+1-k1,2))/h3_f**2+(acof(i,k1,m)*Jacobian_f(k,j,n3_f+1-m)&
                     *((2.d0*mu_f(k,j,n3_f+1-m)+lambda_f(k,j,n3_f+1-m))*XI33_f(k,j,n3_f+1-m)**2+mu_f(k,j,n3_f+1-m)&
                     *(XI13_f(k,j,n3_f+1-m)**2+XI23_f(k,j,n3_f+1-m)**2))*u_f_t(k,j,n3_f+1-k1,3))/h3_f**2
                end do
             end do
          end do
          ! first set equation
          lh_f(k,j,n3_f,1) = lh_f(k,j,n3_f,1) + (u_f_t(k,j,n3_f+1,1)*ghcof(1)*Jacobian_f(k,j,n3_f)*((2.d0*mu_f(k,j,n3_f)&
            +lambda_f(k,j,n3_f))*XI13_f(k,j,n3_f)**2+mu_f(k,j,n3_f)*(XI23_f(k,j,n3_f)**2+XI33_f(k,j,n3_f)**2)))/h3_f**2&
            +(u_f_t(k,j,n3_f+1,2)*ghcof(1)*Jacobian_f(k,j,n3_f)*(mu_f(k,j,n3_f)+lambda_f(k,j,n3_f))*XI13_f(k,j,n3_f)&
            *XI23_f(k,j,n3_f))/h3_f**2+(u_f_t(k,j,n3_f+1,3)*ghcof(1)*Jacobian_f(k,j,n3_f)*(mu_f(k,j,n3_f)+lambda_f(k,j,n3_f))&
            *XI13_f(k,j,n3_f)*XI33_f(k,j,n3_f))/h3_f**2
          ! second set equation
          lh_f(k,j,n3_f,2) = lh_f(k,j,n3_f,2) + (u_f_t(k,j,n3_f+1,1)*ghcof(1)*Jacobian_f(k,j,n3_f)*(mu_f(k,j,n3_f)&
            +lambda_f(k,j,n3_f))*XI13_f(k,j,n3_f)*XI23_f(k,j,n3_f))/h3_f**2+(u_f_t(k,j,n3_f+1,2)*ghcof(1)*Jacobian_f(k,j,n3_f)&
            *((2.d0*mu_f(k,j,n3_f) +lambda_f(k,j,n3_f))*XI23_f(k,j,n3_f)**2+mu_f(k,j,n3_f)*(XI13_f(k,j,n3_f)**2&
            +XI33_f(k,j,n3_f)**2)))/h3_f**2+(u_f_t(k,j,n3_f+1,3)*ghcof(1)*Jacobian_f(k,j,n3_f)*(mu_f(k,j,n3_f)+lambda_f(k,j,n3_f))&
            *XI23_f(k,j,n3_f)*XI33_f(k,j,n3_f))/h3_f**2
          ! third set equation
          lh_f(k,j,n3_f,3) = lh_f(k,j,n3_f,3) + (u_f_t(k,j,n3_f+1,1)*ghcof(1)*Jacobian_f(k,j,n3_f)*(mu_f(k,j,n3_f)&
            +lambda_f(k,j,n3_f))*XI13_f(k,j,n3_f)*XI33_f(k,j,n3_f))/h3_f**2+(u_f_t(k,j,n3_f+1,2)*ghcof(1)*Jacobian_f(k,j,n3_f)&
            *(mu_f(k,j,n3_f)+lambda_f(k,j,n3_f))*XI23_f(k,j,n3_f)*XI33_f(k,j,n3_f))/h3_f**2+(u_f_t(k,j,n3_f+1,3)*ghcof(1)&
            *Jacobian_f(k,j,n3_f)*((2.d0*mu_f(k,j,n3_f)+lambda_f(k,j,n3_f))*XI33_f(k,j,n3_f)**2+mu_f(k,j,n3_f)&
            *(XI13_f(k,j,n3_f)**2+XI23_f(k,j,n3_f)**2)))/h3_f**2
       end do
    end do
    ! mixed derivative 13
    do k = 1,4
       do i = 1,n2_f
          do j = 1,n1_f
             do k1 = 1,6
                ! mixed derivative 13 & 23
                ! first set equation
                lh_f(j,i,k,1) = lh_f(j,i,k,1)+(Jacobian_f(j-2,i,k)*(2.d0*mu_f(j-2,i,k)+lambda_f(j-2,i,k))*XI13_f(j-2,i,k)*bof(k,k1)&
                  *u_f_t(j-2,i,k1,1)/12.d0-Jacobian_f(j-1,i,k)*(2.d0*mu_f(j-1,i,k)+lambda_f(j-1,i,k))*XI13_f(j-1,i,k)*bof(k,k1)&
                  *u_f_t(j-1,i,k1,1)*2.d0/3.d0+Jacobian_f(j+1,i,k)*(2.d0*mu_f(j+1,i,k)+lambda_f(j+1,i,k))*XI13_f(j+1,i,k)*bof(k,k1)&
                  *u_f_t(j+1,i,k1,1)*2.d0/3.d0-Jacobian_f(j+2,i,k)*(2.d0*mu_f(j+2,i,k)+lambda_f(j+2,i,k))*XI13_f(j+2,i,k)*bof(k,k1)&
                  *u_f_t(j+2,i,k1,1)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j-2,i,k)*lambda_f(j-2,i,k)*XI23_f(j-2,i,k)*bof(k,k1)&
                  *u_f_t(j-2,i,k1,2)/12.d0-Jacobian_f(j-1,i,k)*lambda_f(j-1,i,k)*XI23_f(j-1,i,k)*bof(k,k1)*u_f_t(j-1,i,k1,2)&
                  *2.d0/3.d0+Jacobian_f(j+1,i,k)*lambda_f(j+1,i,k)*XI23_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,2)*2.d0/3.d0 &
                  -Jacobian_f(j+2,i,k)*lambda_f(j+2,i,k)*XI23_f(j+2,i,k)*bof(k,k1)*u_f_t(j+2,i,k1,2)/12.d0)/l1/h1_f/h3_f&
                  +(Jacobian_f(j-2,i,k)*lambda_f(j-2,i,k)*XI33_f(j-2,i,k)*bof(k,k1)*u_f_t(j-2,i,k1,3)/12.d0-Jacobian_f(j-1,i,k)&
                  *lambda_f(j-1,i,k)*XI33_f(j-1,i,k)*bof(k,k1)*u_f_t(j-1,i,k1,3)*2.d0/3.d0+Jacobian_f(j+1,i,k)*lambda_f(j+1,i,k)&
                  *XI33_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,3)*2.d0/3.d0-Jacobian_f(j+2,i,k)*lambda_f(j+2,i,k)*XI33_f(j+2,i,k)&
                  *bof(k,k1)*u_f_t(j+2,i,k1,3)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI23_f(j,i-2,k)*bof(k,k1)&
                  *u_f_t(j,i-2,k1,1)/12.d0-Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI23_f(j,i-1,k)*bof(k,k1)*u_f_t(j,i-1,k1,1)*2.d0/3.d0 &
                  +Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI23_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,1)*2.d0/3.d0-Jacobian_f(j,i+2,k)&
                  *mu_f(j,i+2,k)*XI23_f(j,i+2,k)*bof(k,k1)*u_f_t(j,i+2,k1,1)/12.d0)/l2/h2_f/h3_f+(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)&
                  *XI13_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,2)/12.d0-Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI13_f(j,i-1,k)*bof(k,k1)&
                  *u_f_t(j,i-1,k1,2)*2.d0/3.d0+Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI13_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,2)&
                  *2.d0/3.d0-Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI13_f(j,i+2,k)*bof(k,k1)*u_f_t(j,i+2,k1,2)/12.d0)/l2/h2_f/h3_f

                lh_f(j,i,n3_f+1-k,1) = lh_f(j,i,n3_f+1-k,1)+(-Jacobian_f(j-2,i,n3_f+1-k)*(2.d0*mu_f(j-2,i,n3_f+1-k)&
                  +lambda_f(j-2,i,n3_f+1-k))*XI13_f(j-2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,1)/12.d0&
                  +Jacobian_f(j-1,i,n3_f+1-k)*(2.d0*mu_f(j-1,i,n3_f+1-k)+lambda_f(j-1,i,n3_f+1-k))*XI13_f(j-1,i,n3_f+1-k)&
                  *bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,1)*2.d0/3.d0-Jacobian_f(j+1,i,n3_f+1-k)*(2.d0*mu_f(j+1,i,n3_f+1-k)&
                  +lambda_f(j+1,i,n3_f+1-k))*XI13_f(j+1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,1)*2.d0/3.d0 &
                  +Jacobian_f(j+2,i,n3_f+1-k)*(2.d0*mu_f(j+2,i,n3_f+1-k)+lambda_f(j+2,i,n3_f+1-k))*XI13_f(j+2,i,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j+2,i,n3_f+1-k1,1)/12.d0)/l1/h1_f/h3_f+(-Jacobian_f(j-2,i,n3_f+1-k)*lambda_f(j-2,i,n3_f+1-k)&
                  *XI23_f(j-2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,2)/12.d0+Jacobian_f(j-1,i,n3_f+1-k)&
                  *lambda_f(j-1,i,n3_f+1-k)*XI23_f(j-1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,2)*2.d0/3.d0&
                  -Jacobian_f(j+1,i,n3_f+1-k)*lambda_f(j+1,i,n3_f+1-k)*XI23_f(j+1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,2)&
                  *2.d0/3.d0+Jacobian_f(j+2,i,n3_f+1-k)*lambda_f(j+2,i,n3_f+1-k)*XI23_f(j+2,i,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j+2,i,n3_f+1-k1,2)/12.d0)/l1/h1_f/h3_f+(-Jacobian_f(j-2,i,n3_f+1-k)*lambda_f(j-2,i,n3_f+1-k)&
                  *XI33_f(j-2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,3)/12.d0+Jacobian_f(j-1,i,n3_f+1-k)&
                  *lambda_f(j-1,i,n3_f+1-k)*XI33_f(j-1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,3)*2.d0/3.d0 &
                  -Jacobian_f(j+1,i,n3_f+1-k)*lambda_f(j+1,i,n3_f+1-k)*XI33_f(j+1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,3)&
                  *2.d0/3.d0+Jacobian_f(j+2,i,n3_f+1-k)*lambda_f(j+2,i,n3_f+1-k)*XI33_f(j+2,i,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j+2,i,n3_f+1-k1,3)/12.d0)/l1/h1_f/h3_f+(-Jacobian_f(j,i-2,n3_f+1-k)*mu_f(j,i-2,n3_f+1-k)&
                  *XI23_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,1)/12.d0+Jacobian_f(j,i-1,n3_f+1-k)*mu_f(j,i-1,n3_f+1-k)&
                  *XI23_f(j,i-1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-1,n3_f+1-k1,1)*2.d0/3.d0-Jacobian_f(j,i+1,n3_f+1-k)&
                  *mu_f(j,i+1,n3_f+1-k)*XI23_f(j,i+1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,1)*2.d0/3.d0&
                  +Jacobian_f(j,i+2,n3_f+1-k)*mu_f(j,i+2,n3_f+1-k)*XI23_f(j,i+2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i+2,n3_f+1-k1,1)&
                  /12.d0)/l2/h2_f/h3_f+(-Jacobian_f(j,i-2,n3_f+1-k)*mu_f(j,i-2,n3_f+1-k)*XI13_f(j,i-2,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i-2,n3_f+1-k1,2)/12.d0+Jacobian_f(j,i-1,n3_f+1-k)*mu_f(j,i-1,n3_f+1-k)*XI13_f(j,i-1,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i-1,n3_f+1-k1,2)*2.d0/3.d0-Jacobian_f(j,i+1,n3_f+1-k)*mu_f(j,i+1,n3_f+1-k)*XI13_f(j,i+1,n3_f+1-k)&
                  *bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,2)*2.d0/3.d0+Jacobian_f(j,i+2,n3_f+1-k)*mu_f(j,i+2,n3_f+1-k)&
                  *XI13_f(j,i+2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i+2,n3_f+1-k1,2)/12.d0)/l2/h2_f/h3_f
                ! second set equation
                lh_f(j,i,k,2) = lh_f(j,i,k,2)+(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI23_f(j-2,i,k)*bof(k,k1)*u_f_t(j-2,i,k1,1)/12.d0&
                  -Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI23_f(j-1,i,k)*bof(k,k1)*u_f_t(j-1,i,k1,1)*2.d0/3.d0+Jacobian_f(j+1,i,k)&
                  *mu_f(j+1,i,k)*XI23_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,1)*2.d0/3.d0-Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)&
                  *XI23_f(j+2,i,k)*bof(k,k1)*u_f_t(j+2,i,k1,1)/12.d0)/l1/h1_f/h3_f+(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)&
                  *XI13_f(j-2,i,k)*bof(k,k1)*u_f_t(j-2,i,k1,2)/12.d0-Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI13_f(j-1,i,k)*bof(k,k1)&
                  *u_f_t(j-1,i,k1,2)*2.d0/3.d0+Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI13_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,2)&
                  *2.d0/3.d0-Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI13_f(j+2,i,k)*bof(k,k1)*u_f_t(j+2,i,k1,2)/12.d0)/l1/h1_f/h3_f&
                  +(Jacobian_f(j,i-2,k)*lambda_f(j,i-2,k)*XI13_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,1)/12.d0-Jacobian_f(j,i-1,k)&
                  *lambda_f(j,i-1,k)*XI13_f(j,i-1,k)*bof(k,k1)*u_f_t(j,i-1,k1,1)*2.d0/3.d0+Jacobian_f(j,i+1,k)*lambda_f(j,i+1,k)&
                  *XI13_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,1)*2.d0/3.d0-Jacobian_f(j,i+2,k)*lambda_f(j,i+2,k)*XI13_f(j,i+2,k)&
                  *bof(k,k1)*u_f_t(j,i+2,k1,1)/12.d0)/l2/h2_f/h3_f+(Jacobian_f(j,i-2,k)*(2.d0*mu_f(j,i-2,k)+lambda_f(j,i-2,k))&
                  *XI23_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,2)/12.d0-Jacobian_f(j,i-1,k)*(2.d0*mu_f(j,i-1,k)+lambda_f(j,i-1,k))&
                  *XI23_f(j,i-1,k)*bof(k,k1)*u_f_t(j,i-1,k1,2)*2.d0/3.d0+Jacobian_f(j,i+1,k)*(2.d0*mu_f(j,i+1,k)+lambda_f(j,i+1,k))&
                  *XI23_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,2)*2.d0/3.d0-Jacobian_f(j,i+2,k)*(2.d0*mu_f(j,i+2,k)+lambda_f(j,i+2,k))&
                  *XI23_f(j,i+2,k)*bof(k,k1)*u_f_t(j,i+2,k1,2)/12.d0)/l2/h2_f/h3_f+(Jacobian_f(j,i-2,k)*lambda_f(j,i-2,k)&
                  *XI33_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,3)/12.d0-Jacobian_f(j,i-1,k)*lambda_f(j,i-1,k)*XI33_f(j,i-1,k)&
                  *bof(k,k1)*u_f_t(j,i-1,k1,3)*2.d0/3.d0+Jacobian_f(j,i+1,k)*lambda_f(j,i+1,k)*XI33_f(j,i+1,k)*bof(k,k1)&
                  *u_f_t(j,i+1,k1,3)*2.d0/3.d0-Jacobian_f(j,i+2,k)*lambda_f(j,i+2,k)*XI33_f(j,i+2,k)*bof(k,k1)&
                  *u_f_t(j,i+2,k1,3)/12.d0)/l2/h2_f/h3_f

                lh_f(j,i,n3_f+1-k,2) = lh_f(j,i,n3_f+1-k,2)+(-Jacobian_f(j-2,i,n3_f+1-k)*mu_f(j-2,i,n3_f+1-k)&
                  *XI23_f(j-2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,1)/12.d0+Jacobian_f(j-1,i,n3_f+1-k)*mu_f(j-1,i,n3_f+1-k)&
                  *XI23_f(j-1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,1)*2.d0/3.d0-Jacobian_f(j+1,i,n3_f+1-k)&
                  *mu_f(j+1,i,n3_f+1-k)*XI23_f(j+1,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,1)*2.d0/3.d0&
                  +Jacobian_f(j+2,i,n3_f+1-k)*mu_f(j+2,i,n3_f+1-k)*XI23_f(j+2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+2,i,n3_f+1-k1,1)/12.d0)&
                  /l1/h1_f/h3_f+(-Jacobian_f(j-2,i,n3_f+1-k)*mu_f(j-2,i,n3_f+1-k)*XI13_f(j-2,i,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j-2,i,n3_f+1-k1,2)/12.d0+Jacobian_f(j-1,i,n3_f+1-k)*mu_f(j-1,i,n3_f+1-k)*XI13_f(j-1,i,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j-1,i,n3_f+1-k1,2)*2.d0/3.d0-Jacobian_f(j+1,i,n3_f+1-k)*mu_f(j+1,i,n3_f+1-k)*XI13_f(j+1,i,n3_f+1-k)&
                  *bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,2)*2.d0/3.d0+Jacobian_f(j+2,i,n3_f+1-k)*mu_f(j+2,i,n3_f+1-k)&
                  *XI13_f(j+2,i,n3_f+1-k)*bof(k,k1)*u_f_t(j+2,i,n3_f+1-k1,2)/12.d0)/l1/h1_f/h3_f+(-Jacobian_f(j,i-2,n3_f+1-k)&
                  *lambda_f(j,i-2,n3_f+1-k)*XI13_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,1)/12.d0&
                  +Jacobian_f(j,i-1,n3_f+1-k)*lambda_f(j,i-1,n3_f+1-k)*XI13_f(j,i-1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-1,n3_f+1-k1,1)&
                  *2.d0/3.d0-Jacobian_f(j,i+1,n3_f+1-k)*lambda_f(j,i+1,n3_f+1-k)*XI13_f(j,i+1,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i+1,n3_f+1-k1,1)*2.d0/3.d0+Jacobian_f(j,i+2,n3_f+1-k)*lambda_f(j,i+2,n3_f+1-k)*XI13_f(j,i+2,n3_f+1-k)&
                  *bof(k,k1)*u_f_t(j,i+2,n3_f+1-k1,1)/12.d0)/l2/h2_f/h3_f+(-Jacobian_f(j,i-2,n3_f+1-k)*(2.d0*mu_f(j,i-2,n3_f+1-k)&
                  +lambda_f(j,i-2,n3_f+1-k))*XI23_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,2)/12.d0 &
                  +Jacobian_f(j,i-1,n3_f+1-k)*(2.d0*mu_f(j,i-1,n3_f+1-k)+lambda_f(j,i-1,n3_f+1-k))*XI23_f(j,i-1,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i-1,n3_f+1-k1,2)*2.d0/3.d0-Jacobian_f(j,i+1,n3_f+1-k)*(2.d0*mu_f(j,i+1,n3_f+1-k)&
                  +lambda_f(j,i+1,n3_f+1-k))*XI23_f(j,i+1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,2)*2.d0/3.d0&
                  +Jacobian_f(j,i+2,n3_f+1-k)*(2.d0*mu_f(j,i+2,n3_f+1-k)+lambda_f(j,i+2,n3_f+1-k))*XI23_f(j,i+2,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i+2,n3_f+1-k1,2)/12.d0)/l2/h2_f/h3_f+(-Jacobian_f(j,i-2,n3_f+1-k)*lambda_f(j,i-2,n3_f+1-k)&
                  *XI33_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,3)/12.d0+Jacobian_f(j,i-1,n3_f+1-k)&
                  *lambda_f(j,i-1,n3_f+1-k)*XI33_f(j,i-1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-1,n3_f+1-k1,3)*2.d0/3.d0 &
                  -Jacobian_f(j,i+1,n3_f+1-k)*lambda_f(j,i+1,n3_f+1-k)*XI33_f(j,i+1,n3_f+1-k)*bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,3)&
                  *2.d0/3.d0+Jacobian_f(j,i+2,n3_f+1-k)*lambda_f(j,i+2,n3_f+1-k)*XI33_f(j,i+2,n3_f+1-k)*bof(k,k1)&
                  *u_f_t(j,i+2,n3_f+1-k1,3)/12.d0)/l2/h2_f/h3_f
                ! third set equation
                lh_f(j,i,k,3) = lh_f(j,i,k,3) &
                      +(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI33_f(j-2,i,k)*bof(k,k1)*u_f_t(j-2,i,k1,1)/12.d0 &
                      -Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI33_f(j-1,i,k)*bof(k,k1)*u_f_t(j-1,i,k1,1)*2.d0/3.d0 &
                      +Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI33_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI33_f(j+2,i,k)*bof(k,k1)*u_f_t(j+2,i,k1,1)/12.d0)/l1/h1_f/h3_f&
                      +(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI13_f(j-2,i,k)*bof(k,k1)*u_f_t(j-2,i,k1,3)/12.d0 &
                      -Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI13_f(j-1,i,k)*bof(k,k1)*u_f_t(j-1,i,k1,3)*2.d0/3.d0 &
                      +Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI13_f(j+1,i,k)*bof(k,k1)*u_f_t(j+1,i,k1,3)*2.d0/3.d0 &
                      -Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI13_f(j+2,i,k)*bof(k,k1)*u_f_t(j+2,i,k1,3)/12.d0)/l1/h1_f/h3_f&
                      +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI33_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI33_f(j,i-1,k)*bof(k,k1)*u_f_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI33_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI33_f(j,i+2,k)*bof(k,k1)*u_f_t(j,i+2,k1,2)/12.d0)/l2/h2_f/h3_f&
                      +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI23_f(j,i-2,k)*bof(k,k1)*u_f_t(j,i-2,k1,3)/12.d0 &
                      -Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI23_f(j,i-1,k)*bof(k,k1)*u_f_t(j,i-1,k1,3)*2.d0/3.d0 &
                      +Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI23_f(j,i+1,k)*bof(k,k1)*u_f_t(j,i+1,k1,3)*2.d0/3.d0 &
                      -Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI23_f(j,i+2,k)*bof(k,k1)*u_f_t(j,i+2,k1,3)/12.d0)/l2/h2_f/h3_f

                lh_f(j,i,n3_f+1-k,3) = lh_f(j,i,n3_f+1-k,3) +&
                      (-Jacobian_f(j-2,i,n3_f+1-k)*mu_f(j-2,i,n3_f+1-k)*XI33_f(j-2,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,1)/12.d0 &
                      +Jacobian_f(j-1,i,n3_f+1-k)*mu_f(j-1,i,n3_f+1-k)*XI33_f(j-1,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_f(j+1,i,n3_f+1-k)*mu_f(j+1,i,n3_f+1-k)*XI33_f(j+1,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_f(j+2,i,n3_f+1-k)*mu_f(j+2,i,n3_f+1-k)*XI33_f(j+2,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j+2,i,n3_f+1-k1,1)/12.d0)/l1/h1_f/h3_f&
                      +(-Jacobian_f(j-2,i,n3_f+1-k)*mu_f(j-2,i,n3_f+1-k)*XI13_f(j-2,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j-2,i,n3_f+1-k1,3)/12.d0 &
                      +Jacobian_f(j-1,i,n3_f+1-k)*mu_f(j-1,i,n3_f+1-k)*XI13_f(j-1,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j-1,i,n3_f+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_f(j+1,i,n3_f+1-k)*mu_f(j+1,i,n3_f+1-k)*XI13_f(j+1,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j+1,i,n3_f+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_f(j+2,i,n3_f+1-k)*mu_f(j+2,i,n3_f+1-k)*XI13_f(j+2,i,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j+2,i,n3_f+1-k1,3)/12.d0)/l1/h1_f/h3_f&
                      +(-Jacobian_f(j,i-2,n3_f+1-k)*mu_f(j,i-2,n3_f+1-k)&
                      *XI33_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,2)/12.d0 &
                      +Jacobian_f(j,i-1,n3_f+1-k)*mu_f(j,i-1,n3_f+1-k)*XI33_f(j,i-1,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i-1,n3_f+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_f(j,i+1,n3_f+1-k)*mu_f(j,i+1,n3_f+1-k)*XI33_f(j,i+1,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_f(j,i+2,n3_f+1-k)*mu_f(j,i+2,n3_f+1-k)*XI33_f(j,i+2,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i+2,n3_f+1-k1,2)/12.d0)/l2/h2_f/h3_f&
                      +(-Jacobian_f(j,i-2,n3_f+1-k)*mu_f(j,i-2,n3_f+1-k)&
                      *XI23_f(j,i-2,n3_f+1-k)*bof(k,k1)*u_f_t(j,i-2,n3_f+1-k1,3)/12.d0 &
                      +Jacobian_f(j,i-1,n3_f+1-k)*mu_f(j,i-1,n3_f+1-k)*XI23_f(j,i-1,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i-1,n3_f+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_f(j,i+1,n3_f+1-k)*mu_f(j,i+1,n3_f+1-k)*XI23_f(j,i+1,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i+1,n3_f+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_f(j,i+2,n3_f+1-k)*mu_f(j,i+2,n3_f+1-k)*XI23_f(j,i+2,n3_f+1-k)&
                      *bof(k,k1)*u_f_t(j,i+2,n3_f+1-k1,3)/12.d0)/l2/h2_f/h3_f
             end do
          end do
       end do
    end do
    do k = 5,n3_f-4
       do i = 1,n2_f
          do j = 1,n1_f
             ! mixed derivative 13 & 23
             ! first set equation
             lh_f(j,i,k,1) = lh_f(j,i,k,1)+(Jacobian_f(j-2,i,k)*(2.d0*mu_f(j-2,i,k)+lambda_f(j-2,i,k))*XI13_f(j-2,i,k)&
                   *(u_f_t(j-2,i,k-2,1)/12.d0-u_f_t(j-2,i,k-1,1)*2.d0/3.d0&
                   +u_f_t(j-2,i,k+1,1)*2.d0/3.d0-u_f_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*(2.d0*mu_f(j-1,i,k)+lambda_f(j-1,i,k))*XI13_f(j-1,i,k)*(u_f_t(j-1,i,k-2,1)/12.d0&
                   -u_f_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,1)*2.d0/3.d0-u_f_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*(2.d0*mu_f(j+1,i,k)+lambda_f(j+1,i,k))&
                   *XI13_f(j+1,i,k)*(u_f_t(j+1,i,k-2,1)/12.d0-u_f_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,1)*2.d0/3.d0-u_f_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*(2.d0*mu_f(j+2,i,k)+lambda_f(j+2,i,k))*XI13_f(j+2,i,k)&
                   *(u_f_t(j+2,i,k-2,1)/12.d0-u_f_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,1)*2.d0/3.d0-u_f_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j-2,i,k)*lambda_f(j-2,i,k)*XI23_f(j-2,i,k)*(u_f_t(j-2,i,k-2,2)/12.d0&
                   -u_f_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,2)*2.d0/3.d0-u_f_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*lambda_f(j-1,i,k)*XI23_f(j-1,i,k)*(u_f_t(j-1,i,k-2,2)/12.d0&
                   -u_f_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,2)*2.d0/3.d0-u_f_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*lambda_f(j+1,i,k)*XI23_f(j+1,i,k)*(u_f_t(j+1,i,k-2,2)/12.d0&
                   -u_f_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,2)*2.d0/3.d0-u_f_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*lambda_f(j+2,i,k)*XI23_f(j+2,i,k)*(u_f_t(j+2,i,k-2,2)/12.d0&
                   -u_f_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,2)*2.d0/3.d0-u_f_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j-2,i,k)*lambda_f(j-2,i,k)*XI33_f(j-2,i,k)*(u_f_t(j-2,i,k-2,3)/12.d0&
                   -u_f_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,3)*2.d0/3.d0-u_f_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*lambda_f(j-1,i,k)*XI33_f(j-1,i,k)*(u_f_t(j-1,i,k-2,3)/12.d0-u_f_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,3)*2.d0/3.d0-u_f_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*lambda_f(j+1,i,k)*XI33_f(j+1,i,k)*(u_f_t(j+1,i,k-2,3)/12.d0-u_f_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,3)*2.d0/3.d0-u_f_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*lambda_f(j+2,i,k)*XI33_f(j+2,i,k)*(u_f_t(j+2,i,k-2,3)/12.d0-u_f_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,3)*2.d0/3.d0-u_f_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI23_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,1)/12.d0-u_f_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,1)*2.d0/3.d0-u_f_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI23_f(j,i-1,k)*(u_f_t(j,i-1,k-2,1)/12.d0-u_f_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,1)*2.d0/3.d0-u_f_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI23_f(j,i+1,k)*(u_f_t(j,i+1,k-2,1)/12.d0-u_f_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,1)*2.d0/3.d0-u_f_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI23_f(j,i+2,k)*(u_f_t(j,i+2,k-2,1)/12.d0-u_f_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,1)*2.d0/3.d0-u_f_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI13_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,2)/12.d0-u_f_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,2)*2.d0/3.d0-u_f_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI13_f(j,i-1,k)*(u_f_t(j,i-1,k-2,2)/12.d0-u_f_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,2)*2.d0/3.d0-u_f_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI13_f(j,i+1,k)*(u_f_t(j,i+1,k-2,2)/12.d0-u_f_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,2)*2.d0/3.d0-u_f_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI13_f(j,i+2,k)*(u_f_t(j,i+2,k-2,2)/12.d0-u_f_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,2)*2.d0/3.d0-u_f_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_f/h3_f
             ! second set equation
             lh_f(j,i,k,2) = lh_f(j,i,k,2)+(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI23_f(j-2,i,k)*(u_f_t(j-2,i,k-2,1)/12.d0&
                   -u_f_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,1)*2.d0/3.d0-u_f_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI23_f(j-1,i,k)*(u_f_t(j-1,i,k-2,1)/12.d0-u_f_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,1)*2.d0/3.d0-u_f_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI23_f(j+1,i,k)*(u_f_t(j+1,i,k-2,1)/12.d0-u_f_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,1)*2.d0/3.d0-u_f_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI23_f(j+2,i,k)*(u_f_t(j+2,i,k-2,1)/12.d0-u_f_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,1)*2.d0/3.d0-u_f_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI13_f(j-2,i,k)*(u_f_t(j-2,i,k-2,2)/12.d0&
                   -u_f_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,2)*2.d0/3.d0-u_f_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI13_f(j-1,i,k)*(u_f_t(j-1,i,k-2,2)/12.d0-u_f_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,2)*2.d0/3.d0-u_f_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI13_f(j+1,i,k)*(u_f_t(j+1,i,k-2,2)/12.d0-u_f_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,2)*2.d0/3.d0-u_f_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI13_f(j+2,i,k)*(u_f_t(j+2,i,k-2,2)/12.d0-u_f_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,2)*2.d0/3.d0-u_f_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*lambda_f(j,i-2,k)*XI13_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,1)/12.d0-u_f_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,1)*2.d0/3.d0-u_f_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*lambda_f(j,i-1,k)*XI13_f(j,i-1,k)*(u_f_t(j,i-1,k-2,1)/12.d0&
                   -u_f_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,1)*2.d0/3.d0-u_f_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*lambda_f(j,i+1,k)*XI13_f(j,i+1,k)*(u_f_t(j,i+1,k-2,1)/12.d0&
                   -u_f_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,1)*2.d0/3.d0-u_f_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*lambda_f(j,i+2,k)*XI13_f(j,i+2,k)*(u_f_t(j,i+2,k-2,1)/12.d0&
                   -u_f_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,1)*2.d0/3.d0-u_f_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*(2.d0*mu_f(j,i-2,k)+lambda_f(j,i-2,k))&
                   *XI23_f(j,i-2,k)*(u_f_t(j,i-2,k-2,2)/12.d0-u_f_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,2)*2.d0/3.d0-u_f_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*(2.d0*mu_f(j,i-1,k)+lambda_f(j,i-1,k))*XI23_f(j,i-1,k)&
                   *(u_f_t(j,i-1,k-2,2)/12.d0-u_f_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,2)*2.d0/3.d0-u_f_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*(2.d0*mu_f(j,i+1,k)+lambda_f(j,i+1,k))*XI23_f(j,i+1,k)&
                   *(u_f_t(j,i+1,k-2,2)/12.d0-u_f_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,2)*2.d0/3.d0-u_f_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*(2.d0*mu_f(j,i+2,k)+lambda_f(j,i+2,k))*XI23_f(j,i+2,k)&
                   *(u_f_t(j,i+2,k-2,2)/12.d0-u_f_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,2)*2.d0/3.d0-u_f_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*lambda_f(j,i-2,k)*XI33_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,3)/12.d0-u_f_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,3)*2.d0/3.d0-u_f_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*lambda_f(j,i-1,k)*XI33_f(j,i-1,k)*(u_f_t(j,i-1,k-2,3)/12.d0&
                   -u_f_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,3)*2.d0/3.d0-u_f_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*lambda_f(j,i+1,k)*XI33_f(j,i+1,k)*(u_f_t(j,i+1,k-2,3)/12.d0&
                   -u_f_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,3)*2.d0/3.d0-u_f_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*lambda_f(j,i+2,k)*XI33_f(j,i+2,k)*(u_f_t(j,i+2,k-2,3)/12.d0&
                   -u_f_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,3)*2.d0/3.d0-u_f_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_f/h3_f
             ! third set equation
             lh_f(j,i,k,3) = lh_f(j,i,k,3)+(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI33_f(j-2,i,k)*(u_f_t(j-2,i,k-2,1)/12.d0&
                   -u_f_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,1)*2.d0/3.d0-u_f_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI33_f(j-1,i,k)*(u_f_t(j-1,i,k-2,1)/12.d0-u_f_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,1)*2.d0/3.d0-u_f_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI33_f(j+1,i,k)*(u_f_t(j+1,i,k-2,1)/12.d0-u_f_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,1)*2.d0/3.d0-u_f_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI33_f(j+2,i,k)*(u_f_t(j+2,i,k-2,1)/12.d0-u_f_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,1)*2.d0/3.d0-u_f_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j-2,i,k)*mu_f(j-2,i,k)*XI13_f(j-2,i,k)*(u_f_t(j-2,i,k-2,3)/12.d0&
                   -u_f_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j-2,i,k+1,3)*2.d0/3.d0-u_f_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_f(j-1,i,k)*mu_f(j-1,i,k)*XI13_f(j-1,i,k)*(u_f_t(j-1,i,k-2,3)/12.d0-u_f_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j-1,i,k+1,3)*2.d0/3.d0-u_f_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j+1,i,k)*mu_f(j+1,i,k)*XI13_f(j+1,i,k)*(u_f_t(j+1,i,k-2,3)/12.d0-u_f_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j+1,i,k+1,3)*2.d0/3.d0-u_f_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j+2,i,k)*mu_f(j+2,i,k)*XI13_f(j+2,i,k)*(u_f_t(j+2,i,k-2,3)/12.d0-u_f_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j+2,i,k+1,3)*2.d0/3.d0-u_f_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI33_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,2)/12.d0-u_f_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,2)*2.d0/3.d0-u_f_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI33_f(j,i-1,k)*(u_f_t(j,i-1,k-2,2)/12.d0&
                   -u_f_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,2)*2.d0/3.d0-u_f_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI33_f(j,i+1,k)*(u_f_t(j,i+1,k-2,2)/12.d0&
                   -u_f_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,2)*2.d0/3.d0-u_f_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI33_f(j,i+2,k)*(u_f_t(j,i+2,k-2,2)/12.d0&
                   -u_f_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,2)*2.d0/3.d0-u_f_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_f/h3_f&
                   +(Jacobian_f(j,i-2,k)*mu_f(j,i-2,k)*XI23_f(j,i-2,k)&
                   *(u_f_t(j,i-2,k-2,3)/12.d0-u_f_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i-2,k+1,3)*2.d0/3.d0-u_f_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_f(j,i-1,k)*mu_f(j,i-1,k)*XI23_f(j,i-1,k)*(u_f_t(j,i-1,k-2,3)/12.d0&
                   -u_f_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i-1,k+1,3)*2.d0/3.d0-u_f_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_f(j,i+1,k)*mu_f(j,i+1,k)*XI23_f(j,i+1,k)*(u_f_t(j,i+1,k-2,3)/12.d0&
                   -u_f_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i+1,k+1,3)*2.d0/3.d0-u_f_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_f(j,i+2,k)*mu_f(j,i+2,k)*XI23_f(j,i+2,k)*(u_f_t(j,i+2,k-2,3)/12.d0&
                   -u_f_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_f_t(j,i+2,k+1,3)*2.d0/3.d0-u_f_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_f/h3_f
          end do
       end do
    end do
    !
    do i = 5,n3_f-4
       do j = 1,n2_f
          do k = 1,n1_f
             ! mixed derivative 31 & 32
             ! first set equation
             lh_f(k,j,i,1) = lh_f(k,j,i,1)+(Jacobian_f(k,j,i-2)*(2.d0*mu_f(k,j,i-2)&
                     +lambda_f(k,j,i-2))*XI13_f(k,j,i-2)&
                     *(u_f_t(k-2,j,i-2,1)/12.d0-u_f_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,1)*2.d0/3.d0-u_f_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*(2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI13_f(k,j,i-1)&
                     *(u_f_t(k-2,j,i-1,1)/12.d0-u_f_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,1)*2.d0/3.d0-u_f_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*(2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI13_f(k,j,i+1)&
                     *(u_f_t(k-2,j,i+1,1)/12.d0-u_f_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,1)*2.d0/3.d0-u_f_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*(2.d0*mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI13_f(k,j,i+2)&
                     *(u_f_t(k-2,j,i+2,1)/12.d0-u_f_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,1)*2.d0/3.d0-u_f_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_f/h1_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI23_f(k,j,i-2)*(u_f_t(k-2,j,i-2,2)/12.d0&
                     -u_f_t(k-1,j,i-2,2)*2.d0/3.d0+u_f_t(k+1,j,i-2,2)*2.d0/3.d0-u_f_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI23_f(k,j,i-1)*(u_f_t(k-2,j,i-1,2)/12.d0&
                     -u_f_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,2)*2.d0/3.d0-u_f_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI23_f(k,j,i+1)*(u_f_t(k-2,j,i+1,2)/12.d0&
                     -u_f_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,2)*2.d0/3.d0-u_f_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI23_f(k,j,i+2)*(u_f_t(k-2,j,i+2,2)/12.d0&
                     -u_f_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,2)*2.d0/3.d0-u_f_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h3_f/h1_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI33_f(k,j,i-2)*(u_f_t(k-2,j,i-2,3)/12.d0&
                     -u_f_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,3)*2.d0/3.d0-u_f_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI33_f(k,j,i-1)*(u_f_t(k-2,j,i-1,3)/12.d0&
                     -u_f_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,3)*2.d0/3.d0-u_f_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI33_f(k,j,i+1)*(u_f_t(k-2,j,i+1,3)/12.d0&
                     -u_f_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,3)*2.d0/3.d0-u_f_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI33_f(k,j,i+2)*(u_f_t(k-2,j,i+2,3)/12.d0&
                     -u_f_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,3)*2.d0/3.d0-u_f_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h1_f/h3_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI23_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,1)/12.d0-u_f_t(k,j-1,i-2,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,1)*2.d0/3.d0-u_f_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI23_f(k,j,i-1)*(u_f_t(k,j-2,i-1,1)/12.d0&
                     -u_f_t(k,j-1,i-1,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,1)*2.d0/3.d0-u_f_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI23_f(k,j,i+1)*(u_f_t(k,j-2,i+1,1)/12.d0&
                     -u_f_t(k,j-1,i+1,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,1)*2.d0/3.d0-u_f_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI23_f(k,j,i+2)*(u_f_t(k,j-2,i+2,1)/12.d0&
                     -u_f_t(k,j-1,i+2,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,1)*2.d0/3.d0-u_f_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_f/h2_f&
                     +(Jacobian_f(k,j,i-2)*lambda_f(k,j,i-2)*XI13_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,2)/12.d0-u_f_t(k,j-1,i-2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,2)*2.d0/3.d0-u_f_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*lambda_f(k,j,i-1)*XI13_f(k,j,i-1)*(u_f_t(k,j-2,i-1,2)/12.d0&
                     -u_f_t(k,j-1,i-1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,2)*2.d0/3.d0-u_f_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*lambda_f(k,j,i+1)*XI13_f(k,j,i+1)*(u_f_t(k,j-2,i+1,2)/12.d0&
                     -u_f_t(k,j-1,i+1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,2)*2.d0/3.d0-u_f_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*lambda_f(k,j,i+2)*XI13_f(k,j,i+2)*(u_f_t(k,j-2,i+2,2)/12.d0&
                     -u_f_t(k,j-1,i+2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,2)*2.d0/3.d0-u_f_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_f/h2_f
             ! second set equation
             lh_f(k,j,i,2) = lh_f(k,j,i,2)+(Jacobian_f(k,j,i-2)*lambda_f(k,j,i-2)*XI23_f(k,j,i-2)&
                     *(u_f_t(k-2,j,i-2,1)/12.d0-u_f_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,1)*2.d0/3.d0-u_f_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*lambda_f(k,j,i-1)*XI23_f(k,j,i-1)*(u_f_t(k-2,j,i-1,1)/12.d0&
                     -u_f_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,1)*2.d0/3.d0-u_f_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*lambda_f(k,j,i+1)*XI23_f(k,j,i+1)*(u_f_t(k-2,j,i+1,1)/12.d0&
                     -u_f_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,1)*2.d0/3.d0-u_f_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*lambda_f(k,j,i+2)*XI23_f(k,j,i+2)*(u_f_t(k-2,j,i+2,1)/12.d0&
                     -u_f_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,1)*2.d0/3.d0-u_f_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h1_f/h3_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI13_f(k,j,i-2)&
                     *(u_f_t(k-2,j,i-2,2)/12.d0-u_f_t(k-1,j,i-2,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,2)*2.d0/3.d0-u_f_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI13_f(k,j,i-1)*(u_f_t(k-2,j,i-1,2)/12.d0&
                     -u_f_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,2)*2.d0/3.d0-u_f_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI13_f(k,j,i+1)*(u_f_t(k-2,j,i+1,2)/12.d0&
                     -u_f_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,2)*2.d0/3.d0-u_f_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI13_f(k,j,i+2)&
                     *(u_f_t(k-2,j,i+2,2)/12.d0-u_f_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,2)*2.d0/3.d0-u_f_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h1_f/h3_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI13_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,1)/12.d0-u_f_t(k,j-1,i-2,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,1)*2.d0/3.d0-u_f_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI13_f(k,j,i-1)*(u_f_t(k,j-2,i-1,1)/12.d0&
                     -u_f_t(k,j-1,i-1,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,1)*2.d0/3.d0-u_f_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI13_f(k,j,i+1)*(u_f_t(k,j-2,i+1,1)/12.d0&
                     -u_f_t(k,j-1,i+1,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,1)*2.d0/3.d0-u_f_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI13_f(k,j,i+2)*(u_f_t(k,j-2,i+2,1)/12.d0&
                     -u_f_t(k,j-1,i+2,1)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,1)*2.d0/3.d0-u_f_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_f/h2_f&
                     +(Jacobian_f(k,j,i-2)*(2.d0*mu_f(k,j,i-2)+lambda_f(k,j,i-2))&
                     *XI23_f(k,j,i-2)*(u_f_t(k,j-2,i-2,2)/12.d0-u_f_t(k,j-1,i-2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,2)*2.d0/3.d0-u_f_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*(2.d0*mu_f(k,j,i-1)+lambda_f(k,j,i-1))*XI23_f(k,j,i-1)&
                     *(u_f_t(k,j-2,i-1,2)/12.d0-u_f_t(k,j-1,i-1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,2)*2.d0/3.d0-u_f_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*(2.d0*mu_f(k,j,i+1)+lambda_f(k,j,i+1))*XI23_f(k,j,i+1)&
                     *(u_f_t(k,j-2,i+1,2)/12.d0-u_f_t(k,j-1,i+1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,2)*2.d0/3.d0-u_f_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*(2.d0*mu_f(k,j,i+2)+lambda_f(k,j,i+2))*XI23_f(k,j,i+2)&
                     *(u_f_t(k,j-2,i+2,2)/12.d0-u_f_t(k,j-1,i+2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,2)*2.d0/3.d0-u_f_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_f/h2_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI33_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,3)/12.d0-u_f_t(k,j-1,i-2,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,3)*2.d0/3.d0-u_f_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI33_f(k,j,i-1)*(u_f_t(k,j-2,i-1,3)/12.d0&
                     -u_f_t(k,j-1,i-1,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,3)*2.d0/3.d0-u_f_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI33_f(k,j,i+1)*(u_f_t(k,j-2,i+1,3)/12.d0&
                     -u_f_t(k,j-1,i+1,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,3)*2.d0/3.d0-u_f_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI33_f(k,j,i+2)*(u_f_t(k,j-2,i+2,3)/12.d0&
                     -u_f_t(k,j-1,i+2,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,3)*2.d0/3.d0-u_f_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h3_f/h2_f
             ! third set equation
             lh_f(k,j,i,3) = lh_f(k,j,i,3)+(Jacobian_f(k,j,i-2)*lambda_f(k,j,i-2)*XI33_f(k,j,i-2)&
                     *(u_f_t(k-2,j,i-2,1)/12.d0-u_f_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,1)*2.d0/3.d0-u_f_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*lambda_f(k,j,i-1)*XI33_f(k,j,i-1)*(u_f_t(k-2,j,i-1,1)/12.d0&
                     -u_f_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,1)*2.d0/3.d0-u_f_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*lambda_f(k,j,i+1)*XI33_f(k,j,i+1)*(u_f_t(k-2,j,i+1,1)/12.d0&
                     -u_f_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,1)*2.d0/3.d0-u_f_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*lambda_f(k,j,i+2)*XI33_f(k,j,i+2)*(u_f_t(k-2,j,i+2,1)/12.d0&
                     -u_f_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,1)*2.d0/3.d0-u_f_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_f/h1_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI13_f(k,j,i-2)&
                     *(u_f_t(k-2,j,i-2,3)/12.d0-u_f_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-2,3)*2.d0/3.d0-u_f_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI13_f(k,j,i-1)*(u_f_t(k-2,j,i-1,3)/12.d0&
                     -u_f_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i-1,3)*2.d0/3.d0-u_f_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI13_f(k,j,i+1)*(u_f_t(k-2,j,i+1,3)/12.d0&
                     -u_f_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+1,3)*2.d0/3.d0-u_f_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI13_f(k,j,i+2)*(u_f_t(k-2,j,i+2,3)/12.d0&
                     -u_f_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_f_t(k+1,j,i+2,3)*2.d0/3.d0-u_f_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h3_f/h1_f&
                     +(Jacobian_f(k,j,i-2)*lambda_f(k,j,i-2)*XI33_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,2)/12.d0-u_f_t(k,j-1,i-2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,2)*2.d0/3.d0-u_f_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*lambda_f(k,j,i-1)*XI33_f(k,j,i-1)*(u_f_t(k,j-2,i-1,2)/12.d0&
                     -u_f_t(k,j-1,i-1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,2)*2.d0/3.d0-u_f_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*lambda_f(k,j,i+1)*XI33_f(k,j,i+1)&
                     *(u_f_t(k,j-2,i+1,2)/12.d0-u_f_t(k,j-1,i+1,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,2)*2.d0/3.d0-u_f_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*lambda_f(k,j,i+2)*XI33_f(k,j,i+2)*(u_f_t(k,j-2,i+2,2)/12.d0&
                     -u_f_t(k,j-1,i+2,2)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,2)*2.d0/3.d0-u_f_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h2_f/h3_f&
                     +(Jacobian_f(k,j,i-2)*mu_f(k,j,i-2)*XI23_f(k,j,i-2)&
                     *(u_f_t(k,j-2,i-2,3)/12.d0-u_f_t(k,j-1,i-2,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-2,3)*2.d0/3.d0-u_f_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_f(k,j,i-1)*mu_f(k,j,i-1)*XI23_f(k,j,i-1)*(u_f_t(k,j-2,i-1,3)/12.d0&
                     -u_f_t(k,j-1,i-1,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i-1,3)*2.d0/3.d0-u_f_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_f(k,j,i+1)*mu_f(k,j,i+1)*XI23_f(k,j,i+1)*(u_f_t(k,j-2,i+1,3)/12.d0&
                     -u_f_t(k,j-1,i+1,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+1,3)*2.d0/3.d0-u_f_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_f(k,j,i+2)*mu_f(k,j,i+2)*XI23_f(k,j,i+2)*(u_f_t(k,j-2,i+2,3)/12.d0&
                     -u_f_t(k,j-1,i+2,3)*2.d0/3.d0 &
                     +u_f_t(k,j+1,i+2,3)*2.d0/3.d0-u_f_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h2_f/h3_f
          end do
       end do
    end do
    do i = 1,4
       do j = 1,n2_f
          do k = 1,n1_f
             do k1 = 1,6
                ! mixed derivative 31 & 32
                ! first set equation
                lh_f(k,j,i,1) = lh_f(k,j,i,1)+(bof(i,k1)*Jacobian_f(k,j,k1)*(2.d0*mu_f(k,j,k1)+lambda_f(k,j,k1))&
                      *XI13_f(k,j,k1)*(u_f_t(k-2,j,k1,1)/12.d0-u_f_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,1)*2.d0/3.d0-u_f_t(k+2,j,k1,1)/12.d0))/l1/h3_f/h1_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI23_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,2)/12.d0-u_f_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,2)*2.d0/3.d0-u_f_t(k+2,j,k1,2)/12.d0))/l1/h3_f/h1_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI33_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,3)/12.d0-u_f_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,3)*2.d0/3.d0-u_f_t(k+2,j,k1,3)/12.d0))/l1/h1_f/h3_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI23_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,1)/12.d0-u_f_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,1)*2.d0/3.d0-u_f_t(k,j+2,k1,1)/12.d0))/l2/h3_f/h2_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*lambda_f(k,j,k1)*XI13_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,2)/12.d0-u_f_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,2)*2.d0/3.d0-u_f_t(k,j+2,k1,2)/12.d0))/l2/h3_f/h2_f

                lh_f(k,j,n3_f+1-i,1) = lh_f(k,j,n3_f+1-i,1)+(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)*(2.d0&
                   *mu_f(k,j,n3_f+1-k1)+lambda_f(k,j,n3_f+1-k1))*XI13_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,1)/12.d0 &
                   -u_f_t(k-1,j,n3_f+1-k1,1)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,1)*2.d0/3.d0&
                   -u_f_t(k+2,j,n3_f+1-k1,1)/12.d0))/l1/h3_f/h1_f&
                   +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)*mu_f(k,j,n3_f+1-k1)&
                   *XI23_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,2)/12.d0 &
                   -u_f_t(k-1,j,n3_f+1-k1,2)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,2)*2.d0/3.d0&
                   -u_f_t(k+2,j,n3_f+1-k1,2)/12.d0))/l1/h3_f/h1_f&
                   +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)*mu_f(k,j,n3_f+1-k1)&
                   *XI33_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,3)/12.d0 &
                   -u_f_t(k-1,j,n3_f+1-k1,3)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,3)*2.d0/3.d0&
                   -u_f_t(k+2,j,n3_f+1-k1,3)/12.d0))/l1/h1_f/h3_f&
                   +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)*mu_f(k,j,n3_f+1-k1)&
                   *XI23_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,1)/12.d0 &
                   -u_f_t(k,j-1,n3_f+1-k1,1)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,1)*2.d0/3.d0&
                   -u_f_t(k,j+2,n3_f+1-k1,1)/12.d0))/l2/h3_f/h2_f&
                   +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                   *lambda_f(k,j,n3_f+1-k1)*XI13_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,2)/12.d0 &
                   -u_f_t(k,j-1,n3_f+1-k1,2)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,2)*2.d0/3.d0&
                   -u_f_t(k,j+2,n3_f+1-k1,2)/12.d0))/l2/h3_f/h2_f
                ! second set equation
                lh_f(k,j,i,2) = lh_f(k,j,i,2)+(bof(i,k1)*Jacobian_f(k,j,k1)*lambda_f(k,j,k1)*XI23_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,1)/12.d0-u_f_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,1)*2.d0/3.d0-u_f_t(k+2,j,k1,1)/12.d0))/l1/h1_f/h3_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI13_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,2)/12.d0-u_f_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,2)*2.d0/3.d0-u_f_t(k+2,j,k1,2)/12.d0))/l1/h1_f/h3_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI13_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,1)/12.d0-u_f_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,1)*2.d0/3.d0-u_f_t(k,j+2,k1,1)/12.d0))/l2/h3_f/h2_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*(2.d0*mu_f(k,j,k1)&
                      +lambda_f(k,j,k1))*XI23_f(k,j,k1)*(u_f_t(k,j-2,k1,2)/12.d0-u_f_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,2)*2.d0/3.d0-u_f_t(k,j+2,k1,2)/12.d0))/l2/h3_f/h2_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI33_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,3)/12.d0-u_f_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,3)*2.d0/3.d0-u_f_t(k,j+2,k1,3)/12.d0))/l2/h3_f/h2_f

                lh_f(k,j,n3_f+1-i,2) = lh_f(k,j,n3_f+1-i,2)+(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *lambda_f(k,j,n3_f+1-k1)*XI23_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,1)/12.d0 &
                      -u_f_t(k-1,j,n3_f+1-k1,1)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,1)*2.d0/3.d0&
                      -u_f_t(k+2,j,n3_f+1-k1,1)/12.d0))/l1/h1_f/h3_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *mu_f(k,j,n3_f+1-k1)*XI13_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,2)/12.d0 &
                      -u_f_t(k-1,j,n3_f+1-k1,2)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,2)*2.d0/3.d0&
                      -u_f_t(k+2,j,n3_f+1-k1,2)/12.d0))/l1/h1_f/h3_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *mu_f(k,j,n3_f+1-k1)*XI13_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,1)/12.d0 &
                      -u_f_t(k,j-1,n3_f+1-k1,1)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,1)*2.d0/3.d0&
                      -u_f_t(k,j+2,n3_f+1-k1,1)/12.d0))/l2/h3_f/h2_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *(2.d0*mu_f(k,j,n3_f+1-k1)+lambda_f(k,j,n3_f+1-k1))*XI23_f(k,j,n3_f+1-k1)&
                      *(u_f_t(k,j-2,n3_f+1-k1,2)/12.d0-u_f_t(k,j-1,n3_f+1-k1,2)*2.d0/3.d0&
                      +u_f_t(k,j+1,n3_f+1-k1,2)*2.d0/3.d0-u_f_t(k,j+2,n3_f+1-k1,2)/12.d0))/l2/h3_f/h2_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *mu_f(k,j,n3_f+1-k1)*XI33_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,3)/12.d0 &
                      -u_f_t(k,j-1,n3_f+1-k1,3)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,3)*2.d0/3.d0&
                      -u_f_t(k,j+2,n3_f+1-k1,3)/12.d0))/l2/h3_f/h2_f
                ! third set equation
                lh_f(k,j,i,3) = lh_f(k,j,i,3)+(bof(i,k1)*Jacobian_f(k,j,k1)*lambda_f(k,j,k1)*XI33_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,1)/12.d0-u_f_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,1)*2.d0/3.d0-u_f_t(k+2,j,k1,1)/12.d0))/l1/h1_f/h3_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI13_f(k,j,k1)&
                      *(u_f_t(k-2,j,k1,3)/12.d0-u_f_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_f_t(k+1,j,k1,3)*2.d0/3.d0-u_f_t(k+2,j,k1,3)/12.d0))/l1/h3_f/h1_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*lambda_f(k,j,k1)*XI33_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,2)/12.d0-u_f_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,2)*2.d0/3.d0-u_f_t(k,j+2,k1,2)/12.d0))/l2/h2_f/h3_f&
                      +(bof(i,k1)*Jacobian_f(k,j,k1)*mu_f(k,j,k1)*XI23_f(k,j,k1)&
                      *(u_f_t(k,j-2,k1,3)/12.d0-u_f_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_f_t(k,j+1,k1,3)*2.d0/3.d0-u_f_t(k,j+2,k1,3)/12.d0))/l2/h2_f/h3_f

                lh_f(k,j,n3_f+1-i,3) = lh_f(k,j,n3_f+1-i,3)+(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *lambda_f(k,j,n3_f+1-k1)*XI33_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,1)/12.d0 &
                      -u_f_t(k-1,j,n3_f+1-k1,1)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,1)*2.d0/3.d0&
                      -u_f_t(k+2,j,n3_f+1-k1,1)/12.d0))/l1/h3_f/h1_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)*mu_f(k,j,n3_f+1-k1)&
                      *XI13_f(k,j,n3_f+1-k1)*(u_f_t(k-2,j,n3_f+1-k1,3)/12.d0 &
                      -u_f_t(k-1,j,n3_f+1-k1,3)*2.d0/3.d0+u_f_t(k+1,j,n3_f+1-k1,3)*2.d0/3.d0&
                      -u_f_t(k+2,j,n3_f+1-k1,3)/12.d0))/l1/h3_f/h1_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *lambda_f(k,j,n3_f+1-k1)*XI33_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,2)/12.d0 &
                      -u_f_t(k,j-1,n3_f+1-k1,2)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,2)*2.d0/3.d0&
                      -u_f_t(k,j+2,n3_f+1-k1,2)/12.d0))/l2/h2_f/h3_f&
                      +(-bof(i,k1)*Jacobian_f(k,j,n3_f+1-k1)&
                      *mu_f(k,j,n3_f+1-k1)*XI23_f(k,j,n3_f+1-k1)*(u_f_t(k,j-2,n3_f+1-k1,3)/12.d0 &
                      -u_f_t(k,j-1,n3_f+1-k1,3)*2.d0/3.d0+u_f_t(k,j+1,n3_f+1-k1,3)*2.d0/3.d0&
                      -u_f_t(k,j+2,n3_f+1-k1,3)/12.d0))/l2/h2_f/h3_f
             end do
          end do
       end do
    end do
    !
    ! coarse mesh
    lh_c = 0.d0
    do k = 1,n3_c
       do j = 1,n2_c
          do i = 1,n1_c
             ! second derivative 11 & 22 & 12 & 21
             ! first set
             lh_c(i,j,k,1) = lh_c(i,j,k,1)+((-Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/8.d0 &
                      + Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0)*u_c_t(i-2,j,k,1) &
                      +(Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/6.d0 &
                      + Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/2.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/6.d0)*u_c_t(i-1,j,k,1) &
                      +(-Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/24.d0 &
                      - Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))*5.d0/6.d0 &
                      -Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))*3.d0/4.d0 &
                      -Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))*5.d0/6.d0 &
                      -Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/24.d0)*u_c_t(i-0,j,k,1) &
                      +(Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/6.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/2.d0 &
                      + Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/6.d0)*u_c_t(i+1,j,k,1) &
                      +(-Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/6.d0 &
                      - Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/8.d0)*u_c_t(i+2,j,k,1))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/8.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i,j-2,k,1) &
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/6.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0)*u_c_t(i,j-1,k,1) &
                      +(-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/24.d0 - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*5.d0/6.d0 -Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/24.d0)*u_c_t(i,j-0,k,1) &
                      +(Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/2.d0 &
                      + Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/6.d0)*u_c_t(i,j+1,k,1) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/8.d0)*u_c_t(i,j+2,k,1))/h2_c**2/l2**2&
                      +(Jacobian_c(i-2,j,k)*lambda_c(i-2,j,k)*(u_c_t(i-2,j-2,k,2)/12.d0&
                      -u_c_t(i-2,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i-2,j+1,k,2)*2.d0/3.d0 - u_c_t(i-2,j+2,k,2)/12.d0)/12.d0 &
                      - Jacobian_c(i-1,j,k)*lambda_c(i-1,j,k)*(u_c_t(i-1,j-2,k,2)/12.d0-u_c_t(i-1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i-1,j+1,k,2)*2.d0/3.d0 - u_c_t(i-1,j+2,k,2)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i+1,j,k)*lambda_c(i+1,j,k)*(u_c_t(i+1,j-2,k,2)/12.d0-u_c_t(i+1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,2)*2.d0/3.d0 - u_c_t(i+1,j+2,k,2)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i+2,j,k)*lambda_c(i+2,j,k)*(u_c_t(i+2,j-2,k,2)/12.d0-u_c_t(i+2,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+2,j+1,k,2)*2.d0/3.d0 - u_c_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_c/h2_c&
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)*(u_c_t(i-2,j-2,k,2)/12.d0&
                      -u_c_t(i-1,j-2,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j-2,k,2)*2.d0/3.d0 - u_c_t(i+2,j-2,k,2)/12.d0)/12.d0 &
                      - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*(u_c_t(i-2,j-1,k,2)/12.d0-u_c_t(i-1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j-1,k,2)*2.d0/3.d0 - u_c_t(i+2,j-1,k,2)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*(u_c_t(i-2,j+1,k,2)/12.d0-u_c_t(i-1,j+1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,2)*2.d0/3.d0 - u_c_t(i+2,j+1,k,2)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)*(u_c_t(i-2,j+2,k,2)/12.d0-u_c_t(i-1,j+2,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+2,k,2)*2.d0/3.d0 - u_c_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
             ! second set
             lh_c(i,j,k,2) = lh_c(i,j,k,2)+((-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/8.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i-2,j,k,2) &
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/6.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0)*u_c_t(i-1,j,k,2) &
                      +(-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/24.d0 - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*5.d0/6.d0 -Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/24.d0)*u_c_t(i-0,j,k,2) &
                      +(Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/2.d0 &
                      + Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/6.d0)*u_c_t(i+1,j,k,2) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/8.d0)*u_c_t(i+2,j,k,2))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/8.d0 &
                      + Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0)*u_c_t(i,j-2,k,2) &
                      +(Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/6.d0 &
                      + Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/2.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/6.d0)*u_c_t(i,j-1,k,2) &
                      +(-Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/24.d0 &
                      - Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))*5.d0/6.d0 &
                      -Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/24.d0)*u_c_t(i,j-0,k,2) &
                      +(Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/6.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/2.d0 &
                      + Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/6.d0)*u_c_t(i,j+1,k,2) &
                      +(-Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/6.d0 &
                      - Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/8.d0)*u_c_t(i,j+2,k,2))/h2_c**2/l2**2&
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)*(u_c_t(i-2,j-2,k,1)/12.d0&
                      -u_c_t(i-2,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i-2,j+1,k,1)*2.d0/3.d0 - u_c_t(i-2,j+2,k,1)/12.d0)/12.d0 &
                      - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*(u_c_t(i-1,j-2,k,1)/12.d0-u_c_t(i-1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i-1,j+1,k,1)*2.d0/3.d0 - u_c_t(i-1,j+2,k,1)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*(u_c_t(i+1,j-2,k,1)/12.d0-u_c_t(i+1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,1)*2.d0/3.d0 - u_c_t(i+1,j+2,k,1)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)*(u_c_t(i+2,j-2,k,1)/12.d0-u_c_t(i+2,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+2,j+1,k,1)*2.d0/3.d0 - u_c_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_c/h2_c&
                      +(Jacobian_c(i,j-2,k)*lambda_c(i,j-2,k)*(u_c_t(i-2,j-2,k,1)/12.d0&
                      -u_c_t(i-1,j-2,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j-2,k,1)*2.d0/3.d0 - u_c_t(i+2,j-2,k,1)/12.d0)/12.d0 &
                      - Jacobian_c(i,j-1,k)*lambda_c(i,j-1,k)*(u_c_t(i-2,j-1,k,1)/12.d0-u_c_t(i-1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j-1,k,1)*2.d0/3.d0 - u_c_t(i+2,j-1,k,1)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i,j+1,k)*lambda_c(i,j+1,k)*(u_c_t(i-2,j+1,k,1)/12.d0-u_c_t(i-1,j+1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,1)*2.d0/3.d0 - u_c_t(i+2,j+1,k,1)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i,j+2,k)*lambda_c(i,j+2,k)*(u_c_t(i-2,j+2,k,1)/12.d0-u_c_t(i-1,j+2,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+2,k,1)*2.d0/3.d0 - u_c_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
             ! third set
             lh_c(i,j,k,3) = lh_c(i,j,k,3)+((-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/8.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i-2,j,k,3) &
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/6.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0)*u_c_t(i-1,j,k,3) &
                      +(-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/24.d0 - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*5.d0/6.d0 -Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/24.d0)*u_c_t(i-0,j,k,3) &
                      +(Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/2.d0 &
                      + Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/6.d0)*u_c_t(i+1,j,k,3) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/8.d0)*u_c_t(i+2,j,k,3))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/8.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i,j-2,k,3) &
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/6.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0)*u_c_t(i,j-1,k,3) &
                      +(-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/24.d0 - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*5.d0/6.d0 -Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/24.d0)*u_c_t(i,j-0,k,3) &
                      +(Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/2.d0 &
                      + Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/6.d0)*u_c_t(i,j+1,k,3) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/8.d0)*u_c_t(i,j+2,k,3))/h2_c**2/l2**2
          end do
       end do
    end do
    !
    do k = 1,4
       do i = 1,n2_c
          do j = 1,n1_c
             do k1 = 1,6
                ! mixed derivative 13 & 23
                ! first set equation
                lh_c(j,i,k,1) = lh_c(j,i,k,1) &
                      +(Jacobian_c(j-2,i,k)*(2.d0*mu_c(j-2,i,k)+lambda_c(j-2,i,k))*XI13_c(j-2,i,k)*bof(k,k1)&
                      *u_c_t(j-2,i,k1,1)/12.d0-Jacobian_c(j-1,i,k)*(2.d0*mu_c(j-1,i,k)+lambda_c(j-1,i,k))&
                      *XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0+Jacobian_c(j+1,i,k)*(2.d0*mu_c(j+1,i,k)&
                      +lambda_c(j+1,i,k))*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*(2.d0*mu_c(j+2,i,k)&
                      +lambda_c(j+2,i,k))*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI23_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,2)/12.d0 &
                      -Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI23_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI23_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI23_c(j+2,i,k)*bof(k,k1)&
                      *u_c_t(j+2,i,k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI33_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,3)/12.d0 &
                      -Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI33_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI33_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI33_c(j+2,i,k)&
                      *bof(k,k1)*u_c_t(j+2,i,k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,1)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI13_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI13_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI13_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI13_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,1) = lh_c(j,i,n3_c+1-k,1)+ &
                      (-Jacobian_c(j-2,i,n3_c+1-k)*(2.d0*mu_c(j-2,i,n3_c+1-k)+lambda_c(j-2,i,n3_c+1-k))&
                      *XI13_c(j-2,i,n3_c+1-k)*bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*(2.d0*mu_c(j-1,i,n3_c+1-k)+lambda_c(j-1,i,n3_c+1-k))&
                      *XI13_c(j-1,i,n3_c+1-k)*bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*(2.d0*mu_c(j+1,i,n3_c+1-k)+lambda_c(j+1,i,n3_c+1-k))&
                      *XI13_c(j+1,i,n3_c+1-k)*bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*(2.d0*mu_c(j+2,i,n3_c+1-k)+lambda_c(j+2,i,n3_c+1-k))&
                      *XI13_c(j+2,i,n3_c+1-k)*bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*lambda_c(j-2,i,n3_c+1-k)*XI23_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*lambda_c(j-1,i,n3_c+1-k)*XI23_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*lambda_c(j+1,i,n3_c+1-k)*XI23_c(j+1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*lambda_c(j+2,i,n3_c+1-k)*XI23_c(j+2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+2,i,n3_c+1-k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*lambda_c(j-2,i,n3_c+1-k)*XI33_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*lambda_c(j-1,i,n3_c+1-k)*XI33_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*lambda_c(j+1,i,n3_c+1-k)*XI33_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*lambda_c(j+2,i,n3_c+1-k)*XI33_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI23_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI23_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI23_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI13_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI13_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI13_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI13_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c
                ! second set equation
                lh_c(j,i,k,2) = lh_c(j,i,k,2) &
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI23_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,1)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI23_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI23_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI23_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,2)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI13_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,1)/12.d0 &
                      -Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI13_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI13_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI13_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*(2.d0*mu_c(j,i-2,k)+lambda_c(j,i-2,k))&
                      *XI23_c(j,i-2,k)*bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*(2.d0*mu_c(j,i-1,k)+lambda_c(j,i-1,k))*XI23_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*(2.d0*mu_c(j,i+1,k)+lambda_c(j,i+1,k))*XI23_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*(2.d0*mu_c(j,i+2,k)+lambda_c(j,i+2,k))*XI23_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI33_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,3)/12.d0 &
                      -Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI33_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI33_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI33_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,3)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,2) = lh_c(j,i,n3_c+1-k,2)+ &
                      (-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI23_c(j-2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI23_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI23_c(j+1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI23_c(j+2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI13_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI13_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI13_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI13_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*lambda_c(j,i-2,n3_c+1-k)&
                      *XI13_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*lambda_c(j,i-1,n3_c+1-k)*XI13_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*lambda_c(j,i+1,n3_c+1-k)*XI13_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*lambda_c(j,i+2,n3_c+1-k)*XI13_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*(2.d0*mu_c(j,i-2,n3_c+1-k)&
                      +lambda_c(j,i-2,n3_c+1-k))*XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*(2.d0*mu_c(j,i-1,n3_c+1-k)+lambda_c(j,i-1,n3_c+1-k))&
                      *XI23_c(j,i-1,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*(2.d0*mu_c(j,i+1,n3_c+1-k)+lambda_c(j,i+1,n3_c+1-k))&
                      *XI23_c(j,i+1,n3_c+1-k)*bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*(2.d0*mu_c(j,i+2,n3_c+1-k)+lambda_c(j,i+2,n3_c+1-k))&
                      *XI23_c(j,i+2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*lambda_c(j,i-2,n3_c+1-k)&
                      *XI33_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*lambda_c(j,i-1,n3_c+1-k)*XI33_c(j,i-1,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j,i-1,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*lambda_c(j,i+1,n3_c+1-k)*XI33_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*lambda_c(j,i+2,n3_c+1-k)*XI33_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,3)/12.d0)/l2/h2_c/h3_c
                ! third set equation
                lh_c(j,i,k,3) = lh_c(j,i,k,3) &
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI33_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,1)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI33_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI33_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI33_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,3)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI33_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI33_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI33_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI33_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,3)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,3)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,3) = lh_c(j,i,n3_c+1-k,3) +&
                      (-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI33_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI33_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI33_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI33_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI13_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI13_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI13_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI13_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI33_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI33_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI33_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI33_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI23_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI23_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI23_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,3)/12.d0)/l2/h2_c/h3_c
             end do
          end do
       end do
    end do
    !
    do k = 5,n3_c-4
       do i = 1,n2_c
          do j = 1,n1_c
             ! mixed derivative 13 & 23
             ! first set equation
             lh_c(j,i,k,1) = lh_c(j,i,k,1)+(Jacobian_c(j-2,i,k)*(2.d0*mu_c(j-2,i,k)+lambda_c(j-2,i,k))*XI13_c(j-2,i,k)&
                   *(u_c_t(j-2,i,k-2,1)/12.d0-u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*(2.d0*mu_c(j-1,i,k)+lambda_c(j-1,i,k))*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0&
                   -u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*(2.d0*mu_c(j+1,i,k)+lambda_c(j+1,i,k))&
                   *XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*(2.d0*mu_c(j+2,i,k)+lambda_c(j+2,i,k))*XI13_c(j+2,i,k)&
                   *(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI23_c(j-2,i,k)*(u_c_t(j-2,i,k-2,2)/12.d0&
                   -u_c_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,2)*2.d0/3.d0-u_c_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI23_c(j-1,i,k)*(u_c_t(j-1,i,k-2,2)/12.d0&
                   -u_c_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,2)*2.d0/3.d0-u_c_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI23_c(j+1,i,k)*(u_c_t(j+1,i,k-2,2)/12.d0&
                   -u_c_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,2)*2.d0/3.d0-u_c_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI23_c(j+2,i,k)*(u_c_t(j+2,i,k-2,2)/12.d0&
                   -u_c_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,2)*2.d0/3.d0-u_c_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI33_c(j-2,i,k)*(u_c_t(j-2,i,k-2,3)/12.d0&
                   -u_c_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,3)*2.d0/3.d0-u_c_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI33_c(j-1,i,k)*(u_c_t(j-1,i,k-2,3)/12.d0-u_c_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,3)*2.d0/3.d0-u_c_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI33_c(j+1,i,k)*(u_c_t(j+1,i,k-2,3)/12.d0-u_c_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,3)*2.d0/3.d0-u_c_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI33_c(j+2,i,k)*(u_c_t(j+2,i,k-2,3)/12.d0-u_c_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,3)*2.d0/3.d0-u_c_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,1)/12.d0-u_c_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,1)*2.d0/3.d0-u_c_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*(u_c_t(j,i-1,k-2,1)/12.d0-u_c_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,1)*2.d0/3.d0-u_c_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*(u_c_t(j,i+1,k-2,1)/12.d0-u_c_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,1)*2.d0/3.d0-u_c_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*(u_c_t(j,i+2,k-2,1)/12.d0-u_c_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,1)*2.d0/3.d0-u_c_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI13_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI13_c(j,i-1,k)*(u_c_t(j,i-1,k-2,2)/12.d0-u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI13_c(j,i+1,k)*(u_c_t(j,i+1,k-2,2)/12.d0-u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI13_c(j,i+2,k)*(u_c_t(j,i+2,k-2,2)/12.d0-u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c
             ! second set equation
             lh_c(j,i,k,2) = lh_c(j,i,k,2)+(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI23_c(j-2,i,k)*(u_c_t(j-2,i,k-2,1)/12.d0&
                   -u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI23_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0-u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI23_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI23_c(j+2,i,k)*(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*(u_c_t(j-2,i,k-2,2)/12.d0&
                   -u_c_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,2)*2.d0/3.d0-u_c_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,2)/12.d0-u_c_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,2)*2.d0/3.d0-u_c_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,2)/12.d0-u_c_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,2)*2.d0/3.d0-u_c_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*(u_c_t(j+2,i,k-2,2)/12.d0-u_c_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,2)*2.d0/3.d0-u_c_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI13_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,1)/12.d0-u_c_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,1)*2.d0/3.d0-u_c_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI13_c(j,i-1,k)*(u_c_t(j,i-1,k-2,1)/12.d0&
                   -u_c_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,1)*2.d0/3.d0-u_c_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI13_c(j,i+1,k)*(u_c_t(j,i+1,k-2,1)/12.d0&
                   -u_c_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,1)*2.d0/3.d0-u_c_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI13_c(j,i+2,k)*(u_c_t(j,i+2,k-2,1)/12.d0&
                   -u_c_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,1)*2.d0/3.d0-u_c_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*(2.d0*mu_c(j,i-2,k)+lambda_c(j,i-2,k))&
                   *XI23_c(j,i-2,k)*(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*(2.d0*mu_c(j,i-1,k)+lambda_c(j,i-1,k))*XI23_c(j,i-1,k)&
                   *(u_c_t(j,i-1,k-2,2)/12.d0-u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*(2.d0*mu_c(j,i+1,k)+lambda_c(j,i+1,k))*XI23_c(j,i+1,k)&
                   *(u_c_t(j,i+1,k-2,2)/12.d0-u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*(2.d0*mu_c(j,i+2,k)+lambda_c(j,i+2,k))*XI23_c(j,i+2,k)&
                   *(u_c_t(j,i+2,k-2,2)/12.d0-u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI33_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,3)/12.d0-u_c_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,3)*2.d0/3.d0-u_c_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI33_c(j,i-1,k)*(u_c_t(j,i-1,k-2,3)/12.d0&
                   -u_c_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,3)*2.d0/3.d0-u_c_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI33_c(j,i+1,k)*(u_c_t(j,i+1,k-2,3)/12.d0&
                   -u_c_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,3)*2.d0/3.d0-u_c_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI33_c(j,i+2,k)*(u_c_t(j,i+2,k-2,3)/12.d0&
                   -u_c_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,3)*2.d0/3.d0-u_c_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
             ! third set equation
             lh_c(j,i,k,3) = lh_c(j,i,k,3)+(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI33_c(j-2,i,k)*(u_c_t(j-2,i,k-2,1)/12.d0&
                   -u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI33_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0-u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI33_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI33_c(j+2,i,k)*(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*(u_c_t(j-2,i,k-2,3)/12.d0&
                   -u_c_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,3)*2.d0/3.d0-u_c_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,3)/12.d0-u_c_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,3)*2.d0/3.d0-u_c_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,3)/12.d0-u_c_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,3)*2.d0/3.d0-u_c_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*(u_c_t(j+2,i,k-2,3)/12.d0-u_c_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,3)*2.d0/3.d0-u_c_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI33_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI33_c(j,i-1,k)*(u_c_t(j,i-1,k-2,2)/12.d0&
                   -u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI33_c(j,i+1,k)*(u_c_t(j,i+1,k-2,2)/12.d0&
                   -u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI33_c(j,i+2,k)*(u_c_t(j,i+2,k-2,2)/12.d0&
                   -u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,3)/12.d0-u_c_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,3)*2.d0/3.d0-u_c_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*(u_c_t(j,i-1,k-2,3)/12.d0&
                   -u_c_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,3)*2.d0/3.d0-u_c_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*(u_c_t(j,i+1,k-2,3)/12.d0&
                   -u_c_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,3)*2.d0/3.d0-u_c_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*(u_c_t(j,i+2,k-2,3)/12.d0&
                   -u_c_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,3)*2.d0/3.d0-u_c_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
           end do
        end do
    end do
    !
    do i = 5,n3_c-4
       do j = 1,n2_c
          do k = 1,n1_c
             ! mixed derivative 31
             ! first set equation
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+(Jacobian_c(k,j,i-2)*(2.d0*mu_c(k,j,i-2)&
                     +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*(2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)&
                     *(u_c_t(k-2,j,i-1,1)/12.d0-u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*(2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)&
                     *(u_c_t(k-2,j,i+1,1)/12.d0-u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*(2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)&
                     *(u_c_t(k-2,j,i+2,1)/12.d0-u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)*(u_c_t(k-2,j,i-2,2)/12.d0&
                     -u_c_t(k-1,j,i-2,2)*2.d0/3.d0+u_c_t(k+1,j,i-2,2)*2.d0/3.d0-u_c_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k-2,j,i-1,2)/12.d0&
                     -u_c_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,2)*2.d0/3.d0-u_c_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k-2,j,i+1,2)/12.d0&
                     -u_c_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,2)*2.d0/3.d0-u_c_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k-2,j,i+2,2)/12.d0&
                     -u_c_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,2)*2.d0/3.d0-u_c_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI33_c(k,j,i-2)*(u_c_t(k-2,j,i-2,3)/12.d0&
                     -u_c_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,3)*2.d0/3.d0-u_c_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k-2,j,i-1,3)/12.d0&
                     -u_c_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,3)*2.d0/3.d0-u_c_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k-2,j,i+1,3)/12.d0&
                     -u_c_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,3)*2.d0/3.d0-u_c_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k-2,j,i+2,3)/12.d0&
                     -u_c_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,3)*2.d0/3.d0-u_c_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c
             ! second set equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI23_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k-2,j,i-1,1)/12.d0&
                     -u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k-2,j,i+1,1)/12.d0&
                     -u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k-2,j,i+2,1)/12.d0&
                     -u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,2)/12.d0-u_c_t(k-1,j,i-2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,2)*2.d0/3.d0-u_c_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k-2,j,i-1,2)/12.d0&
                     -u_c_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,2)*2.d0/3.d0-u_c_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k-2,j,i+1,2)/12.d0&
                     -u_c_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,2)*2.d0/3.d0-u_c_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)&
                     *(u_c_t(k-2,j,i+2,2)/12.d0-u_c_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,2)*2.d0/3.d0-u_c_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI33_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k-2,j,i-1,1)/12.d0&
                     -u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k-2,j,i+1,1)/12.d0&
                     -u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k-2,j,i+2,1)/12.d0&
                     -u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,3)/12.d0-u_c_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,3)*2.d0/3.d0-u_c_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k-2,j,i-1,3)/12.d0&
                     -u_c_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,3)*2.d0/3.d0-u_c_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k-2,j,i+1,3)/12.d0&
                     -u_c_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,3)*2.d0/3.d0-u_c_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k-2,j,i+2,3)/12.d0&
                     -u_c_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,3)*2.d0/3.d0-u_c_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h3_c/h1_c


             ! 32
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,1)/12.d0-u_c_t(k,j-1,i-2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,1)*2.d0/3.d0-u_c_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k,j-2,i-1,1)/12.d0&
                   -u_c_t(k,j-1,i-1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,1)*2.d0/3.d0-u_c_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k,j-2,i+1,1)/12.d0&
                   -u_c_t(k,j-1,i+1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,1)*2.d0/3.d0-u_c_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k,j-2,i+2,1)/12.d0&
                   -u_c_t(k,j-1,i+2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,1)*2.d0/3.d0-u_c_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI13_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k,j-2,i-1,2)/12.d0&
                   -u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k,j-2,i+1,2)/12.d0&
                   -u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k,j-2,i+2,2)/12.d0&
                   -u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_c/h2_c
             ! second set equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,1)/12.d0-u_c_t(k,j-1,i-2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,1)*2.d0/3.d0-u_c_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k,j-2,i-1,1)/12.d0&
                   -u_c_t(k,j-1,i-1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,1)*2.d0/3.d0-u_c_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k,j-2,i+1,1)/12.d0&
                   -u_c_t(k,j-1,i+1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,1)*2.d0/3.d0-u_c_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k,j-2,i+2,1)/12.d0&
                   -u_c_t(k,j-1,i+2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,1)*2.d0/3.d0-u_c_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*(2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))&
                   *XI23_c(k,j,i-2)*(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*(2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)&
                   *(u_c_t(k,j-2,i-1,2)/12.d0-u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*(2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)&
                   *(u_c_t(k,j-2,i+1,2)/12.d0-u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*(2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)&
                   *(u_c_t(k,j-2,i+2,2)/12.d0-u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI33_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,3)/12.d0-u_c_t(k,j-1,i-2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,3)*2.d0/3.d0-u_c_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k,j-2,i-1,3)/12.d0&
                   -u_c_t(k,j-1,i-1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,3)*2.d0/3.d0-u_c_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k,j-2,i+1,3)/12.d0&
                   -u_c_t(k,j-1,i+1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,3)*2.d0/3.d0-u_c_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k,j-2,i+2,3)/12.d0&
                   -u_c_t(k,j-1,i+2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,3)*2.d0/3.d0-u_c_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h3_c/h2_c
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI33_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k,j-2,i-1,2)/12.d0&
                   -u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI33_c(k,j,i+1)&
                   *(u_c_t(k,j-2,i+1,2)/12.d0-u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k,j-2,i+2,2)/12.d0&
                   -u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,3)/12.d0-u_c_t(k,j-1,i-2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,3)*2.d0/3.d0-u_c_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k,j-2,i-1,3)/12.d0&
                   -u_c_t(k,j-1,i-1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,3)*2.d0/3.d0-u_c_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k,j-2,i+1,3)/12.d0&
                   -u_c_t(k,j-1,i+1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,3)*2.d0/3.d0-u_c_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k,j-2,i+2,3)/12.d0&
                   -u_c_t(k,j-1,i+2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,3)*2.d0/3.d0-u_c_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
          end do
       end do
    end do
    do i = 1,4
       do j = 1,n2_c
          do k = 1,n1_c
             do k1 = 1,6
                ! mixed derivative 31
                ! first set equation
                lh_c(k,j,i,1) = lh_c(k,j,i,1)+(bof(i,k1)*Jacobian_c(k,j,k1)*(2.d0*mu_c(k,j,k1)+lambda_c(k,j,k1))&
                      *XI13_c(k,j,k1)*(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h3_c/h1_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,2)/12.d0-u_c_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,2)*2.d0/3.d0-u_c_t(k+2,j,k1,2)/12.d0))/l1/h3_c/h1_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,3)/12.d0-u_c_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,3)*2.d0/3.d0-u_c_t(k+2,j,k1,3)/12.d0))/l1/h1_c/h3_c

                lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*(2.d0&
                   *mu_c(k,j,n3_c+1-k1)+lambda_c(k,j,n3_c+1-k1))*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h3_c/h1_c&
                   +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                   *XI23_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,2)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,2)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,2)/12.d0))/l1/h3_c/h1_c+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                   *XI33_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,3)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,3)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,3)/12.d0))/l1/h1_c/h3_c
                ! second set equation
                lh_c(k,j,i,2) = lh_c(k,j,i,2)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,2)/12.d0-u_c_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,2)*2.d0/3.d0-u_c_t(k+2,j,k1,2)/12.d0))/l1/h1_c/h3_c

                lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI23_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,2)/12.d0))/l1/h1_c/h3_c
                ! third set equation
                lh_c(k,j,i,3) = lh_c(k,j,i,3)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,3)/12.d0-u_c_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,3)*2.d0/3.d0-u_c_t(k+2,j,k1,3)/12.d0))/l1/h3_c/h1_c

                lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h3_c/h1_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                      *XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,3)/12.d0))/l1/h3_c/h1_c



                !31
                ! first set equation
                lh_c(k,j,i,1) = lh_c(k,j,i,1)+(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,1)/12.d0-u_c_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,1)*2.d0/3.d0-u_c_t(k,j+2,k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h3_c/h2_c

                lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                      *XI23_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h3_c/h2_c
                ! second set equation
                lh_c(k,j,i,2) = lh_c(k,j,i,2)+(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,1)/12.d0-u_c_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,1)*2.d0/3.d0-u_c_t(k,j+2,k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*(2.d0*mu_c(k,j,k1)&
                      +lambda_c(k,j,k1))*XI23_c(k,j,k1)*(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,3)/12.d0-u_c_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,3)*2.d0/3.d0-u_c_t(k,j+2,k1,3)/12.d0))/l2/h3_c/h2_c

                lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *(2.d0*mu_c(k,j,n3_c+1-k1)+lambda_c(k,j,n3_c+1-k1))*XI23_c(k,j,n3_c+1-k1)&
                      *(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0-u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0&
                      +u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0-u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,3)/12.d0))/l2/h3_c/h2_c
                ! third set equation
                lh_c(k,j,i,3) = lh_c(k,j,i,3)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h2_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,3)/12.d0-u_c_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,3)*2.d0/3.d0-u_c_t(k,j+2,k1,3)/12.d0))/l2/h2_c/h3_c

                lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h2_c/h3_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI23_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,3)/12.d0))/l2/h2_c/h3_c
             end do
          end do
       end do
    end do
    !
    do i = 7,n3_c-6
       do j = 1,n2_c
          do k = 1,n1_c
             ! second derivative 33
             ! first set equation
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2 &
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2 &
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/6.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/2.d0&
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2 &
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))*5.d0/6.d0 -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2&
                   +XI33_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0 &
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/8.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/24.d0&
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
             ! second set of equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)*5.d0/6.d0&
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/6.d0 &
                   + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/2.d0 + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI23_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0&
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/6.d0&
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/6.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/2.d0&
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI33_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))*5.d0/6.d0&
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2&
                   +XI23_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/6.d0 + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI33_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI23_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/6.d0 &
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2&
                   +XI23_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
          end do
       end do
    end do
    do j = 1,n2_c
       do k = 1,n1_c
          do i = 1,6
             do k1 = 1,8
                do m = 1,8
                   ! second derivative 33
                   ! first set equation
                   lh_c(k,j,i,1) = lh_c(k,j,i,1) +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)**2+mu_c(k,j,m)*(XI23_c(k,j,m)**2+XI33_c(k,j,m)**2))*u_c_t(k,j,k1,1))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)*XI23_c(k,j,m)*u_c_t(k,j,k1,2))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI23_c(k,j,n3_c+1-m)**2&
                     +XI33_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                   ! second set equation
                   lh_c(k,j,i,2) = lh_c(k,j,i,2) +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI13_c(k,j,m)*XI23_c(k,j,m)*u_c_t(k,j,k1,1))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI23_c(k,j,m)**2+mu_c(k,j,m)*(XI13_c(k,j,m)**2&
                       +XI33_c(k,j,m)**2))*u_c_t(k,j,k1,2))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI23_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI13_c(k,j,n3_c+1-m)**2&
                       +XI33_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                   ! third set equation
                   lh_c(k,j,i,3) = lh_c(k,j,i,3) +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI13_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,1))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI23_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,2))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI33_c(k,j,m)**2+mu_c(k,j,m)*(XI13_c(k,j,m)**2&
                        +XI23_c(k,j,m)**2))*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI33_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI13_c(k,j,n3_c+1-m)**2&
                        +XI23_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                end do
             end do
          end do
          ! first set equation
          lh_c(k,j,1,1) = lh_c(k,j,1,1) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)**2+mu_c(k,j,1)*(XI23_c(k,j,1)**2+XI33_c(k,j,1)**2)))/h3_c**2&
                + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)*XI23_c(k,j,1))/h3_c**2&
                + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)*XI33_c(k,j,1))/h3_c**2

          lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
                + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
                + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! second set equation
          lh_c(k,j,1,2) = lh_c(k,j,1,2) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI13_c(k,j,1)*XI23_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)**2+mu_c(k,j,1)*(XI13_c(k,j,1)**2+XI33_c(k,j,1)**2)))/h3_c**2&
               + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)*XI33_c(k,j,1))/h3_c**2

          lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! third set equation
          lh_c(k,j,1,3) = lh_c(k,j,1,3) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI13_c(k,j,1)*XI33_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)*XI33_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI33_c(k,j,1)**2+mu_c(k,j,1)*(XI13_c(k,j,1)**2+XI23_c(k,j,1)**2)))/h3_c**2

          lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
  !
  end subroutine Update_interior

  subroutine Update_gp(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,u_c,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3,u_f,tv,index) bind(c)
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
    integer,value,intent(in) ::  index
    real(dp) tv
    real(dp) Xgrid_c_1(1-nrg:n1_c+nrg)
    real(dp) Xgrid_c_2(1-nrg:n2_c+nrg)
    real(dp) Xgrid_f_1(1-nrg:n1_f+nrg)
    real(dp) Xgrid_f_2(1-nrg:n2_f+nrg)
    real(dp) Xgrid_c_3(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
    real(dp) Xgrid_f_3(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)
    real(dp) u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
    real(dp) u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
    integer :: i,j,k
    ! Update ghost point values on the left and right domain
    ! fine mesh
    do i=1-nrg,n3_f+nrg
       do j = 1-nrg,n2_f+nrg
          do k = 1-nrg,0
             call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
          do k = n1_f+1,n1_f+nrg
             call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
    end do
    ! coarse mesh
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,0
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
          do k = n1_c+1,n1_c+nrg
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do
    ! Update ghost point values on the front and back domain
    ! fine mesh
    do i=1-nrg,n3_f+nrg
       do j = 1-nrg,0
          do k = 1,n1_f
             call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
       do j = n2_f+1,n2_f+nrg
          do k = 1,n1_f
             call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
    end do
    ! coarse mesh
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,0
          do k = 1,n1_c
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
       do j = n2_c+1,n2_c+nrg
          do k = 1,n1_c
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do
  end subroutine Update_gp

  subroutine Update_traction(traction_rhs,traction_data,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3,&
       u_f,mu_f,lambda_f,Jacobian_f,XI13_f,XI23_f,XI33_f,Sb,tv,index) bind(c)
    ! traction B.C. on the top of the fine mesh
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
    
 

    integer,value,intent(in)  :: index
    real(dp) :: traction_rhs(3)
    real(dp) :: traction_data(1:n1_f,1:n2_f,1:dim)
    real(dp) :: mat_det,tv
    real(dp) Xgrid_f_1(1-nrg:n1_f+nrg)
    real(dp) Xgrid_f_2(1-nrg:n2_f+nrg)
    real(dp) Xgrid_f_3(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)
    real(dp) u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: mu_f,lambda_f,Jacobian_f,XI13_f,XI23_f,XI33_f
    real(dp), dimension (0:4):: Sb
    integer :: i,j,k
    !
    do j = 1,n2_f
       do i=1,n1_f
          call top_normal_data(Xgrid_f_1(i),Xgrid_f_2(j),Xgrid_f_3(i,j,n3_f), &
              l1,l2,tv,traction_data(i,j,:))
       end do
    end do
    !
    do j = 1,n2_f
       do i = 1,n1_f
          traction_rhs = 0.d0
          ! 33
          ! first set equation
          do k = 1,4
             ! first component
             traction_rhs(1) = traction_rhs(1) &
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,1,index)*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)&
                +lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,2,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,3,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)
          end do
          ! 31 & 32
          ! first component
          traction_rhs(1) = traction_rhs(1) + Jacobian_f(i,j,n3_f)*(2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))/l1&
             *XI13_f(i,j,n3_f)*(u_f(i-2,j,n3_f,1,index)/12.d0-u_f(i-1,j,n3_f,1,index)*2.d0/3.d0&
             +u_f(i+1,j,n3_f,1,index)*2.d0/3.d0-u_f(i+2,j,n3_f,1,index)/12.d0)/h1_f + Jacobian_f(i,j,n3_f)&
             *mu_f(i,j,n3_f)/l2*XI23_f(i,j,n3_f)*(u_f(i,j-2,n3_f,1,index)/12.d0&
             -u_f(i,j-1,n3_f,1,index)*2.d0/3.d0+u_f(i,j+1,n3_f,1,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,1,index)/12.d0)/h2_f&
             + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l1*XI23_f(i,j,n3_f)&
             *(u_f(i-2,j,n3_f,2,index)/12.d0-u_f(i-1,j,n3_f,2,index)*2.d0/3.d0+u_f(i+1,j,n3_f,2,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,2,index)/12.d0)/h1_f + Jacobian_f(i,j,n3_f)*lambda_f(i,j,n3_f)/l2*XI13_f(i,j,n3_f)&
             *(u_f(i,j-2,n3_f,2,index)/12.d0-u_f(i,j-1,n3_f,2,index)*2.d0/3.d0+u_f(i,j+1,n3_f,2,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,2,index)/12.d0)/h2_f&
             + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l1*XI33_f(i,j,n3_f)&
             *(u_f(i-2,j,n3_f,3,index)/12.d0-u_f(i-1,j,n3_f,3,index)*2.d0/3.d0+u_f(i+1,j,n3_f,3,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,3,index)/12.d0)/h1_f
          ! scale
          traction_rhs(1) = traction_rhs(1) &
            - sqrt(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2)*Jacobian_f(i,j,n3_f)*traction_data(i,j,1)
          ! 33
          ! second set equation
          do k = 1,4
             ! first component
             traction_rhs(2) = traction_rhs(2) &
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,1,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                *XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,2,index)*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
                *XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,3,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                *XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)
          end do
          ! 31 & 32
          ! first component
          traction_rhs(2) = traction_rhs(2) + Jacobian_f(i,j,n3_f)*lambda_f(i,j,n3_f)/l1&
             *XI23_f(i,j,n3_f)*(u_f(i-2,j,n3_f,1,index)/12.d0&
             -u_f(i-1,j,n3_f,1,index)*2.d0/3.d0+u_f(i+1,j,n3_f,1,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,1,index)/12.d0)/h1_f + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l2*XI13_f(i,j,n3_f)&
             *(u_f(i,j-2,n3_f,1,index)/12.d0-u_f(i,j-1,n3_f,1,index)*2.d0/3.d0+u_f(i,j+1,n3_f,1,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,1,index)/12.d0)/h2_f&
             + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l1*XI13_f(i,j,n3_f)&
             *(u_f(i-2,j,n3_f,2,index)/12.d0&
             -u_f(i-1,j,n3_f,2,index)*2.d0/3.d0+u_f(i+1,j,n3_f,2,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,2,index)/12.d0)/h1_f + Jacobian_f(i,j,n3_f)*(2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))/l2&
             *XI23_f(i,j,n3_f)*(u_f(i,j-2,n3_f,2,index)/12.d0-u_f(i,j-1,n3_f,2,index)*2.d0/3.d0&
             +u_f(i,j+1,n3_f,2,index)*2.d0/3.d0-u_f(i,j+2,n3_f,2,index)/12.d0)/h2_f&
             + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l2*XI33_f(i,j,n3_f)&
             *(u_f(i,j-2,n3_f,3,index)/12.d0-u_f(i,j-1,n3_f,3,index)*2.d0/3.d0+u_f(i,j+1,n3_f,3,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,3,index)/12.d0)/h2_f
          ! scale
          traction_rhs(2) = traction_rhs(2) &
            - sqrt(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2)*Jacobian_f(i,j,n3_f)*traction_data(i,j,2)
          ! 33
          ! third set equation
          do k = 1,4
             ! first component
             traction_rhs(3) = traction_rhs(3) &
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,1,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,2,index)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                +mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)&
                - 1.d0/h3_f*Sb(k)*u_f(i,j,n3_f+1-k,3,index)*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)&
                +lambda_f(i,j,n3_f))*XI33_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))
          end do
          ! 31 & 32
          ! first component
          traction_rhs(3) = traction_rhs(3) + Jacobian_f(i,j,n3_f)*lambda_f(i,j,n3_f)/l1&
             *XI33_f(i,j,n3_f)*(u_f(i-2,j,n3_f,1,index)/12.d0&
             -u_f(i-1,j,n3_f,1,index)*2.d0/3.d0+u_f(i+1,j,n3_f,1,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,1,index)/12.d0)/h1_f&
             + Jacobian_f(i,j,n3_f)*lambda_f(i,j,n3_f)/l2&
             *XI33_f(i,j,n3_f)*(u_f(i,j-2,n3_f,2,index)/12.d0&
             -u_f(i,j-1,n3_f,2,index)*2.d0/3.d0+u_f(i,j+1,n3_f,2,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,2,index)/12.d0)/h2_f&
             + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l1&
             *XI13_f(i,j,n3_f)*(u_f(i-2,j,n3_f,3,index)/12.d0&
             -u_f(i-1,j,n3_f,3,index)*2.d0/3.d0+u_f(i+1,j,n3_f,3,index)*2.d0/3.d0&
             -u_f(i+2,j,n3_f,3,index)/12.d0)/h1_f + Jacobian_f(i,j,n3_f)*mu_f(i,j,n3_f)/l2&
             *XI23_f(i,j,n3_f)*(u_f(i,j-2,n3_f,3,index)/12.d0&
             -u_f(i,j-1,n3_f,3,index)*2.d0/3.d0+u_f(i,j+1,n3_f,3,index)*2.d0/3.d0&
             -u_f(i,j+2,n3_f,3,index)/12.d0)/h2_f
          ! scale
          traction_rhs(3) = traction_rhs(3) &
            - sqrt(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2)*Jacobian_f(i,j,n3_f)*traction_data(i,j,3)

          ! Update ghost point at the traction boundary
          ! compute the determinant of the coefficient matrix
          mat_det = Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
                  *(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
                  *XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)&
                  *((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI33_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
                  *(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))+Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                  *XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                  +mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                  +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)+Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                  +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                  +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                  *XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                  *XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)&
                  +lambda_f(i,j,n3_f))*XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))&
                  *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f) &
                  - Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
                  *(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                  *XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)&
                  *XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI33_f(i,j,n3_f)**2&
                  +mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
                  +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
                  *XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)
          ! update the first component
          u_f(i,j,n3_f+1,1,index) = h3_f/Sb(0)/mat_det* &
             ((Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
             *(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
             *XI33_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))-Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f))*traction_rhs(1) &
             +(Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)&
             *XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI33_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
             *(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2)))*traction_rhs(2)+(Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
             +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
             *XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)&
             *XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
             *XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2)))*traction_rhs(3))
          ! update the second component
          u_f(i,j,n3_f+1,2,index) = h3_f/Sb(0)/mat_det* &
             ((Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)&
             *XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI33_f(i,j,n3_f)**2&
             +mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
             +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f))*traction_rhs(1)+(Jacobian_f(i,j,n3_f)&
             *((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
             *(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
             *XI33_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI23_f(i,j,n3_f)**2))-Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)*Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f))*traction_rhs(2) &
             +(Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
             -Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
             *(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))&
             *XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f))*traction_rhs(3))
          ! update the third component
          u_f(i,j,n3_f+1,3,index) = h3_f/Sb(0)/mat_det* &
             ((Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)&
             *XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI23_f(i,j,n3_f)**2&
             +mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)&
             +mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI33_f(i,j,n3_f))*traction_rhs(1)+(Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)&
             *XI33_f(i,j,n3_f)-Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2&
             +mu_f(i,j,n3_f)*(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI23_f(i,j,n3_f)*XI33_f(i,j,n3_f))*traction_rhs(2) &
             +(Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))*XI13_f(i,j,n3_f)**2+mu_f(i,j,n3_f)&
             *(XI23_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))*Jacobian_f(i,j,n3_f)*((2.d0*mu_f(i,j,n3_f)+lambda_f(i,j,n3_f))&
             *XI23_f(i,j,n3_f)**2+mu_f(i,j,n3_f)*(XI13_f(i,j,n3_f)**2+XI33_f(i,j,n3_f)**2))-Jacobian_f(i,j,n3_f)&
             *(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f)&
             *Jacobian_f(i,j,n3_f)*(lambda_f(i,j,n3_f)+mu_f(i,j,n3_f))*XI13_f(i,j,n3_f)*XI23_f(i,j,n3_f))*traction_rhs(3))
       end do
    end do
  end subroutine Update_traction




  
  subroutine Update_Dirichlet_BC(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,u_c,tv,index) bind(c)

    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none

    integer,value,intent(in) ::  index
    real(dp) Xgrid_c_1(1-nrg:n1_c+nrg),tv
    real(dp) Xgrid_c_2(1-nrg:n2_c+nrg)
    real(dp) Xgrid_c_3(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)
    real(dp) u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
    integer :: i,j,k
    ! fine mesh
    do j=1-nrg,n2_f+nrg
       do k=1-nrg,n1_f+nrg
          i = 1
          !call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(i),tv, &
                    !u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          i = n3_f
          !call exact_solution(Xgrid_f_1(k),Xgrid_f_2(j),Xgrid_f_3(k,j,i),tv, &
          !          u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
       end do
    end do
    !coarse mesh
    do j=1-nrg,n2_c+nrg
       do k=1-nrg,n1_c+nrg
          i = 1
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          i = n3_c
          !call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,3),tv, &
                   !u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
       end do
    end do
  end subroutine Update_Dirichlet_BC

  subroutine print_uf(u_f) bind(c)
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
    !real(dp) :: u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
    real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
    open (unit = 7, file = "compare.dat")
    write(7,*)u_f
    close(7)
  end subroutine print_uf
  
subroutine print_f(f) bind(c)
    use iso_fortran_env
    use iso_c_binding
    use problemsetup_new_3d
    implicit none
    !real(dp) :: u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4)
    !real(dp) :: u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4)
    real(dp) ::  f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim)
    open (unit = 7, file = "compare.dat")
    write(7,*)f
    close(7)
  end subroutine print_f
