module SBP_operator
  use iso_fortran_env
  use problemsetup_new_3d
  !use problemsetup_3d
  implicit none

contains
  subroutine VARCOEFFS4(acof, ghcof, Sb)
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

  subroutine dx_46(bop)
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

  subroutine varcoeff_NoGhost( acof_no_gp, ghcof_no_gp, sbop_no_gp )
    !use problemsetup_new_3d, only : dp
    use problemsetup_new_3d, only : dp
    !use SBP_operator
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

  subroutine interpolation(P)
    real(dp) :: P(-1:2)
    P(-1) = -1.d0/16.d0
    P(0) = 9.d0/16.d0
    P(1) = 9.d0/16.d0
    P(2) = -1.d0/16.d0
  end subroutine interpolation

  subroutine restriction(Rop)
    real(dp) :: Rop(-4:2)
    Rop(-4) = -1.d0/32.d0
    Rop(-3) = 0.d0
    Rop(-2) = 9.d0/32.d0
    Rop(-1) = 1.d0/2.d0
    Rop(0) = 9.d0/32.d0
    Rop(1) = 0.d0
    Rop(2) = -1.d0/32.d0
  end subroutine restriction

  subroutine interpolation_restriction(P,Rop,RPop)
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

  subroutine central_difference(ux_cof,uxx_cof)
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


end module SBP_operator
