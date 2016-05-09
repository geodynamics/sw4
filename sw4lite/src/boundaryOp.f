!  SW4 LICENSE
! # ----------------------------------------------------------------------
! # SW4 - Seismic Waves, 4th order
! # ----------------------------------------------------------------------
! # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
! # Produced at the Lawrence Livermore National Laboratory. 
! # 
! # Written by:
! # N. Anders Petersson (petersson1@llnl.gov)
! # Bjorn Sjogreen      (sjogreen2@llnl.gov)
! # 
! # LLNL-CODE-643337 
! # 
! # All rights reserved. 
! # 
! # This file is part of SW4, Version: 1.0
! # 
! # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
! # 
! # This program is free software; you can redistribute it and/or modify
! # it under the terms of the GNU General Public License (as published by
! # the Free Software Foundation) version 2, dated June 1991. 
! # 
! # This program is distributed in the hope that it will be useful, but
! # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
! # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
! # conditions of the GNU General Public License for more details. 
! # 
! # You should have received a copy of the GNU General Public License
! # along with this program; if not, write to the Free Software
! # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
c-----------------------------------------------------------------------
      subroutine WAVEPROPBOP_4( iop, iop2, bop, bop2, gh2, h, s )

***********************************************************************
***
*** SBP operator of order 2/4 for approximating 1st derivative. 5pt stencil.
***
***********************************************************************

      implicit none
      real*8 bop(4,6), bop2(4,6), gh2, h(4), s(0:4), iop(5), iop2(5)

*** Norm 
      h(1) = 17d0/48
      h(2) = 59d0/48
      h(3) = 43d0/48
      h(4) = 49d0/48

*** First derivative interior, 4th order
      iop(1) =  1d0/12
      iop(2) = -2d0/3
      iop(3) =  0d0
      iop(4) =  2d0/3
      iop(5) = -1d0/12

*** First derivative bop, 2nd order
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

*** Second derivative interior, 4th order
      iop2(1) = -1d0/12
      iop2(2) =  4d0/3
      iop2(3) = -5d0/2
      iop2(4) =  4d0/3
      iop2(5) = -1d0/12

*** Second derivative bop, 2nd order
      gh2 = 12d0/17
      bop2(1,1) = -14.D0/17.D0
      bop2(1,2) = -13.D0/17.D0
      bop2(1,3) = 20.D0/17.D0
      bop2(1,4) = -5.D0/17.D0
      bop2(1,5) = 0
      bop2(1,6) = 0

      bop2(2,1) = 1
      bop2(2,2) = -2
      bop2(2,3) = 1
      bop2(2,4) = 0
      bop2(2,5) = 0
      bop2(2,6) = 0

      bop2(3,1) = -4.D0/43.D0
      bop2(3,2) = 59.D0/43.D0
      bop2(3,3) = -110.D0/43.D0
      bop2(3,4) = 59.D0/43.D0
      bop2(3,5) = -4.D0/43.D0
      bop2(3,6) = 0

      bop2(4,1) = -1.D0/49.D0
      bop2(4,2) = 0
      bop2(4,3) = 59.D0/49.D0
      bop2(4,4) = -118.D0/49.D0
      bop2(4,5) = 64d0/49
      bop2(4,6) = -4d0/49

*** boundary derivative, 4th order
      s(0) = -1d0/4
      s(1) = -5d0/6
      s(2) =  3d0/2
      s(3) = -1d0/2
      s(4) =  1d0/12

      end

c-----------------------------------------------------------------------
      subroutine WAVEPROPBOP_6( iop, iop2, bop, bop2, gh2, h, s )

***********************************************************************
***
*** SBP operator of order 3/6 for approximating 1st derivative. 7pt stencil
***
***********************************************************************

      implicit none
      real*8 bop(6,9), bop2(6,9), gh2, h(6), s(0:6), iop(7), iop2(7)

*** Norm 
      h(1) = 13649.D0/43200.D0
      h(2) = 12013.D0/8640.D0
      h(3) = 2711.D0/4320.D0
      h(4) = 5359.D0/4320.D0
      h(5) = 7877.D0/8640.D0
      h(6) = 43801.D0/43200.D0

*** First derivative interior, 6th order
      iop(1) = -1d0/60
      iop(2) =  3d0/20
      iop(3) = -3d0/4
      iop(4) =  0d0
      iop(5) = -iop(3)
      iop(6) = -iop(2)
      iop(7) = -iop(1)

*** First derivative bop, 3nd order
      bop(1,1) = -21600.D0/13649.D0
      bop(1,2) = 2089388.D0/1023675.D0
      bop(1,3) = -352679.D0/2047350.D0
      bop(1,4) = -137949.D0/341225.D0
      bop(1,5) = 150971.D0/2047350.D0
      bop(1,6) = 45313.D0/1023675.D0
      bop(1,7) = 0
      bop(1,8) = 0
      bop(1,9) = 0
      bop(2,1) = -2089388.D0/4504875.D0
      bop(2,2) = 0
      bop(2,3) = 91517.D0/300325.D0
      bop(2,4) = 403421.D0/1801950.D0
      bop(2,5) = -12887.D0/300325.D0
      bop(2,6) = -65743.D0/3003250.D0
      bop(2,7) = 0
      bop(2,8) = 0
      bop(2,9) = 0
      bop(3,1) = 352679.D0/4066500.D0
      bop(3,2) = -91517.D0/135550.D0
      bop(3,3) = 0
      bop(3,4) = 139001.D0/203325.D0
      bop(3,5) = -27193.D0/271100.D0
      bop(3,6) = 3451.D0/677750.D0
      bop(3,7) = 0
      bop(3,8) = 0
      bop(3,9) = 0
      bop(4,1) = 137949.D0/1339750.D0
      bop(4,2) = -403421.D0/1607700.D0
      bop(4,3) = -139001.D0/401925.D0
      bop(4,4) = 0
      bop(4,5) = 433951.D0/803850.D0
      bop(4,6) = -478079.D0/8038500.D0
      bop(4,7) = 72.D0/5359.D0
      bop(4,8) = 0
      bop(4,9) = 0
      bop(5,1) = -150971.D0/5907750.D0
      bop(5,2) = 12887.D0/196925.D0
      bop(5,3) = 27193.D0/393850.D0
      bop(5,4) = -433951.D0/590775.D0
      bop(5,5) = 0
      bop(5,6) = 759996.D0/984625.D0
      bop(5,7) = -1296.D0/7877.D0
      bop(5,8) = 144.D0/7877.D0
      bop(5,9) = 0
      bop(6,1) = -45313.D0/3285075.D0
      bop(6,2) = 65743.D0/2190050.D0
      bop(6,3) = -3451.D0/1095025.D0
      bop(6,4) = 478079.D0/6570150.D0
      bop(6,5) = -759996.D0/1095025.D0
      bop(6,6) = 0
      bop(6,7) = 32400.D0/43801.D0
      bop(6,8) = -6480.D0/43801.D0
      bop(6,9) = 720.D0/43801.D0


*** Second derivative interior, 6th order
      iop2(1) = 1d0/90 
      iop2(2) =-3d0/20
      iop2(3) = 3d0/2
      iop2(4) = -49d0/18
      iop2(5) = iop2(3)
      iop2(6) = iop2(2)
      iop2(7) = iop2(1)

*** Second derivative bop, 3rd order
      gh2 = 7200d0/13649
      bop2(1,1) = 10490.D0/40947.D0
      bop2(1,2) = -178907.D0/54596.D0
      bop2(1,3) = 163609.D0/40947.D0
      bop2(1,4) = -147397.D0/81894.D0
      bop2(1,5) = 3747.D0/13649.D0
      bop2(1,6) = 3755.D0/163788.D0
      bop2(1,7) = 0
      bop2(1,8) = 0
      bop2(1,9) = 0
      bop2(2,1) = 6173.D0/5860.D0
      bop2(2,2) = -2066.D0/879.D0
      bop2(2,3) = 3283.D0/1758.D0
      bop2(2,4) = -303.D0/293.D0
      bop2(2,5) = 2111.D0/3516.D0
      bop2(2,6) = -601.D0/4395.D0
      bop2(2,7) = 0
      bop2(2,8) = 0
      bop2(2,9) = 0
      bop2(3,1) = -52391.D0/81330.D0
      bop2(3,2) = 134603.D0/32532.D0
      bop2(3,3) = -21982.D0/2711.D0
      bop2(3,4) = 112915.D0/16266.D0
      bop2(3,5) = -46969.D0/16266.D0
      bop2(3,6) = 30409.D0/54220.D0
      bop2(3,7) = 0
      bop2(3,8) = 0
      bop2(3,9) = 0
      bop2(4,1) = 68603.D0/321540.D0
      bop2(4,2) = -12423.D0/10718.D0
      bop2(4,3) = 112915.D0/32154.D0
      bop2(4,4) = -75934.D0/16077.D0
      bop2(4,5) = 53369.D0/21436.D0
      bop2(4,6) = -54899.D0/160770.D0
      bop2(4,7) = 48.D0/5359.D0
      bop2(4,8) = 0
      bop2(4,9) = 0
      bop2(5,1) = -7053.D0/39385.D0
      bop2(5,2) = 86551.D0/94524.D0
      bop2(5,3) = -46969.D0/23631.D0
      bop2(5,4) = 53369.D0/15754.D0
      bop2(5,5) = -87904.D0/23631.D0
      bop2(5,6) = 820271.D0/472620.D0
      bop2(5,7) = -1296.D0/7877.D0
      bop2(5,8) = 96.D0/7877.D0
      bop2(5,9) = 0
      bop2(6,1) = 21035.D0/525612.D0
      bop2(6,2) = -24641.D0/131403.D0
      bop2(6,3) = 30409.D0/87602.D0
      bop2(6,4) = -54899.D0/131403.D0
      bop2(6,5) = 820271.D0/525612.D0
      bop2(6,6) = -117600.D0/43801.D0
      bop2(6,7) = 64800.D0/43801.D0
      bop2(6,8) = -6480.D0/43801.D0
      bop2(6,9) = 480.D0/43801.D0

      s(0) = -1d0/6
      s(1) = -77d0/60
      s(2) = 5d0/2
      s(3) = -5d0/3
      s(4) = 5d0/6
      s(5) = -1d0/4
      s(6) = 1d0/30

      end

      subroutine GETBOP4S2( bop, hnorm )
      implicit none
      real*8 bop(4,6), hnorm(4)
*** First derivative bop, 2nd order
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
      hnorm(1) = 17d0/48
      hnorm(2) = 59d0/48
      hnorm(3) = 43d0/48
      hnorm(4) = 49d0/48
      end

c-----------------------------------------------------------------------
      subroutine BOPEXT4TH( bop, bope )
      implicit none
      integer i, j
      real*8 bop(4,6), bope(6,8), d4a, d4b
      do j=1,8
         do i=1,6
            bope(i,j) = 0
         enddo
      enddo
      do j=1,6
         do i=1,4
            bope(i,j) = bop(i,j)
         enddo
      enddo
      d4a = 2d0/3
      d4b = -1d0/12
      bope(5,3) = -d4b
      bope(5,4) = -d4a
      bope(5,6) =  d4a
      bope(5,7) =  d4b
      bope(6,4) = -d4b
      bope(6,5) = -d4a
      bope(6,7) =  d4a
      bope(6,8) =  d4b
      end
c-----------------------------------------------------------------------
      subroutine VARCOEFFS4( acof, ghcof )
      implicit none
      real*8 acof(6,8,8), ghcof(6)
*** acofs(i,j,k) is coefficient of a(k) in stencil coefficient (i,j)
*** ghcof is coefficient of ghost point, a(1)*ghcof*u(0) in stencil at i=1.
      ghcof(1) = 12d0/17
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
*** 129 non-zero out of 384.
      end

