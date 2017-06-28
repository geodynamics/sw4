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
      subroutine TESTSRC( f, if, il, jf, jl, kf, kl, nk, wind, zmin, h,
     *                    kx, ky, kz, mom ) bind(c)
      implicit none
      integer if, il, jf, jl, kf, kl, nk, kx(3), ky(3), kz(3), wind(6)
      integer i, j, k
      real*8 f(3,if:il,jf:jl,kf:kl), zmin, h, mom(3), x, y, z
      real*8 h3, normfact, wgh(4), x1, x2, x3, y1, y2, y3, z1, z2, z3
      wgh(1) = 17d0/48
      wgh(2) = 59d0/48
      wgh(3) = 43d0/48
      wgh(4) = 49d0/48
      h3 = h*h*h
!$OMP PARALLEL PRIVATE(k,j,i,z,y,x,normfact,x1,x2,x3,y1,y2,y3,z1,z2,z3)
!$OMP DO REDUCTION(+:mom)
      do k=wind(5),wind(6)
         z = zmin + (k-1)*h
         do j=wind(3),wind(4)
            y = (j-1)*h
            do i=wind(1),wind(2)
               normfact = h3
               x = (i-1)*h
               if( k.le.4 )then
                  normfact = normfact*wgh(k)
               endif
               if( k.ge.nk-3 .and. k.le.nk )then
                  normfact = normfact*wgh(nk-k+1)
               endif
               if( kx(1).eq.0 )then
                  x1 = 1
               else
                  x1 = x**kx(1)
               endif
               if( kx(2).eq.0 )then
                  x2 = 1
               else
                  x2 = x**kx(2)
               endif
               if( kx(3).eq.0 )then
                  x3 = 1
               else
                  x3 = x**kx(3)
               endif
               if( ky(1).eq.0 )then
                  y1 = 1
               else
                  y1 = y**ky(1)
               endif
               if( ky(2).eq.0 )then
                  y2 = 1
               else
                  y2 = y**ky(2)
               endif
               if( ky(3).eq.0 )then
                  y3 = 1
               else
                  y3 = y**ky(3)
               endif
               if( kz(1).eq.0 )then
                  z1 = 1
               else
                  z1 = z**kz(1)
               endif
               if( kz(2).eq.0 )then
                  z2 = 1
               else
                  z2 = z**kz(2)
               endif
               if( kz(3).eq.0 )then
                  z3 = 1
               else
                  z3 = z**kz(3)
               endif
               mom(1) = mom(1) + x1*y1*z1*f(1,i,j,k)*normfact
               mom(2) = mom(2) + x2*y2*z2*f(2,i,j,k)*normfact
               mom(3) = mom(3) + x3*y3*z3*f(3,i,j,k)*normfact
            enddo
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end
