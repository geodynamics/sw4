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
      subroutine VELSUM( is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2,
     *                   mu, lambda, rho, cp, cs, npts ) bind(c)
      implicit none
      integer is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2
      integer i, j, k
      real*8 mu(is:ie,js:je,ks:ke), lambda(is:ie,js:je,ks:ke)
      real*8 rho(is:ie,js:je,ks:ke), cp, cs, npts
      cp = 0
      cs = 0
      npts = DBLE(i2-i1+1)*DBLE(j2-j1+1)*DBLE(k2-k1+1)
!$OMP PARALLEL PRIVATE( i, j, k )
!$OMP DO REDUCTION(+:cp,cs)      
      do k=k1,k2
         do j=j1,j2
            do i=i1,i2
               cp = cp + SQRT( (2*mu(i,j,k)+lambda(i,j,k))/rho(i,j,k) )
               cs = cs + SQRT( mu(i,j,k)/rho(i,j,k) )
            enddo
         enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      end
