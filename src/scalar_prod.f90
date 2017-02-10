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
subroutine scalar_prod( is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2, onesided, &
     u, w, strx, stry, strz, sc_prod ) bind(c)
  implicit none
  integer, value:: is, ie, js, je, ks, ke, i1, i2, j1, j2, k1, k2
  integer onesided(6)
  real*8 u(3,is:ie,js:je,ks:ke), w(3,is:ie,js:je,ks:ke)
  real*8 strx(is:ie), stry(js:je), strz(ks:ke)
  real*8, intent(out):: sc_prod

! local variables
  real*8 normwgh(4), term, normfact
  integer i, j, k
  normwgh(1) = 17d0/48
  normwgh(2) = 59d0/48
  normwgh(3) = 43d0/48
  normwgh(4) = 49d0/48
! initialize  
  sc_prod = 0
  do k=k1,k2
     do j=j1,j2
        do i=i1,i2
! NOTE: the scalar product is scaled by the stretching
           term =(u(1,i,j,k)*w(1,i,j,k) + u(2,i,j,k)*w(2,i,j,k) + u(3,i,j,k)*w(3,i,j,k))/(strx(i)*stry(j)*strz(k))
           normfact = 1               
           if( k.le.4 .and. onesided(5).eq.1 )then
              normfact = normwgh(k)
           endif
           if( k.ge.k2-3 .and. onesided(6).eq.1 )then
              normfact = normwgh(k2-k+1)
           endif
           sc_prod = sc_prod + normfact*term
        enddo
     enddo
  enddo
end subroutine scalar_prod
