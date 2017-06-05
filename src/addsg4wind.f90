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
!-----------------------------------------------------------------------
! Adds 4th order artificial disssipation for super-grid damping layers
!
! Windowed version used with mesh refinement
!
!-----------------------------------------------------------------------
subroutine addsg4wind( dt, h, up, u, um, rho, dcx, dcy, dcz, strx, stry, strz, cox, coy, coz, &
     ifirst, ilast, jfirst, jlast, kfirst, klast, beta, kupb, kupe, kwindb, kwinde ) bind(c)

!***********************************************************************
!*** Version with correct density scaling and supergrid stretching.
!*** cox, coy, coz are corner factors that reduce the damping near edges and corners
!***
!***********************************************************************

  implicit none
  real*8 dt, h
  real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
! up can have different size in k-index
  real*8 up(3,ifirst:ilast,jfirst:jlast,kupb:kupe)
  real*8  rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 dcx(ifirst:ilast), strx(ifirst:ilast), cox(ifirst:ilast)
  real*8 dcy(jfirst:jlast), stry(jfirst:jlast), coy(jfirst:jlast)
  real*8 dcz(kfirst:klast), strz(kfirst:klast), coz(kfirst:klast)
  integer, value:: ifirst, ilast, jfirst, jlast, kfirst, klast, kupb, kupe, kwindb, kwinde
  real*8,value:: beta

! time stepping stability condition on beta?

! this routine uses un-divided differences in x and t
  real*8 coeff, irho

  integer i, j, k, c;
  
  if( beta .eq. 0d0 ) return;

  coeff = beta
! beta is the supergrid damping coefficient as entered in the input file
!
! add in the SG damping
!
  ! There are enough ghost points to always use the interior formula
  ! But only in (i,j) because the k-window may be near a refinement bndry
!
! the corner tapering is applied by replacing
! strx -> strx*coy(j)*coz(k)
! stry -> stry*cox(i)*coz(k)
! strz -> strz*cox(i)*coy(j)
!
!  do k=kfirst+2,klast-2
  do k=kwindb,kwinde
     do j=jfirst+2,jlast-2
        do i=ifirst+2, ilast-2
           irho = 1/rho(i,j,k)
           do c=1,3
              up(c,i,j,k) = up(c,i,j,k) - irho*coeff*( &
! x-differences
                   strx(i)*coy(j)*coz(k)*( rho(i+1,j,k)*dcx(i+1)*( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k)) &
                   -2*rho(i,j,k)*dcx(i)  * ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k)) &
                   +rho(i-1,j,k)*dcx(i-1)*( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) &
                   -rho(i+1,j,k)*dcx(i+1)*(um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) &
                   +2*rho(i,j,k)*dcx(i)*(um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k))&
                   -rho(i-1,j,k)*dcx(i-1)* (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) &
! y-differences
                   +stry(j)*cox(i)*coz(k)*( rho(i,j+1,k)*dcy(j+1)* ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) &
                   -2*rho(i,j,k)*dcy(j)  * ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k)) &
                   +rho(i,j-1,k)*dcy(j-1)* ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) &
                   -rho(i,j+1,k)*dcy(j+1)* (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) &
                   +2*rho(i,j,k)*dcy(j)  * (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) &
                   -rho(i,j-1,k)*dcy(j-1) * (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) )  &
! z-differences (not used at a refinement interface)
                   ! + strz(k)*cox(i)*coy(j)*( rho(i,j,k+1)*dcz(k+1)* ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k) ) & 
                   ! - 2*rho(i,j,k)*dcz(k)  * ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1)) +rho(i,j,k-1)*dcz(k-1)* &
                   ! ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) -rho(i,j,k+1)*dcz(k+1) * &
                   ! (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) +2*rho(i,j,k)*dcz(k)  * &
                   ! (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) -rho(i,j,k-1)*dcz(k-1)* &
                   ! (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) &
                   )
           enddo ! do c
        enddo ! do i
     enddo ! do j
  enddo ! do k
end subroutine addsg4wind
