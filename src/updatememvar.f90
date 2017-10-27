!-----------------------------------------------------------------------
subroutine memvar_pred_fort( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, alm, u, omega, dt, domain ) bind(c)

  !***********************************************************************
  !*** 
  !*** domain = 0 --> Entire domain
  !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  !***
  !***********************************************************************      

  ! AP: this routine implementes the 2nd order predictor step for evolving the memory variables
  implicit none

  integer, value:: ifirst, ilast, jfirst, jlast, kfirst, klast, domain

  real*8 ::alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 ::alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 :: u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8, value:: omega, dt

! local variables
  integer :: i, j, k, c, k1, k2
  real*8 ::dto, icp, cm

  dto = dt*omega

  icp = 1/( 1d0/2 + 1/(2*dto) )
  cm = 1d0/2 - 1/(2*dto)

  if( domain.eq.0 )then
     k1 = kfirst
     k2 = klast
  elseif( domain.eq.1 )then
     k1 = kfirst
     k2 = kfirst+2
  elseif( domain.eq.2 )then
     k1 = klast-2
     k2 = klast
  endif
  do k=k1,k2
     do j=jfirst,jlast
        do i=ifirst,ilast
           do c=1,3
              alp(c,i,j,k) = icp*(-cm*alm(c,i,j,k) + u(c,i,j,k) )
           enddo
        enddo
     enddo
  enddo
end subroutine memvar_pred_fort

!-----------------------------------------------------------------------
subroutine memvar_corr_fort( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, alm, up, u, um, omega, dt, domain ) bind(c)

  !***********************************************************************
  !*** 
  !*** domain = 0 --> Entire domain
  !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  !***
  !***********************************************************************      

  ! AP Apr. 3, 2017: corrector step for updating memory variables
  ! AP June 14, 2017: make corrector step independent of predictor step to simplify
  ! the mesh refinement algorithm
  implicit none
  real*8 i6
  parameter( i6=1d0/6 )
  integer, value:: ifirst, ilast, jfirst, jlast, kfirst, klast, domain
  integer i, j, k, c
  integer k1, k2

  real*8, value:: omega, dt
  real*8 dto, icp, cm, co
  real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast) 
  real*8 alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

  dto = dt*omega

  icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
  cm = 1/(2*dto) + dto/4  - 1d0/2 - dto*dto/12 

  if( domain.eq.0 )then
     k1 = kfirst
     k2 = klast
  elseif( domain.eq.1 )then
     k1 = kfirst
     k2 = kfirst+2
  elseif( domain.eq.2 )then
     k1 = klast-2
     k2 = klast
  endif
  do k=k1,k2
     do j=jfirst,jlast
        do i=ifirst,ilast
           do c=1,3
              ! Note that alp is ASSIGNED by this formula
              alp(c,i,j,k) = icp*( cm*alm(c,i,j,k) + u(c,i,j,k) + i6* ( dto*dto*u(c,i,j,k) + &
                   dto*(up(c,i,j,k)-um(c,i,j,k)) + (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k)) ) )
           enddo
        enddo
     enddo
  enddo
end subroutine memvar_corr_fort

!-----------------------------------------------------------------------
subroutine memvar_corr_fort_wind( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, d1b, d1e, d2b, d2e, d3b, d3e, &
     alm, up, u, um, omega, dt, domain ) bind(c)

  !***********************************************************************
  !*** 
  !*** domain = 0 --> Entire domain
  !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  !***
  !***********************************************************************      

  ! AP Apr. 3, 2017: corrector step for updating memory variables
  ! AP June 14, 2017: make corrector step independent of predictor step to simplify
  ! the mesh refinement algorithm
  implicit none
  real*8 i6
  parameter( i6=1d0/6 )
  integer, value:: ifirst, ilast, jfirst, jlast, kfirst, klast, domain
  integer, value:: d1b, d1e, d2b, d2e, d3b, d3e;
  integer i, j, k, c
  integer k1, k2

  real*8, value:: omega, dt
  real*8 dto, icp, cm, co
  real*8 up(3,d1b:d1e, d2b:d2e, d3b:d3e)
  real*8  u(3,d1b:d1e, d2b:d2e, d3b:d3e);
  real*8 um(3,d1b:d1e, d2b:d2e, d3b:d3e);
  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast) ! different sizes here
  real*8 alm(3,d1b:d1e, d2b:d2e, d3b:d3e)

  dto = dt*omega

  icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
  cm = 1/(2*dto) + dto/4  - 1d0/2 - dto*dto/12 

  if( domain.eq.0 )then
     k1 = kfirst
     k2 = klast
  elseif( domain.eq.1 )then
     k1 = kfirst
     k2 = kfirst+2
  elseif( domain.eq.2 )then
     k1 = klast-2
     k2 = klast
  endif
  do k=k1,k2
     do j=jfirst,jlast
        do i=ifirst,ilast
           do c=1,3
              ! Note that alp is ASSIGNED by this formula
              alp(c,i,j,k) = icp*( cm*alm(c,i,j,k) + u(c,i,j,k) + i6* ( dto*dto*u(c,i,j,k) + &
                   dto*(up(c,i,j,k)-um(c,i,j,k)) + (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k)) ) )
           enddo
        enddo
     enddo
  enddo
end subroutine memvar_corr_fort_wind

! !-----------------------------------------------------------------------
! subroutine memvar_corr_fort( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, alm, up, u, um, omega, dt, domain ) bind(c)

!   !***********************************************************************
!   !*** 
!   !*** domain = 0 --> Entire domain
!   !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
!   !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
!   !***
!   !***********************************************************************      

!   ! AP Apr. 3, 2017: corrector step for updating memory variables
!   ! This routine must be called with the predictor value in 'up'
!   implicit none
!   real*8 i6
!   parameter( i6=1d0/6 )
!   integer, value:: ifirst, ilast, jfirst, jlast, kfirst, klast, domain
!   integer i, j, k, c
!   integer k1, k2

!   real*8, value:: omega, dt
!   real*8 dto, icp, cm, co
!   real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!   real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!   real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!   real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!   real*8 alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

!   dto = dt*omega

!   icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
!   cm =  dto/4 - dto*dto/12
!   co  = 1d0/2 + 1/(2*dto)

!   if( domain.eq.0 )then
!      k1 = kfirst
!      k2 = klast
!   elseif( domain.eq.1 )then
!      k1 = kfirst
!      k2 = kfirst+2
!   elseif( domain.eq.2 )then
!      k1 = klast-2
!      k2 = klast
!   endif
!   do k=k1,k2
!      do j=jfirst,jlast
!         do i=ifirst,ilast
!            do c=1,3
!               ! Note that alp is over-written by this formula
!               alp(c,i,j,k) = icp*( co*alp(c,i,j,k) + cm*alm(c,i,j,k) + i6* ( dto*dto*u(c,i,j,k) + &
!                    dto*(up(c,i,j,k)-um(c,i,j,k)) + (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k)) ) )
!            enddo
!         enddo
!      enddo
!   enddo
! end subroutine memvar_corr_fort

!-----------------------------------------------------------------------
!!$subroutine updatememvar( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, alm, up, u, um, omega, dt, domain, pred ) bind(c)
!!$
!!$  !***********************************************************************
!!$  !*** 
!!$  !*** domain = 0 --> Entire domain
!!$  !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
!!$  !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
!!$  !***
!!$  !***********************************************************************      
!!$
!!$  implicit none
!!$  real*8 i6
!!$  parameter( i6=1d0/6 )
!!$  integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
!!$  integer domain, k1, k2
!!$  !c AP Nov 14, 2016: experimenting with PC formulation for memory variables
!!$  integer pred
!!$  real*8 pcoeff
!!$
!!$  real*8 omega, dt, dto, icp, cm, cof3
!!$  real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!!$  real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!!$  real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!!$  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!!$  real*8 alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
!!$
!!$  dto = dt*omega
!!$  if (pred==1) then
!!$     icp = 1/( 1d0/2 + 1/(2*dto) )
!!$     cm = 1d0/2 - 1/(2*dto)
!!$     pcoeff = 0
!!$  else
!!$     icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
!!$     cm  =     1d0/2 - 1/(2*dto) - dto/4 + dto*dto/12
!!$     pcoeff = 1
!!$  endif
!!$  if( domain.eq.0 )then
!!$     k1 = kfirst
!!$     k2 = klast
!!$  elseif( domain.eq.1 )then
!!$     k1 = kfirst
!!$     k2 = kfirst+2
!!$  elseif( domain.eq.2 )then
!!$     k1 = klast-2
!!$     k2 = klast
!!$  endif
!!$  do k=k1,k2
!!$     do j=jfirst,jlast
!!$        do i=ifirst,ilast
!!$           do c=1,3
!!$              alp(c,i,j,k) = icp*(-cm*alm(c,i,j,k) + u(c,i,j,k) + pcoeff*i6* ( dto*dto*u(c,i,j,k) + &
!!$                   dto*(up(c,i,j,k)-um(c,i,j,k)) + (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k)) ) )
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo
!!$end subroutine updatememvar
