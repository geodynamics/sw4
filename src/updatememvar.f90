!-----------------------------------------------------------------------
subroutine updatememvar( ifirst, ilast, jfirst, jlast, kfirst, klast, alp, alm, up, u, um, omega, dt, domain, pred ) bind(c)

  !***********************************************************************
  !*** 
  !*** domain = 0 --> Entire domain
  !*** domain = 1 --> Only upper (k=1) boundary ghost points + boundary point
  !*** domain = 2 --> Only lower (k=N) boundary ghost points + boundary point
  !***
  !***********************************************************************      

  implicit none
  real*8 i6
  parameter( i6=1d0/6 )
  integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
  integer domain, k1, k2
  !c AP Nov 14, 2016: experimenting with PC formulation for memory variables
  integer pred
  real*8 pcoeff

  real*8 omega, dt, dto, icp, cm, cof3
  real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 alp(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
  real*8 alm(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

  dto = dt*omega
  if (pred==1) then
     icp = 1/( 1d0/2 + 1/(2*dto) )
     cm = 1d0/2 - 1/(2*dto)
     pcoeff = 0
  else
     icp = 1/( 1d0/2 + 1/(2*dto) + dto/4 + dto*dto/12 )
     cm  =     1d0/2 - 1/(2*dto) - dto/4 + dto*dto/12
     pcoeff = 1
  endif
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
              alp(c,i,j,k) = icp*(-cm*alm(c,i,j,k) + u(c,i,j,k) + pcoeff*i6* ( dto*dto*u(c,i,j,k) + &
                   dto*(up(c,i,j,k)-um(c,i,j,k)) + (up(c,i,j,k)-2*u(c,i,j,k)+um(c,i,j,k)) ) )
           enddo
        enddo
     enddo
  enddo
end subroutine updatememvar
