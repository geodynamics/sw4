module problemsetup_new_3d
  use iso_fortran_env
  implicit none

  integer, parameter :: dp=real64
  real(dp),parameter :: pi = dacos(-1.d0)
  real(dp),parameter :: tn = 0.5d0
  integer, parameter :: nrg = 5

  ! parameters for generating meshes
  real(dp),parameter :: l1 = 2.d0*pi, l2 = 2.d0*pi, l3 = 2.d0*pi ! space interval
  real(dp),parameter :: int_pos = 0.5d0 ! This is the position in r3 where the interface is located.

  integer,parameter :: n1_c = 49, n2_c = 49! number of grid points in direction-1,2 in the coarse domain
  real(dp),parameter :: h1phy_c = l1/(n1_c-1), h1phy_f = h1phy_c*0.5d0  ! mesh size in physical space, x
  real(dp),parameter :: h2phy_c = l2/(n2_c-1), h2phy_f = h2phy_c*0.5d0  ! mesh size in physical space, y
  integer,parameter :: n3_c = ceiling(int_pos*l3/h1phy_c)+1 ! number of grid points in direction-3
  integer,parameter :: n1_f = n1_c*2-1,n2_f = n2_c*2-1,n3_f = ceiling((1.d0-int_pos)*l3/(h1phy_f))+1
  real(dp),parameter :: h1_c = 1.d0/(n1_c-1), h2_c = 1.d0/(n2_c-1), h3_c = 1.d0/(n3_c-1)
  real(dp),parameter :: h1_f = 1.d0/(n1_f-1), h2_f = 1.d0/(n2_f-1), h3_f = 1.d0/(n3_f-1)


  integer,parameter :: dim = 3
  real(dp),parameter :: amp = 0.2d0, peak = 0.04d0

contains

  subroutine generate_grid(X1_c,X2_c,X3_c,X1_f,X2_f,X3_f)
    integer i,j,k
    real(dp),dimension(1-nrg:n1_c+nrg) :: X1_c
    real(dp),dimension(1-nrg:n2_c+nrg) :: X2_c
    real(dp),dimension(1-nrg:n1_f+nrg) :: X1_f
    real(dp),dimension(1-nrg:n2_f+nrg) :: X2_f
    real(dp),dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: X3_c
    real(dp),dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: X3_f

    ! coarse domain
    do i = 1-nrg,n1_c+nrg
       X1_c(i) = dble(i-1)/(n1_c-1)*l1
    end do
    do i = 1-nrg,n2_c+nrg
       X2_c(i) = dble(i-1)/(n2_c-1)*l2
    end do
    do i = 1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,n1_c+nrg
             X3_c(k,j,i) = dble(i-1)/(n3_c-1)*interface_cf(dble(k-1)/(n1_c-1),dble(j-1)/(n2_c-1)) &
               + (1.d0-dble(i-1)/(n3_c-1))*bottom(dble(k-1)/(n1_c-1),dble(j-1)/(n2_c-1))
          end do
       end do
    end do

    ! fine domain
    do i = 1-nrg,n1_f+nrg
       X1_f(i) = dble(i-1)/(n1_f-1)*l1
    end do
    do i = 1-nrg,n2_f+nrg
       X2_f(i) = dble(i-1)/(n2_f-1)*l2
    end do
    do i = 1-nrg,n3_f+nrg
       do j = 1-nrg,n2_f+nrg
          do k = 1-nrg,n1_f+nrg
             X3_f(k,j,i) = dble(i-1)/(n3_f-1)*top(dble(k-1)/(n1_f-1),dble(j-1)/(n2_f-1)) &
               + (1.d0-dble(i-1)/(n3_f-1))*interface_cf(dble(k-1)/(n1_f-1),dble(j-1)/(n2_f-1))
          end do
       end do
    end do
    !
  end subroutine generate_grid

  subroutine metric_derivative(r1,r2,r3,xi13,xi23,xi33,J,flag)
    real(dp) :: r1,r2,r3,x3r1,x3r2,x3r3
    real(dp) :: xi13,xi23,xi33,J
    integer :: flag
    ! forward derivatives
    select case (flag)
      case (0) ! coarse
       x3r1 = r3*interface_cfx(r1,r2) + (1.d0-r3)*bottomx(r1,r2)
       x3r2 = r3*interface_cfy(r1,r2) + (1.d0-r3)*bottomy(r1,r2)
       x3r3 = interface_cf(r1,r2) - bottom(r1,r2)
      case (1) ! fine
       x3r1 = r3*topx(r1,r2) + (1.d0-r3)*interface_cfx(r1,r2)
       x3r2 = r3*topy(r1,r2) + (1.d0-r3)*interface_cfy(r1,r2)
       x3r3 = top(r1,r2) - interface_cf(r1,r2)
    end select
       J = l1*l2*x3r3
      ! backward derivative
       xi13 = -l2*x3r1/J
       xi23 = -l1*x3r2/J
       xi33 = l1*l2/J
  !
  end subroutine metric_derivative

  ! describing the top geometry
  real(dp) function top(x,y)
    real(dp) x,y
      top = l3+amp*exp(-(x-0.5d0)**2/peak)+amp*exp(-(y-0.5d0)**2/peak)
  end function top

  real(dp) function topx(x,y)
    real(dp) x,y
    topx = -amp*exp(-(x-0.5d0)**2/peak)*(2.d0*x-1.d0)/peak
  end function topx

  real(dp) function topy(x,y)
    real(dp) x,y
    topy = -amp*exp(-(y-0.5d0)**2/peak)*(2.d0*y-1.d0)/peak
  end function topy

  ! describe the bottom geometry
  real(dp) function bottom(x,y)
    real(dp) x,y
    bottom = amp*exp(-(x-0.6d0)**2/peak)+amp*exp(-(y-0.6d0)**2/peak)
  end function bottom

  real(dp) function bottomx(x,y)
    real(dp) x,y
    bottomx = -amp*exp(-(x-0.6d0)**2/peak)*(2.d0*x-1.2d0)/peak
  end function bottomx

  real(dp) function bottomy(x,y)
    real(dp) x,y
    bottomy = -amp*exp(-(y-0.6d0)**2/peak)*(2.d0*y-1.2d0)/peak
  end function bottomy

  ! describe the interface of geometry
  real(dp) function interface_cf(x,y)
    real(dp) x,y
    interface_cf = int_pos*l3+amp*sin(4.d0*pi*x)+amp*cos(4.d0*pi*y)
  end function interface_cf

  real(dp) function interface_cfx(x,y)
    real(dp) x,y
    interface_cfx = amp*4.d0*pi*cos(4.d0*pi*x)
  end function interface_cfx

  real(dp) function interface_cfy(x,y)
    real(dp) x,y
    interface_cfy = -amp*4.d0*pi*sin(4.d0*pi*y)
  end function interface_cfy

  ! exact solution
  ! We want the time-term to be the same in u1,u2 and u3, so that we don't need to call this function
  ! in every time step. We just need to scale the initial condition.
  subroutine exact_solution(x1,x2,x3,t,u1,u2,u3,flag)
    implicit none
    real(dp) :: x1, x2, x3, t, u1, u2, u3
    integer flag
    u1 = cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    u2 = sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    u3 = sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
  end subroutine exact_solution

  subroutine forcing(x1,x2,x3,t,f1,f2,f3)
    real(dp) :: t, x1, x2, x3, f1, f2, f3
    real(dp) :: d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3
    real(dp) :: d1d1u1, d2d2u1, d3d3u1, d1d1u2, d2d2u2, d3d3u2, d1d1u3, d2d2u3, d3d3u3
    real(dp) :: d1d2u1, d1d3u1, d2d3u1, d1d2u2, d1d3u2, d2d3u2, d1d2u3, d1d3u3, d2d3u3
    real(dp) :: d2d1u1, d3d1u1, d3d2u1, d2d1u2, d3d1u2, d3d2u2, d2d1u3, d3d1u3, d3d2u3
    real(dp) :: mu, lambda, rho, d1mu, d2mu, d3mu, d1lambda, d2lambda, d3lambda

    ! material property
    mu = 3.d0 + 1.d0*sin(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*sin(x3)
    lambda = 21.d0+ 1.d0*cos(x1+0.1d0)*cos(x2+0.1d0)*(sin(3.d0*x3)**2)
    rho = 2.d0 + 1.d0*sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3-0.2d0)

    ! space derivatives
    d1u1 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2u1 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3u1 = cos(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d1u2 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2u2 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3u2 = sin(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d1u3 = cos(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d2u3 = sin(x1+0.2d0)*cos(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d3u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*sin(x3+0.2d0)*sin(t)

    d1d1u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2d2u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3d3u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d1d1u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2d2u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3d3u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d1d1u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d2d2u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d3d3u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)

    d1d2u1 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2d1u1 = d1d2u1
    d1d3u1 = -sin(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d3d1u1 = d1d3u1
    d2d3u1 = cos(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d3d2u1 = d2d3u1
    d1d2u2 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2d1u2 = d1d2u2
    d1d3u2 = cos(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d3d1u2 = d1d3u2
    d2d3u2 = -sin(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d3d2u2 = d2d3u2
    d1d2u3 = cos(x1+0.2d0)*cos(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d2d1u3 = d1d2u3
    d1d3u3 = -cos(x1+0.2d0)*sin(x2+0.2d0)*sin(x3+0.2d0)*sin(t)
    d3d1u3 = d1d3u3
    d2d3u3 = -sin(x1+0.2d0)*cos(x2+0.2d0)*sin(x3+0.2d0)*sin(t)
    d3d2u3 = d2d3u3

    d1mu = 3.d0*cos(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*sin(x3)
    d2mu = 3.d0*sin(3.d0*x1+0.1d0)*cos(3.d0*x2+0.1d0)*sin(x3)
    d3mu = sin(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*cos(x3)
    d1lambda = -sin(x1+0.1d0)*cos(x2+0.1d0)*(sin(3.d0*x3)**2)
    d2lambda = -cos(x1+0.1d0)*sin(x2+0.1d0)*(sin(3.d0*x3)**2)
    d3lambda = cos(x1+0.1d0)*cos(x2+0.1d0)*2.d0*sin(3.d0*x3)*cos(3.d0*x3)*3.d0

    f1 = rho*(-cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(cos(t*t)*4.d0*t*t+sin(t*t)*2.d0)) &
         - (2.d0*d1mu*d1u1+d1lambda*d1u1+(2.d0*mu+lambda)*d1d1u1+d1lambda*d2u2+lambda*d1d2u2+d1lambda*d3u3 &
            +lambda*d1d3u3+d2mu*d1u2+mu*d2d1u2+d2mu*d2u1+mu*d2d2u1+d3mu*d1u3+mu*d3d1u3+d3mu*d3u1+mu*d3d3u1)

    f2 = rho*(-sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(cos(t*t)*4.d0*t*t+sin(t*t)*2.d0)) &
         - (d1mu*d1u2+mu*d1d1u2+d1mu*d2u1+mu*d1d2u1+d2lambda*d1u1+lambda*d2d1u1+2.d0*d2mu*d2u2+d2lambda*d2u2 &
            +(2.d0*mu+lambda)*d2d2u2+d2lambda*d3u3+lambda*d2d3u3+d3mu*d2u3+mu*d3d2u3+d3mu*d3u2+mu*d3d3u2)

    f3 = rho*(-sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)) &
         - (d1mu*d1u3+mu*d1d1u3+d1mu*d3u1+mu*d1d3u1+d2mu*d2u3+mu*d2d2u3+d2mu*d3u2+mu*d2d3u2+d3lambda*d1u1 &
            +lambda*d3d1u1+d3lambda*d2u2+lambda*d3d2u2+2.d0*d3mu*d3u3+d3lambda*d3u3+(2.d0*mu+lambda)*d3d3u3)

  end subroutine forcing

  subroutine forcing_tt(x1,x2,x3,t,f1tt,f2tt,f3tt)
    real(dp) :: x1, x2, x3, t, f1tt, f2tt, f3tt
    real(dp) :: d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3
    real(dp) :: d1d1u1, d2d2u1, d3d3u1, d1d1u2, d2d2u2, d3d3u2, d1d1u3, d2d2u3, d3d3u3
    real(dp) :: d1d2u1, d1d3u1, d2d3u1, d1d2u2, d1d3u2, d2d3u2, d1d2u3, d1d3u3, d2d3u3
    real(dp) :: d2d1u1, d3d1u1, d3d2u1, d2d1u2, d3d1u2, d3d2u2, d2d1u3, d3d1u3, d3d2u3
    real(dp) :: mu, lambda, rho, d1mu, d2mu, d3mu, d1lambda, d2lambda, d3lambda

    ! material property
    mu = 3.d0 + 1.d0*sin(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*sin(x3)
    lambda = 21.d0+ 1.d0*cos(x1+0.1d0)*cos(x2+0.1d0)*(sin(3.d0*x3)**2)
    rho = 2.d0 + 1.d0*sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3-0.2d0)

    ! second time derivatives to space derivatives
    d1u1 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2u1 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3u1 = cos(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d1u2 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2u2 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3u2 = sin(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d1u3 = cos(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))
    d2u3 = sin(x1+0.2d0)*cos(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))
    d3u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*sin(x3+0.2d0)*(-sin(t))

    d1d1u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2d2u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d3u1 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d1d1u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2d2u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d3u2 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d1d1u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))
    d2d2u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))
    d3d3u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))

    d1d2u1 = -sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2d1u1 = d1d2u1
    d1d3u1 = -sin(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d1u1 = d1d3u1
    d2d3u1 = cos(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d2u1 = d2d3u1
    d1d2u2 = -cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d2d1u2 = d1d2u2
    d1d3u2 = cos(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d1u2 = d1d3u2
    d2d3u2 = -sin(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*(-cos(t*t)*4.d0*t*t-sin(t*t)*2.d0)
    d3d2u2 = d2d3u2
    d1d2u3 = cos(x1+0.2d0)*cos(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))
    d2d1u3 = d1d2u3
    d1d3u3 = -cos(x1+0.2d0)*sin(x2+0.2d0)*sin(x3+0.2d0)*(-sin(t))
    d3d1u3 = d1d3u3
    d2d3u3 = -sin(x1+0.2d0)*cos(x2+0.2d0)*sin(x3+0.2d0)*(-sin(t))
    d3d2u3 = d2d3u3

    d1mu = 3.d0*cos(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*sin(x3)
    d2mu = 3.d0*sin(3.d0*x1+0.1d0)*cos(3.d0*x2+0.1d0)*sin(x3)
    d3mu = sin(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*cos(x3)
    d1lambda = -sin(x1+0.1d0)*cos(x2+0.1d0)*(sin(3.d0*x3)**2)
    d2lambda = -cos(x1+0.1d0)*sin(x2+0.1d0)*(sin(3.d0*x3)**2)
    d3lambda = cos(x1+0.1d0)*cos(x2+0.1d0)*2.d0*sin(3.d0*x3)*cos(3.d0*x3)*3.d0

    f1tt = rho*(-cos(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*(cos(t*t)*(12.d0-16.d0*(t**4))-sin(t*t)*48.d0*t*t)) &
         - (2.d0*d1mu*d1u1+d1lambda*d1u1+(2.d0*mu+lambda)*d1d1u1+d1lambda*d2u2+lambda*d1d2u2+d1lambda*d3u3 &
            +lambda*d1d3u3+d2mu*d1u2+mu*d2d1u2+d2mu*d2u1+mu*d2d2u1+d3mu*d1u3+mu*d3d1u3+d3mu*d3u1+mu*d3d3u1)

    f2tt = rho*(-sin(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*(cos(t*t)*(12.d0-16.d0*(t**4))-sin(t*t)*48.d0*t*t)) &
         - (d1mu*d1u2+mu*d1d1u2+d1mu*d2u1+mu*d1d2u1+d2lambda*d1u1+lambda*d2d1u1+2.d0*d2mu*d2u2+d2lambda*d2u2 &
            +(2.d0*mu+lambda)*d2d2u2+d2lambda*d3u3+lambda*d2d3u3+d3mu*d2u3+mu*d3d2u3+d3mu*d3u2+mu*d3d3u2)

    f3tt = rho*(-sin(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*(-sin(t))) &
         - (d1mu*d1u3+mu*d1d1u3+d1mu*d3u1+mu*d1d3u1+d2mu*d2u3+mu*d2d2u3+d2mu*d3u2+mu*d2d3u2+d3lambda*d1u1 &
            +lambda*d3d1u1+d3lambda*d2u2+lambda*d3d2u2+2.d0*d3mu*d3u3+d3lambda*d3u3+(2.d0*mu+lambda)*d3d3u3)

  end subroutine forcing_tt

  subroutine top_normal_data(x1,x2,x3,l1,l2,t,g)
    real(dp) :: l1, l2, x1, x2, x3, t
    real(dp), dimension(1:dim,1:dim) :: St
    real(dp), dimension(1:dim) :: n,g
    real(dp) :: d1u1, d2u1, d3u1, d1u2, d2u2, d3u2, d1u3, d2u3, d3u3, mu, lambda

    ! material property
    mu = 3.d0 + 1.d0*sin(3.d0*x1+0.1d0)*sin(3.d0*x2+0.1d0)*sin(x3)
    lambda = 21.d0+ 1.d0*cos(x1+0.1d0)*cos(x2+0.1d0)*(sin(3.d0*x3)**2)

    ! space derivatives
    d1u1 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2u1 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3u1 = cos(x1+0.3d0)*sin(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d1u2 = cos(x1+0.3d0)*cos(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d2u2 = -sin(x1+0.3d0)*sin(x2+0.3d0)*sin(x3+0.2d0)*cos(t*t)
    d3u2 = sin(x1+0.3d0)*cos(x2+0.3d0)*cos(x3+0.2d0)*cos(t*t)
    d1u3 = cos(x1+0.2d0)*sin(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d2u3 = sin(x1+0.2d0)*cos(x2+0.2d0)*cos(x3+0.2d0)*sin(t)
    d3u3 = -sin(x1+0.2d0)*sin(x2+0.2d0)*sin(x3+0.2d0)*sin(t)

    St(1,1) = (2.d0*mu+lambda)*d1u1 + lambda*d2u2 + lambda*d3u3

    St(1,2) = mu*d1u2 + mu*d2u1

    St(1,3) = mu*d1u3 + mu*d3u1

    St(2,1) = mu*d1u2 + mu*d2u1

    St(2,2) = lambda*d1u1 + (2.d0*mu+lambda)*d2u2 + lambda*d3u3

    St(2,3) = mu*d2u3 + mu*d3u2

    St(3,1) = mu*d1u3 + mu*d3u1

    St(3,2) = mu*d2u3 + mu*d3u2

    St(3,3) = lambda*d1u1 + lambda*d2u2 + (2.d0*mu+lambda)*d3u3

    n(1) = amp*exp(-(x1/l1-0.5d0)**2/peak)*(2.d0*x1/l1-1.d0)/peak/l1
    n(2) = amp*exp(-(x2/l2-0.5d0)**2/peak)*(2.d0*x2/l2-1.d0)/peak/l2
    n(3) = 1.d0
    n = n/sqrt(n(1)**2+n(2)**2+n(3)**2)

    g(1) = St(1,1)*n(1)+St(1,2)*n(2)+St(1,3)*n(3)
    g(2) = St(2,1)*n(1)+St(2,2)*n(2)+St(2,3)*n(3)
    g(3) = St(3,1)*n(1)+St(3,2)*n(2)+St(3,3)*n(3)

  end subroutine top_normal_data

end module problemsetup_new_3d
