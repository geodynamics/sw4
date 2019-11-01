module problemsetup_new_3d
  use iso_fortran_env
  implicit none

  integer, parameter :: dp = real64
  real(dp),parameter :: pi = dacos(-1.d0)
  real(dp),parameter :: tn = 5.d0
  integer, parameter :: nrg = 5

  ! parameters for generating meshes
  real(dp),parameter :: l1 = 2.d0*pi, l2 = 2.d0*pi, l3 = 2.d0*pi ! space interval
  real(dp),parameter :: int_pos = 0.5d0 ! This is the position in r3 where the interface is located.
  integer,parameter :: n1_c = 25, n2_c = 25 ! number of grid points in direction-1 in the coarse domain
  real(dp),parameter :: h1phy_c = l1/(n1_c-1), h1phy_f = h1phy_c*0.5d0  ! mesh size in physical space, x
  real(dp),parameter :: h2phy_c = l2/(n2_c-1), h2phy_f = h2phy_c*0.5d0  ! mesh size in physical space, y
  integer,parameter :: n3_c = ceiling(int_pos*l3/h1phy_c)+1 ! number of grid points in direction-3
  integer,parameter :: n1_f = n1_c*2-1,n2_f = n1_c*2-1,n3_f = ceiling((1.d0-int_pos)*l3/(h1phy_f))+1
  real(dp),parameter :: h1_c = 1.d0/(n1_c-1), h2_c = 1.d0/(n2_c-1), h3_c = 1.d0/(n3_c-1)
  real(dp),parameter :: h1_f = 1.d0/(n1_f-1), h2_f = 1.d0/(n2_f-1), h3_f = 1.d0/(n3_f-1)


  integer,parameter :: dim = 3
  real(dp),parameter :: amp = 0.2d0, peak = 0.04d0

contains

  subroutine generate_grid(R1_c,R2_c,R3_c,X1_c,X2_c,X3_c,R1_f,R2_f,R3_f,X1_f,X2_f,X3_f)
    integer i,j,k
    integer, parameter :: m = 6 ! controls the smoothness of the mapping
    real(dp),dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)::R1_c,R2_c,R3_c,X1_c,X2_c,X3_c,Tau_c,Tau_ci
    real(dp),dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)::R1_f,R2_f,R3_f,X1_f,X2_f,X3_f,Tau_f,Tau_fi

    ! coarse domain, mesh in parameter space
    do i = 1-nrg,n1_c+nrg
       R1_c(i,:,:) = dble(i-1)/(n1_c-1)
    end do

    do i = 1-nrg,n2_c+nrg
       R2_c(:,i,:) = dble(i-1)/(n2_c-1)
    end do

    do i = 1-nrg,n3_c+nrg
       R3_c(:,:,i) = dble(i-1)/(n3_c-1)
    end do

    ! bottom geometry
    do i = 1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,n1_c+nrg
             Tau_c(k,j,i) = bottom(R1_c(k,j,i),R2_c(k,j,i))
          end do
       end do
    end do

    ! interface geometry
    do i = 1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,n1_c+nrg
             Tau_ci(k,j,i) = interface_cf(R1_c(k,j,i),R2_c(k,j,i))
          end do
       end do
    end do

    ! scale the mesh
    X1_c = R1_c*l1
    X2_c = R2_c*l2

    ! X3
    X3_c(:,:,:) = R3_c*Tau_ci(:,:,:) + (1.d0-R3_c)*Tau_c(:,:,:)

    ! fine domain
    do i = 1-nrg,n1_f+nrg
       R1_f(i,:,:) = dble(i-1)/(n1_f-1)
    end do

    do i = 1-nrg,n2_f+nrg
       R2_f(:,i,:) = dble(i-1)/(n2_f-1)
    end do

    do i = 1-nrg,n3_f+nrg
       R3_f(:,:,i) = dble(i-1)/(n3_f-1)
    end do

    ! top geometry
    do i = 1-nrg,n3_f+nrg
       do j = 1-nrg,n2_f+nrg
          do k = 1-nrg,n1_f+nrg
             Tau_f(k,j,i) = top(R1_f(k,j,i),R2_f(k,j,i))
          end do
       end do
    end do

    ! interface geometry
    do i = 1-nrg,n3_f+nrg
       do j = 1-nrg,n2_f+nrg
          do k = 1-nrg,n1_f+nrg
             Tau_fi(k,j,i) = interface_cf(R1_f(k,j,i),R2_f(k,j,i))
          end do
       end do
    end do

    ! scale the mesh
    X1_f = R1_f*l1
    X2_f = R2_f*l2

    ! X3
    X3_f(:,:,:) = R3_f*Tau_f(:,:,:) + (1.d0-R3_f)*Tau_fi(:,:,:)

    !
  end subroutine generate_grid

  subroutine metric_derivative(r1,r2,r3,x1r1,x1r2,x1r3,x2r1,x2r2,x2r3,x3r1,x3r2,x3r3, &
                              xi11,xi21,xi31,xi12,xi22,xi32,xi13,xi23,xi33,J,flag)
    integer, parameter :: m =6 ! controls the smoothness of the mapping
    real(dp) :: r1,r2,r3,x1r1,x1r2,x1r3,x2r1,x2r2,x2r3,x3r1,x3r2,x3r3
    real(dp) :: xi11,xi21,xi31,xi12,xi22,xi32,xi13,xi23,xi33,J
    integer :: flag
    ! forward derivatives
    select case (flag)
      case (0) ! coarse
       x1r1 = l1
       x1r2 = 0.d0
       x1r3 = 0.d0
       x2r1 = 0.d0
       x2r2 = l2
       x2r3 = 0.d0
       x3r1 = r3*interface_cfx(r1,r2) + (1.d0-r3)*bottomx(r1,r2)
       x3r2 = r3*interface_cfy(r1,r2) + (1.d0-r3)*bottomy(r1,r2)
       x3r3 = interface_cf(r1,r2) - bottom(r1,r2)
      case (1) ! fine
       x1r1 = l1
       x1r2 = 0.d0
       x1r3 = 0.d0
       x2r1 = 0.d0
       x2r2 = l2
       x2r3 = 0.d0
       x3r1 = r3*topx(r1,r2) + (1.d0-r3)*interface_cfx(r1,r2)
       x3r2 = r3*topy(r1,r2) + (1.d0-r3)*interface_cfy(r1,r2)
       x3r3 = top(r1,r2) - interface_cf(r1,r2)
    end select
       J = x1r1*x2r2*x3r3+x3r1*x1r2*x2r3+x1r3*x2r1*x3r2 &
          -x3r1*x2r2*x1r3-x1r1*x3r2*x2r3-x3r3*x2r1*x1r2
      ! backward derivative
       xi11 = (x2r2*x3r3-x3r2*x2r3)/J
       xi21 = (x1r3*x3r2-x1r2*x3r3)/J
       xi31 = (x1r2*x2r3-x1r3*x2r2)/J
       xi12 = (x2r3*x3r1-x2r1*x3r3)/J
       xi22 = (x1r1*x3r3-x1r3*x3r1)/J
       xi32 = (x1r3*x2r1-x1r1*x2r3)/J
       xi13 = (x2r1*x3r2-x2r2*x3r1)/J
       xi23 = (x1r2*x3r1-x1r1*x3r2)/J
       xi33 = (x1r1*x2r2-x1r2*x2r1)/J
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
    u1 = 0.d0
    u2 = 0.d0
    u3 = 0.d0
  end subroutine exact_solution

  subroutine initial_solution(l1,l2,l3,x1,x2,x3,u1,u2,u3)
    implicit none
    real(dp) :: l1,l2,l3,x1,x2,x3,u1,u2,u3
    !
    u1 = exp(-((x1-l1/2)**2)/0.1d0)*exp(-((x2-l2/2)**2)/0.1d0)*exp(-((x3-l3/2)**2)/0.1d0)
    u2 = exp(-((x1-l1/2)**2)/0.2d0)*exp(-((x2-l2/2)**2)/0.2d0)*exp(-((x3-l3/2)**2)/0.2d0)
    u3 = exp(-((x1-l1/2)**2)/0.1d0)*exp(-((x2-l2/2)**2)/0.2d0)*exp(-((x3-l3/2)**2)/0.2d0)
    !
  end subroutine initial_solution

end module problemsetup_new_3d
