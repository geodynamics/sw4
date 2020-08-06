program Elastic_3D_sbp4

  use problemsetup_new_3d
  use SBP_operator

  implicit none

  integer i,j,k,i1,j1,k1,l,m,time_index

  ! Physical domain X
  real(dp), dimension (:), allocatable :: Xgrid_c_1, Xgrid_c_2
  real(dp), dimension (:,:,:), allocatable :: Xgrid_c_3

  ! forcing function and spatial discretization
  real(dp), dimension (:,:,:,:), allocatable :: force_c, force_tt_c, lh_c

  ! Metric derivatives
  real(dp), dimension (:,:,:), allocatable :: XI13_c, XI23_c, XI33_c, Jacobian_c

  ! material parameters
  real(dp), dimension (:,:,:), allocatable :: rho_c, mu_c, lambda_c

  ! Solution
  real(dp), dimension (:,:,:,:,:), allocatable :: u_c

  ! Error
  real(dp), dimension (1:3) :: N6
  real(dp) :: l2_err = 0.d0
  real(dp), dimension (:,:,:), allocatable :: err_c

  ! traction B.C. on the top
  real(dp), dimension (1:dim):: traction_rhs = 0.d0
  real(dp), dimension (:,:,:), allocatable :: traction_data

  ! time
  real(dp):: tv
  integer :: cr,c1,c2,nt
  real(dp) :: rate
  real(dp):: cfl, dt0, dt

  ! stencils
  real(dp):: int_cof
  real(dp), dimension (-1:2) :: P
  real(dp):: uxx_cof(-2:2), ux_cof(-2:2)
  real(dp):: acof(6,8,8), ghcof(6), bof(4,6)
  real(dp):: ghcof_no_gp(6), sbop_no_gp(0:5);
  real(dp):: acof_no_gp(6,8,8);
  real(dp), dimension (0:4):: Sb
  real(dp), dimension (-4:2):: Rop
  real(dp), dimension (-3:3):: RPop

  ! allocate memory for Physical domain X
  allocate(Xgrid_c_1(1-nrg:n1_c+nrg))
  allocate(Xgrid_c_2(1-nrg:n2_c+nrg))
  allocate(Xgrid_c_3(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))

  ! allocate memory for forcing function and spatial discretization
  allocate(force_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim))
  allocate(force_tt_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim))
  allocate(lh_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim))

  ! allocate memory for metric derivatives
  allocate(XI13_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(XI23_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(XI33_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(Jacobian_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))

  ! allocate memory for material parameters
  allocate(rho_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(mu_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(lambda_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))

  ! allocate memory for solutions
  allocate(u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4))

  ! allocate memory for error
  allocate(err_c(1:n1_c,1:n2_c,1:n3_c))

  ! allocate memory for top traction data
  allocate(traction_data(1:n1_c,1:n2_c,1:dim))

  ! Difference stencils and interpolations
  call interpolation_restriction(P,Rop,RPop)
  call central_difference(ux_cof,uxx_cof)
  call varcoeff_NoGhost(acof_no_gp, ghcof_no_gp, sbop_no_gp )
  call VARCOEFFS4( acof, ghcof, Sb )
  call dx_46(bof)

  ! metric derivatives and coefficients of the equation
  call  equation_cof(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3, &
       XI13_c,XI23_c,XI33_c,Jacobian_c,rho_c, &
       mu_c,lambda_c)

  ! estimate time step
  call kappa(tv)
  cfl=sqrt(3.d0/1.40/tv)
  dt0 = min(h1_c,h2_c,h3_c)*cfl*0.9
  nt = ceiling(tn/dt0)+1
  nt = 2326
  dt = tn/nt

  print *, nt
  
  ! exact solution at t = -dt, 0
  ! set initial value to be zero everywhere
  do k=1-nrg,n3_c+nrg
     do i=1-nrg,n2_c+nrg
        do j=1-nrg,n1_c+nrg
           call exact_solution(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k), &
                -dt,u_c(j,i,k,1,1),u_c(j,i,k,2,1),u_c(j,i,k,3,1),0)
           call exact_solution(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k), &
                0.d0,u_c(j,i,k,1,2),u_c(j,i,k,2,2),u_c(j,i,k,3,2),0)
        end do
     end do
  end do

  CALL system_clock(count_rate=cr)
  CALL SYSTEM_CLOCK(c1)

  rate = real(cr)
  tv = 0.d0
  ! time stepping
  do time_index = 1, nt
     ! Predictor step

     ! Evaluate the difference operators
     call Update_interior(u_c(:,:,:,:,2))

     ! Compute forcing functions
     do k=1-nrg,n3_c+nrg
        do i=1-nrg,n2_c+nrg
           do j=1-nrg,n1_c+nrg
              call forcing(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k),tv, &
                       force_c(j,i,k,1),force_c(j,i,k,2),force_c(j,i,k,3))
           end do
        end do
     end do

     ! Update the solution in the predictor step
     do k = 1,n3_c
        do j = 1,n2_c
           do i = 1,n1_c
              u_c(i,j,k,1,3) = 2.d0*u_c(i,j,k,1,2) - u_c(i,j,k,1,1) &
                   + dt**2*(lh_c(i,j,k,1) + Jacobian_c(i,j,k)*force_c(i,j,k,1))/rho_c(i,j,k)
              u_c(i,j,k,2,3) = 2.d0*u_c(i,j,k,2,2) - u_c(i,j,k,2,1) &
                   + dt**2*(lh_c(i,j,k,2) + Jacobian_c(i,j,k)*force_c(i,j,k,2))/rho_c(i,j,k)
              u_c(i,j,k,3,3) = 2.d0*u_c(i,j,k,3,2) - u_c(i,j,k,3,1) &
                   + dt**2*(lh_c(i,j,k,3) + Jacobian_c(i,j,k)*force_c(i,j,k,3))/rho_c(i,j,k)
           end do
        end do
     end do

     ! Update time
     tv = tv + dt

     ! Update ghost points outside the left and right boundaries
     ! The argument '3' means time level star
     ! '1', '2' and '4' mean time level n-1, n and n+1, respectively.
     call  Update_gp(3)

     ! Update ghost point values for the traction boundary
     call  Update_traction(3)

     ! Update Dirichlet boundary condition
     call Update_Dirichlet_BC(3)

     ! Corrector step

     ! Evaluate the difference operators
     call Update_interior((u_c(:,:,:,:,3) - 2.d0*u_c(:,:,:,:,2) + u_c(:,:,:,:,1))/dt**2)

     ! Compute the second time derivative of the forcing fucntions
     do k=1-nrg,n3_c+nrg
        do i=1-nrg,n2_c+nrg
           do j=1-nrg,n1_c+nrg
              call forcing_tt(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k), &
                       tv-dt,force_tt_c(j,i,k,1),force_tt_c(j,i,k,2),force_tt_c(j,i,k,3))
           end do
        end do
     end do

     ! Evaluate the difference operators
      do k = 1,n3_c
         do j = 1,n2_c
            do i = 1,n1_c
               u_c(i,j,k,1,4) = u_c(i,j,k,1,3) + dt**4/12.d0*(lh_c(i,j,k,1) &
                   + Jacobian_c(i,j,k)*force_tt_c(i,j,k,1))/rho_c(i,j,k)
               u_c(i,j,k,2,4) = u_c(i,j,k,2,3) + dt**4/12.d0*(lh_c(i,j,k,2) &
                   + Jacobian_c(i,j,k)*force_tt_c(i,j,k,2))/rho_c(i,j,k)
               u_c(i,j,k,3,4) = u_c(i,j,k,3,3) + dt**4/12.d0*(lh_c(i,j,k,3) &
                   + Jacobian_c(i,j,k)*force_tt_c(i,j,k,3))/rho_c(i,j,k)
            end do
         end do
      end do

     ! Update ghost point values outside the left and right boundaries
     call  Update_gp(4)

     ! Update ghost point values for the traction boundary
     call  Update_traction(4)

     ! Update Dirichlet boundary condition
     call  Update_Dirichlet_BC(4)

     ! Update solutions
     u_c(:,:,:,:,1) = u_c(:,:,:,:,2)
     u_c(:,:,:,:,2) = u_c(:,:,:,:,4)
     if (time_index .eq. (nt/5)) then
       print *, tv
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_1.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_1.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_1.txt')
     end if
     if (time_index .eq. 2*(nt/5)) then
       print *, tv
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_2.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_2.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_2.txt')
     end if
     if (time_index .eq. 3*(nt/5)) then
       print *, tv
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_3.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_3.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_3.txt')
     end if
     if (time_index .eq. 4*(nt/5)) then
       print *, tv
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_4.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_4.txt')
       call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_4.txt')
     end if

  end do ! end of time loop

  CALL SYSTEM_CLOCK(c2)
  ! exact solution at final time are stored in u_c(:,:,:,3)

  ! coarse mesh
  !do k=1,n3_c
  !   do i=1,n2_c
  !      do j=1,n1_c
  !         call exact_solution(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k), &
  !                  tv,u_c(j,i,k,1,3),u_c(j,i,k,2,3),u_c(j,i,k,3,3),0)
  !         err_c(j,i,k) = max(abs(u_c(j,i,k,1,4)-u_c(j,i,k,1,3)),abs(u_c(j,i,k,2,4)-u_c(j,i,k,2,3)), &
  !                  abs(u_c(j,i,k,3,4)-u_c(j,i,k,3,3)))
  !         l2_err = l2_err + h1_c*h2_c*h3_c*((u_c(j,i,k,1,4)-u_c(j,i,k,1,3))**2 &
  !            +(u_c(j,i,k,2,4)-u_c(j,i,k,2,3))**2+(u_c(j,i,k,3,4)-u_c(j,i,k,3,3))**2)
  !      end do
  !   end do
  !end do

  !l2_err = sqrt(l2_err)
  N6(1) = n1_c
  N6(2) = n2_c
  N6(3) = n3_c

  call print_array_to_file(1,1,3,N6,'N6.txt')
  call print_array_to_file(1,1,n1_c,Xgrid_c_1(1:n1_c),'X1c.txt')
  call print_array_to_file(1,1,n2_c,Xgrid_c_2(1:n2_c),'X2c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c_3(1:n1_c,1:n2_c,1:n3_c),'X3c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_5.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_5.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_5.txt')

  !call print_array_to_file(n1_c,n2_c,n3_c,err_c(1:n1_c,1:n2_c,1:n3_c),'err_c.txt')

  !write(*,"(A20,6I5)") 'No. of grid points', n1_c,n2_c,n3_c
  !write(*,"(A20,5ES25.15E3)") 'errors',maxval(err_c), l2_err
  !write(*,"(A20,5ES25.15E3)") 'errors/16',maxval(err_c)/16.d0, l2_err/16.d0
  write(*,"(A20,5ES25.15E3)") 'computational time', (c2-c1)/rate
  write(*,"(A20,5ES25.15E3)") 'cfl', cfl


contains

  subroutine equation_cof(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3, &
       XI13_c,XI23_c,XI33_c,Jacobian_c,rho_c,&
       mu_c,lambda_c)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: rho_c,mu_c,lambda_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: Xgrid_c_3
    real(dp), dimension (1-nrg:n1_c+nrg) :: Xgrid_c_1
    real(dp), dimension (1-nrg:n2_c+nrg) :: Xgrid_c_2
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI13_c,XI23_c,XI33_c,Jacobian_c

    call generate_grid(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3)

    ! compute metric derivatives, coarse
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             call metric_derivative(dble(k-1)/(n1_c-1),dble(j-1)/(n2_c-1),dble(i-1)/(n3_c-1),&
                 XI13_c(k,j,i),XI23_c(k,j,i),XI33_c(k,j,i),Jacobian_c(k,j,i))
          end do
       end do
    end do

    ! variable coefficients
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             mu_c(k,j,i) = (1000.d0**2)*1500.d0!3.d0 + sin(3.d0*Xgrid_c_1(k)+0.1d0)*sin(3.d0*Xgrid_c_2(j)+0.1d0)*sin(Xgrid_c_3(k,j,i))
             lambda_c(k,j,i) = (2000.d0**2)*1500.d0-2.d0*(1000.d0**2)*1500.d0 !21.d0+ cos(Xgrid_c_1(k)+0.1d0)*cos(Xgrid_c_2(j)+0.1d0)*sin(3.d0*Xgrid_c_3(k,j,i))**2
             rho_c(k,j,i) = 1500.d0!2.d0 + sin(Xgrid_c_1(k)+0.3d0)*sin(Xgrid_c_2(j)+0.3d0)*sin(Xgrid_c_3(k,j,i)-0.2d0)
          end do
       end do
    end do

    rho_c = rho_c*Jacobian_c

  end subroutine  equation_cof
  !
  subroutine kappa(kappa2)
    real(dp), dimension (3) :: kappa1
    real(dp), dimension (3,3) :: mat_t
    real(dp), dimension (8) :: work
    integer :: info
    real(dp) :: kappa2
    kappa2 = 0.d0
    do k = 1, n3_c
       do j = 1, n2_c
          do i = 1,n1_c
            mat_t(1,1) = (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))/l1**2
            mat_t(1,2) = 0.d0
            mat_t(1,3) = (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))/l1*XI13_c(i,j,k)
            mat_t(2,1) = 0.d0
            mat_t(2,2) = (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))/l2**2
            mat_t(2,3) = (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))/l2*XI23_c(i,j,k)
            mat_t(3,1) = mat_t(1,3)
            mat_t(3,2) = mat_t(2,3)
            mat_t(3,3) = (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))*XI13_c(i,j,k)*XI13_c(i,j,k) &
                       + (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))*XI23_c(i,j,k)*XI23_c(i,j,k) &
                       + (4.d0*mu_c(i,j,k)+lambda_c(i,j,k))*XI33_c(i,j,k)*XI33_c(i,j,k)
            mat_t = Jacobian_c(i,j,k)*mat_t/rho_c(i,j,k)
             call dsyev('N','U',3,mat_t,3,kappa1,work,8,info)
             if (kappa1(3) > kappa2) then
                kappa2 = kappa1(3)
             end if
          end do
       end do
    end do
  end subroutine kappa
  !
  subroutine Update_interior(u_c_t)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t
    integer :: k1

    ! Difference operators in the interior of the domains

    lh_c = 0.d0
    do k = 1,n3_c
       do j = 1,n2_c
          do i = 1,n1_c
             ! second derivative 11 & 22 & 12 & 21
             ! first set
             lh_c(i,j,k,1) = lh_c(i,j,k,1)+((-Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/8.d0 &
                      + Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0)*u_c_t(i-2,j,k,1) &
                      +(Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/6.d0 &
                      + Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/2.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/6.d0)*u_c_t(i-1,j,k,1) &
                      +(-Jacobian_c(i-2,j,k)*(2.d0*mu_c(i-2,j,k)+lambda_c(i-2,j,k))/24.d0 &
                      - Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))*5.d0/6.d0 &
                      -Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))*3.d0/4.d0 &
                      -Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))*5.d0/6.d0 &
                      -Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/24.d0)*u_c_t(i-0,j,k,1) &
                      +(Jacobian_c(i-1,j,k)*(2.d0*mu_c(i-1,j,k)+lambda_c(i-1,j,k))/6.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/2.d0 &
                      + Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/6.d0)*u_c_t(i+1,j,k,1) &
                      +(-Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0 &
                      + Jacobian_c(i+1,j,k)*(2.d0*mu_c(i+1,j,k)+lambda_c(i+1,j,k))/6.d0 &
                      - Jacobian_c(i+2,j,k)*(2.d0*mu_c(i+2,j,k)+lambda_c(i+2,j,k))/8.d0)*u_c_t(i+2,j,k,1))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/8.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i,j-2,k,1) &
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/6.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0)*u_c_t(i,j-1,k,1) &
                      +(-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/24.d0 - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*5.d0/6.d0 -Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/24.d0)*u_c_t(i,j-0,k,1) &
                      +(Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/2.d0 &
                      + Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/6.d0)*u_c_t(i,j+1,k,1) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/8.d0)*u_c_t(i,j+2,k,1))/h2_c**2/l2**2&
                      +(Jacobian_c(i-2,j,k)*lambda_c(i-2,j,k)*(u_c_t(i-2,j-2,k,2)/12.d0&
                      -u_c_t(i-2,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i-2,j+1,k,2)*2.d0/3.d0 - u_c_t(i-2,j+2,k,2)/12.d0)/12.d0 &
                      - Jacobian_c(i-1,j,k)*lambda_c(i-1,j,k)*(u_c_t(i-1,j-2,k,2)/12.d0-u_c_t(i-1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i-1,j+1,k,2)*2.d0/3.d0 - u_c_t(i-1,j+2,k,2)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i+1,j,k)*lambda_c(i+1,j,k)*(u_c_t(i+1,j-2,k,2)/12.d0-u_c_t(i+1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,2)*2.d0/3.d0 - u_c_t(i+1,j+2,k,2)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i+2,j,k)*lambda_c(i+2,j,k)*(u_c_t(i+2,j-2,k,2)/12.d0-u_c_t(i+2,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+2,j+1,k,2)*2.d0/3.d0 - u_c_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_c/h2_c&
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)*(u_c_t(i-2,j-2,k,2)/12.d0&
                      -u_c_t(i-1,j-2,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j-2,k,2)*2.d0/3.d0 - u_c_t(i+2,j-2,k,2)/12.d0)/12.d0 &
                      - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*(u_c_t(i-2,j-1,k,2)/12.d0-u_c_t(i-1,j-1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j-1,k,2)*2.d0/3.d0 - u_c_t(i+2,j-1,k,2)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*(u_c_t(i-2,j+1,k,2)/12.d0-u_c_t(i-1,j+1,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,2)*2.d0/3.d0 - u_c_t(i+2,j+1,k,2)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)*(u_c_t(i-2,j+2,k,2)/12.d0-u_c_t(i-1,j+2,k,2)*2.d0/3.d0 &
                      + u_c_t(i+1,j+2,k,2)*2.d0/3.d0 - u_c_t(i+2,j+2,k,2)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
             ! second set
             lh_c(i,j,k,2) = lh_c(i,j,k,2)+((-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/8.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i-2,j,k,2) &
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/6.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0)*u_c_t(i-1,j,k,2) &
                      +(-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/24.d0 - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*5.d0/6.d0 -Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/24.d0)*u_c_t(i-0,j,k,2) &
                      +(Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/2.d0 &
                      + Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/6.d0)*u_c_t(i+1,j,k,2) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/8.d0)*u_c_t(i+2,j,k,2))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/8.d0 &
                      + Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0)*u_c_t(i,j-2,k,2) &
                      +(Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/6.d0 &
                      + Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/2.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/6.d0)*u_c_t(i,j-1,k,2) &
                      +(-Jacobian_c(i,j-2,k)*(2.d0*mu_c(i,j-2,k)+lambda_c(i,j-2,k))/24.d0 &
                      - Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))*5.d0/6.d0 &
                      -Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/24.d0)*u_c_t(i,j-0,k,2) &
                      +(Jacobian_c(i,j-1,k)*(2.d0*mu_c(i,j-1,k)+lambda_c(i,j-1,k))/6.d0 &
                      + Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/2.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/2.d0 &
                      + Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/6.d0)*u_c_t(i,j+1,k,2) &
                      +(-Jacobian_c(i,j,k)*(2.d0*mu_c(i,j,k)+lambda_c(i,j,k))/8.d0 &
                      + Jacobian_c(i,j+1,k)*(2.d0*mu_c(i,j+1,k)+lambda_c(i,j+1,k))/6.d0 &
                      - Jacobian_c(i,j+2,k)*(2.d0*mu_c(i,j+2,k)+lambda_c(i,j+2,k))/8.d0)*u_c_t(i,j+2,k,2))/h2_c**2/l2**2&
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)*(u_c_t(i-2,j-2,k,1)/12.d0&
                      -u_c_t(i-2,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i-2,j+1,k,1)*2.d0/3.d0 - u_c_t(i-2,j+2,k,1)/12.d0)/12.d0 &
                      - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*(u_c_t(i-1,j-2,k,1)/12.d0-u_c_t(i-1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i-1,j+1,k,1)*2.d0/3.d0 - u_c_t(i-1,j+2,k,1)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*(u_c_t(i+1,j-2,k,1)/12.d0-u_c_t(i+1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,1)*2.d0/3.d0 - u_c_t(i+1,j+2,k,1)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)*(u_c_t(i+2,j-2,k,1)/12.d0-u_c_t(i+2,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+2,j+1,k,1)*2.d0/3.d0 - u_c_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_c/h2_c&
                      +(Jacobian_c(i,j-2,k)*lambda_c(i,j-2,k)*(u_c_t(i-2,j-2,k,1)/12.d0&
                      -u_c_t(i-1,j-2,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j-2,k,1)*2.d0/3.d0 - u_c_t(i+2,j-2,k,1)/12.d0)/12.d0 &
                      - Jacobian_c(i,j-1,k)*lambda_c(i,j-1,k)*(u_c_t(i-2,j-1,k,1)/12.d0-u_c_t(i-1,j-1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j-1,k,1)*2.d0/3.d0 - u_c_t(i+2,j-1,k,1)/12.d0)*2.d0/3.d0 &
                      + Jacobian_c(i,j+1,k)*lambda_c(i,j+1,k)*(u_c_t(i-2,j+1,k,1)/12.d0-u_c_t(i-1,j+1,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+1,k,1)*2.d0/3.d0 - u_c_t(i+2,j+1,k,1)/12.d0)*2.d0/3.d0 &
                      - Jacobian_c(i,j+2,k)*lambda_c(i,j+2,k)*(u_c_t(i-2,j+2,k,1)/12.d0-u_c_t(i-1,j+2,k,1)*2.d0/3.d0 &
                      + u_c_t(i+1,j+2,k,1)*2.d0/3.d0 - u_c_t(i+2,j+2,k,1)/12.d0)/12.d0)/l1/l2/h1_c/h2_c
             ! third set
             lh_c(i,j,k,3) = lh_c(i,j,k,3)+((-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/8.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i-2,j,k,3) &
                      +(Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/6.d0 + Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0)*u_c_t(i-1,j,k,3) &
                      +(-Jacobian_c(i-2,j,k)*mu_c(i-2,j,k)/24.d0 - Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)*5.d0/6.d0 -Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/24.d0)*u_c_t(i-0,j,k,3) &
                      +(Jacobian_c(i-1,j,k)*mu_c(i-1,j,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/2.d0 &
                      + Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/6.d0)*u_c_t(i+1,j,k,3) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i+1,j,k)*mu_c(i+1,j,k)/6.d0 &
                      - Jacobian_c(i+2,j,k)*mu_c(i+2,j,k)/8.d0)*u_c_t(i+2,j,k,3))/h1_c**2/l1**2&
                      +((-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/8.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0)*u_c_t(i,j-2,k,3) &
                      +(Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/6.d0 + Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/2.d0 &
                      + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0)*u_c_t(i,j-1,k,3) &
                      +(-Jacobian_c(i,j-2,k)*mu_c(i,j-2,k)/24.d0 - Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)*5.d0/6.d0 &
                      - Jacobian_c(i,j,k)*mu_c(i,j,k)*3.d0/4.d0 &
                      - Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)*5.d0/6.d0 -Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/24.d0)*u_c_t(i,j-0,k,3) &
                      +(Jacobian_c(i,j-1,k)*mu_c(i,j-1,k)/6.d0 + Jacobian_c(i,j,k)*mu_c(i,j,k)/2.d0 &
                      + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/2.d0 &
                      + Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/6.d0)*u_c_t(i,j+1,k,3) &
                      +(-Jacobian_c(i,j,k)*mu_c(i,j,k)/8.d0 + Jacobian_c(i,j+1,k)*mu_c(i,j+1,k)/6.d0 &
                      - Jacobian_c(i,j+2,k)*mu_c(i,j+2,k)/8.d0)*u_c_t(i,j+2,k,3))/h2_c**2/l2**2
          end do
       end do
    end do
    !
    do k = 1,4
       do i = 1,n2_c
          do j = 1,n1_c
             do k1 = 1,6
                ! mixed derivative 13 & 23
                ! first set equation
                lh_c(j,i,k,1) = lh_c(j,i,k,1) &
                      +(Jacobian_c(j-2,i,k)*(2.d0*mu_c(j-2,i,k)+lambda_c(j-2,i,k))*XI13_c(j-2,i,k)*bof(k,k1)&
                      *u_c_t(j-2,i,k1,1)/12.d0-Jacobian_c(j-1,i,k)*(2.d0*mu_c(j-1,i,k)+lambda_c(j-1,i,k))&
                      *XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0+Jacobian_c(j+1,i,k)*(2.d0*mu_c(j+1,i,k)&
                      +lambda_c(j+1,i,k))*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*(2.d0*mu_c(j+2,i,k)&
                      +lambda_c(j+2,i,k))*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI23_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,2)/12.d0 &
                      -Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI23_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI23_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI23_c(j+2,i,k)*bof(k,k1)&
                      *u_c_t(j+2,i,k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI33_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,3)/12.d0 &
                      -Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI33_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI33_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI33_c(j+2,i,k)&
                      *bof(k,k1)*u_c_t(j+2,i,k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,1)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI13_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI13_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI13_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI13_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,1) = lh_c(j,i,n3_c+1-k,1)+ &
                      (-Jacobian_c(j-2,i,n3_c+1-k)*(2.d0*mu_c(j-2,i,n3_c+1-k)+lambda_c(j-2,i,n3_c+1-k))&
                      *XI13_c(j-2,i,n3_c+1-k)*bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*(2.d0*mu_c(j-1,i,n3_c+1-k)+lambda_c(j-1,i,n3_c+1-k))&
                      *XI13_c(j-1,i,n3_c+1-k)*bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*(2.d0*mu_c(j+1,i,n3_c+1-k)+lambda_c(j+1,i,n3_c+1-k))&
                      *XI13_c(j+1,i,n3_c+1-k)*bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*(2.d0*mu_c(j+2,i,n3_c+1-k)+lambda_c(j+2,i,n3_c+1-k))&
                      *XI13_c(j+2,i,n3_c+1-k)*bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*lambda_c(j-2,i,n3_c+1-k)*XI23_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*lambda_c(j-1,i,n3_c+1-k)*XI23_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*lambda_c(j+1,i,n3_c+1-k)*XI23_c(j+1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*lambda_c(j+2,i,n3_c+1-k)*XI23_c(j+2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+2,i,n3_c+1-k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*lambda_c(j-2,i,n3_c+1-k)*XI33_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*lambda_c(j-1,i,n3_c+1-k)*XI33_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*lambda_c(j+1,i,n3_c+1-k)*XI33_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*lambda_c(j+2,i,n3_c+1-k)*XI33_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI23_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI23_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI23_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI13_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI13_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI13_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI13_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c
                ! second set equation
                lh_c(j,i,k,2) = lh_c(j,i,k,2) &
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI23_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,1)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI23_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI23_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI23_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,2)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI13_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,1)/12.d0 &
                      -Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI13_c(j,i-1,k)*bof(k,k1)*u_c_t(j,i-1,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI13_c(j,i+1,k)*bof(k,k1)*u_c_t(j,i+1,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI13_c(j,i+2,k)*bof(k,k1)*u_c_t(j,i+2,k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*(2.d0*mu_c(j,i-2,k)+lambda_c(j,i-2,k))&
                      *XI23_c(j,i-2,k)*bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*(2.d0*mu_c(j,i-1,k)+lambda_c(j,i-1,k))*XI23_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*(2.d0*mu_c(j,i+1,k)+lambda_c(j,i+1,k))*XI23_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*(2.d0*mu_c(j,i+2,k)+lambda_c(j,i+2,k))*XI23_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI33_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,3)/12.d0 &
                      -Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI33_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI33_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI33_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,3)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,2) = lh_c(j,i,n3_c+1-k,2)+ &
                      (-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI23_c(j-2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI23_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI23_c(j+1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI23_c(j+2,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI13_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI13_c(j-1,i,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j-1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI13_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI13_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,2)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*lambda_c(j,i-2,n3_c+1-k)&
                      *XI13_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*lambda_c(j,i-1,n3_c+1-k)*XI13_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*lambda_c(j,i+1,n3_c+1-k)*XI13_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*lambda_c(j,i+2,n3_c+1-k)*XI13_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,1)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*(2.d0*mu_c(j,i-2,n3_c+1-k)&
                      +lambda_c(j,i-2,n3_c+1-k))*XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*(2.d0*mu_c(j,i-1,n3_c+1-k)+lambda_c(j,i-1,n3_c+1-k))&
                      *XI23_c(j,i-1,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*(2.d0*mu_c(j,i+1,n3_c+1-k)+lambda_c(j,i+1,n3_c+1-k))&
                      *XI23_c(j,i+1,n3_c+1-k)*bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*(2.d0*mu_c(j,i+2,n3_c+1-k)+lambda_c(j,i+2,n3_c+1-k))&
                      *XI23_c(j,i+2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*lambda_c(j,i-2,n3_c+1-k)&
                      *XI33_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*lambda_c(j,i-1,n3_c+1-k)*XI33_c(j,i-1,n3_c+1-k)*bof(k,k1)&
                      *u_c_t(j,i-1,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*lambda_c(j,i+1,n3_c+1-k)*XI33_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*lambda_c(j,i+2,n3_c+1-k)*XI33_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,3)/12.d0)/l2/h2_c/h3_c
                ! third set equation
                lh_c(j,i,k,3) = lh_c(j,i,k,3) &
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI33_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,1)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI33_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI33_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI33_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*bof(k,k1)*u_c_t(j-2,i,k1,3)/12.d0 &
                      -Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*bof(k,k1)*u_c_t(j-1,i,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*bof(k,k1)*u_c_t(j+1,i,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*bof(k,k1)*u_c_t(j+2,i,k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI33_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,2)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI33_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI33_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI33_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                      *bof(k,k1)*u_c_t(j,i-2,k1,3)/12.d0 &
                      -Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)&
                      *bof(k,k1)*u_c_t(j,i-1,k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)&
                      *bof(k,k1)*u_c_t(j,i+1,k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)&
                      *bof(k,k1)*u_c_t(j,i+2,k1,3)/12.d0)/l2/h2_c/h3_c

                lh_c(j,i,n3_c+1-k,3) = lh_c(j,i,n3_c+1-k,3) +&
                      (-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI33_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,1)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI33_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI33_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,1)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI33_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,1)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j-2,i,n3_c+1-k)*mu_c(j-2,i,n3_c+1-k)*XI13_c(j-2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-2,i,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j-1,i,n3_c+1-k)*mu_c(j-1,i,n3_c+1-k)*XI13_c(j-1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j-1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j+1,i,n3_c+1-k)*mu_c(j+1,i,n3_c+1-k)*XI13_c(j+1,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+1,i,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j+2,i,n3_c+1-k)*mu_c(j+2,i,n3_c+1-k)*XI13_c(j+2,i,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j+2,i,n3_c+1-k1,3)/12.d0)/l1/h1_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI33_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,2)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI33_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,2)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI33_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,2)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI33_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,2)/12.d0)/l2/h2_c/h3_c&
                      +(-Jacobian_c(j,i-2,n3_c+1-k)*mu_c(j,i-2,n3_c+1-k)&
                      *XI23_c(j,i-2,n3_c+1-k)*bof(k,k1)*u_c_t(j,i-2,n3_c+1-k1,3)/12.d0 &
                      +Jacobian_c(j,i-1,n3_c+1-k)*mu_c(j,i-1,n3_c+1-k)*XI23_c(j,i-1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i-1,n3_c+1-k1,3)*2.d0/3.d0 &
                      -Jacobian_c(j,i+1,n3_c+1-k)*mu_c(j,i+1,n3_c+1-k)*XI23_c(j,i+1,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+1,n3_c+1-k1,3)*2.d0/3.d0 &
                      +Jacobian_c(j,i+2,n3_c+1-k)*mu_c(j,i+2,n3_c+1-k)*XI23_c(j,i+2,n3_c+1-k)&
                      *bof(k,k1)*u_c_t(j,i+2,n3_c+1-k1,3)/12.d0)/l2/h2_c/h3_c
             end do
          end do
       end do
    end do
    !
    do k = 5,n3_c-4
       do i = 1,n2_c
          do j = 1,n1_c
             ! mixed derivative 13 & 23
             ! first set equation
             lh_c(j,i,k,1) = lh_c(j,i,k,1)+(Jacobian_c(j-2,i,k)*(2.d0*mu_c(j-2,i,k)+lambda_c(j-2,i,k))*XI13_c(j-2,i,k)&
                   *(u_c_t(j-2,i,k-2,1)/12.d0-u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*(2.d0*mu_c(j-1,i,k)+lambda_c(j-1,i,k))*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0&
                   -u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*(2.d0*mu_c(j+1,i,k)+lambda_c(j+1,i,k))&
                   *XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*(2.d0*mu_c(j+2,i,k)+lambda_c(j+2,i,k))*XI13_c(j+2,i,k)&
                   *(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI23_c(j-2,i,k)*(u_c_t(j-2,i,k-2,2)/12.d0&
                   -u_c_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,2)*2.d0/3.d0-u_c_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI23_c(j-1,i,k)*(u_c_t(j-1,i,k-2,2)/12.d0&
                   -u_c_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,2)*2.d0/3.d0-u_c_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI23_c(j+1,i,k)*(u_c_t(j+1,i,k-2,2)/12.d0&
                   -u_c_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,2)*2.d0/3.d0-u_c_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI23_c(j+2,i,k)*(u_c_t(j+2,i,k-2,2)/12.d0&
                   -u_c_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,2)*2.d0/3.d0-u_c_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*lambda_c(j-2,i,k)*XI33_c(j-2,i,k)*(u_c_t(j-2,i,k-2,3)/12.d0&
                   -u_c_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,3)*2.d0/3.d0-u_c_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*lambda_c(j-1,i,k)*XI33_c(j-1,i,k)*(u_c_t(j-1,i,k-2,3)/12.d0-u_c_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,3)*2.d0/3.d0-u_c_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*lambda_c(j+1,i,k)*XI33_c(j+1,i,k)*(u_c_t(j+1,i,k-2,3)/12.d0-u_c_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,3)*2.d0/3.d0-u_c_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*lambda_c(j+2,i,k)*XI33_c(j+2,i,k)*(u_c_t(j+2,i,k-2,3)/12.d0-u_c_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,3)*2.d0/3.d0-u_c_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,1)/12.d0-u_c_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,1)*2.d0/3.d0-u_c_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*(u_c_t(j,i-1,k-2,1)/12.d0-u_c_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,1)*2.d0/3.d0-u_c_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*(u_c_t(j,i+1,k-2,1)/12.d0-u_c_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,1)*2.d0/3.d0-u_c_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*(u_c_t(j,i+2,k-2,1)/12.d0-u_c_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,1)*2.d0/3.d0-u_c_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI13_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI13_c(j,i-1,k)*(u_c_t(j,i-1,k-2,2)/12.d0-u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI13_c(j,i+1,k)*(u_c_t(j,i+1,k-2,2)/12.d0-u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI13_c(j,i+2,k)*(u_c_t(j,i+2,k-2,2)/12.d0-u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c
             ! second set equation
             lh_c(j,i,k,2) = lh_c(j,i,k,2)+(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI23_c(j-2,i,k)*(u_c_t(j-2,i,k-2,1)/12.d0&
                   -u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI23_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0-u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI23_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI23_c(j+2,i,k)*(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*(u_c_t(j-2,i,k-2,2)/12.d0&
                   -u_c_t(j-2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,2)*2.d0/3.d0-u_c_t(j-2,i,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,2)/12.d0-u_c_t(j-1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,2)*2.d0/3.d0-u_c_t(j-1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,2)/12.d0-u_c_t(j+1,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,2)*2.d0/3.d0-u_c_t(j+1,i,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*(u_c_t(j+2,i,k-2,2)/12.d0-u_c_t(j+2,i,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,2)*2.d0/3.d0-u_c_t(j+2,i,k+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI13_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,1)/12.d0-u_c_t(j,i-2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,1)*2.d0/3.d0-u_c_t(j,i-2,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI13_c(j,i-1,k)*(u_c_t(j,i-1,k-2,1)/12.d0&
                   -u_c_t(j,i-1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,1)*2.d0/3.d0-u_c_t(j,i-1,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI13_c(j,i+1,k)*(u_c_t(j,i+1,k-2,1)/12.d0&
                   -u_c_t(j,i+1,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,1)*2.d0/3.d0-u_c_t(j,i+1,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI13_c(j,i+2,k)*(u_c_t(j,i+2,k-2,1)/12.d0&
                   -u_c_t(j,i+2,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,1)*2.d0/3.d0-u_c_t(j,i+2,k+2,1)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*(2.d0*mu_c(j,i-2,k)+lambda_c(j,i-2,k))&
                   *XI23_c(j,i-2,k)*(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*(2.d0*mu_c(j,i-1,k)+lambda_c(j,i-1,k))*XI23_c(j,i-1,k)&
                   *(u_c_t(j,i-1,k-2,2)/12.d0-u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*(2.d0*mu_c(j,i+1,k)+lambda_c(j,i+1,k))*XI23_c(j,i+1,k)&
                   *(u_c_t(j,i+1,k-2,2)/12.d0-u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*(2.d0*mu_c(j,i+2,k)+lambda_c(j,i+2,k))*XI23_c(j,i+2,k)&
                   *(u_c_t(j,i+2,k-2,2)/12.d0-u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*lambda_c(j,i-2,k)*XI33_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,3)/12.d0-u_c_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,3)*2.d0/3.d0-u_c_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*lambda_c(j,i-1,k)*XI33_c(j,i-1,k)*(u_c_t(j,i-1,k-2,3)/12.d0&
                   -u_c_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,3)*2.d0/3.d0-u_c_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*lambda_c(j,i+1,k)*XI33_c(j,i+1,k)*(u_c_t(j,i+1,k-2,3)/12.d0&
                   -u_c_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,3)*2.d0/3.d0-u_c_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*lambda_c(j,i+2,k)*XI33_c(j,i+2,k)*(u_c_t(j,i+2,k-2,3)/12.d0&
                   -u_c_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,3)*2.d0/3.d0-u_c_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
             ! third set equation
             lh_c(j,i,k,3) = lh_c(j,i,k,3)+(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI33_c(j-2,i,k)*(u_c_t(j-2,i,k-2,1)/12.d0&
                   -u_c_t(j-2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,1)*2.d0/3.d0-u_c_t(j-2,i,k+2,1)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI33_c(j-1,i,k)*(u_c_t(j-1,i,k-2,1)/12.d0-u_c_t(j-1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,1)*2.d0/3.d0-u_c_t(j-1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI33_c(j+1,i,k)*(u_c_t(j+1,i,k-2,1)/12.d0-u_c_t(j+1,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,1)*2.d0/3.d0-u_c_t(j+1,i,k+2,1)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI33_c(j+2,i,k)*(u_c_t(j+2,i,k-2,1)/12.d0-u_c_t(j+2,i,k-1,1)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,1)*2.d0/3.d0-u_c_t(j+2,i,k+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j-2,i,k)*mu_c(j-2,i,k)*XI13_c(j-2,i,k)*(u_c_t(j-2,i,k-2,3)/12.d0&
                   -u_c_t(j-2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-2,i,k+1,3)*2.d0/3.d0-u_c_t(j-2,i,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j-1,i,k)*mu_c(j-1,i,k)*XI13_c(j-1,i,k)*(u_c_t(j-1,i,k-2,3)/12.d0-u_c_t(j-1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j-1,i,k+1,3)*2.d0/3.d0-u_c_t(j-1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j+1,i,k)*mu_c(j+1,i,k)*XI13_c(j+1,i,k)*(u_c_t(j+1,i,k-2,3)/12.d0-u_c_t(j+1,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+1,i,k+1,3)*2.d0/3.d0-u_c_t(j+1,i,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j+2,i,k)*mu_c(j+2,i,k)*XI13_c(j+2,i,k)*(u_c_t(j+2,i,k-2,3)/12.d0-u_c_t(j+2,i,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j+2,i,k+1,3)*2.d0/3.d0-u_c_t(j+2,i,k+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI33_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,2)/12.d0-u_c_t(j,i-2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,2)*2.d0/3.d0-u_c_t(j,i-2,k+2,2)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI33_c(j,i-1,k)*(u_c_t(j,i-1,k-2,2)/12.d0&
                   -u_c_t(j,i-1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,2)*2.d0/3.d0-u_c_t(j,i-1,k+2,2)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI33_c(j,i+1,k)*(u_c_t(j,i+1,k-2,2)/12.d0&
                   -u_c_t(j,i+1,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,2)*2.d0/3.d0-u_c_t(j,i+1,k+2,2)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI33_c(j,i+2,k)*(u_c_t(j,i+2,k-2,2)/12.d0&
                   -u_c_t(j,i+2,k-1,2)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,2)*2.d0/3.d0-u_c_t(j,i+2,k+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(j,i-2,k)*mu_c(j,i-2,k)*XI23_c(j,i-2,k)&
                   *(u_c_t(j,i-2,k-2,3)/12.d0-u_c_t(j,i-2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-2,k+1,3)*2.d0/3.d0-u_c_t(j,i-2,k+2,3)/12.d0)/12.d0 &
                   - Jacobian_c(j,i-1,k)*mu_c(j,i-1,k)*XI23_c(j,i-1,k)*(u_c_t(j,i-1,k-2,3)/12.d0&
                   -u_c_t(j,i-1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i-1,k+1,3)*2.d0/3.d0-u_c_t(j,i-1,k+2,3)/12.d0)*2.d0/3.d0 &
                   + Jacobian_c(j,i+1,k)*mu_c(j,i+1,k)*XI23_c(j,i+1,k)*(u_c_t(j,i+1,k-2,3)/12.d0&
                   -u_c_t(j,i+1,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+1,k+1,3)*2.d0/3.d0-u_c_t(j,i+1,k+2,3)/12.d0)*2.d0/3.d0 &
                   - Jacobian_c(j,i+2,k)*mu_c(j,i+2,k)*XI23_c(j,i+2,k)*(u_c_t(j,i+2,k-2,3)/12.d0&
                   -u_c_t(j,i+2,k-1,3)*2.d0/3.d0 &
                   + u_c_t(j,i+2,k+1,3)*2.d0/3.d0-u_c_t(j,i+2,k+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
           end do
        end do
    end do
    !
    do i = 5,n3_c-4
       do j = 1,n2_c
          do k = 1,n1_c
             ! mixed derivative 31
             ! first set equation
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+(Jacobian_c(k,j,i-2)*(2.d0*mu_c(k,j,i-2)&
                     +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*(2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)&
                     *(u_c_t(k-2,j,i-1,1)/12.d0-u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*(2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)&
                     *(u_c_t(k-2,j,i+1,1)/12.d0-u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*(2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)&
                     *(u_c_t(k-2,j,i+2,1)/12.d0-u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)*(u_c_t(k-2,j,i-2,2)/12.d0&
                     -u_c_t(k-1,j,i-2,2)*2.d0/3.d0+u_c_t(k+1,j,i-2,2)*2.d0/3.d0-u_c_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k-2,j,i-1,2)/12.d0&
                     -u_c_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,2)*2.d0/3.d0-u_c_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k-2,j,i+1,2)/12.d0&
                     -u_c_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,2)*2.d0/3.d0-u_c_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k-2,j,i+2,2)/12.d0&
                     -u_c_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,2)*2.d0/3.d0-u_c_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI33_c(k,j,i-2)*(u_c_t(k-2,j,i-2,3)/12.d0&
                     -u_c_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,3)*2.d0/3.d0-u_c_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k-2,j,i-1,3)/12.d0&
                     -u_c_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,3)*2.d0/3.d0-u_c_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k-2,j,i+1,3)/12.d0&
                     -u_c_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,3)*2.d0/3.d0-u_c_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k-2,j,i+2,3)/12.d0&
                     -u_c_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,3)*2.d0/3.d0-u_c_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h1_c/h3_c
             ! second set equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI23_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k-2,j,i-1,1)/12.d0&
                     -u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k-2,j,i+1,1)/12.d0&
                     -u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k-2,j,i+2,1)/12.d0&
                     -u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h1_c/h3_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,2)/12.d0-u_c_t(k-1,j,i-2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,2)*2.d0/3.d0-u_c_t(k+2,j,i-2,2)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k-2,j,i-1,2)/12.d0&
                     -u_c_t(k-1,j,i-1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,2)*2.d0/3.d0-u_c_t(k+2,j,i-1,2)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k-2,j,i+1,2)/12.d0&
                     -u_c_t(k-1,j,i+1,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,2)*2.d0/3.d0-u_c_t(k+2,j,i+1,2)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)&
                     *(u_c_t(k-2,j,i+2,2)/12.d0-u_c_t(k-1,j,i+2,2)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,2)*2.d0/3.d0-u_c_t(k+2,j,i+2,2)/12.d0)/12.d0)/l1/h1_c/h3_c
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI33_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,1)/12.d0-u_c_t(k-1,j,i-2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,1)*2.d0/3.d0-u_c_t(k+2,j,i-2,1)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k-2,j,i-1,1)/12.d0&
                     -u_c_t(k-1,j,i-1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,1)*2.d0/3.d0-u_c_t(k+2,j,i-1,1)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k-2,j,i+1,1)/12.d0&
                     -u_c_t(k-1,j,i+1,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,1)*2.d0/3.d0-u_c_t(k+2,j,i+1,1)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k-2,j,i+2,1)/12.d0&
                     -u_c_t(k-1,j,i+2,1)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,1)*2.d0/3.d0-u_c_t(k+2,j,i+2,1)/12.d0)/12.d0)/l1/h3_c/h1_c&
                     +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                     *(u_c_t(k-2,j,i-2,3)/12.d0-u_c_t(k-1,j,i-2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-2,3)*2.d0/3.d0-u_c_t(k+2,j,i-2,3)/12.d0)/12.d0 &
                     -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k-2,j,i-1,3)/12.d0&
                     -u_c_t(k-1,j,i-1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i-1,3)*2.d0/3.d0-u_c_t(k+2,j,i-1,3)/12.d0)*2.d0/3.d0 &
                     +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k-2,j,i+1,3)/12.d0&
                     -u_c_t(k-1,j,i+1,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+1,3)*2.d0/3.d0-u_c_t(k+2,j,i+1,3)/12.d0)*2.d0/3.d0 &
                     -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k-2,j,i+2,3)/12.d0&
                     -u_c_t(k-1,j,i+2,3)*2.d0/3.d0 &
                     +u_c_t(k+1,j,i+2,3)*2.d0/3.d0-u_c_t(k+2,j,i+2,3)/12.d0)/12.d0)/l1/h3_c/h1_c


             ! 32
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,1)/12.d0-u_c_t(k,j-1,i-2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,1)*2.d0/3.d0-u_c_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k,j-2,i-1,1)/12.d0&
                   -u_c_t(k,j-1,i-1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,1)*2.d0/3.d0-u_c_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k,j-2,i+1,1)/12.d0&
                   -u_c_t(k,j-1,i+1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,1)*2.d0/3.d0-u_c_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k,j-2,i+2,1)/12.d0&
                   -u_c_t(k,j-1,i+2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,1)*2.d0/3.d0-u_c_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI13_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k,j-2,i-1,2)/12.d0&
                   -u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k,j-2,i+1,2)/12.d0&
                   -u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k,j-2,i+2,2)/12.d0&
                   -u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_c/h2_c
             ! second set equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI13_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,1)/12.d0-u_c_t(k,j-1,i-2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,1)*2.d0/3.d0-u_c_t(k,j+2,i-2,1)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI13_c(k,j,i-1)*(u_c_t(k,j-2,i-1,1)/12.d0&
                   -u_c_t(k,j-1,i-1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,1)*2.d0/3.d0-u_c_t(k,j+2,i-1,1)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI13_c(k,j,i+1)*(u_c_t(k,j-2,i+1,1)/12.d0&
                   -u_c_t(k,j-1,i+1,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,1)*2.d0/3.d0-u_c_t(k,j+2,i+1,1)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI13_c(k,j,i+2)*(u_c_t(k,j-2,i+2,1)/12.d0&
                   -u_c_t(k,j-1,i+2,1)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,1)*2.d0/3.d0-u_c_t(k,j+2,i+2,1)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*(2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))&
                   *XI23_c(k,j,i-2)*(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*(2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)&
                   *(u_c_t(k,j-2,i-1,2)/12.d0-u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*(2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)&
                   *(u_c_t(k,j-2,i+1,2)/12.d0-u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*(2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)&
                   *(u_c_t(k,j-2,i+2,2)/12.d0-u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h3_c/h2_c&
                   +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI33_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,3)/12.d0-u_c_t(k,j-1,i-2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,3)*2.d0/3.d0-u_c_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k,j-2,i-1,3)/12.d0&
                   -u_c_t(k,j-1,i-1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,3)*2.d0/3.d0-u_c_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI33_c(k,j,i+1)*(u_c_t(k,j-2,i+1,3)/12.d0&
                   -u_c_t(k,j-1,i+1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,3)*2.d0/3.d0-u_c_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k,j-2,i+2,3)/12.d0&
                   -u_c_t(k,j-1,i+2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,3)*2.d0/3.d0-u_c_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h3_c/h2_c
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+(Jacobian_c(k,j,i-2)*lambda_c(k,j,i-2)*XI33_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,2)/12.d0-u_c_t(k,j-1,i-2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,2)*2.d0/3.d0-u_c_t(k,j+2,i-2,2)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*lambda_c(k,j,i-1)*XI33_c(k,j,i-1)*(u_c_t(k,j-2,i-1,2)/12.d0&
                   -u_c_t(k,j-1,i-1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,2)*2.d0/3.d0-u_c_t(k,j+2,i-1,2)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*lambda_c(k,j,i+1)*XI33_c(k,j,i+1)&
                   *(u_c_t(k,j-2,i+1,2)/12.d0-u_c_t(k,j-1,i+1,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,2)*2.d0/3.d0-u_c_t(k,j+2,i+1,2)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*lambda_c(k,j,i+2)*XI33_c(k,j,i+2)*(u_c_t(k,j-2,i+2,2)/12.d0&
                   -u_c_t(k,j-1,i+2,2)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,2)*2.d0/3.d0-u_c_t(k,j+2,i+2,2)/12.d0)/12.d0)/l2/h2_c/h3_c&
                   +(Jacobian_c(k,j,i-2)*mu_c(k,j,i-2)*XI23_c(k,j,i-2)&
                   *(u_c_t(k,j-2,i-2,3)/12.d0-u_c_t(k,j-1,i-2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-2,3)*2.d0/3.d0-u_c_t(k,j+2,i-2,3)/12.d0)/12.d0 &
                   -Jacobian_c(k,j,i-1)*mu_c(k,j,i-1)*XI23_c(k,j,i-1)*(u_c_t(k,j-2,i-1,3)/12.d0&
                   -u_c_t(k,j-1,i-1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i-1,3)*2.d0/3.d0-u_c_t(k,j+2,i-1,3)/12.d0)*2.d0/3.d0 &
                   +Jacobian_c(k,j,i+1)*mu_c(k,j,i+1)*XI23_c(k,j,i+1)*(u_c_t(k,j-2,i+1,3)/12.d0&
                   -u_c_t(k,j-1,i+1,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+1,3)*2.d0/3.d0-u_c_t(k,j+2,i+1,3)/12.d0)*2.d0/3.d0 &
                   -Jacobian_c(k,j,i+2)*mu_c(k,j,i+2)*XI23_c(k,j,i+2)*(u_c_t(k,j-2,i+2,3)/12.d0&
                   -u_c_t(k,j-1,i+2,3)*2.d0/3.d0 &
                   +u_c_t(k,j+1,i+2,3)*2.d0/3.d0-u_c_t(k,j+2,i+2,3)/12.d0)/12.d0)/l2/h2_c/h3_c
          end do
       end do
    end do
    do i = 1,4
       do j = 1,n2_c
          do k = 1,n1_c
             do k1 = 1,6
                ! mixed derivative 31
                ! first set equation
                lh_c(k,j,i,1) = lh_c(k,j,i,1)+(bof(i,k1)*Jacobian_c(k,j,k1)*(2.d0*mu_c(k,j,k1)+lambda_c(k,j,k1))&
                      *XI13_c(k,j,k1)*(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h3_c/h1_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,2)/12.d0-u_c_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,2)*2.d0/3.d0-u_c_t(k+2,j,k1,2)/12.d0))/l1/h3_c/h1_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,3)/12.d0-u_c_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,3)*2.d0/3.d0-u_c_t(k+2,j,k1,3)/12.d0))/l1/h1_c/h3_c

                lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*(2.d0&
                   *mu_c(k,j,n3_c+1-k1)+lambda_c(k,j,n3_c+1-k1))*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h3_c/h1_c&
                   +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                   *XI23_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,2)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,2)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,2)/12.d0))/l1/h3_c/h1_c+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                   *XI33_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,3)/12.d0 &
                   -u_c_t(k-1,j,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,3)*2.d0/3.d0&
                   -u_c_t(k+2,j,n3_c+1-k1,3)/12.d0))/l1/h1_c/h3_c
                ! second set equation
                lh_c(k,j,i,2) = lh_c(k,j,i,2)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,2)/12.d0-u_c_t(k-1,j,k1,2)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,2)*2.d0/3.d0-u_c_t(k+2,j,k1,2)/12.d0))/l1/h1_c/h3_c

                lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI23_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,2)/12.d0))/l1/h1_c/h3_c
                ! third set equation
                lh_c(k,j,i,3) = lh_c(k,j,i,3)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,1)/12.d0-u_c_t(k-1,j,k1,1)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,1)*2.d0/3.d0-u_c_t(k+2,j,k1,1)/12.d0))/l1/h1_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k-2,j,k1,3)/12.d0-u_c_t(k-1,j,k1,3)*2.d0/3.d0 &
                      +u_c_t(k+1,j,k1,3)*2.d0/3.d0-u_c_t(k+2,j,k1,3)/12.d0))/l1/h3_c/h1_c

                lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,1)/12.d0))/l1/h3_c/h1_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                      *XI13_c(k,j,n3_c+1-k1)*(u_c_t(k-2,j,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k-1,j,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k+1,j,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k+2,j,n3_c+1-k1,3)/12.d0))/l1/h3_c/h1_c



                !31
                ! first set equation
                lh_c(k,j,i,1) = lh_c(k,j,i,1)+(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,1)/12.d0-u_c_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,1)*2.d0/3.d0-u_c_t(k,j+2,k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h3_c/h2_c

                lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)*mu_c(k,j,n3_c+1-k1)&
                      *XI23_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h3_c/h2_c
                ! second set equation
                lh_c(k,j,i,2) = lh_c(k,j,i,2)+(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI13_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,1)/12.d0-u_c_t(k,j-1,k1,1)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,1)*2.d0/3.d0-u_c_t(k,j+2,k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*(2.d0*mu_c(k,j,k1)&
                      +lambda_c(k,j,k1))*XI23_c(k,j,k1)*(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h3_c/h2_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,3)/12.d0-u_c_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,3)*2.d0/3.d0-u_c_t(k,j+2,k1,3)/12.d0))/l2/h3_c/h2_c

                lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI13_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,1)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,1)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,1)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,1)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *(2.d0*mu_c(k,j,n3_c+1-k1)+lambda_c(k,j,n3_c+1-k1))*XI23_c(k,j,n3_c+1-k1)&
                      *(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0-u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0&
                      +u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0-u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h3_c/h2_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,3)/12.d0))/l2/h3_c/h2_c
                ! third set equation
                lh_c(k,j,i,3) = lh_c(k,j,i,3)+(bof(i,k1)*Jacobian_c(k,j,k1)*lambda_c(k,j,k1)*XI33_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,2)/12.d0-u_c_t(k,j-1,k1,2)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,2)*2.d0/3.d0-u_c_t(k,j+2,k1,2)/12.d0))/l2/h2_c/h3_c&
                      +(bof(i,k1)*Jacobian_c(k,j,k1)*mu_c(k,j,k1)*XI23_c(k,j,k1)&
                      *(u_c_t(k,j-2,k1,3)/12.d0-u_c_t(k,j-1,k1,3)*2.d0/3.d0 &
                      +u_c_t(k,j+1,k1,3)*2.d0/3.d0-u_c_t(k,j+2,k1,3)/12.d0))/l2/h2_c/h3_c

                lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3)+(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *lambda_c(k,j,n3_c+1-k1)*XI33_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,2)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,2)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,2)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,2)/12.d0))/l2/h2_c/h3_c&
                      +(-bof(i,k1)*Jacobian_c(k,j,n3_c+1-k1)&
                      *mu_c(k,j,n3_c+1-k1)*XI23_c(k,j,n3_c+1-k1)*(u_c_t(k,j-2,n3_c+1-k1,3)/12.d0 &
                      -u_c_t(k,j-1,n3_c+1-k1,3)*2.d0/3.d0+u_c_t(k,j+1,n3_c+1-k1,3)*2.d0/3.d0&
                      -u_c_t(k,j+2,n3_c+1-k1,3)/12.d0))/l2/h2_c/h3_c
             end do
          end do
       end do
    end do
    !
    do i = 7,n3_c-6
       do j = 1,n2_c
          do k = 1,n1_c
             ! second derivative 33
             ! first set equation
             lh_c(k,j,i,1) = lh_c(k,j,i,1)+((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2 &
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2 &
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/6.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/2.d0&
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)**2 &
                   +mu_c(k,j,i-2)*(XI23_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))*5.d0/6.d0 -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2&
                   +XI33_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI23_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI23_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI13_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI23_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0 &
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI23_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/8.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/24.d0&
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
             ! second set of equation
             lh_c(k,j,i,2) = lh_c(k,j,i,2)+((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI23_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)*5.d0/6.d0&
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI23_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI23_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI23_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI23_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/6.d0 &
                   + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/2.d0 + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI23_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI33_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI33_c(k,j,i-1)**2))/6.d0 &
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/2.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI33_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI23_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI33_c(k,j,i+1)**2))/6.d0&
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI33_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
             ! third set equation
             lh_c(k,j,i,3) = lh_c(k,j,i,3)+((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0&
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,1) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,1) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI13_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,1) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI13_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,1) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI13_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI13_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI13_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,1))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)&
                   +lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/8.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0)*u_c_t(k,j,i-2,2) &
                   +(Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/6.d0 &
                   + Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/2.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0)*u_c_t(k,j,i-1,2) &
                   +(-Jacobian_c(k,j,i-2)*(mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI23_c(k,j,i-2)*XI33_c(k,j,i-2)/24.d0 &
                   - Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)*5.d0/6.d0 &
                   - Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/24.d0)*u_c_t(k,j,i,2) &
                   +(Jacobian_c(k,j,i-1)*(mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI23_c(k,j,i-1)*XI33_c(k,j,i-1)/6.d0 &
                   + Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/2.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/2.d0 &
                   + Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/6.d0)*u_c_t(k,j,i+1,2) &
                   +(-Jacobian_c(k,j,i)*(mu_c(k,j,i)+lambda_c(k,j,i))*XI23_c(k,j,i)*XI33_c(k,j,i)/8.d0 &
                   + Jacobian_c(k,j,i+1)*(mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI23_c(k,j,i+1)*XI33_c(k,j,i+1)/6.d0 &
                   - Jacobian_c(k,j,i+2)*(mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI23_c(k,j,i+2)*XI33_c(k,j,i+2)/8.d0)*u_c_t(k,j,i+2,2))/h3_c**2&
                   +((-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/8.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/6.d0&
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/8.d0)*u_c_t(k,j,i-2,3) &
                   +(Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/6.d0 + Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/2.d0&
                   + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI33_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/6.d0)*u_c_t(k,j,i-1,3) &
                   +(-Jacobian_c(k,j,i-2)*((2.d0*mu_c(k,j,i-2)+lambda_c(k,j,i-2))*XI33_c(k,j,i-2)**2&
                   +mu_c(k,j,i-2)*(XI13_c(k,j,i-2)**2+XI23_c(k,j,i-2)**2))/24.d0 - Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)&
                   +lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2+mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))*5.d0/6.d0&
                   - Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))*3.d0/4.d0 &
                   - Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)+lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2&
                   +mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))*5.d0/6.d0 &
                   -Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2&
                   +XI23_c(k,j,i+2)**2))/24.d0)*u_c_t(k,j,i,3) &
                   +(Jacobian_c(k,j,i-1)*((2.d0*mu_c(k,j,i-1)+lambda_c(k,j,i-1))*XI33_c(k,j,i-1)**2&
                   +mu_c(k,j,i-1)*(XI13_c(k,j,i-1)**2+XI23_c(k,j,i-1)**2))/6.d0 + Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)&
                   +lambda_c(k,j,i))*XI33_c(k,j,i)**2+mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/2.d0 &
                   + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/2.d0 &
                   + Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)+lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2&
                   +mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2+XI23_c(k,j,i+2)**2))/6.d0)*u_c_t(k,j,i+1,3) &
                   +(-Jacobian_c(k,j,i)*((2.d0*mu_c(k,j,i)+lambda_c(k,j,i))*XI33_c(k,j,i)**2&
                   +mu_c(k,j,i)*(XI13_c(k,j,i)**2+XI23_c(k,j,i)**2))/8.d0 + Jacobian_c(k,j,i+1)*((2.d0*mu_c(k,j,i+1)&
                   +lambda_c(k,j,i+1))*XI33_c(k,j,i+1)**2+mu_c(k,j,i+1)*(XI13_c(k,j,i+1)**2+XI23_c(k,j,i+1)**2))/6.d0 &
                   - Jacobian_c(k,j,i+2)*((2.d0*mu_c(k,j,i+2)&
                   +lambda_c(k,j,i+2))*XI33_c(k,j,i+2)**2+mu_c(k,j,i+2)*(XI13_c(k,j,i+2)**2&
                   +XI23_c(k,j,i+2)**2))/8.d0)*u_c_t(k,j,i+2,3))/h3_c**2
          end do
       end do
    end do
    do j = 1,n2_c
       do k = 1,n1_c
          do i = 1,6
             do k1 = 1,8
                do m = 1,8
                   ! second derivative 33
                   ! first set equation
                   lh_c(k,j,i,1) = lh_c(k,j,i,1) +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)**2+mu_c(k,j,m)*(XI23_c(k,j,m)**2+XI33_c(k,j,m)**2))*u_c_t(k,j,k1,1))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)*XI23_c(k,j,m)*u_c_t(k,j,k1,2))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                     +lambda_c(k,j,m))*XI13_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,1) = lh_c(k,j,n3_c+1-i,1) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI23_c(k,j,n3_c+1-m)**2&
                     +XI33_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                     +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                     +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                   ! second set equation
                   lh_c(k,j,i,2) = lh_c(k,j,i,2) +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI13_c(k,j,m)*XI23_c(k,j,m)*u_c_t(k,j,k1,1))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI23_c(k,j,m)**2+mu_c(k,j,m)*(XI13_c(k,j,m)**2&
                       +XI33_c(k,j,m)**2))*u_c_t(k,j,k1,2))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                       +lambda_c(k,j,m))*XI23_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,2) = lh_c(k,j,n3_c+1-i,2) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI23_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI13_c(k,j,n3_c+1-m)**2&
                       +XI33_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                       +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                       +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                   ! third set equation
                   lh_c(k,j,i,3) = lh_c(k,j,i,3) +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI13_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,1))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,m)*(mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI23_c(k,j,m)*XI33_c(k,j,m)*u_c_t(k,j,k1,2))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,m)*((2.d0*mu_c(k,j,m)&
                        +lambda_c(k,j,m))*XI33_c(k,j,m)**2+mu_c(k,j,m)*(XI13_c(k,j,m)**2&
                        +XI23_c(k,j,m)**2))*u_c_t(k,j,k1,3))/h3_c**2

                   lh_c(k,j,n3_c+1-i,3) = lh_c(k,j,n3_c+1-i,3) +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI13_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,1))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*(mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI23_c(k,j,n3_c+1-m)*XI33_c(k,j,n3_c+1-m)*u_c_t(k,j,n3_c+1-k1,2))/h3_c**2&
                        +(acof(i,k1,m)*Jacobian_c(k,j,n3_c+1-m)*((2.d0*mu_c(k,j,n3_c+1-m)&
                        +lambda_c(k,j,n3_c+1-m))*XI33_c(k,j,n3_c+1-m)**2+mu_c(k,j,n3_c+1-m)*(XI13_c(k,j,n3_c+1-m)**2&
                        +XI23_c(k,j,n3_c+1-m)**2))*u_c_t(k,j,n3_c+1-k1,3))/h3_c**2
                end do
             end do
          end do
          ! first set equation
          lh_c(k,j,1,1) = lh_c(k,j,1,1) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)**2+mu_c(k,j,1)*(XI23_c(k,j,1)**2+XI33_c(k,j,1)**2)))/h3_c**2&
                + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)*XI23_c(k,j,1))/h3_c**2&
                + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
                +lambda_c(k,j,1))*XI13_c(k,j,1)*XI33_c(k,j,1))/h3_c**2

          lh_c(k,j,n3_c,1) = lh_c(k,j,n3_c,1) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI23_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
                + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
                + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
                +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! second set equation
          lh_c(k,j,1,2) = lh_c(k,j,1,2) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI13_c(k,j,1)*XI23_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)**2+mu_c(k,j,1)*(XI13_c(k,j,1)**2+XI33_c(k,j,1)**2)))/h3_c**2&
               + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)*XI33_c(k,j,1))/h3_c**2

          lh_c(k,j,n3_c,2) = lh_c(k,j,n3_c,2) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI23_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI33_c(k,j,n3_c)**2)))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2
          ! third set equation
          lh_c(k,j,1,3) = lh_c(k,j,1,3) + (u_c_t(k,j,0,1)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI13_c(k,j,1)*XI33_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,2)*ghcof(1)*Jacobian_c(k,j,1)*(mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI23_c(k,j,1)*XI33_c(k,j,1))/h3_c**2&
               + (u_c_t(k,j,0,3)*ghcof(1)*Jacobian_c(k,j,1)*((2.d0*mu_c(k,j,1)&
               +lambda_c(k,j,1))*XI33_c(k,j,1)**2+mu_c(k,j,1)*(XI13_c(k,j,1)**2+XI23_c(k,j,1)**2)))/h3_c**2

          lh_c(k,j,n3_c,3) = lh_c(k,j,n3_c,3) + (u_c_t(k,j,n3_c+1,1)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI13_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,2)*ghcof(1)*Jacobian_c(k,j,n3_c)*(mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI23_c(k,j,n3_c)*XI33_c(k,j,n3_c))/h3_c**2&
               + (u_c_t(k,j,n3_c+1,3)*ghcof(1)*Jacobian_c(k,j,n3_c)*((2.d0*mu_c(k,j,n3_c)&
               +lambda_c(k,j,n3_c))*XI33_c(k,j,n3_c)**2+mu_c(k,j,n3_c)*(XI13_c(k,j,n3_c)**2+XI23_c(k,j,n3_c)**2)))/h3_c**2
       end do
    end do
  !
  end subroutine Update_interior

  subroutine Update_gp(index)
    integer index
    ! Update ghost point values on the left and right domain
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,0
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
          do k = n1_c+1,n1_c+nrg
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do
    ! Update ghost point values on the front and back domain
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,0
          do k = 1,n1_c
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
       do j = n2_c+1,n2_c+nrg
          do k = 1,n1_c
             call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do
  end subroutine Update_gp

  !
  subroutine Update_traction(index)
   ! traction B.C. on the top of the fine mesh

   integer :: index
   real(dp) :: mat_det
   !
   do j = 1,n2_c
      do i=1,n1_c
         call top_normal_data(Xgrid_c_1(i),Xgrid_c_2(j),Xgrid_c_3(i,j,n3_c), &
             l1,l2,tv,traction_data(i,j,:))
      end do
   end do
   !
   do j = 1,n2_c
      do i = 1,n1_c
         traction_rhs = 0.d0
         ! 33
         ! first set equation
         do k = 1,4
            ! first component
            traction_rhs(1) = traction_rhs(1) &
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,1,index)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)&
               +lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,2,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
               +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,3,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
               +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)
         end do
         ! 31 & 32
         ! first component
         traction_rhs(1) = traction_rhs(1) + Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/l1&
            *XI13_c(i,j,n3_c)*(u_c(i-2,j,n3_c,1,index)/12.d0-u_c(i-1,j,n3_c,1,index)*2.d0/3.d0&
            +u_c(i+1,j,n3_c,1,index)*2.d0/3.d0-u_c(i+2,j,n3_c,1,index)/12.d0)/h1_c + Jacobian_c(i,j,n3_c)&
            *mu_c(i,j,n3_c)/l2*XI23_c(i,j,n3_c)*(u_c(i,j-2,n3_c,1,index)/12.d0&
            -u_c(i,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i,j+1,n3_c,1,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,1,index)/12.d0)/h2_c&
            + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l1*XI23_c(i,j,n3_c)&
            *(u_c(i-2,j,n3_c,2,index)/12.d0-u_c(i-1,j,n3_c,2,index)*2.d0/3.d0+u_c(i+1,j,n3_c,2,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,2,index)/12.d0)/h1_c + Jacobian_c(i,j,n3_c)*lambda_c(i,j,n3_c)/l2*XI13_c(i,j,n3_c)&
            *(u_c(i,j-2,n3_c,2,index)/12.d0-u_c(i,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i,j+1,n3_c,2,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,2,index)/12.d0)/h2_c&
            + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l1*XI33_c(i,j,n3_c)&
            *(u_c(i-2,j,n3_c,3,index)/12.d0-u_c(i-1,j,n3_c,3,index)*2.d0/3.d0+u_c(i+1,j,n3_c,3,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,3,index)/12.d0)/h1_c
         ! scale
         traction_rhs(1) = traction_rhs(1) &
           - sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)*Jacobian_c(i,j,n3_c)*traction_data(i,j,1)
         ! 33
         ! second set equation
         do k = 1,4
            ! first component
            traction_rhs(2) = traction_rhs(2) &
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,1,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
               *XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,2,index)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
               *XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,3,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
               *XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)
         end do
         ! 31 & 32
         ! first component
         traction_rhs(2) = traction_rhs(2) + Jacobian_c(i,j,n3_c)*lambda_c(i,j,n3_c)/l1&
            *XI23_c(i,j,n3_c)*(u_c(i-2,j,n3_c,1,index)/12.d0&
            -u_c(i-1,j,n3_c,1,index)*2.d0/3.d0+u_c(i+1,j,n3_c,1,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,1,index)/12.d0)/h1_c + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l2*XI13_c(i,j,n3_c)&
            *(u_c(i,j-2,n3_c,1,index)/12.d0-u_c(i,j-1,n3_c,1,index)*2.d0/3.d0+u_c(i,j+1,n3_c,1,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,1,index)/12.d0)/h2_c&
            + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l1*XI13_c(i,j,n3_c)&
            *(u_c(i-2,j,n3_c,2,index)/12.d0&
            -u_c(i-1,j,n3_c,2,index)*2.d0/3.d0+u_c(i+1,j,n3_c,2,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,2,index)/12.d0)/h1_c + Jacobian_c(i,j,n3_c)*(2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/l2&
            *XI23_c(i,j,n3_c)*(u_c(i,j-2,n3_c,2,index)/12.d0-u_c(i,j-1,n3_c,2,index)*2.d0/3.d0&
            +u_c(i,j+1,n3_c,2,index)*2.d0/3.d0-u_c(i,j+2,n3_c,2,index)/12.d0)/h2_c&
            + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l2*XI33_c(i,j,n3_c)&
            *(u_c(i,j-2,n3_c,3,index)/12.d0-u_c(i,j-1,n3_c,3,index)*2.d0/3.d0+u_c(i,j+1,n3_c,3,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,3,index)/12.d0)/h2_c
         ! scale
         traction_rhs(2) = traction_rhs(2) &
           - sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)*Jacobian_c(i,j,n3_c)*traction_data(i,j,2)
         ! 33
         ! third set equation
         do k = 1,4
            ! first component
            traction_rhs(3) = traction_rhs(3) &
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,1,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
               +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,2,index)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
               +mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)&
               - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,3,index)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)&
               +lambda_c(i,j,n3_c))*XI33_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))
         end do
         ! 31 & 32
         ! first component
         traction_rhs(3) = traction_rhs(3) + Jacobian_c(i,j,n3_c)*lambda_c(i,j,n3_c)/l1&
            *XI33_c(i,j,n3_c)*(u_c(i-2,j,n3_c,1,index)/12.d0&
            -u_c(i-1,j,n3_c,1,index)*2.d0/3.d0+u_c(i+1,j,n3_c,1,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,1,index)/12.d0)/h1_c&
            + Jacobian_c(i,j,n3_c)*lambda_c(i,j,n3_c)/l2&
            *XI33_c(i,j,n3_c)*(u_c(i,j-2,n3_c,2,index)/12.d0&
            -u_c(i,j-1,n3_c,2,index)*2.d0/3.d0+u_c(i,j+1,n3_c,2,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,2,index)/12.d0)/h2_c&
            + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l1&
            *XI13_c(i,j,n3_c)*(u_c(i-2,j,n3_c,3,index)/12.d0&
            -u_c(i-1,j,n3_c,3,index)*2.d0/3.d0+u_c(i+1,j,n3_c,3,index)*2.d0/3.d0&
            -u_c(i+2,j,n3_c,3,index)/12.d0)/h1_c + Jacobian_c(i,j,n3_c)*mu_c(i,j,n3_c)/l2&
            *XI23_c(i,j,n3_c)*(u_c(i,j-2,n3_c,3,index)/12.d0&
            -u_c(i,j-1,n3_c,3,index)*2.d0/3.d0+u_c(i,j+1,n3_c,3,index)*2.d0/3.d0&
            -u_c(i,j+2,n3_c,3,index)/12.d0)/h2_c
         ! scale
         traction_rhs(3) = traction_rhs(3) &
           - sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)*Jacobian_c(i,j,n3_c)*traction_data(i,j,3)

         ! Update ghost point at the traction boundary
         ! compute the determinant of the coefficient matrix
         mat_det = Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
                 *(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
                 *XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)&
                 *((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI33_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
                 *(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))+Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
                 *XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
                 +mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
                 +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)+Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
                 +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
                 +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
                 *XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
                 *XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)&
                 +lambda_c(i,j,n3_c))*XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))&
                 *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c) &
                 - Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
                 *(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
                 *XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)&
                 *XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI33_c(i,j,n3_c)**2&
                 +mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
                 +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
                 *XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)
         ! update the first component
         u_c(i,j,n3_c+1,1,index) = h3_c/Sb(0)/mat_det* &
            ((Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
            *(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
            *XI33_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))-Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c))*traction_rhs(1) &
            +(Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)&
            *XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI33_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
            *(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2)))*traction_rhs(2)+(Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
            +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
            *XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)&
            *XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
            *XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)))*traction_rhs(3))
         ! update the second component
         u_c(i,j,n3_c+1,2,index) = h3_c/Sb(0)/mat_det* &
            ((Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)&
            *XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI33_c(i,j,n3_c)**2&
            +mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
            +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c))*traction_rhs(1)+(Jacobian_c(i,j,n3_c)&
            *((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
            *(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
            *XI33_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2))-Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)*Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c))*traction_rhs(2) &
            +(Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
            -Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
            *(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))&
            *XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c))*traction_rhs(3))
         ! update the third component
         u_c(i,j,n3_c+1,3,index) = h3_c/Sb(0)/mat_det* &
            ((Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)&
            *XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI23_c(i,j,n3_c)**2&
            +mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)&
            +mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c))*traction_rhs(1)+(Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)&
            *XI33_c(i,j,n3_c)-Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2&
            +mu_c(i,j,n3_c)*(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c))*traction_rhs(2) &
            +(Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))*XI13_c(i,j,n3_c)**2+mu_c(i,j,n3_c)&
            *(XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))&
            *XI23_c(i,j,n3_c)**2+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2))-Jacobian_c(i,j,n3_c)&
            *(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c)&
            *Jacobian_c(i,j,n3_c)*(lambda_c(i,j,n3_c)+mu_c(i,j,n3_c))*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c))*traction_rhs(3))
      end do
   end do
 end subroutine Update_traction

  subroutine Update_Dirichlet_BC(index)
    integer index

    do j=1-nrg,n2_c+nrg
       do k=1-nrg,n1_c+nrg
          i = 1
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          i = n3_c
          !call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,3),tv, &
                   !u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
       end do
    end do

    !!! left and right
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          k = 1
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          k = n1_c
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
       end do
    end do

    !!!! front and back
    do i=1-nrg,n3_c+nrg
       do k=1-nrg,n1_c+nrg
          j = 1
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          j = n2_c
          call exact_solution(Xgrid_c_1(k),Xgrid_c_2(j),Xgrid_c_3(k,j,i),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
       end do
    end do

  end subroutine Update_Dirichlet_BC

end program Elastic_3D_sbp4

subroutine print_array_to_file(n1,n2,n3,A,file_name)
  use problemsetup_new_3d, only: dp
  integer n1,n2,n3
  real(dp), dimension (1:n1,1:n2,1:n3) :: A
  character(len=*) :: file_name
  open (unit = 7, file = file_name)
  do i = 1,n3
     do j = 1,n2
        do k = 1,n1
           write(7,"(ES15.5E3)") A(k,j,i)
        end do
     end do
  end do
  close(7)
end  subroutine print_array_to_file
