program energy_test_cartesian

  use problemsetup_new_3d
  use SBP_operator

  implicit none

  integer i,j,k,l,m,i1,j1,time_index,iset,jj

  ! Reference domain R, Physical domain X
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) ::  Xgrid_c,Rgrid_c

  ! forcing function and spatial discretization
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) :: force_c,force_tt_c,lh_c

  ! Metric derivatives
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: X1R1_c,X1R2_c,X1R3_c,X2R1_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI11_c,XI21_c,XI31_c,XI12_c,XI22_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c

  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: Jacobian_c_3

  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N11_c,N12_c,N13_c,N21_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N22_c,N23_c,N31_c,N32_c,N33_c

  ! material parameters
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: rho_c

  ! Difference operators
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: G1_c,G2_c,G3_c,D12_c,D13_c,D21_c,D23_c,D31_c,D32_c
  real(dp), dimension (-1:n1_c+2,1:n2_c,1:n3_c) :: D3_1c,D2_1c
  real(dp), dimension (1:n1_c,-1:n2_c+2,1:n3_c) :: D3_2c,D1_2c
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: D1_3c,D2_3c
  !! interface and top boundary
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1:4) :: int_temp_c,int_temp_c_1,int_temp_c_2

  ! Solution
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4) :: u_c, u_c_2, uold
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c,1:dim) :: u_c_n

  ! Error
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: err_c
  real(dp), dimension (1:3) :: N6
  real(dp) :: l2_err = 0.d0

  ! discrete energy
  real(dp), dimension (1:10000) :: energy_discrete
  real(dp), dimension (1:10000) :: times
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c,3) :: Lh_new
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c,3) :: Lh_old
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c,3) :: Lh_total


  ! traction B.C. on the top
  real(dp), dimension (1:n1_c,1:n2_c,1:dim):: traction_data
  real(dp), dimension (1:dim):: traction_rhs = 0.d0

  ! time
  real(dp):: tv
  integer :: cr,c1,c2,nt
  real(dp) :: rate
  real(dp):: cfl, dt0, dt

  ! stencils
  !real(dp):: int_cof
  real(dp), dimension (-1:2) :: P
  real(dp):: uxx_cof(-2:2), ux_cof(-2:2)
  real(dp):: acof(6,8,8), ghcof(6), bof(4,6)
  real(dp):: ghcof_no_gp(6), sbop_no_gp(0:5);
  real(dp):: acof_no_gp(6,8,8);
  real(dp), dimension (0:4):: Sb
  real(dp), dimension (-4:2):: Rop
  real(dp), dimension (-3:3):: RPop

  ! Difference stencils and interpolations
  call interpolation_restriction(P,Rop,RPop)
  call central_difference(ux_cof,uxx_cof)
  call varcoeff_NoGhost(acof_no_gp, ghcof_no_gp, sbop_no_gp )
  call VARCOEFFS4( acof, ghcof, Sb )
  call dx_46(bof)

  ! metric derivatives and coefficients of the equation
  call  equation_cof(Xgrid_c,X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c, &
       XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c,Jacobian_c_3,rho_c, &
       N11_c,N12_c,N13_c,N21_c,N22_c,N23_c,N31_c,N32_c,N33_c)

  ! estimate time step
  call kappa(tv)
  cfl=sqrt(3.d0/1.40/tv)
  dt0 = min(h1_c,h2_c,h3_c)*cfl*0.9
  nt = ceiling(tn/dt0)+1
  dt = tn/nt

  print *, nt
  ! exact solution at t = -dt, 0, dt
  do k=1-nrg,n3_c+nrg
     do i=1-nrg,n2_c+nrg
        do j=1-nrg,n1_c+nrg
           call initial_solution(l1,l2,l3,Xgrid_c(j,i,k,1),Xgrid_c(j,i,k,2),Xgrid_c(j,i,k,3), &
                u_c(j,i,k,1,1),u_c(j,i,k,2,1),u_c(j,i,k,3,1))
           call initial_solution(l1,l2,l3,Xgrid_c(j,i,k,1),Xgrid_c(j,i,k,2),Xgrid_c(j,i,k,3), &
                u_c(j,i,k,1,2),u_c(j,i,k,2,2),u_c(j,i,k,3,2))
        end do
     end do
  end do

  CALL system_clock(count_rate=cr)
  CALL SYSTEM_CLOCK(c1)

  rate = real(cr)
  tv = 0.d0
  ! make sure the initial condition satisfy the free surface boundarycondition
  ! Update ghost point values for the traction boundary
  !call Update_traction(1)
  !call Update_traction(2)

  ! time stepping
  do time_index = 1, nt
     ! Predictor step

     ! Evaluate the difference operators
     call Update_interior(u_c(:,:,:,:,2))

     ! Update the solution in the predictor step
     u_c(1:n1_c,1:n2_c,1:n3_c,:,3) = 2.d0*u_c(1:n1_c,1:n2_c,1:n3_c,:,2) - u_c(1:n1_c,1:n2_c,1:n3_c,:,1) &
          + dt**2*lh_c(1:n1_c,1:n2_c,1:n3_c,:)/rho_c(1:n1_c,1:n2_c,1:n3_c,:)

     ! Update time
     tv = tv + dt
     times(time_index) = tv

     ! Update ghost points outside the left and right boundaries
     ! The argument '3' means time level star
     ! '1', '2' and '4' mean time level n-1, n and n+1, respectively.
     call  Update_gp(3)

     ! Update ghost point values for the traction boundary
     !call  Update_traction(3)

     ! Update Dirichlet boundary condition
     call  Update_Dirichlet_BC(3)

     ! Corrector step

     ! Evaluate the difference operators
     call Update_interior((u_c(:,:,:,:,3) - 2.d0*u_c(:,:,:,:,2) + u_c(:,:,:,:,1))/dt**2)

     ! Evaluate the difference operators
     u_c(1:n1_c,1:n2_c,1:n3_c,:,4) = u_c(1:n1_c,1:n2_c,1:n3_c,:,3) &
          + dt**4/12.d0*lh_c(1:n1_c,1:n2_c,1:n3_c,:)/rho_c(1:n1_c,1:n2_c,1:n3_c,:)


     ! Update ghost point values outside the left and right boundaries
     call  Update_gp(4)

     ! Update ghost point values for the traction boundary
     !call  Update_traction(4)

     ! Update Dirichlet boundary condition
     call  Update_Dirichlet_BC(4)

     ! Update solutions
     u_c(:,:,:,:,1) = u_c(:,:,:,:,2)
     u_c(:,:,:,:,2) = u_c(:,:,:,:,4)

     ! compute the discrete energy
     call discrete_energy(u_c(:,:,:,:,1),u_c(:,:,:,:,2),dt,energy_discrete(time_index))
  end do ! end of time loop

  N6(1) = n1_c
  N6(2) = n2_c
  N6(3) = n3_c

  call print_array_to_file(1,1,3,N6,'N6.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,1),'X1c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,2),'X2c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,3),'X3c.txt')

  call print_array_to_file(1,1,nt,energy_discrete,'energy.txt')
  call print_array_to_file(1,1,nt,times,'times.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'u1_1.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'u2_1.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'u3_1.txt')

  CALL SYSTEM_CLOCK(c2)


contains

  subroutine equation_cof(Xgrid_c,X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c, &
       XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c,Jacobian_c_3,rho_c, &
       N11_c,N12_c,N13_c,N21_c,N22_c,N23_c,N31_c,N32_c,N33_c)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N11_c,N12_c,N13_c,N21_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N22_c,N23_c,N31_c,N32_c,N33_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: mu_c,lambda_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) :: Xgrid_c,Rgrid_c,Jacobian_c_3,rho_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)::X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c, &
         X3R1_c,X3R2_c,X3R3_c,XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c


    call generate_grid(Rgrid_c(:,:,:,1),Rgrid_c(:,:,:,2),Rgrid_c(:,:,:,3), &
         Xgrid_c(:,:,:,1),Xgrid_c(:,:,:,2),Xgrid_c(:,:,:,3))

    ! compute metric derivatives, coarse
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             call  metric_derivative(Rgrid_c(k,j,i,1),Rgrid_c(k,j,i,2),Rgrid_c(k,j,i,3),X1R1_c(k,j,i),X1R2_c(k,j,i),X1R3_c(k,j,i), &
                 X2R1_c(k,j,i),X2R2_c(k,j,i),X2R3_c(k,j,i),X3R1_c(k,j,i),X3R2_c(k,j,i),X3R3_c(k,j,i), &
                 XI11_c(k,j,i),XI21_c(k,j,i),XI31_c(k,j,i),XI12_c(k,j,i),XI22_c(k,j,i),XI32_c(k,j,i), &
                 XI13_c(k,j,i),XI23_c(k,j,i),XI33_c(k,j,i),Jacobian_c(k,j,i))
          end do
       end do
    end do

    ! variable coefficients
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             mu_c(k,j,i) = 3.d0 + sin(3.d0*Xgrid_c(k,j,i,1)+0.1d0)*sin(3.d0*Xgrid_c(k,j,i,2)+0.1d0)*sin(Xgrid_c(k,j,i,3))
             lambda_c(k,j,i) = 21.d0+ cos(Xgrid_c(k,j,i,1)+0.1d0)*cos(Xgrid_c(k,j,i,2)+0.1d0)*sin(3.d0*Xgrid_c(k,j,i,3))**2
             rho_c(k,j,i,:) = 2.d0 + sin(Xgrid_c(k,j,i,1)+0.3d0)*sin(Xgrid_c(k,j,i,2)+0.3d0)*sin(Xgrid_c(k,j,i,3)-0.2d0)
          end do
       end do
    end do

    Jacobian_c_3(:,:,:,1) = Jacobian_c
    Jacobian_c_3(:,:,:,2) = Jacobian_c
    Jacobian_c_3(:,:,:,3) = Jacobian_c

    rho_c = rho_c*Jacobian_c_3


    ! coarse
    N11_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI11_c*XI11_c+mu_c*(XI21_c*XI21_c+XI31_c*XI31_c))
    N11_c(:,:,:,1,2) = Jacobian_c*XI11_c*XI21_c*(mu_c+lambda_c)
    N11_c(:,:,:,1,3) = Jacobian_c*XI11_c*XI31_c*(mu_c+lambda_c)
    N11_c(:,:,:,2,1) = N11_c(:,:,:,1,2)
    N11_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI21_c*XI21_c+mu_c*(XI11_c*XI11_c+XI31_c*XI31_c))
    N11_c(:,:,:,2,3) = Jacobian_c*XI21_c*XI31_c*(mu_c+lambda_c)
    N11_c(:,:,:,3,1) = N11_c(:,:,:,1,3)
    N11_c(:,:,:,3,2) = N11_c(:,:,:,2,3)
    N11_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI31_c*XI31_c+mu_c*(XI11_c*XI11_c+XI21_c*XI21_c))

    N22_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI12_c*XI12_c+mu_c*XI22_c*XI22_c+mu_c*XI32_c*XI32_c)
    N22_c(:,:,:,1,2) = Jacobian_c*XI12_c*XI22_c*(mu_c+lambda_c)
    N22_c(:,:,:,1,3) = Jacobian_c*XI12_c*XI32_c*(mu_c+lambda_c)
    N22_c(:,:,:,2,1) = N22_c(:,:,:,1,2)
    N22_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI22_c*XI22_c+mu_c*XI12_c*XI12_c+mu_c*XI32_c*XI32_c)
    N22_c(:,:,:,2,3) = Jacobian_c*XI22_c*XI32_c*(mu_c+lambda_c)
    N22_c(:,:,:,3,1) = N22_c(:,:,:,1,3)
    N22_c(:,:,:,3,2) = N22_c(:,:,:,2,3)
    N22_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI32_c*XI32_c+mu_c*XI12_c*XI12_c+mu_c*XI22_c*XI22_c)

    N33_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI13_c*XI13_c+mu_c*(XI23_c*XI23_c+XI33_c*XI33_c))
    N33_c(:,:,:,1,2) = Jacobian_c*XI13_c*XI23_c*(mu_c+lambda_c)
    N33_c(:,:,:,1,3) = Jacobian_c*XI13_c*XI33_c*(mu_c+lambda_c)
    N33_c(:,:,:,2,1) = N33_c(:,:,:,1,2)
    N33_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI23_c*XI23_c+XI13_c*XI13_c*mu_c+XI33_c*XI33_c*mu_c)
    N33_c(:,:,:,2,3) = Jacobian_c*XI23_c*XI33_c*(mu_c+lambda_c)
    N33_c(:,:,:,3,1) = N33_c(:,:,:,1,3)
    N33_c(:,:,:,3,2) = N33_c(:,:,:,2,3)
    N33_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI33_c*XI33_c+mu_c*(XI13_c*XI13_c+XI23_c*XI23_c))

    N12_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI11_c*XI12_c+mu_c*XI21_c*XI22_c+mu_c*XI31_c*XI32_c)
    N12_c(:,:,:,1,2) = Jacobian_c*(XI11_c*XI22_c*lambda_c+XI21_c*XI12_c*mu_c)
    N12_c(:,:,:,1,3) = Jacobian_c*(XI11_c*XI32_c*lambda_c+XI31_c*XI12_c*mu_c)
    N12_c(:,:,:,2,1) = Jacobian_c*(XI11_c*XI22_c*mu_c+XI21_c*XI12_c*lambda_c)
    N12_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI21_c*XI22_c+mu_c*XI11_c*XI12_c+mu_c*XI31_c*XI32_c)
    N12_c(:,:,:,2,3) = Jacobian_c*(XI21_c*XI32_c*lambda_c+XI31_c*XI22_c*mu_c)
    N12_c(:,:,:,3,1) = Jacobian_c*(XI11_c*XI32_c*mu_c+XI31_c*XI12_c*lambda_c)
    N12_c(:,:,:,3,2) = Jacobian_c*(XI21_c*XI32_c*mu_c+XI31_c*XI22_c*lambda_c)
    N12_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI31_c*XI32_c+mu_c*XI11_c*XI12_c+mu_c*XI21_c*XI22_c)

    N21_c(:,:,:,1,1) = N12_c(:,:,:,1,1)
    N21_c(:,:,:,1,2) = N12_c(:,:,:,2,1)
    N21_c(:,:,:,1,3) = N12_c(:,:,:,3,1)
    N21_c(:,:,:,2,1) = N12_c(:,:,:,1,2)
    N21_c(:,:,:,2,2) = N12_c(:,:,:,2,2)
    N21_c(:,:,:,2,3) = N12_c(:,:,:,3,2)
    N21_c(:,:,:,3,1) = N12_c(:,:,:,1,3)
    N21_c(:,:,:,3,2) = N12_c(:,:,:,2,3)
    N21_c(:,:,:,3,3) = N12_c(:,:,:,3,3)

    N13_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI11_c*XI13_c+XI21_c*XI23_c*mu_c+XI31_c*XI33_c*mu_c)
    N13_c(:,:,:,1,2) = Jacobian_c*(XI11_c*XI23_c*lambda_c+XI21_c*XI13_c*mu_c)
    N13_c(:,:,:,1,3) = Jacobian_c*(XI11_c*XI33_c*lambda_c+XI31_c*XI13_c*mu_c)
    N13_c(:,:,:,2,1) = Jacobian_c*(XI11_c*XI23_c*mu_c+XI21_c*XI13_c*lambda_c)
    N13_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI21_c*XI23_c+mu_c*XI11_c*XI13_c+mu_c*XI31_c*XI33_c)
    N13_c(:,:,:,2,3) = Jacobian_c*(XI21_c*XI33_c*lambda_c+XI31_c*XI23_c*mu_c)
    N13_c(:,:,:,3,1) = Jacobian_c*(XI11_c*XI33_c*mu_c+XI31_c*XI13_c*lambda_c)
    N13_c(:,:,:,3,2) = Jacobian_c*(XI21_c*XI33_c*mu_c+XI31_c*XI23_c*lambda_c)
    N13_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI31_c*XI33_c+mu_c*XI11_c*XI13_c+mu_c*XI21_c*XI23_c)

    N31_c(:,:,:,1,1) = N13_c(:,:,:,1,1)
    N31_c(:,:,:,1,2) = N13_c(:,:,:,2,1)
    N31_c(:,:,:,1,3) = N13_c(:,:,:,3,1)
    N31_c(:,:,:,2,1) = N13_c(:,:,:,1,2)
    N31_c(:,:,:,2,2) = N13_c(:,:,:,2,2)
    N31_c(:,:,:,2,3) = N13_c(:,:,:,3,2)
    N31_c(:,:,:,3,1) = N13_c(:,:,:,1,3)
    N31_c(:,:,:,3,2) = N13_c(:,:,:,2,3)
    N31_c(:,:,:,3,3) = N13_c(:,:,:,3,3)

    N23_c(:,:,:,1,1) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI12_c*XI13_c+mu_c*XI22_c*XI23_c+mu_c*XI32_c*XI33_c)
    N23_c(:,:,:,1,2) = Jacobian_c*(XI12_c*XI23_c*lambda_c+XI22_c*XI13_c*mu_c)
    N23_c(:,:,:,1,3) = Jacobian_c*(XI12_c*XI33_c*lambda_c+XI32_c*XI13_c*mu_c)
    N23_c(:,:,:,2,1) = Jacobian_c*(XI12_c*XI23_c*mu_c+XI22_c*XI13_c*lambda_c)
    N23_c(:,:,:,2,2) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI22_c*XI23_c+mu_c*XI12_c*XI13_c+mu_c*XI32_c*XI33_c)
    N23_c(:,:,:,2,3) = Jacobian_c*(XI22_c*XI33_c*lambda_c+XI32_c*XI23_c*mu_c)
    N23_c(:,:,:,3,1) = Jacobian_c*(XI12_c*XI33_c*mu_c+XI32_c*XI13_c*lambda_c)
    N23_c(:,:,:,3,2) = Jacobian_c*(XI22_c*XI33_c*mu_c+XI32_c*XI23_c*lambda_c)
    N23_c(:,:,:,3,3) = Jacobian_c*((2.d0*mu_c+lambda_c)*XI32_c*XI33_c+mu_c*XI12_c*XI13_c+mu_c*XI22_c*XI23_c)

    N32_c(:,:,:,1,1) = N23_c(:,:,:,1,1)
    N32_c(:,:,:,1,2) = N23_c(:,:,:,2,1)
    N32_c(:,:,:,1,3) = N23_c(:,:,:,3,1)
    N32_c(:,:,:,2,1) = N23_c(:,:,:,1,2)
    N32_c(:,:,:,2,2) = N23_c(:,:,:,2,2)
    N32_c(:,:,:,2,3) = N23_c(:,:,:,3,2)
    N32_c(:,:,:,3,1) = N23_c(:,:,:,1,3)
    N32_c(:,:,:,3,2) = N23_c(:,:,:,2,3)
    N32_c(:,:,:,3,3) = N23_c(:,:,:,3,3)

  end subroutine  equation_cof


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
             mat_t(1,1) = (N11_c(i,j,k,1,1)+N11_c(i,j,k,2,2)+N11_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(1,2) = (N12_c(i,j,k,1,1)+N12_c(i,j,k,2,2)+N12_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(1,3) = (N13_c(i,j,k,1,1)+N13_c(i,j,k,2,2)+N13_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(2,1) = (N21_c(i,j,k,1,1)+N21_c(i,j,k,2,2)+N21_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(2,2) = (N22_c(i,j,k,1,1)+N22_c(i,j,k,2,2)+N22_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(2,3) = (N23_c(i,j,k,1,1)+N23_c(i,j,k,2,2)+N23_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(3,1) = (N31_c(i,j,k,1,1)+N31_c(i,j,k,2,2)+N31_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(3,2) = (N32_c(i,j,k,1,1)+N32_c(i,j,k,2,2)+N32_c(i,j,k,3,3))/rho_c(i,j,k,1)
             mat_t(3,3) = (N33_c(i,j,k,1,1)+N33_c(i,j,k,2,2)+N33_c(i,j,k,3,3))/rho_c(i,j,k,1)
             call dsyev('N','U',3,mat_t,3,kappa1,work,8,info)
             if (kappa1(3) > kappa2) then
                kappa2 = kappa1(3)
             end if
          end do
       end do
    end do
  end subroutine kappa

  subroutine Update_interior(u_c_t)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t

    ! Difference operators in the interior of the domains
    lh_c = 0.d0
    do iset = 1,3
      do jj = 1,3
         call FD_r1r1(h1_c,n1_c,n2_c,n3_c,N11_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G1_c)
         call FD_r2r2(h2_c,n1_c,n2_c,n3_c,N22_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G2_c)
         call FD_r3r3(h3_c,n1_c,n2_c,n3_c,N33_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G3_c,acof,ghcof)

         call FD_r3(h3_c,-1,n1_c+2,1,n2_c,1,n3_c,u_c_t(-1:n1_c+2,1:n2_c,1:n3_c,jj),D3_1c,bof)
         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,N13_c(-1:n1_c+2,1:n2_c,1:n3_c,iset,jj)*D3_1c,D31_c)

         call FD_r2(h2_c,-1,n1_c+2,-1,n2_c+2,1,n3_c,u_c_t(-1:n1_c+2,-1:n2_c+2,1:n3_c,jj),D2_1c)
         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,N12_c(-1:n1_c+2,1:n2_c,1:n3_c,iset,jj)*D2_1c,D21_c)

         call FD_r1(h1_c,-1,n1_c+2,-1,n2_c+2,1,n3_c,u_c_t(-1:n1_c+2,-1:n2_c+2,1:n3_c,jj),D1_2c)
         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,N21_c(1:n1_c,-1:n2_c+2,1:n3_c,iset,jj)*D1_2c,D12_c)

         call FD_r3(h3_c,1,n1_c,-1,n2_c+2,1,n3_c,u_c_t(1:n1_c,-1:n2_c+2,1:n3_c,jj),D3_2c,bof)
         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,N23_c(1:n1_c,-1:n2_c+2,1:n3_c,iset,jj)*D3_2c,D32_c)

         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,u_c_t(-1:n1_c+2,1:n2_c,1:n3_c,jj),D1_3c)
         call FD_r3(h3_c,1,n1_c,1,n2_c,1,n3_c,N31_c(1:n1_c,1:n2_c,1:n3_c,iset,jj)*D1_3c,D13_c,bof)

         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,u_c_t(1:n1_c,-1:n2_c+2,1:n3_c,jj),D2_3c)
         call FD_r3(h3_c,1,n1_c,1,n2_c,1,n3_c,N32_c(1:n1_c,1:n2_c,1:n3_c,iset,jj)*D2_3c,D23_c,bof)

         lh_c(1:n1_c,1:n2_c,1:n3_c,iset) = lh_c(1:n1_c,1:n2_c,1:n3_c,iset)+G1_c+G2_c+G3_c+D31_c+D21_c+D12_c+D32_c+D13_c+D23_c

      end do
    end do

  end subroutine Update_interior

  subroutine Update_traction(index)
    ! traction B.C. on the top of the fine mesh

    integer :: index
    real(dp) :: mat_det

    int_temp_c_1 = 0.d0
    int_temp_c_2 = 0.d0

    do j = 1,n2_c
       do i=1,n1_c
          call top_normal_data(traction_data(i,j,:))
       end do
    end do

    do jj = 1,3
       call FD_r1(h1_c,-1,n1_c+2,1,n2_c,n3_c,n3_c,u_c(-1:n1_c+2,1:n2_c,n3_c,jj,index),int_temp_c_1(1:n1_c,1:n2_c,jj))
       call FD_r2(h2_c,1,n1_c,-1,n2_c+2,n3_c,n3_c,u_c(1:n1_c,-1:n2_c+2,n3_c,jj,index),int_temp_c_2(1:n1_c,1:n2_c,jj))
    end do

    do j = 1,n2_c
       do i = 1,n1_c
          traction_rhs = 0.d0
          do iset = 1,3
             do k = 1,4
                do jj = 1,3
                   traction_rhs(iset) = traction_rhs(iset) &
                       - 1.d0/h3_c*Sb(k)*u_c(i,j,n3_c+1-k,jj,index)*N33_c(i,j,n3_c,iset,jj)
                end do
             end do
             do jj = 1,3
                traction_rhs(iset) = traction_rhs(iset) + N31_c(i,j,n3_c,iset,jj)*int_temp_c_1(i,j,jj) &
                                     + N32_c(i,j,n3_c,iset,jj)*int_temp_c_2(i,j,jj)
             end do
             traction_rhs(iset) = traction_rhs(iset) &
               - sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)*Jacobian_c(i,j,n3_c)*traction_data(i,j,iset)
          end do

          ! Update ghost point at the traction boundary
          mat_det = N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,2,2)*N33_c(i,j,n3_c,3,3) &
                  +N33_c(i,j,n3_c,2,1)*N33_c(i,j,n3_c,3,2)*N33_c(i,j,n3_c,1,3) &
                  + N33_c(i,j,n3_c,3,1)*N33_c(i,j,n3_c,1,2)*N33_c(i,j,n3_c,2,3) &
                  -N33_c(i,j,n3_c,1,3)*N33_c(i,j,n3_c,2,2)*N33_c(i,j,n3_c,3,1) &
                  - N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,3,2)*N33_c(i,j,n3_c,2,3) &
                  -N33_c(i,j,n3_c,3,3)*N33_c(i,j,n3_c,2,1)*N33_c(i,j,n3_c,1,2)

          u_c(i,j,n3_c+1,1,index) = h3_c/Sb(0)/mat_det* &
              ((N33_c(i,j,n3_c,2,2)*N33_c(i,j,n3_c,3,3)-N33_c(i,j,n3_c,2,3)*N33_c(i,j,n3_c,3,2))*traction_rhs(1) &
             +(N33_c(i,j,n3_c,1,3)*N33_c(i,j,n3_c,3,2)-N33_c(i,j,n3_c,1,2)*N33_c(i,j,n3_c,3,3))*traction_rhs(2) &
             +(N33_c(i,j,n3_c,1,2)*N33_c(i,j,n3_c,2,3)-N33_c(i,j,n3_c,1,3)*N33_c(i,j,n3_c,2,2))*traction_rhs(3))
          u_c(i,j,n3_c+1,2,index) = h3_c/Sb(0)/mat_det* &
              ((N33_c(i,j,n3_c,2,3)*N33_c(i,j,n3_c,3,1)-N33_c(i,j,n3_c,3,3)*N33_c(i,j,n3_c,2,1))*traction_rhs(1) &
             +(N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,3,3)-N33_c(i,j,n3_c,1,3)*N33_c(i,j,n3_c,3,1))*traction_rhs(2) &
             +(N33_c(i,j,n3_c,1,3)*N33_c(i,j,n3_c,2,1)-N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,2,3))*traction_rhs(3))
          u_c(i,j,n3_c+1,3,index) = h3_c/Sb(0)/mat_det* &
              ((N33_c(i,j,n3_c,2,1)*N33_c(i,j,n3_c,3,2)-N33_c(i,j,n3_c,2,2)*N33_c(i,j,n3_c,3,1))*traction_rhs(1) &
             +(N33_c(i,j,n3_c,1,2)*N33_c(i,j,n3_c,3,1)-N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,3,2))*traction_rhs(2) &
             +(N33_c(i,j,n3_c,1,1)*N33_c(i,j,n3_c,2,2)-N33_c(i,j,n3_c,1,2)*N33_c(i,j,n3_c,2,1))*traction_rhs(3))
       end do
    end do
  end subroutine Update_traction

  subroutine Update_gp(index)
    integer index
    ! Update ghost point values on the left and right domain
    ! coarse fine
    do i=0,n3_c+1
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,0
             u_c(k,j,i,1,index) = 0.d0
             u_c(k,j,i,2,index) = 0.d0
             u_c(k,j,i,3,index) = 0.d0
          end do
          do k = n1_c+1,n1_c+nrg
             u_c(k,j,i,1,index) = 0.d0
             u_c(k,j,i,2,index) = 0.d0
             u_c(k,j,i,3,index) = 0.d0
          end do
       end do
    end do
    ! Update ghost point values on the front and back domain
    ! coarse fine
    do i=0,n3_c+1
       do j = 1-nrg,0
          do k = 1,n1_c
             u_c(k,j,i,1,index) = 0.d0
             u_c(k,j,i,2,index) = 0.d0
             u_c(k,j,i,3,index) = 0.d0
          end do
       end do
       do j = n2_c+1,n2_c+nrg
          do k = 1,n1_c
             u_c(k,j,i,1,index) = 0.d0
             u_c(k,j,i,2,index) = 0.d0
             u_c(k,j,i,3,index) = 0.d0
          end do
       end do
    end do

  end subroutine Update_gp

  subroutine Update_Dirichlet_BC(index)
    integer index

    do j=1,n2_c
       do k=1,n1_c
          i = 1
          u_c(k,j,i,1,index) = 0.d0
          u_c(k,j,i,2,index) = 0.d0
          u_c(k,j,i,3,index) = 0.d0
          i = n3_c
          u_c(k,j,i,1,index) = 0.d0
          u_c(k,j,i,2,index) = 0.d0
          u_c(k,j,i,3,index) = 0.d0
       end do
    end do
  end subroutine Update_Dirichlet_BC
  !
  subroutine Lh_interior(u_c_t)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t

    ! Difference operators in the interior of the domains
    lh_total = 0.d0
    do iset = 1,3
      do jj = 1,3
         call FD_r1r1(h1_c,n1_c,n2_c,n3_c,N11_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G1_c)
         call FD_r2r2(h2_c,n1_c,n2_c,n3_c,N22_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G2_c)
         call FD_r3r3(h3_c,n1_c,n2_c,n3_c,N33_c(:,:,:,iset,jj),u_c_t(:,:,:,jj),1,n1_c,1,n2_c,1,n3_c,G3_c,acof,ghcof)

         call FD_r3(h3_c,-1,n1_c+2,1,n2_c,1,n3_c,u_c_t(-1:n1_c+2,1:n2_c,1:n3_c,jj),D3_1c,bof)
         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,N13_c(-1:n1_c+2,1:n2_c,1:n3_c,iset,jj)*D3_1c,D31_c)

         call FD_r2(h2_c,-1,n1_c+2,-1,n2_c+2,1,n3_c,u_c_t(-1:n1_c+2,-1:n2_c+2,1:n3_c,jj),D2_1c)
         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,N12_c(-1:n1_c+2,1:n2_c,1:n3_c,iset,jj)*D2_1c,D21_c)

         call FD_r1(h1_c,-1,n1_c+2,-1,n2_c+2,1,n3_c,u_c_t(-1:n1_c+2,-1:n2_c+2,1:n3_c,jj),D1_2c)
         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,N21_c(1:n1_c,-1:n2_c+2,1:n3_c,iset,jj)*D1_2c,D12_c)

         call FD_r3(h3_c,1,n1_c,-1,n2_c+2,1,n3_c,u_c_t(1:n1_c,-1:n2_c+2,1:n3_c,jj),D3_2c,bof)
         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,N23_c(1:n1_c,-1:n2_c+2,1:n3_c,iset,jj)*D3_2c,D32_c)

         call FD_r1(h1_c,-1,n1_c+2,1,n2_c,1,n3_c,u_c_t(-1:n1_c+2,1:n2_c,1:n3_c,jj),D1_3c)
         call FD_r3(h3_c,1,n1_c,1,n2_c,1,n3_c,N31_c(1:n1_c,1:n2_c,1:n3_c,iset,jj)*D1_3c,D13_c,bof)

         call FD_r2(h2_c,1,n1_c,-1,n2_c+2,1,n3_c,u_c_t(1:n1_c,-1:n2_c+2,1:n3_c,jj),D2_3c)
         call FD_r3(h3_c,1,n1_c,1,n2_c,1,n3_c,N32_c(1:n1_c,1:n2_c,1:n3_c,iset,jj)*D2_3c,D23_c,bof)

         lh_total(1:n1_c,1:n2_c,1:n3_c,iset) = &
             lh_total(1:n1_c,1:n2_c,1:n3_c,iset)+G1_c+G2_c+G3_c+D31_c+D21_c+D12_c+D32_c+D13_c+D23_c

      end do
    end do

  end subroutine Lh_interior

  ! compute the discrete energy
  subroutine discrete_energy(uc_old,uc_new,dt,energy_num)
    real(dp), dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: uc_old,uc_new
    real(dp), dimension(1:n1_c,1:n2_c,1:n3_c) :: Dc1,Dc2
    real(dp), dimension(1:n1_c,1:n2_c,1:n3_c) :: energy_temp_c,energy_num_temp_c
    real(dp) :: dt,energy_num

    integer :: l1,l2

    ! coarse domain
    energy_num_temp_c = 0.d0
    ! term w.r.t time derivative
    do l1 = 1,3
      energy_num_temp_c = energy_num_temp_c &
        + rho_c(1:n1_c,1:n2_c,1:n3_c,1)* &
        (uc_new(1:n1_c,1:n2_c,1:n3_c,l1)-uc_old(1:n1_c,1:n2_c,1:n3_c,l1))**2/(dt**2)
    end do
    energy_num_temp_c = energy_num_temp_c*Jacobian_c(1:n1_c,1:n2_c,1:n3_c)
    !
    call Lh_interior(uc_old)
    Lh_old = Lh_total
    ! term w.r.t S_h(unew,uold) = -(unew,L_h uold)
    energy_num_temp_c = energy_num_temp_c - uc_new(1:n1_c,1:n2_c,1:n3_c,1)*Lh_old(:,:,:,1) &
       - uc_new(1:n1_c,1:n2_c,1:n3_c,2)*Lh_old(:,:,:,2) &
       - uc_new(1:n1_c,1:n2_c,1:n3_c,3)*Lh_old(:,:,:,3)
    ! term w.r.t (L_h unew,L_h uold)
    call Lh_interior(uc_new)
    Lh_new = Lh_total
    energy_num_temp_c = energy_num_temp_c-dt*dt*(Lh_new(:,:,:,1)*Lh_old(:,:,:,1) &
      +Lh_new(:,:,:,2)*Lh_old(:,:,:,2)+Lh_new(:,:,:,3)*Lh_old(:,:,:,3)) &
      /Jacobian_c(1:n1_c,1:n2_c,1:n3_c)/12.d0
    ! add grid points with corresponding weights
    energy_num = 0.d0
    do k = 5,n3_c-4
      do j = 1,n2_c
        do i = 1,n1_c
          energy_num = energy_num + h1_c*h2_c*h3_c*energy_num_temp_c(i,j,k)
        end do
      end do
    end do
    !
    do j = 1,n2_c
      do i = 1,n1_c
        energy_num = energy_num + h1_c*h2_c*h3_c*17.d0/48.d0*energy_num_temp_c(i,j,1) &
          + h1_c*h2_c*h3_c*17.d0/48.d0*energy_num_temp_c(i,j,n3_c)
      end do
    end do
    !
    do j = 1,n2_c
      do i = 1,n1_c
        energy_num = energy_num + h1_c*h2_c*h3_c*59.d0/48.d0*energy_num_temp_c(i,j,2) &
          + h1_c*h2_c*h3_c*59.d0/48.d0*energy_num_temp_c(i,j,n3_c-1)
      end do
    end do
    !
    do j = 1,n2_c
      do i = 1,n1_c
        energy_num = energy_num + h1_c*h2_c*h3_c*43.d0/48.d0*energy_num_temp_c(i,j,3) &
          + h1_c*h2_c*h3_c*43.d0/48.d0*energy_num_temp_c(i,j,n3_c-2)
      end do
    end do
    !
    do j = 1,n2_c
      do i = 1,n1_c
        energy_num = energy_num + h1_c*h2_c*h3_c*49.d0/48.d0*energy_num_temp_c(i,j,4) &
          + h1_c*h2_c*h3_c*49.d0/48.d0*energy_num_temp_c(i,j,n3_c-3)
      end do
    end do
    !
  end subroutine discrete_energy
end program energy_test_cartesian


subroutine FD_r1r1(h1,n1,n2,n3,N,u,m11,m12,m21,m22,m31,m32,Gu)
  use problemsetup_new_3d, only: dp, nrg
  integer i,j,k
  integer n1,n2,n3
  integer m11,m12,m21,m22,m31,m32 ! dimension of Gu is (m11:m12,m21:m22,m31:m32)
  real(dp), dimension (m11:m12,m21:m22,m31:m32) :: Gu
  real(dp), dimension (1-nrg:n1+nrg,1-nrg:n2+nrg,1-nrg:n3+nrg) :: u, N
  real(dp) :: h1
  Gu = 0.d0
  do i = m31,m32
     do j = m21,m22
        do k = m11,m12
           Gu(k,j,i) = (-N(k-2,j,i)/8.d0 + N(k-1,j,i)/6.d0 - N(k,j,i)/8.d0)*u(k-2,j,i) &
             +(N(k-2,j,i)/6.d0 + N(k-1,j,i)/2.d0 + N(k,j,i)/2.d0 + N(k+1,j,i)/6.d0)*u(k-1,j,i) &
             +(-N(k-2,j,i)/24.d0 - N(k-1,j,i)*5.d0/6.d0 - N(k,j,i)*3.d0/4.d0 &
             - N(k+1,j,i)*5.d0/6.d0 -N(k+2,j,i)/24.d0)*u(k-0,j,i) &
             +(N(k-1,j,i)/6.d0 + N(k,j,i)/2.d0 + N(k+1,j,i)/2.d0 + N(k+2,j,i)/6.d0)*u(k+1,j,i) &
             +(-N(k,j,i)/8.d0 + N(k+1,j,i)/6.d0 - N(k+2,j,i)/8.d0)*u(k+2,j,i)
        end do
     end do
  end do
  Gu = Gu/(h1**2)
end subroutine FD_r1r1

subroutine FD_r2r2(h2,n1,n2,n3,N,u,m11,m12,m21,m22,m31,m32,Gu)
  use problemsetup_new_3d, only: dp, nrg
  integer i,j,k
  integer n1,n2,n3
  integer m11,m12,m21,m22,m31,m32 ! dimension of Gu is (m11:m12,m21:m22,m31:m32)
  real(dp), dimension (m11:m12,m21:m22,m31:m32) :: Gu
  real(dp), dimension (1-nrg:n1+nrg,1-nrg:n2+nrg,1-nrg:n3+nrg) :: u, N
  real(dp) :: h2
  Gu = 0.d0
  do i = m31,m32
     do j = m21,m22
        do k = m11,m12
           Gu(k,j,i) = (-N(k,j-2,i)/8.d0 + N(k,j-1,i)/6.d0 - N(k,j,i)/8.d0)*u(k,j-2,i) &
             +(N(k,j-2,i)/6.d0 + N(k,j-1,i)/2.d0 + N(k,j,i)/2.d0 + N(k,j+1,i)/6.d0)*u(k,j-1,i) &
             +(-N(k,j-2,i)/24.d0 - N(k,j-1,i)*5.d0/6.d0 - N(k,j,i)*3.d0/4.d0 &
             - N(k,j+1,i)*5.d0/6.d0 -N(k,j+2,i)/24.d0)*u(k,j-0,i) &
             +(N(k,j-1,i)/6.d0 + N(k,j,i)/2.d0 + N(k,j+1,i)/2.d0 + N(k,j+2,i)/6.d0)*u(k,j+1,i) &
             +(-N(k,j,i)/8.d0 + N(k,j+1,i)/6.d0 - N(k,j+2,i)/8.d0)*u(k,j+2,i)
        end do
     end do
  end do
  Gu = Gu/(h2**2)
end subroutine FD_r2r2

subroutine FD_r3r3(h3,n1,n2,n3,N,u,m11,m12,m21,m22,m31,m32,Gu,acof,ghcof)
  use problemsetup_new_3d, only: dp,nrg
  integer i,j,k,k1,m,n1,n2,n3
  integer m11,m12,m21,m22,m31,m32
  real(dp), dimension (m11:m12,m21:m22,m31:m32) :: Gu
  real(dp), dimension (1-nrg:n1+nrg,1-nrg:n2+nrg,1-nrg:n3+nrg) :: u, N
  real(dp):: acof(6,8,8), ghcof(6), h3
  Gu = 0.d0
  do i = m31+6,m32-6
     do j = m21,m22
        do k = m11,m12
           Gu(k,j,i) = (-N(k,j,i-2)/8.d0 + N(k,j,i-1)/6.d0 - N(k,j,i)/8.d0)*u(k,j,i-2) &
             +(N(k,j,i-2)/6.d0 + N(k,j,i-1)/2.d0 + N(k,j,i)/2.d0 + N(k,j,i+1)/6.d0)*u(k,j,i-1) &
             +(-N(k,j,i-2)/24.d0 - N(k,j,i-1)*5.d0/6.d0 - N(k,j,i)*3.d0/4.d0 &
             - N(k,j,i+1)*5.d0/6.d0 -N(k,j,i+2)/24.d0)*u(k,j,i) &
             +(N(k,j,i-1)/6.d0 + N(k,j,i)/2.d0 + N(k,j,i+1)/2.d0 + N(k,j,i+2)/6.d0)*u(k,j,i+1) &
             +(-N(k,j,i)/8.d0 + N(k,j,i+1)/6.d0 - N(k,j,i+2)/8.d0)*u(k,j,i+2)
         end do
     end do
  end do

  do j = m21,m22
     do k = m11,m12
        do i = 1,6
           do k1 = 1,8
              do m = 1,8
                 Gu(k,j,m31-1+i) = Gu(k,j,m31-1+i) +acof(i,k1,m)*N(k,j,m31-1+m)*u(k,j,m31-1+k1)
                 Gu(k,j,m32+1-i) = Gu(k,j,m32+1-i) +acof(i,k1,m)*N(k,j,m32+1-m)*u(k,j,m32+1-k1)
              end do
          end do
        end do
        Gu(k,j,m31) =  Gu(k,j,m31) + u(k,j,0)*ghcof(1)*N(k,j,1)
        Gu(k,j,m32) =  Gu(k,j,m32) + u(k,j,m32+1)*ghcof(1)*N(k,j,m32)
     end do
  end do
  Gu = Gu/(h3**2)
end subroutine FD_r3r3

subroutine FD_r1(h1,n11,n12,n21,n22,n31,n32,u,D1u)
  use problemsetup_new_3d, only: dp
  integer i,j,k
  integer n11,n12,n21,n22,n31,n32 ! Dimension of u
  real(dp), dimension (n11+2:n12-2,n21:n22,n31:n32) :: D1u
  real(dp), dimension (n11:n12,n21:n22,n31:n32) :: u
  real(dp) :: bof(4,6),h1
  D1u = 0.d0
  do i = n31,n32
     do j = n21,n22
        do k = n11+2,n12-2
           D1u(k,j,i) = u(k-2,j,i)/12.d0 - u(k-1,j,i)*2.d0/3.d0 &
             + u(k+1,j,i)*2.d0/3.d0 - u(k+2,j,i)/12.d0
        end do
     end do
  end do
  D1u = D1u/h1
end subroutine FD_r1

subroutine FD_r2(h2,n11,n12,n21,n22,n31,n32,u,D2u)
  use problemsetup_new_3d, only: dp
  integer i,j,k
  integer n11,n12,n21,n22,n31,n32 ! Dimension of u
  real(dp), dimension (n11:n12,n21+2:n22-2,n31:n32) :: D2u
  real(dp), dimension (n11:n12,n21:n22,n31:n32) :: u
  real(dp) :: bof(4,6),h2
  D2u = 0.d0
  do i = n31,n32
     do j = n21+2,n22-2
        do k = n11,n12
           D2u(k,j,i) = u(k,j-2,i)/12.d0 - u(k,j-1,i)*2.d0/3.d0 &
             + u(k,j+1,i)*2.d0/3.d0 - u(k,j+2,i)/12.d0
        end do
     end do
  end do
  D2u = D2u/h2
end subroutine FD_r2

subroutine FD_r3(h3,n11,n12,n21,n22,n31,n32,u,D3u,bof)
  use problemsetup_new_3d, only: dp
  integer i,j,k,k1
  integer n11,n12,n21,n22,n31,n32
  real(dp), dimension (n11:n12,n21:n22,n31:n32) :: D3u,u
  real(dp) :: bof(4,6),h3
  D3u = 0.d0
  do k = n11,n12
     do j = n21,n22
        do i = n31+4,n32-4
           D3u(k,j,i) = u(k,j,i-2)/12.d0 - u(k,j,i-1)*2.d0/3.d0 &
             + u(k,j,i+1)*2.d0/3.d0 - u(k,j,i+2)/12.d0
        end do

        do i = 1,4
           do k1 = 1,6
              D3u(k,j,n31-1+i) = D3u(k,j,n31-1+i) + bof(i,k1)*u(k,j,n31-1+k1)
              D3u(k,j,n32+1-i) = D3u(k,j,n32+1-i) - bof(i,k1)*u(k,j,n32+1-k1) ! Change sign !!!
           end do
        end do
     end do
  end do
  D3u = D3u/h3
end subroutine FD_r3

subroutine FD_r3r3_no_gp(n1,n2,n3,h3,N,u,Gu,acof_no_gp,acof,ghcof)
  use problemsetup_new_3d, only: dp,nrg
  integer i,j,k,k1,m,n1,n2,n3
  real(dp), dimension (1:n1,1:n2,1:n3) :: Gu
  real(dp), dimension (1-nrg:n1+nrg,1-nrg:n2+nrg,1-nrg:n3+nrg) :: u, N
  real(dp):: acof_no_gp(6,8,8),h3
  real(dp):: acof(6,8,8), ghcof(6)
  Gu = 0.d0
  do i = 7,n3-6
     do j = 1,n2
        do k = 1,n1
           Gu(k,j,i) = (-N(k,j,i-2)/8.d0 + N(k,j,i-1)/6.d0 - N(k,j,i)/8.d0)*u(k,j,i-2) &
             +(N(k,j,i-2)/6.d0 + N(k,j,i-1)/2.d0 + N(k,j,i)/2.d0 + N(k,j,i+1)/6.d0)*u(k,j,i-1) &
             +(-N(k,j,i-2)/24.d0 - N(k,j,i-1)*5.d0/6.d0 - N(k,j,i)*3.d0/4.d0 &
             - N(k,j,i+1)*5.d0/6.d0 -N(k,j,i+2)/24.d0)*u(k,j,i) &
             +(N(k,j,i-1)/6.d0 + N(k,j,i)/2.d0 + N(k,j,i+1)/2.d0 + N(k,j,i+2)/6.d0)*u(k,j,i+1) &
             +(-N(k,j,i)/8.d0 + N(k,j,i+1)/6.d0 - N(k,j,i+2)/8.d0)*u(k,j,i+2)
        end do
     end do
  end do

  do j = 1,n2
     do k = 1,n1
        do i = 1,6
           do k1 = 1,8
              do m = 1,8
                 Gu(k,j,i) = Gu(k,j,i) +acof_no_gp(i,k1,m)*N(k,j,m)*u(k,j,k1)
                 Gu(k,j,n3+1-i) = Gu(k,j,n3+1-i) +acof(i,k1,m)*N(k,j,n3+1-m)*u(k,j,n3+1-k1)
              end do
           end do
        end do
        Gu(k,j,n3) = Gu(k,j,n3) + u(k,j,n3+1)*ghcof(1)*N(k,j,n3)
     end do
  end do
  Gu = Gu/(h3**2)
end subroutine FD_r3r3_no_gp

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
