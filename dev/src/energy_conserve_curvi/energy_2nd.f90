program energy_2nd

  use problemsetup_new_3d
  use SBP_operator

  implicit none

  integer i,j,k,l,m,i1,j1,time_index,iset,jj

  ! Reference domain R, Physical domain X
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) ::  Xgrid_c,Rgrid_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim) ::  Xgrid_f,Rgrid_f

  ! forcing function and spatial discretization
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) :: force_c,force_tt_c,lh_c,eta_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim) :: force_f,force_tt_f,lh_f,eta_f

  ! Metric derivatives
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: X1R1_c,X1R2_c,X1R3_c,X2R1_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI11_c,XI21_c,XI31_c,XI12_c,XI22_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: X1R1_f,X1R2_f,X1R3_f,X2R1_f
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: X2R2_f,X2R3_f,X3R1_f,X3R2_f,X3R3_f
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: XI11_f,XI21_f,XI31_f,XI12_f,XI22_f
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: XI32_f,XI13_f,XI23_f,XI33_f,Jacobian_f

  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: Jacobian_c_3
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: Jacobian_f_3

  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N11_c,N12_c,N13_c,N21_c
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N22_c,N23_c,N31_c,N32_c,N33_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:dim) :: N11_f,N12_f,N13_f,N21_f
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:dim) :: N22_f,N23_f,N31_f,N32_f,N33_f

  ! material parameters
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: rho_f
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: rho_c

  ! Difference operators
  real(dp), dimension (1:n1_f,1:n2_f,1:n3_f) :: G1_f,G2_f,G3_f,D12_f,D13_f,D21_f,D23_f,D31_f,D32_f
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: G1_c,G2_c,G3_c,D12_c,D13_c,D21_c,D23_c,D31_c,D32_c
  real(dp), dimension (-1:n1_f+2,1:n2_f,1:n3_f) :: D3_1f,D2_1f
  real(dp), dimension (1:n1_f,-1:n2_f+2,1:n3_f) :: D3_2f,D1_2f
  real(dp), dimension (-1:n1_c+2,1:n2_c,1:n3_c) :: D3_1c,D2_1c
  real(dp), dimension (1:n1_c,-1:n2_c+2,1:n3_c) :: D3_2c,D1_2c
  real(dp), dimension (1:n1_f,1:n2_f,1:n3_f) :: D1_3f,D2_3f
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: D1_3c,D2_3c
  !! interface and top boundary
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1:4) :: int_temp_c
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1:4) :: int_temp_f,int_temp_f_1,int_temp_f_2

  ! Solution
  real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4) :: u_c, u_c_2, uold
  real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4) :: u_f, u_f_2
  real(dp), dimension (1:n1_f,1:n2_f,1:n3_f,1:dim) :: u_f_n
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c,1:dim) :: u_c_n

  !
  real(dp), dimension(1:1000) :: energy_discrete, Times, energy_bdry

  ! Error
  real(dp), dimension (1:n1_f,1:n2_f,1:n3_f) :: err_f
  real(dp), dimension (1:n1_c,1:n2_c,1:n3_c) :: err_c
  real(dp), dimension (1:6) :: N6
  real(dp) :: l2_err = 0.d0


  ! System of linear equations
  real(dp), dimension (1:n1_c*n2_c*3,1:n1_c*n2_c*3) :: Mass
  real(dp), dimension (1:n1_c*n2_c*3) :: Vass, IPIV
  real(dp), dimension (1:3,1:3,1:n1_c,1:n2_c) :: DJI_matrix, DJI_matrix1
  real(dp), dimension (1:3*n1_c*n2_c) :: LHS_off,LHS,residual,a,a1
  real(dp), dimension (1:3) :: IPIV_DJ,Vass1
  real(dp) :: tol = 1e-7
  integer INFO, i0, j0
  real(dp), dimension (-2:n1_f+3,-2:n2_f+3) :: Mass_f1 = 0.d0

  ! traction B.C. on the top
  real(dp), dimension (1:n1_f,1:n2_f,1:dim):: traction_data
  real(dp), dimension (1:dim):: traction_rhs = 0.d0

  ! time
  real(dp):: tv
  integer :: cr,c1,c2,nt
  real(dp) :: rate
  real(dp):: cfl, dt0, dt

  !!!!! test
  real(dp),dimension(1:n1_c,1:n2_c,1:dim) :: u_test
  real(dp) :: err_test
  real(dp),dimension(1:n1_c,1:n2_c) :: err_array
  !!!!!!

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

  ! Difference stencils and interpolations
  call interpolation_restriction(P,Rop,RPop)
  call central_difference(ux_cof,uxx_cof)
  call varcoeff_NoGhost(acof_no_gp, ghcof_no_gp, sbop_no_gp )
  call VARCOEFFS4( acof, ghcof, Sb )
  call dx_46(bof)

  ! metric derivatives and coefficients of the equation
  call  equation_cof(Xgrid_c,Xgrid_f,X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c, &
       XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c,Jacobian_c_3,rho_c,rho_f, &
       X1R1_f,X1R2_f,X1R3_f,X2R1_f,X2R2_f,X2R3_f,X3R1_f,X3R2_f,X3R3_f, &
       XI11_f,XI21_f,XI31_f,XI12_f,XI22_f,XI32_f,XI13_f,XI23_f,XI33_f,Jacobian_f,Jacobian_f_3,N11_c,N12_c, &
       N13_c,N21_c,N22_c,N23_c,N31_c,N32_c,N33_c,N11_f,N12_f,N13_f,N21_f,N22_f,N23_f,N31_f,N32_f,N33_f)

  ! estimate time step
  call kappa(tv)
  cfl=sqrt(3.d0/1.40/tv)
  dt0 = min(h1_f,h2_f,h3_f)*cfl*0.9
  nt = ceiling(tn/dt0)+1
  dt = tn/nt

  print *, nt
  ! exact solution at t = -dt, 0
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
  do k=1-nrg,n3_f+nrg
     do i=1-nrg,n2_f+nrg
        do j=1-nrg,n1_f+nrg
           call initial_solution(l1,l2,l3,Xgrid_f(j,i,k,1),Xgrid_f(j,i,k,2),Xgrid_f(j,i,k,3), &
                u_f(j,i,k,1,1),u_f(j,i,k,2,1),u_f(j,i,k,3,1))
           call initial_solution(l1,l2,l3,Xgrid_f(j,i,k,1),Xgrid_f(j,i,k,2),Xgrid_f(j,i,k,3), &
                u_f(j,i,k,1,2),u_f(j,i,k,2,2),u_f(j,i,k,3,2))
        end do
     end do
  end do

  call Update_gp(1)
  call Update_gp(2)
  call Update_Dirichlet_BC(1)
  call Update_Dirichlet_BC(2)

  ! Construct the system matrix for computing ghost points values on the interface
  ! We have 3*n1_c*n2_c equations in 3D
  ! There are three sets of equations, one for the first component of u, one for the
  ! second component and another for the third component of u
  ! Each set consists of n1_c*n2_c equations

  call Interface_system(Mass)

  ! lu factorization
  call dgetrf(3*n1_c*n2_c,3*n1_c*n2_c,Mass,3*n1_c*n2_c,IPIV,INFO)
  if (INFO .ne. 0) then
     write(*,"(A20,I5)") 'LU fails, INFO=', INFO
  end if

  ! Before the time loop, we make the initial conditions compatible with interface conditions
  call Injection(1)

  call Interface_RHS(1)
  ! Solve system of linear equations to get values on ghost points
  call dgetrs('N',3*n1_c*n2_c,1,Mass,3*n1_c*n2_c,IPIV,Vass,3*n1_c*n2_c,INFO)
  if (INFO .ne. 0) then
     write(*,"(A20)") 'Solving fails'
  end if

  ! Update ghost points values for the interface
  do j = 1,n2_c
     do i = 1,n1_c
        u_c(i,j,n3_c+1,1,1) = Vass((j-1)*3*n1_c+3*(i-1)+1)
        u_c(i,j,n3_c+1,2,1) = Vass((j-1)*3*n1_c+3*(i-1)+2)
        u_c(i,j,n3_c+1,3,1) = Vass((j-1)*3*n1_c+3*(i-1)+3)
     end do
  end do

  call Injection(2)

  call Interface_RHS(2)
  ! Solve system of linear equations to get values on ghost points
  call dgetrs('N',3*n1_c*n2_c,1,Mass,3*n1_c*n2_c,IPIV,Vass,3*n1_c*n2_c,INFO)
  if (INFO .ne. 0) then
     write(*,"(A20)") 'Solving fails'
  end if

  ! Update ghost point values for the interface
  do j = 1,n2_c
     do i = 1,n1_c
        u_c(i,j,n3_c+1,1,2) = Vass((j-1)*3*n1_c+3*(i-1)+1)
        u_c(i,j,n3_c+1,2,2) = Vass((j-1)*3*n1_c+3*(i-1)+2)
        u_c(i,j,n3_c+1,3,2) = Vass((j-1)*3*n1_c+3*(i-1)+3)
     end do
  end do

  tv = 0.d0

  ! time stepping
  do time_index = 1, nt

    ! Evaluate the difference operators
    call Update_interior(u_c(:,:,:,:,2),u_f(:,:,:,:,2))
    call compute_eta(u_c(:,:,:,:,2),u_f(:,:,:,:,2))

     ! Update the solution  (it seems that we don't add the stability term eta)
     u_f(1:n1_f,1:n2_f,2:n3_f,:,4) = 2.d0*u_f(1:n1_f,1:n2_f,2:n3_f,:,2) - u_f(1:n1_f,1:n2_f,2:n3_f,:,1) &
          + dt**2*(lh_f(1:n1_f,1:n2_f,2:n3_f,:))/rho_f(1:n1_f,1:n2_f,2:n3_f,:)
     !
     u_f(1:n1_f,1:n2_f,1,:,4) = 2.d0*u_f(1:n1_f,1:n2_f,1,:,2) - u_f(1:n1_f,1:n2_f,1,:,1) &
          + dt**2*(lh_f(1:n1_f,1:n2_f,1,:)+eta_f(1:n1_f,1:n2_f,1,:))/rho_f(1:n1_f,1:n2_f,1,:)

     u_c(1:n1_c,1:n2_c,1:n3_c,:,4) = 2.d0*u_c(1:n1_c,1:n2_c,1:n3_c,:,2) - u_c(1:n1_c,1:n2_c,1:n3_c,:,1) &
          + dt**2*(lh_c(1:n1_c,1:n2_c,1:n3_c,:))/rho_c(1:n1_c,1:n2_c,1:n3_c,:)

     ! Update time
     tv = tv + dt
     Times(time_index) = tv

     ! Update ghost points outside the left and right boundaries
     ! The argument '3' means time level star
     ! '1', '2' and '4' mean time level n-1, n and n+1, respectively.
     call Update_gp(4)

     ! Update ghost point values for the traction boundary
     !call  Update_traction(3)

     ! Update Dirichlet boundary condition
     call Update_Dirichlet_BC(4)

     ! Injection at the interface
     call Injection(4)

     ! Build the right hand side vector for the system of linear equations for the interface
     call Interface_RHS(4)

     ! Solve system of linear equations to get values on ghost points
     call dgetrs('N',3*n1_c*n2_c,1,Mass,3*n1_c*n2_c,IPIV,Vass,3*n1_c*n2_c,INFO)
     if (INFO .ne. 0) then
        write(*,"(A20)") 'Solving fails'
     end if

     ! Update ghost point values for the interface
     do j = 1,n2_c
        do i = 1,n1_c
           u_c(i,j,n3_c+1,1,4) = Vass((j-1)*3*n1_c+3*(i-1)+1)
           u_c(i,j,n3_c+1,2,4) = Vass((j-1)*3*n1_c+3*(i-1)+2)
           u_c(i,j,n3_c+1,3,4) = Vass((j-1)*3*n1_c+3*(i-1)+3)
        end do
     end do

     ! compute energy
     call discrete_energy(u_f(:,:,:,:,1),u_f(:,:,:,:,2),u_f(:,:,:,:,3),u_f(:,:,:,:,4),&
        u_c(:,:,:,:,1),u_c(:,:,:,:,2),u_c(:,:,:,:,3),u_c(:,:,:,:,4),dt,&
        energy_discrete(time_index),energy_bdry(time_index))

     ! Update solutions
     u_f(:,:,:,:,1) = u_f(:,:,:,:,2)
     u_f(:,:,:,:,2) = u_f(:,:,:,:,4)
     u_c(:,:,:,:,1) = u_c(:,:,:,:,2)
     u_c(:,:,:,:,2) = u_c(:,:,:,:,4)

  end do ! end of time loop

  call print_array_to_file(1,1,nt,Times,'times.txt')
  call print_array_to_file(1,1,nt,energy_discrete,'energy.txt')
  call print_array_to_file(1,1,nt,energy_bdry,'energy_bdry.txt')

  N6(1) = n1_c
  N6(2) = n2_c
  N6(3) = n3_c
  N6(4) = n1_f
  N6(5) = n2_f
  N6(6) = n3_f

  call print_array_to_file(1,1,6,N6,'N6.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,1),'X1c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,2),'X2c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,Xgrid_c(1:n1_c,1:n2_c,1:n3_c,3),'X3c.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,1,4),'uc1.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,2,4),'uc2.txt')
  call print_array_to_file(n1_c,n2_c,n3_c,u_c(1:n1_c,1:n2_c,1:n3_c,3,4),'uc3.txt')

  call print_array_to_file(n1_f,n2_f,n3_f,Xgrid_f(1:n1_f,1:n2_f,1:n3_f,1),'X1f.txt')
  call print_array_to_file(n1_f,n2_f,n3_f,Xgrid_f(1:n1_f,1:n2_f,1:n3_f,2),'X2f.txt')
  call print_array_to_file(n1_f,n2_f,n3_f,Xgrid_f(1:n1_f,1:n2_f,1:n3_f,3),'X3f.txt')
  call print_array_to_file(n1_f,n2_f,n3_f,u_f(1:n1_f,1:n2_f,1:n3_f,1,4),'uf1.txt')
  call print_array_to_file(n1_f,n2_f,n3_f,u_f(1:n1_f,1:n2_f,1:n3_f,2,4),'uf2.txt')
  call print_array_to_file(n1_f,n2_f,n3_f,u_f(1:n1_f,1:n2_f,1:n3_f,3,4),'uf3.txt')

contains

  subroutine equation_cof(Xgrid_c,Xgrid_f,X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c,X3R1_c,X3R2_c,X3R3_c, &
       XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c,Jacobian_c_3,rho_c,rho_f, &
       X1R1_f,X1R2_f,X1R3_f,X2R1_f,X2R2_f,X2R3_f,X3R1_f,X3R2_f,X3R3_f, &
       XI11_f,XI21_f,XI31_f,XI12_f,XI22_f,XI32_f,XI13_f,XI23_f,XI33_f,Jacobian_f,Jacobian_f_3,N11_c,N12_c, &
       N13_c,N21_c,N22_c,N23_c,N31_c,N32_c,N33_c,N11_f,N12_f,N13_f,N21_f,N22_f,N23_f,N31_f,N32_f,N33_f)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N11_c,N12_c,N13_c,N21_c
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:dim) :: N22_c,N23_c,N31_c,N32_c,N33_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:dim) :: N11_f,N12_f,N13_f,N21_f
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:dim) :: N22_f,N23_f,N31_f,N32_f,N33_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: mu_c,lambda_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: mu_f,lambda_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,dim) :: Xgrid_c,Rgrid_c,Jacobian_c_3,rho_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,dim) :: Xgrid_f,Rgrid_f,Jacobian_f_3,rho_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg)::X1R1_c,X1R2_c,X1R3_c,X2R1_c,X2R2_c,X2R3_c, &
         X3R1_c,X3R2_c,X3R3_c,XI11_c,XI21_c,XI31_c,XI12_c,XI22_c,XI32_c,XI13_c,XI23_c,XI33_c,Jacobian_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg)::X1R1_f,X1R2_f,X1R3_f,X2R1_f,X2R2_f,X2R3_f, &
         X3R1_f,X3R2_f,X3R3_f,XI11_f,XI21_f,XI31_f,XI12_f,XI22_f,XI32_f,XI13_f,XI23_f,XI33_f,Jacobian_f


    call generate_grid(Rgrid_c(:,:,:,1),Rgrid_c(:,:,:,2),Rgrid_c(:,:,:,3), &
         Xgrid_c(:,:,:,1),Xgrid_c(:,:,:,2),Xgrid_c(:,:,:,3), &
         Rgrid_f(:,:,:,1),Rgrid_f(:,:,:,2),Rgrid_f(:,:,:,3), &
         Xgrid_f(:,:,:,1),Xgrid_f(:,:,:,2),Xgrid_f(:,:,:,3))

    ! compute metric derivatives, coarse
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             call  metric_derivative(Rgrid_c(k,j,i,1),Rgrid_c(k,j,i,2),Rgrid_c(k,j,i,3),X1R1_c(k,j,i),X1R2_c(k,j,i),X1R3_c(k,j,i), &
                 X2R1_c(k,j,i),X2R2_c(k,j,i),X2R3_c(k,j,i),X3R1_c(k,j,i),X3R2_c(k,j,i),X3R3_c(k,j,i), &
                 XI11_c(k,j,i),XI21_c(k,j,i),XI31_c(k,j,i),XI12_c(k,j,i),XI22_c(k,j,i),XI32_c(k,j,i), &
                 XI13_c(k,j,i),XI23_c(k,j,i),XI33_c(k,j,i),Jacobian_c(k,j,i),0)
          end do
       end do
    end do
    ! compute metric derivatives, fine
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
             call  metric_derivative(Rgrid_f(k,j,i,1),Rgrid_f(k,j,i,2),Rgrid_f(k,j,i,3),X1R1_f(k,j,i),X1R2_f(k,j,i),X1R3_f(k,j,i), &
                 X2R1_f(k,j,i),X2R2_f(k,j,i),X2R3_f(k,j,i),X3R1_f(k,j,i),X3R2_f(k,j,i),X3R3_f(k,j,i), &
                 XI11_f(k,j,i),XI21_f(k,j,i),XI31_f(k,j,i),XI12_f(k,j,i),XI22_f(k,j,i),XI32_f(k,j,i), &
                 XI13_f(k,j,i),XI23_f(k,j,i),XI33_f(k,j,i),Jacobian_f(k,j,i),1)
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
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
             mu_f(k,j,i) = 3.d0 + sin(3.d0*Xgrid_f(k,j,i,1)+0.1d0)*sin(3.d0*Xgrid_f(k,j,i,2)+0.1d0)*sin(Xgrid_f(k,j,i,3))
             lambda_f(k,j,i) = 21.d0+ cos(Xgrid_f(k,j,i,1)+0.1d0)*cos(Xgrid_f(k,j,i,2)+0.1d0)*sin(3.d0*Xgrid_f(k,j,i,3))**2
             rho_f(k,j,i,:) = 2.d0 + sin(Xgrid_f(k,j,i,1)+0.3d0)*sin(Xgrid_f(k,j,i,2)+0.3d0)*sin(Xgrid_f(k,j,i,3)-0.2d0)
          end do
       end do
    end do

    Jacobian_c_3(:,:,:,1) = Jacobian_c
    Jacobian_c_3(:,:,:,2) = Jacobian_c
    Jacobian_c_3(:,:,:,3) = Jacobian_c
    Jacobian_f_3(:,:,:,1) = Jacobian_f
    Jacobian_f_3(:,:,:,2) = Jacobian_f
    Jacobian_f_3(:,:,:,3) = Jacobian_f
    rho_c = rho_c*Jacobian_c_3
    rho_f = rho_f*Jacobian_f_3


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

    ! fine
    N11_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI11_f*XI11_f+mu_f*(XI21_f*XI21_f+XI31_f*XI31_f))
    N11_f(:,:,:,1,2) = Jacobian_f*XI11_f*XI21_f*(mu_f+lambda_f)
    N11_f(:,:,:,1,3) = Jacobian_f*XI11_f*XI31_f*(mu_f+lambda_f)
    N11_f(:,:,:,2,1) = N11_f(:,:,:,1,2)
    N11_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI21_f*XI21_f+mu_f*(XI11_f*XI11_f+XI31_f*XI31_f))
    N11_f(:,:,:,2,3) = Jacobian_f*XI21_f*XI31_f*(mu_f+lambda_f)
    N11_f(:,:,:,3,1) = N11_f(:,:,:,1,3)
    N11_f(:,:,:,3,2) = N11_f(:,:,:,2,3)
    N11_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI31_f*XI31_f+mu_f*(XI11_f*XI11_f+XI21_f*XI21_f))

    N22_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI12_f*XI12_f+mu_f*XI22_f*XI22_f+mu_f*XI32_f*XI32_f)
    N22_f(:,:,:,1,2) = Jacobian_f*XI12_f*XI22_f*(mu_f+lambda_f)
    N22_f(:,:,:,1,3) = Jacobian_f*XI12_f*XI32_f*(mu_f+lambda_f)
    N22_f(:,:,:,2,1) = N22_f(:,:,:,1,2)
    N22_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI22_f*XI22_f+mu_f*XI12_f*XI12_f+mu_f*XI32_f*XI32_f)
    N22_f(:,:,:,2,3) = Jacobian_f*XI22_f*XI32_f*(mu_f+lambda_f)
    N22_f(:,:,:,3,1) = N22_f(:,:,:,1,3)
    N22_f(:,:,:,3,2) = N22_f(:,:,:,2,3)
    N22_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI32_f*XI32_f+mu_f*XI12_f*XI12_f+mu_f*XI22_f*XI22_f)

    N33_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI13_f*XI13_f+mu_f*(XI23_f*XI23_f+XI33_f*XI33_f))
    N33_f(:,:,:,1,2) = Jacobian_f*XI13_f*XI23_f*(mu_f+lambda_f)
    N33_f(:,:,:,1,3) = Jacobian_f*XI13_f*XI33_f*(mu_f+lambda_f)
    N33_f(:,:,:,2,1) = N33_f(:,:,:,1,2)
    N33_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI23_f*XI23_f+XI13_f*XI13_f*mu_f+XI33_f*XI33_f*mu_f)
    N33_f(:,:,:,2,3) = Jacobian_f*XI23_f*XI33_f*(mu_f+lambda_f)
    N33_f(:,:,:,3,1) = N33_f(:,:,:,1,3)
    N33_f(:,:,:,3,2) = N33_f(:,:,:,2,3)
    N33_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI33_f*XI33_f+mu_f*(XI13_f*XI13_f+XI23_f*XI23_f))

    N12_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI11_f*XI12_f+mu_f*XI21_f*XI22_f+mu_f*XI31_f*XI32_f)
    N12_f(:,:,:,1,2) = Jacobian_f*(XI11_f*XI22_f*lambda_f+XI21_f*XI12_f*mu_f)
    N12_f(:,:,:,1,3) = Jacobian_f*(XI11_f*XI32_f*lambda_f+XI31_f*XI12_f*mu_f)
    N12_f(:,:,:,2,1) = Jacobian_f*(XI11_f*XI22_f*mu_f+XI21_f*XI12_f*lambda_f)
    N12_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI21_f*XI22_f+mu_f*XI11_f*XI12_f+mu_f*XI31_f*XI32_f)
    N12_f(:,:,:,2,3) = Jacobian_f*(XI21_f*XI32_f*lambda_f+XI31_f*XI22_f*mu_f)
    N12_f(:,:,:,3,1) = Jacobian_f*(XI11_f*XI32_f*mu_f+XI31_f*XI12_f*lambda_f)
    N12_f(:,:,:,3,2) = Jacobian_f*(XI21_f*XI32_f*mu_f+XI31_f*XI22_f*lambda_f)
    N12_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI31_f*XI32_f+mu_f*XI11_f*XI12_f+mu_f*XI21_f*XI22_f)

    N21_f(:,:,:,1,1) = N12_f(:,:,:,1,1)
    N21_f(:,:,:,1,2) = N12_f(:,:,:,2,1)
    N21_f(:,:,:,1,3) = N12_f(:,:,:,3,1)
    N21_f(:,:,:,2,1) = N12_f(:,:,:,1,2)
    N21_f(:,:,:,2,2) = N12_f(:,:,:,2,2)
    N21_f(:,:,:,2,3) = N12_f(:,:,:,3,2)
    N21_f(:,:,:,3,1) = N12_f(:,:,:,1,3)
    N21_f(:,:,:,3,2) = N12_f(:,:,:,2,3)
    N21_f(:,:,:,3,3) = N12_f(:,:,:,3,3)

    N13_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI11_f*XI13_f+XI21_f*XI23_f*mu_f+XI31_f*XI33_f*mu_f)
    N13_f(:,:,:,1,2) = Jacobian_f*(XI11_f*XI23_f*lambda_f+XI21_f*XI13_f*mu_f)
    N13_f(:,:,:,1,3) = Jacobian_f*(XI11_f*XI33_f*lambda_f+XI31_f*XI13_f*mu_f)
    N13_f(:,:,:,2,1) = Jacobian_f*(XI11_f*XI23_f*mu_f+XI21_f*XI13_f*lambda_f)
    N13_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI21_f*XI23_f+mu_f*XI11_f*XI13_f+mu_f*XI31_f*XI33_f)
    N13_f(:,:,:,2,3) = Jacobian_f*(XI21_f*XI33_f*lambda_f+XI31_f*XI23_f*mu_f)
    N13_f(:,:,:,3,1) = Jacobian_f*(XI11_f*XI33_f*mu_f+XI31_f*XI13_f*lambda_f)
    N13_f(:,:,:,3,2) = Jacobian_f*(XI21_f*XI33_f*mu_f+XI31_f*XI23_f*lambda_f)
    N13_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI31_f*XI33_f+mu_f*XI11_f*XI13_f+mu_f*XI21_f*XI23_f)

    N31_f(:,:,:,1,1) = N13_f(:,:,:,1,1)
    N31_f(:,:,:,1,2) = N13_f(:,:,:,2,1)
    N31_f(:,:,:,1,3) = N13_f(:,:,:,3,1)
    N31_f(:,:,:,2,1) = N13_f(:,:,:,1,2)
    N31_f(:,:,:,2,2) = N13_f(:,:,:,2,2)
    N31_f(:,:,:,2,3) = N13_f(:,:,:,3,2)
    N31_f(:,:,:,3,1) = N13_f(:,:,:,1,3)
    N31_f(:,:,:,3,2) = N13_f(:,:,:,2,3)
    N31_f(:,:,:,3,3) = N13_f(:,:,:,3,3)

    N23_f(:,:,:,1,1) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI12_f*XI13_f+mu_f*XI22_f*XI23_f+mu_f*XI32_f*XI33_f)
    N23_f(:,:,:,1,2) = Jacobian_f*(XI12_f*XI23_f*lambda_f+XI22_f*XI13_f*mu_f)
    N23_f(:,:,:,1,3) = Jacobian_f*(XI12_f*XI33_f*lambda_f+XI32_f*XI13_f*mu_f)
    N23_f(:,:,:,2,1) = Jacobian_f*(XI12_f*XI23_f*mu_f+XI22_f*XI13_f*lambda_f)
    N23_f(:,:,:,2,2) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI22_f*XI23_f+mu_f*XI12_f*XI13_f+mu_f*XI32_f*XI33_f)
    N23_f(:,:,:,2,3) = Jacobian_f*(XI22_f*XI33_f*lambda_f+XI32_f*XI23_f*mu_f)
    N23_f(:,:,:,3,1) = Jacobian_f*(XI12_f*XI33_f*mu_f+XI32_f*XI13_f*lambda_f)
    N23_f(:,:,:,3,2) = Jacobian_f*(XI22_f*XI33_f*mu_f+XI32_f*XI23_f*lambda_f)
    N23_f(:,:,:,3,3) = Jacobian_f*((2.d0*mu_f+lambda_f)*XI32_f*XI33_f+mu_f*XI12_f*XI13_f+mu_f*XI22_f*XI23_f)

    N32_f(:,:,:,1,1) = N23_f(:,:,:,1,1)
    N32_f(:,:,:,1,2) = N23_f(:,:,:,2,1)
    N32_f(:,:,:,1,3) = N23_f(:,:,:,3,1)
    N32_f(:,:,:,2,1) = N23_f(:,:,:,1,2)
    N32_f(:,:,:,2,2) = N23_f(:,:,:,2,2)
    N32_f(:,:,:,2,3) = N23_f(:,:,:,3,2)
    N32_f(:,:,:,3,1) = N23_f(:,:,:,1,3)
    N32_f(:,:,:,3,2) = N23_f(:,:,:,2,3)
    N32_f(:,:,:,3,3) = N23_f(:,:,:,3,3)
  end subroutine  equation_cof

  subroutine kappa(kappa2)
    real(dp), dimension (3) :: kappa1
    real(dp), dimension (3,3) :: mat_t
    real(dp), dimension (8) :: work
    integer :: info
    real(dp) :: kappa2
    kappa2 = 0.d0
    do k = 1, n3_f
       do j = 1, n2_f
          do i = 1,n1_f
             mat_t(1,1) = (N11_f(i,j,k,1,1)+N11_f(i,j,k,2,2)+N11_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(1,2) = (N12_f(i,j,k,1,1)+N12_f(i,j,k,2,2)+N12_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(1,3) = (N13_f(i,j,k,1,1)+N13_f(i,j,k,2,2)+N13_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(2,1) = (N21_f(i,j,k,1,1)+N21_f(i,j,k,2,2)+N21_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(2,2) = (N22_f(i,j,k,1,1)+N22_f(i,j,k,2,2)+N22_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(2,3) = (N23_f(i,j,k,1,1)+N23_f(i,j,k,2,2)+N23_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(3,1) = (N31_f(i,j,k,1,1)+N31_f(i,j,k,2,2)+N31_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(3,2) = (N32_f(i,j,k,1,1)+N32_f(i,j,k,2,2)+N32_f(i,j,k,3,3))/rho_f(i,j,k,1)
             mat_t(3,3) = (N33_f(i,j,k,1,1)+N33_f(i,j,k,2,2)+N33_f(i,j,k,3,3))/rho_f(i,j,k,1)
             call dsyev('N','U',3,mat_t,3,kappa1,work,8,info)
             if (kappa1(3) > kappa2) then
                kappa2 = kappa1(3)
             end if
          end do
       end do
    end do
  end subroutine kappa

  subroutine Interface_system(Mass)

    ! System of linear equations
    real(dp), dimension (1:n1_c*n2_c*3,1:n1_c*n2_c*3) :: Mass
    real(dp), dimension (-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3) :: Mass_p = 0.d0
    real(dp), dimension (1:n1_c,1:n2_c,1:n1_c,1:n2_c) :: Mass_r = 0.d0
    real(dp) :: int_cof
    real(dp), dimension (1:n1_c,1:n2_c) :: int_cof_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg) :: int_cof_f

    Mass =0.d0

    int_cof = 17.d0/48.d0*h3_f*ghcof(1)/h3_c**2

    do j = 1,n2_c
       do i = 1,n1_c
          int_cof_c(i,j) = Jacobian_c(i,j,n3_c)*sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)
       end do
    end do

    do j = 1-nrg,n2_f+nrg
       do i = 1-nrg,n1_f+nrg
          int_cof_f(i,j) = Jacobian_f(i,j,1)*sqrt(XI13_f(i,j,1)**2+XI23_f(i,j,1)**2+XI33_f(i,j,1)**2)
       end do
    end do

    do iset = 1,3
       do jj = 1,3
         !
          do k = -1,n2_c+1
             do i = -1,n1_c+1
               if (i .gt. -1 .and. k .gt. -1) then
                   Mass_p(2*i-1,2*k-1,i,k) = N33_c(i,k,n3_c,iset,jj)/rho_c(i,k,n3_c,iset)
               end if
               if (i .gt. -1 .and. k .ge. -1) then
                  do j = -1,2
                     Mass_p(2*i-1,2*k,i,k+j) = Mass_p(2*i-1,2*k,i,k+j) &
                        + P(j)*N33_c(i,k+j,n3_c,iset,jj)/rho_c(i,k+j,n3_c,iset)
                  end do
               end if
               if (i .ge. -1 .and. k .gt. -1) then
                  do j = -1,2
                     Mass_p(2*i,2*k-1,i+j,k) = Mass_p(2*i,2*k-1,i+j,k) &
                        + P(j)*N33_c(i+j,k,n3_c,iset,jj)/rho_c(i+j,k,n3_c,iset)
                  end do
               end if
                do l = -1,2
                   do j = -1,2
                      Mass_p(2*i,2*k,i+j,k+l) = Mass_p(2*i,2*k,i+j,k+l) &
                         +P(l)*(P(j)*N33_c(i+j,k+l,n3_c,iset,jj)/rho_c(i+j,k+l,n3_c,iset))
                  end do
              end do
           end do
        end do

        Mass_p(:,:,-2:0,:) = 0.d0
        Mass_p(:,:,n1_c+1:n1_c+3,:) = 0.d0
        Mass_p(:,:,:,-2:0) = 0.d0
        Mass_p(:,:,:,n2_c+1:n2_c+3) = 0.d0

        do k = 1,n2_c
           do i = 1,n1_c
              do l = -4,2
                 do j = -4,2
                      Mass_r(i,k,1:n1_c,1:n2_c) = Mass_r(i,k,1:n1_c,1:n2_c) &
                      + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1,iset)*Mass_p(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                      /int_cof_f(2*i+j,2*k+l))
                 end do
              end do
           end do
        end do
        do l = 1,n2_c
           do k = 1,n1_c
              do j = 1,n2_c
                 do i = 1,n1_c
                    Mass((l-1)*3*n1_c+(k-1)*3+iset,(j-1)*3*n1_c+(i-1)*3+jj) = Mass_r(k,l,i,j)*int_cof
                 end do
              end do
           end do
        end do
        do j = 1, n2_c
           do i = 1, n1_c
              Mass((j-1)*3*n1_c+(i-1)*3+iset,(j-1)*3*n1_c+(i-1)*3+jj) = &
                  Mass((j-1)*3*n1_c+(i-1)*3+iset,(j-1)*3*n1_c+(i-1)*3+jj) &
                  -Sb(0)*N33_c(i,j,n3_c,iset,jj)/h3_c/int_cof_c(i,j)
           end do
        end do
        Mass_p = 0.d0
        Mass_r = 0.d0
      end do
      Mass_p = 0.d0
      Mass_r = 0.d0
    end do

  end subroutine Interface_system

  subroutine Interface_RHS(index)
    integer index
    real(dp), dimension (1:n1_c,1:n2_c) :: int_cof_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg) :: int_cof_f

    do j = 1,n2_c
       do i = 1,n1_c
          int_cof_c(i,j) = Jacobian_c(i,j,n3_c)*sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)
       end do
    end do

    do j = 1-nrg,n2_f+nrg
       do i = 1-nrg,n1_f+nrg
          int_cof_f(i,j) = Jacobian_f(i,j,1)*sqrt(XI13_f(i,j,1)**2+XI23_f(i,j,1)**2+XI33_f(i,j,1)**2)
       end do
    end do

    Vass = 0.d0
    !
    do iset = 1,3
       ! term 1
       do k = 1,n2_c
          do i = 1,n1_c
             do j = 1,4
                do jj = 1,3
                   Vass((k-1)*3*n1_c+(i-1)*3+iset) = Vass((k-1)*3*n1_c+(i-1)*3+iset) &
                     + N33_c(i,k,n3_c,iset,jj)*Sb(j)*u_c(i,k,n3_c+1-j,jj,index)/h3_c/int_cof_c(i,k)
                end do
             end do
             do j = -2,2
                do jj = 1,3
                   Vass((k-1)*3*n1_c+(i-1)*3+iset) = Vass((k-1)*3*n1_c+(i-1)*3+iset) &
                     - N31_c(i,k,n3_c,iset,jj)*ux_cof(j)*u_c(i+j,k,n3_c,jj,index)/h1_c/int_cof_c(i,k) &
                     - N32_c(i,k,n3_c,iset,jj)*ux_cof(j)*u_c(i,k+j,n3_c,jj,index)/h2_c/int_cof_c(i,k)
                end do
             end do
          end do
       end do

       ! term 2
       ! interior
       ! second derivatives 11
       lh_c = 0.d0
       int_temp_c = 0.d0
       do jj = 1,3
          call FD_r1r1(h1_c,n1_c,n2_c,n3_c,N11_c(:,:,:,iset,jj),u_c(:,:,:,jj,index), &
                -2,n1_c+3,-2,n2_c+3,n3_c,n3_c,int_temp_c(-2:n1_c+3,-2:n2_c+3,1))
          lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)+int_temp_c(-2:n1_c+3,-2:n2_c+3,1)
       end do

       ! second derivatives 22
       int_temp_c = 0.d0
       do jj = 1,3
          call FD_r2r2(h2_c,n1_c,n2_c,n3_c,N22_c(:,:,:,iset,jj),u_c(:,:,:,jj,index), &
                -2,n1_c+3,-2,n2_c+3,n3_c,n3_c,int_temp_c(-2:n1_c+3,-2:n2_c+3,1))
          lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)+int_temp_c(-2:n1_c+3,-2:n2_c+3,1)
       end do

       ! second derivatives 33
       do k = 1,8
          do m = 1,8
             do jj = 1,3
                lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) =lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                 +acof(1,k,m)*N33_c(-2:n1_c+3,-2:n2_c+3,n3_c+1-m,iset,jj)*u_c(-2:n1_c+3,-2:n2_c+3,n3_c+1-k,jj,index)/h3_c**2
             end do
          end do
       end do
       ! ghost points
       do jj = 1,3
          lh_c(-2:0,-2:n2_c+3,n3_c,iset) =lh_c(-2:0,-2:n2_c+3,n3_c,iset) &
            + N33_c(-2:0,-2:n2_c+3,n3_c,iset,jj)*ghcof(1)*u_c(-2:0,-2:n2_c+3,n3_c+1,jj,index)/h3_c**2
          lh_c(n1_c+1:n1_c+3,-2:n2_c+3,n3_c,iset) = lh_c(n1_c+1:n1_c+3,-2:n2_c+3,n3_c,iset) &
            + N33_c(n1_c+1:n1_c+3,-2:n2_c+3,n3_c,iset,jj)*ghcof(1)*u_c(n1_c+1:n1_c+3,-2:n2_c+3,n3_c+1,jj,index)/h3_c**2
          lh_c(1:n1_c,-2:0,n3_c,iset) =lh_c(1:n1_c,-2:0,n3_c,iset) &
            + N33_c(1:n1_c,-2:0,n3_c,iset,jj)*ghcof(1)*u_c(1:n1_c,-2:0,n3_c+1,jj,index)/h3_c**2
          lh_c(1:n1_c,n2_c+1:n2_c+3,n3_c,iset) = lh_c(1:n1_c,n2_c+1:n2_c+3,n3_c,iset) &
            + N33_c(1:n1_c,n2_c+1:n2_c+3,n3_c,iset,jj)*ghcof(1)*u_c(1:n1_c,n2_c+1:n2_c+3,n3_c+1,jj,index)/h3_c**2
       end do

       ! mixed derivatives 12
       do jj = 1,3
          int_temp_c = 0.d0
          do i=-2,2
             int_temp_c(-4:n1_c+5,-2:n2_c+3,1)=int_temp_c(-4:n1_c+5,-2:n2_c+3,1) &
                       +ux_cof(i)*u_c(-4:n1_c+5,-2+i:n2_c+3+i,n3_c,jj,index)/h2_c
          end do
          do i=-2,2
             lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                +ux_cof(i)*N12_c(-2+i:n1_c+3+i,-2:n2_c+3,n3_c,iset,jj)*int_temp_c(-2+i:n1_c+3+i,-2:n2_c+3,1)/h1_c
          end do
       end do

       ! mixed derivatives 13
       do jj = 1,3
          int_temp_c = 0.d0
          do i =1,4
             int_temp_c(-4:n1_c+5,-2:n2_c+3,1)=int_temp_c(-4:n1_c+5,-2:n2_c+3,1) &
                        -u_c(-4:n1_c+5,-2:n2_c+3,n3_c+1-i,jj,index)*bof(1,i)/h3_c
          end do
          do i =-2,2
             lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                +ux_cof(i)*N13_c(-2+i:n1_c+3+i,-2:n2_c+3,n3_c,iset,jj)*int_temp_c(-2+i:n1_c+3+i,-2:n2_c+3,1)/h1_c
          end do
       end do

       ! mixed derivatives 21
       do jj = 1,3
          int_temp_c = 0.d0
          do i=-2,2
             int_temp_c(-2:n1_c+3,-4:n2_c+5,1)=int_temp_c(-2:n1_c+3,-4:n2_c+5,1) &
                       +ux_cof(i)*u_c(-2+i:n1_c+3+i,-4:n2_c+5,n3_c,jj,index)/h1_c
          end do
          do i=-2,2
             lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                +ux_cof(i)*N21_c(-2:n1_c+3,-2+i:n2_c+3+i,n3_c,iset,jj)*int_temp_c(-2:n1_c+3,-2+i:n2_c+3+i,1)/h2_c
          end do
      end do

      ! mixed derivatives 23
      do jj = 1,3
         int_temp_c = 0.d0
         do i =1,4
            int_temp_c(-2:n1_c+3,-4:n2_c+5,1)=int_temp_c(-2:n1_c+3,-4:n2_c+5,1) &
                        -u_c(-2:n1_c+3,-4:n2_c+5,n3_c+1-i,jj,index)*bof(1,i)/h3_c
         end do
         do i =-2,2
            lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                +ux_cof(i)*N23_c(-2:n1_c+3,-2+i:n2_c+3+i,n3_c,iset,jj)*int_temp_c(-2:n1_c+3,-2+i:n2_c+3+i,1)/h2_c
         end do
      end do

      ! mixed derivatives 31
      do jj = 1,3
         int_temp_c = 0.d0
         do i = -2,2
            int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4)=int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4) &
               +ux_cof(i)*u_c(-2+i:n1_c+3+i,-2:n2_c+3,n3_c:n3_c-3:-1,jj,index)/h1_c
         end do
         do i = 1,4
            lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
               -N31_c(-2:n1_c+3,-2:n2_c+3,n3_c-i+1,iset,jj)*int_temp_c(-2:n1_c+3,-2:n2_c+3,i)*bof(1,i)/h3_c
         end do
      end do

      ! mixed derivatives 32
      do jj = 1,3
         int_temp_c = 0.d0
         do i = -2,2
            int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4)=int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4) &
               +ux_cof(i)*u_c(-2:n1_c+3,-2+i:n2_c+3+i,n3_c:n3_c-3:-1,jj,index)/h2_c
         end do
         do i = 1,4
            lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=lh_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
                -N32_c(-2:n1_c+3,-2:n2_c+3,n3_c-i+1,iset,jj)*int_temp_c(-2:n1_c+3,-2:n2_c+3,i)*bof(1,i)/h3_c
         end do
      end do

      ! project
      Mass_f1 = 0.d0
      do k = -1,n2_c+1
         do i = -1,n1_c+1
            if (k .gt. -1 .and. i .gt. -1) then
                Mass_f1(2*i-1,2*k-1)=lh_c(i,k,n3_c,iset)/rho_c(i,k,n3_c,1)
            end if
            if (k .ge. -1 .and. i .gt. -1) then
               do j = -1,2
                  Mass_f1(2*i-1,2*k)=Mass_f1(2*i-1,2*k)+P(j)*lh_c(i,k+j,n3_c,iset)/rho_c(i,j+k,n3_c,1)
               end do
            end if
            if (k .gt. -1 .and. i .ge. -1) then
               do j = -1,2
                  Mass_f1(2*i,2*k-1)=Mass_f1(2*i,2*k-1)+P(j)*lh_c(i+j,k,n3_c,iset)/rho_c(i+j,k,n3_c,1)
               end do
            end if
            do j = -1,2
               do l = -1,2
                   Mass_f1(2*i,2*k)=Mass_f1(2*i,2*k)+P(j)*(P(l)*lh_c(i+l,k+j,n3_c,iset)/rho_c(i+l,k+j,n3_c,1))
                end do
             end do
           end do
       end do
      ! restrict
      do k = 1,n2_c
         do i = 1,n1_c
            do j = -4,2
               do l = -4,2
                  Vass((k-1)*3*n1_c+(i-1)*3+iset)=Vass((k-1)*3*n1_c+(i-1)*3+iset) &
                  -17.d0/48.d0*h3_f*Rop(j)*(Rop(l)*rho_f(2*i+l,2*k+j,1,1)*Mass_f1(2*i+l,2*k+j)*1.d0/int_cof_f(2*i+l,2*k+j))
               end do
            end do
         end do
      end do

      ! term 3
      ! second derivatives 11
      lh_f = 0.d0
      int_temp_f = 0.d0
      do jj = 1,3
         call FD_r1r1(h1_f,n1_f,n2_f,n3_f,N11_f(:,:,:,iset,jj),u_f(:,:,:,jj,index), &
                      -2,n1_f+3,-2,n2_f+3,1,1,int_temp_f(-2:n1_f+3,-2:n2_f+3,1))
         lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) = lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)+int_temp_f(-2:n1_f+3,-2:n2_f+3,1)
      end do

      ! second derivatives 22
      int_temp_f = 0.d0
      do jj = 1,3
         call FD_r2r2(h2_f,n1_f,n2_f,n3_f,N22_f(:,:,:,iset,jj),u_f(:,:,:,jj,index), &
                      -2,n1_f+3,-2,n2_f+3,1,1,int_temp_f(-2:n1_f+3,-2:n2_f+3,1))
         lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) = lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)+int_temp_f(-2:n1_f+3,-2:n2_f+3,1)
      end do

      ! second derivatives 33
      do k = 1,8
         do m = 1,8
            do jj = 1,3
               lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) =lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
                 +acof_no_gp(1,k,m)*N33_f(-2:n1_f+3,-2:n2_f+3,m,iset,jj)*u_f(-2:n1_f+3,-2:n2_f+3,k,jj,index)/h3_f**2
            end do
         end do
      end do

      ! mixed derivatives 12
      do jj = 1,3
         int_temp_f = 0.d0
         do i=-2,2
            int_temp_f(-4:n1_f+5,-2:n2_f+3,1)=int_temp_f(-4:n1_f+5,-2:n2_f+3,1) &
                          +ux_cof(i)*u_f(-4:n1_f+5,-2+i:n2_f+3+i,1,jj,index)/h2_f
         end do
         do i=-2,2
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
                +ux_cof(i)*N12_f(-2+i:n1_f+3+i,-2:n2_f+3,1,iset,jj)*int_temp_f(-2+i:n1_f+3+i,-2:n2_f+3,1)/h1_f
         end do
      end do

      ! mixed derivatives 13
      do jj = 1,3
         int_temp_f = 0.d0
         do i =1,4
            int_temp_f(-4:n1_f+5,-2:n2_f+3,1)=int_temp_f(-4:n1_f+5,-2:n2_f+3,1) &
                               +u_f(-4:n1_f+5,-2:n2_f+3,i,jj,index)*bof(1,i)/h3_f
         end do
         do i =-2,2
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
               +ux_cof(i)*N13_f(-2+i:n1_f+3+i,-2:n2_f+3,1,iset,jj)*int_temp_f(-2+i:n1_f+3+i,-2:n2_f+3,1)/h1_f
         end do
      end do

      ! mixed derivatives 21
      do jj = 1,3
         int_temp_f = 0.d0
         do i=-2,2
            int_temp_f(-2:n1_f+3,-4:n2_f+5,1)=int_temp_f(-2:n1_f+3,-4:n2_f+5,1) &
                          +ux_cof(i)*u_f(-2+i:n1_f+3+i,-4:n2_f+5,1,jj,index)/h1_f
         end do
         do i=-2,2
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
               +ux_cof(i)*N21_f(-2:n1_f+3,-2+i:n2_f+3+i,1,iset,jj)*int_temp_f(-2:n1_f+3,-2+i:n2_f+3+i,1)/h2_f
        end do
      end do

      ! mixed derivatives 23
      do jj = 1,3
         int_temp_f = 0.d0
         do i =1,4
            int_temp_f(-2:n1_f+3,-4:n2_f+5,1)=int_temp_f(-2:n1_f+3,-4:n2_f+5,1) &
                               +u_f(-2:n1_f+3,-4:n2_f+5,i,jj,index)*bof(1,i)/h3_f
         end do
         do i =-2,2
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
                +ux_cof(i)*N23_f(-2:n1_f+3,-2+i:n2_f+3+i,1,iset,jj)*int_temp_f(-2:n1_f+3,-2+i:n2_f+3+i,1)/h2_f
         end do
      end do

      ! mixed derivatives 31
      do jj = 1,3
         int_temp_f = 0.d0
         do i = -2,2
            int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4)=int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4) &
                        +ux_cof(i)*u_f(-2+i:n1_f+3+i,-2:n2_f+3,1:4,jj,index)/h1_f
         end do
         do i = 1,4
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
               +N31_f(-2:n1_f+3,-2:n2_f+3,i,iset,jj)*int_temp_f(-2:n1_f+3,-2:n2_f+3,i)*bof(1,i)/h3_f
         end do
      end do

      ! mixed derivatives 32
      do jj = 1,3
         int_temp_f = 0.d0
         do i = -2,2
            int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4)=int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4) &
                        +ux_cof(i)*u_f(-2:n1_f+3,-2+i:n2_f+3+i,1:4,jj,index)/h2_f
         end do
         do i = 1,4
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
                +N32_f(-2:n1_f+3,-2:n2_f+3,i,iset,jj)*int_temp_f(-2:n1_f+3,-2:n2_f+3,i)*bof(1,i)/h3_f
         end do
      end do

      ! Scale L_f1 !!!!!!!
      lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)=lh_f(-2:n1_f+3,-2:n2_f+3,1,iset)*17.d0/48.d0*h3_f

      ! Add another three terms that need to be restricted
      do j = 1,5
         do jj = 1,3
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) = lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
              + N33_f(-2:n1_f+3,-2:n2_f+3,1,iset,jj)*sbop_no_gp(j)*u_f(-2:n1_f+3,-2:n2_f+3,j,jj,index)/h3_f
          end do
      end do
      do j = -2,2
         do jj = 1,3
            lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) = lh_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
              + (N31_f(-2:n1_f+3,-2:n2_f+3,1,iset,jj)*ux_cof(j)*u_f(-2+j:n1_f+3+j,-2:n2_f+3,1,jj,index)/h1_f &
              + N32_f(-2:n1_f+3,-2:n2_f+3,1,iset,jj)*ux_cof(j)*u_f(-2:n1_f+3,-2+j:n2_f+3+j,1,jj,index)/h2_f)
         end do
      end do

      ! now restrict it to the coarse grid
      do k = 1,n2_c
         do i = 1,n1_c
            do j = -4,2
               do l = -4,2
                  Vass((k-1)*3*n1_c+(i-1)*3+iset) = Vass((k-1)*3*n1_c+(i-1)*3+iset) &
                        + Rop(j)*(Rop(l)*lh_f(2*i+l,2*k+j,1,iset)/int_cof_f(2*i+l,2*k+j))
               end do
            end do
         end do
      end do
      !
   end do

  end subroutine Interface_RHS

  subroutine compute_eta(u_c_t,u_f_t)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: u_f_t
    !
    do iset = 1,3

    eta_c = 0.d0
    int_temp_c = 0.d0

    ! first derivative 11
    do jj = 1,3
       call FD_r1r1(h1_c,n1_c,n2_c,n3_c,N11_c(:,:,:,iset,jj),u_c_t(:,:,:,jj), &
             -2,n1_c+3,-2,n2_c+3,n3_c,n3_c,int_temp_c(-2:n1_c+3,-2:n2_c+3,1))
       eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)+int_temp_c(-2:n1_c+3,-2:n2_c+3,1)
    end do

    ! second derivatives 22
    int_temp_c = 0.d0
    do jj = 1,3
       call FD_r2r2(h2_c,n1_c,n2_c,n3_c,N22_c(:,:,:,iset,jj),u_c_t(:,:,:,jj), &
             -2,n1_c+3,-2,n2_c+3,n3_c,n3_c,int_temp_c(-2:n1_c+3,-2:n2_c+3,1))
       eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)+int_temp_c(-2:n1_c+3,-2:n2_c+3,1)
    end do

    ! second derivatives 33
    do k = 1,8
       do m = 1,8
          do jj = 1,3
             eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) =eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
              +acof(1,k,m)*N33_c(-2:n1_c+3,-2:n2_c+3,n3_c+1-m,iset,jj)*u_c_t(-2:n1_c+3,-2:n2_c+3,n3_c+1-k,jj)/h3_c**2
          end do
       end do
    end do
    do jj = 1,3
       eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) =eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
         + N33_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset,jj)*ghcof(1)*u_c_t(-2:n1_c+3,-2:n2_c+3,n3_c+1,jj)/h3_c**2
    end do

    ! mixed derivatives 12
    do jj = 1,3
       int_temp_c = 0.d0
       do i=-2,2
          int_temp_c(-4:n1_c+5,-2:n2_c+3,1)=int_temp_c(-4:n1_c+5,-2:n2_c+3,1) &
                    +ux_cof(i)*u_c_t(-4:n1_c+5,-2+i:n2_c+3+i,n3_c,jj)/h2_c
       end do
       do i=-2,2
          eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
             +ux_cof(i)*N12_c(-2+i:n1_c+3+i,-2:n2_c+3,n3_c,iset,jj)*int_temp_c(-2+i:n1_c+3+i,-2:n2_c+3,1)/h1_c
       end do
    end do

    ! mixed derivatives 13
    do jj = 1,3
       int_temp_c = 0.d0
       do i =1,4
          int_temp_c(-4:n1_c+5,-2:n2_c+3,1)=int_temp_c(-4:n1_c+5,-2:n2_c+3,1) &
                     -u_c_t(-4:n1_c+5,-2:n2_c+3,n3_c+1-i,jj)*bof(1,i)/h3_c
       end do
       do i =-2,2
          eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
             +ux_cof(i)*N13_c(-2+i:n1_c+3+i,-2:n2_c+3,n3_c,iset,jj)*int_temp_c(-2+i:n1_c+3+i,-2:n2_c+3,1)/h1_c
       end do
    end do

    ! mixed derivatives 21
    do jj = 1,3
       int_temp_c = 0.d0
       do i=-2,2
          int_temp_c(-2:n1_c+3,-4:n2_c+5,1)=int_temp_c(-2:n1_c+3,-4:n2_c+5,1) &
                    +ux_cof(i)*u_c_t(-2+i:n1_c+3+i,-4:n2_c+5,n3_c,jj)/h1_c
       end do
       do i=-2,2
          eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
             +ux_cof(i)*N21_c(-2:n1_c+3,-2+i:n2_c+3+i,n3_c,iset,jj)*int_temp_c(-2:n1_c+3,-2+i:n2_c+3+i,1)/h2_c
       end do
   end do

   ! mixed derivatives 23
   do jj = 1,3
      int_temp_c = 0.d0
      do i =1,4
         int_temp_c(-2:n1_c+3,-4:n2_c+5,1)=int_temp_c(-2:n1_c+3,-4:n2_c+5,1) &
                     -u_c_t(-2:n1_c+3,-4:n2_c+5,n3_c+1-i,jj)*bof(1,i)/h3_c
      end do
      do i =-2,2
         eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
             +ux_cof(i)*N23_c(-2:n1_c+3,-2+i:n2_c+3+i,n3_c,iset,jj)*int_temp_c(-2:n1_c+3,-2+i:n2_c+3+i,1)/h2_c
      end do
   end do

   ! mixed derivatives 31
   do jj = 1,3
      int_temp_c = 0.d0
      do i = -2,2
         int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4)=int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4) &
            +ux_cof(i)*u_c_t(-2+i:n1_c+3+i,-2:n2_c+3,n3_c:n3_c-3:-1,jj)/h1_c
      end do
      do i = 1,4
         eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
            -N31_c(-2:n1_c+3,-2:n2_c+3,n3_c-i+1,iset,jj)*int_temp_c(-2:n1_c+3,-2:n2_c+3,i)*bof(1,i)/h3_c
      end do
   end do

   ! mixed derivatives 32
   do jj = 1,3
      int_temp_c = 0.d0
      do i = -2,2
         int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4)=int_temp_c(-2:n1_c+3,-2:n2_c+3,1:4) &
            +ux_cof(i)*u_c_t(-2:n1_c+3,-2+i:n2_c+3+i,n3_c:n3_c-3:-1,jj)/h2_c
      end do
      do i = 1,4
         eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset)=eta_c(-2:n1_c+3,-2:n2_c+3,n3_c,iset) &
             -N32_c(-2:n1_c+3,-2:n2_c+3,n3_c-i+1,iset,jj)*int_temp_c(-2:n1_c+3,-2:n2_c+3,i)*bof(1,i)/h3_c
      end do
   end do

   ! project
   Mass_f1 = 0.d0
   do k = -1,n2_c+1
      do i = -1,n1_c+1
         if (k .gt. -1 .and. i .gt. -1) then
             Mass_f1(2*i-1,2*k-1)=eta_c(i,k,n3_c,iset)/rho_c(i,k,n3_c,1)
         end if
         if (k .ge. -1 .and. i .gt. -1) then
            do j = -1,2
               Mass_f1(2*i-1,2*k)=Mass_f1(2*i-1,2*k)+P(j)*eta_c(i,k+j,n3_c,iset)/rho_c(i,j+k,n3_c,1)
            end do
         end if
         if (k .gt. -1 .and. i .ge. -1) then
            do j = -1,2
               Mass_f1(2*i,2*k-1)=Mass_f1(2*i,2*k-1)+P(j)*eta_c(i+j,k,n3_c,iset)/rho_c(i+j,k,n3_c,1)
            end do
         end if
         do j = -1,2
            do l = -1,2
                Mass_f1(2*i,2*k)=Mass_f1(2*i,2*k)+P(j)*(P(l)*eta_c(i+l,k+j,n3_c,iset)/rho_c(i+l,k+j,n3_c,1))
             end do
          end do
        end do
    end do
    !
    eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) = -Mass_f1(-2:n1_f+3,-2:n2_f+3)*rho_f(-2:n1_f+3,-2:n2_f+3,1,1)

    ! second derivatives 11
    !eta_f = 0.d0
    int_temp_f = 0.d0
    do jj = 1,3
       call FD_r1r1(h1_f,n1_f,n2_f,n3_f,N11_f(:,:,:,iset,jj),u_f_t(:,:,:,jj), &
                    -2,n1_f+3,-2,n2_f+3,1,1,int_temp_f(-2:n1_f+3,-2:n2_f+3,1))
       eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) = eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)+int_temp_f(-2:n1_f+3,-2:n2_f+3,1)
    end do

    ! second derivatives 22
    int_temp_f = 0.d0
    do jj = 1,3
       call FD_r2r2(h2_f,n1_f,n2_f,n3_f,N22_f(:,:,:,iset,jj),u_f_t(:,:,:,jj), &
                    -2,n1_f+3,-2,n2_f+3,1,1,int_temp_f(-2:n1_f+3,-2:n2_f+3,1))
       eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) = eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)+int_temp_f(-2:n1_f+3,-2:n2_f+3,1)
    end do

    ! second derivatives 33
    do k = 1,8
       do m = 1,8
          do jj = 1,3
             eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) =eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
               +acof_no_gp(1,k,m)*N33_f(-2:n1_f+3,-2:n2_f+3,m,iset,jj)*u_f_t(-2:n1_f+3,-2:n2_f+3,k,jj)/h3_f**2
          end do
       end do
    end do

    ! mixed derivatives 12
    do jj = 1,3
       int_temp_f = 0.d0
       do i=-2,2
          int_temp_f(-4:n1_f+5,-2:n2_f+3,1)=int_temp_f(-4:n1_f+5,-2:n2_f+3,1) &
                        +ux_cof(i)*u_f_t(-4:n1_f+5,-2+i:n2_f+3+i,1,jj)/h2_f
       end do
       do i=-2,2
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
              +ux_cof(i)*N12_f(-2+i:n1_f+3+i,-2:n2_f+3,1,iset,jj)*int_temp_f(-2+i:n1_f+3+i,-2:n2_f+3,1)/h1_f
       end do
    end do

    ! mixed derivatives 13
    do jj = 1,3
       int_temp_f = 0.d0
       do i =1,4
          int_temp_f(-4:n1_f+5,-2:n2_f+3,1)=int_temp_f(-4:n1_f+5,-2:n2_f+3,1) &
                             +u_f_t(-4:n1_f+5,-2:n2_f+3,i,jj)*bof(1,i)/h3_f
       end do
       do i =-2,2
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
             +ux_cof(i)*N13_f(-2+i:n1_f+3+i,-2:n2_f+3,1,iset,jj)*int_temp_f(-2+i:n1_f+3+i,-2:n2_f+3,1)/h1_f
       end do
    end do

    ! mixed derivatives 21
    do jj = 1,3
       int_temp_f = 0.d0
       do i=-2,2
          int_temp_f(-2:n1_f+3,-4:n2_f+5,1)=int_temp_f(-2:n1_f+3,-4:n2_f+5,1) &
                        +ux_cof(i)*u_f_t(-2+i:n1_f+3+i,-4:n2_f+5,1,jj)/h1_f
       end do
       do i=-2,2
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
             +ux_cof(i)*N21_f(-2:n1_f+3,-2+i:n2_f+3+i,1,iset,jj)*int_temp_f(-2:n1_f+3,-2+i:n2_f+3+i,1)/h2_f
      end do
    end do

    ! mixed derivatives 23
    do jj = 1,3
       int_temp_f = 0.d0
       do i =1,4
          int_temp_f(-2:n1_f+3,-4:n2_f+5,1)=int_temp_f(-2:n1_f+3,-4:n2_f+5,1) &
                             +u_f_t(-2:n1_f+3,-4:n2_f+5,i,jj)*bof(1,i)/h3_f
       end do
       do i =-2,2
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
              +ux_cof(i)*N23_f(-2:n1_f+3,-2+i:n2_f+3+i,1,iset,jj)*int_temp_f(-2:n1_f+3,-2+i:n2_f+3+i,1)/h2_f
       end do
    end do

    ! mixed derivatives 31
    do jj = 1,3
       int_temp_f = 0.d0
       do i = -2,2
          int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4)=int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4) &
                      +ux_cof(i)*u_f_t(-2+i:n1_f+3+i,-2:n2_f+3,1:4,jj)/h1_f
       end do
       do i = 1,4
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
             +N31_f(-2:n1_f+3,-2:n2_f+3,i,iset,jj)*int_temp_f(-2:n1_f+3,-2:n2_f+3,i)*bof(1,i)/h3_f
       end do
    end do

    ! mixed derivatives 32
    do jj = 1,3
       int_temp_f = 0.d0
       do i = -2,2
          int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4)=int_temp_f(-2:n1_f+3,-2:n2_f+3,1:4) &
                      +ux_cof(i)*u_f_t(-2:n1_f+3,-2+i:n2_f+3+i,1:4,jj)/h2_f
       end do
       do i = 1,4
          eta_f(-2:n1_f+3,-2:n2_f+3,1,iset)=eta_f(-2:n1_f+3,-2:n2_f+3,1,iset) &
              +N32_f(-2:n1_f+3,-2:n2_f+3,i,iset,jj)*int_temp_f(-2:n1_f+3,-2:n2_f+3,i)*bof(1,i)/h3_f
       end do
    end do
    !
    end do ! iset
    eta_f = -eta_f   !!! eta does not divided by Jacobian_f as in paper
  end subroutine compute_eta

  subroutine Update_interior(u_c_t,u_f_t)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: u_c_t
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: u_f_t

    ! Difference operators in the interior of the domains
    ! fine mesh
    lh_f = 0.d0
    do iset = 1,3
      do jj = 1,3
         call FD_r1r1(h1_f,n1_f,n2_f,n3_f,N11_f(:,:,:,iset,jj),u_f_t(:,:,:,jj),1,n1_f,1,n2_f,1,n3_f,G1_f)
         call FD_r2r2(h2_f,n1_f,n2_f,n3_f,N22_f(:,:,:,iset,jj),u_f_t(:,:,:,jj),1,n1_f,1,n2_f,1,n3_f,G2_f)
         call FD_r3r3_no_gp(n1_f,n2_f,n3_f,h3_f,N33_f(:,:,:,iset,jj),u_f_t(:,:,:,jj),G3_f,acof_no_gp,acof,ghcof)

         call FD_r3(h3_f,-1,n1_f+2,1,n2_f,1,n3_f,u_f_t(-1:n1_f+2,1:n2_f,1:n3_f,jj),D3_1f,bof)
         call FD_r1(h1_f,-1,n1_f+2,1,n2_f,1,n3_f,N13_f(-1:n1_f+2,1:n2_f,1:n3_f,iset,jj)*D3_1f,D31_f)

         call FD_r2(h2_f,-1,n1_f+2,-1,n2_f+2,1,n3_f,u_f_t(-1:n1_f+2,-1:n2_f+2,1:n3_f,jj),D2_1f)
         call FD_r1(h1_f,-1,n1_f+2,1,n2_f,1,n3_f,N12_f(-1:n1_f+2,1:n2_f,1:n3_f,iset,jj)*D2_1f,D21_f)

         call FD_r1(h1_f,-1,n1_f+2,-1,n2_f+2,1,n3_f,u_f_t(-1:n1_f+2,-1:n2_f+2,1:n3_f,jj),D1_2f)
         call FD_r2(h2_f,1,n1_f,-1,n2_f+2,1,n3_f,N21_f(1:n1_f,-1:n2_f+2,1:n3_f,iset,jj)*D1_2f,D12_f)

         call FD_r3(h3_f,1,n1_f,-1,n2_f+2,1,n3_f,u_f_t(1:n1_f,-1:n2_f+2,1:n3_f,jj),D3_2f,bof)
         call FD_r2(h2_f,1,n1_f,-1,n2_f+2,1,n3_f,N23_f(1:n1_f,-1:n2_f+2,1:n3_f,iset,jj)*D3_2f,D32_f)

         call FD_r1(h1_f,-1,n1_f+2,1,n2_f,1,n3_f,u_f_t(-1:n1_f+2,1:n2_f,1:n3_f,jj),D1_3f)
         call FD_r3(h3_f,1,n1_f,1,n2_f,1,n3_f,N31_f(1:n1_f,1:n2_f,1:n3_f,iset,jj)*D1_3f,D13_f,bof)

         call FD_r2(h2_f,1,n1_f,-1,n2_f+2,1,n3_f,u_f_t(1:n1_f,-1:n2_f+2,1:n3_f,jj),D2_3f)
         call FD_r3(h3_f,1,n1_f,1,n2_f,1,n3_f,N32_f(1:n1_f,1:n2_f,1:n3_f,iset,jj)*D2_3f,D23_f,bof)

         lh_f(1:n1_f,1:n2_f,1:n3_f,iset) = lh_f(1:n1_f,1:n2_f,1:n3_f,iset)+G1_f+G2_f+G3_f+D31_f+D21_f+D12_f+D32_f+D13_f+D23_f
      end do
    end do


    ! coarse mesh
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

  subroutine Update_gp(index)
    integer index

    ! Update ghost point values on the left and right domain
    ! fine mesh
    do i=1-nrg,n3_f+nrg
       do j = 1-nrg,n2_f+nrg
          do k = 1-nrg,1
             call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
          do k = n1_f,n1_f+nrg
             call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
    end do
    ! coarse fine
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,n2_c+nrg
          do k = 1-nrg,1
             call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,i,3),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
          do k = n1_c,n1_c+nrg
             call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,i,3),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do
    ! Update ghost point values on the front and back domain
    ! fine mesh
    do i=1-nrg,n3_f+nrg
       do j = 1-nrg,1
          do k = 1-nrg,n1_f+nrg
             call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
       do j = n2_f,n2_f+nrg
          do k = 1-nrg,n1_f+nrg
             call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                     u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
    end do
    ! coarse fine
    do i=1-nrg,n3_c+nrg
       do j = 1-nrg,1
          do k = 1-nrg,n1_c+nrg
             call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,i,3),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
       do j = n2_c,n2_c+nrg
          do k = 1-nrg,n1_c+nrg
             call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,i,3),tv, &
                     u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
       end do
    end do

  end subroutine Update_gp

  subroutine Update_Dirichlet_BC(index)
    integer index
    ! fine mesh
    do j=1-nrg,n2_f+nrg
       do k=1-nrg,n1_f+nrg
          !i = 1
          !call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                    !u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          do i = n3_f,n3_f+nrg
             call exact_solution(Xgrid_f(k,j,i,1),Xgrid_f(k,j,i,2),Xgrid_f(k,j,i,3),tv, &
                    u_f(k,j,i,1,index),u_f(k,j,i,2,index),u_f(k,j,i,3,index),1)
          end do
       end do
    end do
    !coarse mesh
    do j=1-nrg,n2_c+nrg
       do k=1-nrg,n1_c+nrg
          do i = 1-nrg,1
             call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,i,3),tv, &
                   u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
          end do
          !i = n3_c
          !call exact_solution(Xgrid_c(k,j,i,1),Xgrid_c(k,j,i,2),Xgrid_c(k,j,3),tv, &
                   !u_c(k,j,i,1,index),u_c(k,j,i,2,index),u_c(k,j,i,3,index),0)
       end do
    end do
  end subroutine Update_Dirichlet_BC

  subroutine Injection(index)
    integer index
    ! Injection at the interface
    do j = 1, n2_c
       do i = 1, n1_c
          u_f(2*i-1,2*j-1,1,:,index) = u_c(i,j,n3_c,:,index)
          u_f(2*i,2*j-1,1,:,index) = P(-1)*u_c(i-1,j,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
                                   +P(1)*u_c(i+1,j,n3_c,:,index)+P(2)*u_c(i+2,j,n3_c,:,index)
          u_f(2*i-1,2*j,1,:,index) = P(-1)*u_c(i,j-1,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
                                   +P(1)*u_c(i,j+1,n3_c,:,index)+P(2)*u_c(i,j+2,n3_c,:,index)
          u_f(2*i,2*j,1,:,index) = P(-1)*(P(-1)*u_c(i-1,j-1,n3_c,:,index)+P(0)*u_c(i,j-1,n3_c,:,index)&
                                           +P(1)*u_c(i+1,j-1,n3_c,:,index)+P(2)*u_c(i+2,j-1,n3_c,:,index))&
                                     +P(0)*(P(-1)*u_c(i-1,j,n3_c,:,index)+P(0)*u_c(i,j,n3_c,:,index)&
                                           +P(1)*u_c(i+1,j,n3_c,:,index)+P(2)*u_c(i+2,j,n3_c,:,index))&
                                     +P(1)*(P(-1)*u_c(i-1,j+1,n3_c,:,index)+P(0)*u_c(i,j+1,n3_c,:,index)&
                                           +P(1)*u_c(i+1,j+1,n3_c,:,index)+P(2)*u_c(i+2,j+1,n3_c,:,index))&
                                     +P(2)*(P(-1)*u_c(i-1,j+2,n3_c,:,index)+P(0)*u_c(i,j+2,n3_c,:,index)&
                                           +P(1)*u_c(i+1,j+2,n3_c,:,index)+P(2)*u_c(i+2,j+2,n3_c,:,index))
       end do
    end do
  end subroutine Injection

  subroutine discrete_energy(uf_oldold,uf_old,uf_star,uf_new,uc_oldold,uc_old,uc_star,uc_new,dt,energy_num,energy_num_bdry)
    real(dp), dimension(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim) :: uc_oldold,uc_old,uc_new,uc_star
    real(dp), dimension(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim) :: uf_oldold,uf_old,uf_new,uf_star
    real(dp), dimension(1:n1_c,1:n2_c,1:n3_c) :: energy_num_temp_c
    real(dp), dimension(1:n1_f,1:n2_f,1:n3_f) :: energy_num_temp_f
    real(dp), dimension(1:n1_c,1:n2_c,1:n3_c,1:3) :: Lh_old_c
    real(dp), dimension(1:n1_f,1:n2_f,1:n3_f,1:3) :: Lh_old_f
    real(dp), dimension(1:n1_c,1:n2_c,1:3) :: boundary_c
    real(dp), dimension(1:n1_f,1:n2_f,1:3) :: boundary_f
    real(dp) :: dt,energy_num,energy_num_bdry

    ! coarse domain
    energy_num_temp_c = 0.d0
    ! term w.r.t time derivative
    energy_num_temp_c = energy_num_temp_c + rho_c(1:n1_c,1:n2_c,1:n3_c,1) &
        * ((uc_new(1:n1_c,1:n2_c,1:n3_c,1)-uc_old(1:n1_c,1:n2_c,1:n3_c,1))**2 &
        + (uc_new(1:n1_c,1:n2_c,1:n3_c,2)-uc_old(1:n1_c,1:n2_c,1:n3_c,2))**2 &
        + (uc_new(1:n1_c,1:n2_c,1:n3_c,3)-uc_old(1:n1_c,1:n2_c,1:n3_c,3))**2)/(dt**2)

	  Lh_old_c = rho_c(1:n1_c,1:n2_c,1:n3_c,:)* &
	      (uc_new(1:n1_c,1:n2_c,1:n3_c,:) - 2.d0*uc_old(1:n1_c,1:n2_c,1:n3_c,:)&
        +uc_oldold(1:n1_c,1:n2_c,1:n3_c,:))/dt**2

	  ! term w.r.t S_h(unew,uold) = -(unew,L_h uold)
	  energy_num_temp_c = energy_num_temp_c - (uc_new(1:n1_c,1:n2_c,1:n3_c,1)*Lh_old_c(:,:,:,1) &
       + uc_new(1:n1_c,1:n2_c,1:n3_c,2)*Lh_old_c(:,:,:,2) &
       + uc_new(1:n1_c,1:n2_c,1:n3_c,3)*Lh_old_c(:,:,:,3))

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
        energy_num = energy_num + h1_c*h2_c*h3_c &
          * (17.d0/48.d0*(energy_num_temp_c(i,j,1) + energy_num_temp_c(i,j,n3_c)) &
          + 59.d0/48.d0*(energy_num_temp_c(i,j,2) + energy_num_temp_c(i,j,n3_c-1)) &
          + 43.d0/48.d0*(energy_num_temp_c(i,j,3) + energy_num_temp_c(i,j,n3_c-2)) &
          + 49.d0/48.d0*(energy_num_temp_c(i,j,4) + energy_num_temp_c(i,j,n3_c-3)))
      end do
    end do
    !
    ! interface integral
    energy_num_bdry = 0.d0
    boundary_c = 0.d0
    do l = 1,n2_c
       do i = 1,n1_c
          do i1 = 1,3
             do k = 1,3
                do j = -2,2
                   boundary_c(i,l,i1) = boundary_c(i,l,i1) &
                    + N31_c(i+j,l,n3_c,i1,k)*ux_cof(j)*uc_old(i+j,l,n3_c,k)/h1_c &
                    + N32_c(i,l+j,n3_c,i1,k)*ux_cof(j)*uc_old(i,l+j,n3_c,k)/h1_c
                end do
                do j = 0,4
                   boundary_c(i,l,i1) = boundary_c(i,l,i1) + N33_c(i,l,n3_c+1-j,i1,k)*Sb(j)*uc_old(i,l,n3_c+1-j,k)/h1_c
                end do
             end do
             energy_num_bdry = energy_num_bdry + h1_c*h2_c*uc_new(i,l,n3_c,i1)*boundary_c(i,l,i1)
          end do
       end do
    end do
    !
    ! fine domain
    energy_num_temp_f = 0.d0
    ! term w.r.t time derivative
    energy_num_temp_f = energy_num_temp_f + rho_f(1:n1_f,1:n2_f,1:n3_f,1) &
        * ((uf_new(1:n1_f,1:n2_f,1:n3_f,1)-uf_old(1:n1_f,1:n2_f,1:n3_f,1))**2 &
        + (uf_new(1:n1_f,1:n2_f,1:n3_f,2)-uf_old(1:n1_f,1:n2_f,1:n3_f,2))**2 &
        + (uf_new(1:n1_f,1:n2_f,1:n3_f,3)-uf_old(1:n1_f,1:n2_f,1:n3_f,3))**2)/(dt**2)

    Lh_old_f = rho_f(1:n1_f,1:n2_f,1:n3_f,:)* &
         (uf_new(1:n1_f,1:n2_f,1:n3_f,:) - 2.d0*uf_old(1:n1_f,1:n2_f,1:n3_f,:)&
         +uf_oldold(1:n1_f,1:n2_f,1:n3_f,:))/dt**2

     !term w.r.t S_h(unew,uold) = -(unew,L_h uold)
     energy_num_temp_f = energy_num_temp_f - (uf_new(1:n1_f,1:n2_f,1:n3_f,1)*Lh_old_f(:,:,:,1) &
       + uf_new(1:n1_f,1:n2_f,1:n3_f,2)*Lh_old_f(:,:,:,2) &
       + uf_new(1:n1_f,1:n2_f,1:n3_f,3)*Lh_old_f(:,:,:,3))


    ! add grid points with corresponding weights
    do k = 5,n3_f-4
      do j = 1,n2_f
        do i = 1,n1_f
          energy_num = energy_num + h1_f*h2_f*h3_f*energy_num_temp_f(i,j,k)
        end do
      end do
    end do
    !
    do j = 1,n2_f
      do i = 1,n1_f
        energy_num = energy_num + h1_f*h2_f*h3_f &
          * (17.d0/48.d0*(energy_num_temp_f(i,j,1)+ energy_num_temp_f(i,j,n3_f)) &
          + 59.d0/48.d0*(energy_num_temp_f(i,j,2) + energy_num_temp_f(i,j,n3_f-1)) &
          + 43.d0/48.d0*(energy_num_temp_f(i,j,3) + energy_num_temp_f(i,j,n3_f-2)) &
          + 49.d0/48.d0*(energy_num_temp_f(i,j,4) + energy_num_temp_f(i,j,n3_f-3)))
      end do
    end do
    ! interface integral
    boundary_f = 0.d0
    do l = 1,n2_f
       do i = 1,n1_f
          do i1 = 1,3
             do k = 1,3
                do j = -2,2
                   boundary_f(i,l,i1) = boundary_f(i,l,i1) &
                    + N31_f(i+j,l,1,i1,k)*ux_cof(j)*uf_old(i+j,l,1,k)/h1_f &
                    + N32_f(i,l+j,1,i1,k)*ux_cof(j)*uf_old(i,l+j,1,k)/h2_f
                end do
                do j = 0,5
                   boundary_f(i,l,i1) = boundary_f(i,l,i1) + N33_f(i,l,j,i1,k)*sbop_no_gp(j)*uf_old(i,l,j,k)/h3_f
                end do
             end do
             energy_num_bdry = energy_num_bdry - h1_f*h2_f*uf_new(i,l,1,i1)*boundary_f(i,l,i1)
          end do
       end do
    end do
    ! term w.r.t eta
    !energy_num_bdry = 0.d0
    do l = 1,n2_f
       do i = 1,n1_f
          energy_num_bdry = energy_num_bdry + h3_f*17.d0/48.d0*h2_f*h1_f &
          *(uf_new(i,l,1,1)*eta_f(i,l,1,1)+uf_new(i,l,1,2)*eta_f(i,l,1,2)+uf_new(i,l,1,3)*eta_f(i,l,1,3))
      end do
    end do
    ! add terms together
    energy_num = energy_num !+ energy_num_bdry
  end subroutine discrete_energy

end program energy_2nd

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
           write(7,"(ES25.15E3)") A(k,j,i)
        end do
     end do
  end do
  close(7)
end  subroutine print_array_to_file
