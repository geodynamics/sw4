program condition_number

  use problemsetup_new_3d
  use SBP_operator

  implicit none

  integer i,j,k,i1,j1,k1,l,m,time_index

  ! Physical domain X
  real(dp), dimension (:), allocatable :: Xgrid_c_1, Xgrid_f_1, Xgrid_c_2, Xgrid_f_2
  real(dp), dimension (:,:,:), allocatable :: Xgrid_c_3, Xgrid_f_3

  ! Metric derivatives
  real(dp), dimension (:,:,:), allocatable :: XI13_c, XI23_c, XI33_c, Jacobian_c
  real(dp), dimension (:,:,:), allocatable :: XI13_f, XI23_f, XI33_f, Jacobian_f

  ! material parameters
  real(dp), dimension (:,:,:), allocatable :: rho_c, mu_c, lambda_c
  real(dp), dimension (:,:,:), allocatable :: rho_f, mu_f, lambda_f

  ! Solution
  real(dp), dimension (:,:,:,:,:), allocatable :: u_c, u_f

  ! Error
  real(dp), dimension (1:6) :: N6

  ! interface system
  real(dp), dimension (:,:), allocatable :: Mass,Mass0,Mass_pre0,Mass_temp
  real(dp), dimension (:), allocatable :: pdir,pdir1,LHS
  real(dp), dimension (:,:,:,:), allocatable :: Mass_pre
  integer, dimension (:,:,:), allocatable :: IPIV_pre
  integer, dimension (:), allocatable :: IPIV, iwork
  real(dp), dimension (:), allocatable :: work
  integer :: INFO
  real(dp) :: anorm,colsum,rcon

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
  allocate(Xgrid_f_1(1-nrg:n1_f+nrg))
  allocate(Xgrid_f_2(1-nrg:n2_f+nrg))
  allocate(Xgrid_c_3(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(Xgrid_f_3(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))

  ! allocate memory for metric derivatives
  allocate(XI13_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(XI23_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(XI33_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(Jacobian_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(XI13_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))
  allocate(XI23_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))
  allocate(XI33_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))
  allocate(Jacobian_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))

  ! allocate memory for material parameters
  allocate(rho_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(mu_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(lambda_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg))
  allocate(rho_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))
  allocate(mu_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))
  allocate(lambda_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg))

  ! allocate memory for solutions
  allocate(u_c(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg,1:dim,1:4))
  allocate(u_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg,1:dim,1:4))

  ! allocate memory for interface linear system
  allocate(Mass(1:3*n1_c*n2_c,1:3*n1_c*n2_c))
  allocate(Mass0(1:3*n1_c*n2_c,1:3*n1_c*n2_c))
  allocate(Mass_temp(1:3*n1_c*n2_c,1:3*n1_c*n2_c))
  allocate(pdir(1:n1_c*n2_c*3))
  allocate(pdir1(1:n1_c*n2_c*3))
  allocate(LHS(1:n1_c*n2_c*3))
  allocate(Mass_pre(1:3,1:3,1:n1_c,1:n2_c))
  allocate(Mass_pre0(1:3*n1_c*n2_c,1:3*n1_c*n2_c))
  allocate(IPIV_pre(1:3,1:n1_c,1:n2_c))
  allocate(IPIV(1:3*n1_c*n2_c))
  allocate(work(1:4*3*n1_c*n2_c))
  allocate(iwork(3*n1_c*n2_c))


  ! Difference stencils and interpolations
  call interpolation_restriction(P,Rop,RPop)
  call central_difference(ux_cof,uxx_cof)
  call varcoeff_NoGhost(acof_no_gp, ghcof_no_gp, sbop_no_gp )
  call VARCOEFFS4( acof, ghcof, Sb )
  call dx_46(bof)

  ! metric derivatives and coefficients of the equation
  call  equation_cof(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3, &
       XI13_c,XI23_c,XI33_c,Jacobian_c,rho_c,rho_f, &
       XI13_f,XI23_f,XI33_f,Jacobian_f,&
       mu_c,mu_f,lambda_c,lambda_f)

  !do k=1-nrg,n3_c+nrg
  !   do i=1-nrg,n2_c+nrg
  !      do j=1-nrg,n1_c+nrg
  !         call exact_solution(Xgrid_c_1(j),Xgrid_c_2(i),Xgrid_c_3(j,i,k), &
  !              0.d0,u_c(j,i,k,1,1),u_c(j,i,k,2,1),u_c(j,i,k,3,1),0)
  !      end do
  !   end do
  !end do

  ! Construct the system matrix for computing ghost points values on the interface
  ! We have 3*n1_c*n2_c equations in 3D
  ! There are three sets of equations, one for the first component of u, one for the
  ! second component and another for the third component of u
  ! Each set consists of n1_c*n2_c equations

  ! interface system
  call Interface_system(Mass)
  Mass0 = Mass
  !
  !pdir1 = 0.d0
  !do j = 1,n2_c
  !  do i = 1,n1_c
  !    do l = 1,n2_c
  !      do k = 1,n1_c
  !        pdir1((j-1)*3*n1_c+(i-1)*3+1) = pdir1((j-1)*3*n1_c+(i-1)*3+1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+1,(l-1)*3*n1_c+(k-1)*3+1)*u_c(k,l,n3_c+1,1,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+1,(l-1)*3*n1_c+(k-1)*3+2)*u_c(k,l,n3_c+1,2,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+1,(l-1)*3*n1_c+(k-1)*3+3)*u_c(k,l,n3_c+1,3,1)
          !
  !        pdir1((j-1)*3*n1_c+(i-1)*3+2) = pdir1((j-1)*3*n1_c+(i-1)*3+2) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+2,(l-1)*3*n1_c+(k-1)*3+1)*u_c(k,l,n3_c+1,1,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+2,(l-1)*3*n1_c+(k-1)*3+2)*u_c(k,l,n3_c+1,2,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+2,(l-1)*3*n1_c+(k-1)*3+3)*u_c(k,l,n3_c+1,3,1)
          !
  !        pdir1((j-1)*3*n1_c+(i-1)*3+3) = pdir1((j-1)*3*n1_c+(i-1)*3+3) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+3,(l-1)*3*n1_c+(k-1)*3+1)*u_c(k,l,n3_c+1,1,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+3,(l-1)*3*n1_c+(k-1)*3+2)*u_c(k,l,n3_c+1,2,1) &
  !          + Mass((j-1)*3*n1_c+(i-1)*3+3,(l-1)*3*n1_c+(k-1)*3+3)*u_c(k,l,n3_c+1,3,1)
  !      end do
  !    end do
  !  end do
  !end do

  ! construct preconditioner
  call Interface_preconditioner(Mass_pre)
  !
  do j = 1,n2_c
    do i = 1,n1_c
      Mass_pre0((j-1)*3*n1_c+(i-1)*3+1:(j-1)*3*n1_c+(i-1)*3+3, &
         (j-1)*3*n1_c+(i-1)*3+1:(j-1)*3*n1_c+(i-1)*3+3) = Mass_pre(:,:,i,j)
    end do
  end do
  !
  do j = 1,n2_c
     do i = 1,n1_c
        call dgetrf(3,3,Mass_pre(:,:,i,j),3,IPIV_pre(:,i,j),INFO)
        if (INFO .ne. 0) then
           write(*,"(A20,I5)") 'LU fails at (i,j) equals', i,j, 'INFO=', INFO
        end if
     end do
  end do

  !do j = 1,n2_c
  !   do i = 1,n1_c
  !      pdir((j-1)*3*n1_c+(i-1)*3+1) = u_c(i,j,n3_c+1,1,1)
  !      pdir((j-1)*3*n1_c+(i-1)*3+2) = u_c(i,j,n3_c+1,2,1)
  !      pdir((j-1)*3*n1_c+(i-1)*3+3) = u_c(i,j,n3_c+1,3,1)
  !   end do
  !end do
  !call Interface_LHS_cg(pdir)
  !print *, pdir1-LHS

  ! pre-conditioned matrix * coefficient matrx
  do j = 1,n2_c
    do i = 1,n1_c
      do l = 1,n2_c
        do k = 1,n1_c
          call dgetrs('N',3,1,Mass_pre(:,:,i,j),3,IPIV_pre(:,i,j),&
           Mass((j-1)*3*n1_c+3*(i-1)+1:(j-1)*3*n1_c+3*(i-1)+3,(l-1)*3*n1_c+3*(k-1)+1),3,INFO)
          if (INFO .ne. 0) then
           write(*,"(A20)") 'Solving fails at (i,j,k,l,ele) = ', i,j,k,l,1
           stop
          end if
          !
          call dgetrs('N',3,1,Mass_pre(:,:,i,j),3,IPIV_pre(:,i,j),&
           Mass((j-1)*3*n1_c+3*(i-1)+1:(j-1)*3*n1_c+3*(i-1)+3,(l-1)*3*n1_c+3*(k-1)+2),3,INFO)
          if (INFO .ne. 0) then
           write(*,"(A20)") 'Solving fails at (i,j,k,l,ele) = ', i,j,k,l,2
           stop
          end if
          !
          call dgetrs('N',3,1,Mass_pre(:,:,i,j),3,IPIV_pre(:,i,j),&
           Mass((j-1)*3*n1_c+3*(i-1)+1:(j-1)*3*n1_c+3*(i-1)+3,(l-1)*3*n1_c+3*(k-1)+3),3,INFO)
          if (INFO .ne. 0) then
           write(*,"(A20)") 'Solving fails at (i,j,k,l,ele) = ', i,j,k,l,3
           stop
          end if
        end do
      end do
    end do
  end do
  !
  ! compute the condition number of the matrix (original matrix, block Jacobian, preconditioned origianl matrix)
  ! one norm of the matrix
  Mass_temp = Mass0
  anorm = 0.d0
  do j = 1,3*n1_c*n2_c
    colsum = 0.d0
    do i = 1,3*n1_c*n2_c
      colsum = colsum + abs(Mass_temp(i,j))
    end do
    anorm = max(anorm,colsum)
  end do
  call dgetrf(3*n1_c*n2_c,3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,IPIV,INFO)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgetrf fails at Mass'
    stop
  end if
  call dgecon('1',3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,anorm,rcon,work,iwork,info)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgecon fails at Mass'
    stop
  end if
  print *, 'condition number of the original coefficient matrix = ',1.d0/rcon
  !
  Mass_temp = Mass_pre0
  anorm = 0.d0
  do j = 1,3*n1_c*n2_c
    colsum = 0.d0
    do i = 1,3*n1_c*n2_c
      colsum = colsum + abs(Mass_temp(i,j))
    end do
    anorm = max(anorm,colsum)
  end do
  call dgetrf(3*n1_c*n2_c,3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,IPIV,INFO)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgetrf fails at Mass'
    stop
  end if
  call dgecon('1',3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,anorm,rcon,work,iwork,info)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgecon fails at Mass'
    stop
  end if
  print *, 'condition number of the block jacobian matrix = ',1.d0/rcon
  !
  Mass_temp = Mass
  anorm = 0.d0
  do j = 1,3*n1_c*n2_c
    colsum = 0.d0
    do i = 1,3*n1_c*n2_c
      colsum = colsum + abs(Mass_temp(i,j))
    end do
    anorm = max(anorm,colsum)
  end do
  call dgetrf(3*n1_c*n2_c,3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,IPIV,INFO)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgetrf fails at Mass'
    stop
  end if
  call dgecon('1',3*n1_c*n2_c,Mass_temp,3*n1_c*n2_c,anorm,rcon,work,iwork,info)
  if (INFO .ne. 0) then
    write(*,"(A20)") 'dgecon fails at Mass'
    stop
  end if
  print *, 'condition number of the original coefficient matrix after preconditioner',1.d0/rcon

contains

  subroutine equation_cof(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3, &
       XI13_c,XI23_c,XI33_c,Jacobian_c,rho_c,rho_f,XI13_f,XI23_f,XI33_f,Jacobian_f,&
       mu_c,mu_f,lambda_c,lambda_f)
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: rho_c,mu_c,lambda_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: rho_f,mu_f,lambda_f
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: Xgrid_c_3
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: Xgrid_f_3
    real(dp), dimension (1-nrg:n1_c+nrg) :: Xgrid_c_1
    real(dp), dimension (1-nrg:n1_f+nrg) :: Xgrid_f_1
    real(dp), dimension (1-nrg:n2_c+nrg) :: Xgrid_c_2
    real(dp), dimension (1-nrg:n2_f+nrg) :: Xgrid_f_2
    real(dp), dimension (1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1-nrg:n3_c+nrg) :: XI13_c,XI23_c,XI33_c,Jacobian_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg,1-nrg:n3_f+nrg) :: XI13_f,XI23_f,XI33_f,Jacobian_f


    call generate_grid(Xgrid_c_1,Xgrid_c_2,Xgrid_c_3,Xgrid_f_1,Xgrid_f_2,Xgrid_f_3)


    ! compute metric derivatives, coarse
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             call metric_derivative(dble(k-1)/(n1_c-1),dble(j-1)/(n2_c-1),dble(i-1)/(n3_c-1),&
                 XI13_c(k,j,i),XI23_c(k,j,i),XI33_c(k,j,i),Jacobian_c(k,j,i),0)
          end do
       end do
    end do
    ! compute metric derivatives, fine
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
              call metric_derivative(dble(k-1)/(n1_f-1),dble(j-1)/(n2_f-1),dble(i-1)/(n3_f-1),&
                  XI13_f(k,j,i),XI23_f(k,j,i),XI33_f(k,j,i),Jacobian_f(k,j,i),1)
          end do
       end do
    end do
    ! variable coefficients
    do i=1-nrg,n3_c+nrg
       do j=1-nrg,n2_c+nrg
          do k=1-nrg,n1_c+nrg
             mu_c(k,j,i) = 3.d0 + sin(3.d0*Xgrid_c_1(k)+0.1d0)*sin(3.d0*Xgrid_c_2(j)+0.1d0)*sin(Xgrid_c_3(k,j,i))
             lambda_c(k,j,i) = 21.d0+ cos(Xgrid_c_1(k)+0.1d0)*cos(Xgrid_c_2(j)+0.1d0)*sin(3.d0*Xgrid_c_3(k,j,i))**2
             rho_c(k,j,i) = 2.d0 + sin(Xgrid_c_1(k)+0.3d0)*sin(Xgrid_c_2(j)+0.3d0)*sin(Xgrid_c_3(k,j,i)-0.2d0)
          end do
       end do
    end do
    do i=1-nrg,n3_f+nrg
       do j=1-nrg,n2_f+nrg
          do k=1-nrg,n1_f+nrg
             mu_f(k,j,i) = 3.d0 + sin(3.d0*Xgrid_f_1(k)+0.1d0)*sin(3.d0*Xgrid_f_2(j)+0.1d0)*sin(Xgrid_f_3(k,j,i))
             lambda_f(k,j,i) = 21.d0+ cos(Xgrid_f_1(k)+0.1d0)*cos(Xgrid_f_2(j)+0.1d0)*sin(3.d0*Xgrid_f_3(k,j,i))**2
             rho_f(k,j,i) = 2.d0 + sin(Xgrid_f_1(k)+0.3d0)*sin(Xgrid_f_2(j)+0.3d0)*sin(Xgrid_f_3(k,j,i)-0.2d0)
          end do
       end do
    end do

    rho_c = rho_c*Jacobian_c
    rho_f = rho_f*Jacobian_f

  end subroutine  equation_cof
  !
  subroutine Interface_preconditioner(Mass_pre)
    !
    real(dp), dimension (1:3,1:3,1:n1_c,1:n2_c) :: Mass_pre
    real(dp) :: int_cof
    real(dp), dimension (:,:), allocatable :: int_cof_c, int_cof_f
    !
    allocate(int_cof_c(1:n1_c,1:n2_c))
    allocate(int_cof_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg))
    !
    do j = 1,n2_c
       do i = 1,n1_c
          int_cof_c(i,j) = Jacobian_c(i,j,n3_c)*sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)
       end do
    end do
    !
    do j = 1-nrg,n2_f+nrg
       do i = 1-nrg,n1_f+nrg
          int_cof_f(i,j) = Jacobian_f(i,j,1)*sqrt(XI13_f(i,j,1)**2+XI23_f(i,j,1)**2+XI33_f(i,j,1)**2)
       end do
    end do
    !
    int_cof = 17.d0/48.d0*h3_f*ghcof(1)/h3_c**2
    !
    Mass_pre = 0.d0
    do l = 1,n2_c
       do k = 1,n1_c
          !
          do j = -4,2,2
             do i = -4,2,2
                ! first set equation w.r.t the first component
                Mass_pre(1,1,k,l) = Mass_pre(1,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                  * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
                  *XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))&
                  /rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                  ! first set equation w.r.t the second component
                Mass_pre(1,2,k,l) = Mass_pre(1,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                  * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
                  *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                  ! first set equation w.r.t the third component
                Mass_pre(1,3,k,l) = Mass_pre(1,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                  * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
                  *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                  ! second set equation w.r.t the first component
                Mass_pre(2,1,k,l) = Mass_pre(2,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                  * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
                  +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                ! second set equation w.r.t the second component
                Mass_pre(2,2,k,l) = Mass_pre(2,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                   * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
                   +lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
                   +XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                ! second set equation w.r.t the third component
                Mass_pre(2,3,k,l) = Mass_pre(2,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1) &
                   * P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
                   +mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l+j)
                ! third set equation w.r.t the first component
                Mass_pre(3,1,k,l) = Mass_pre(3,1,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
                   /int_cof_f(2*k+i,2*l+j)
                ! third set equation w.r.t the second component
                Mass_pre(3,2,k,l) = Mass_pre(3,2,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
                   /int_cof_f(2*k+i,2*l+j)
                ! third set equation w.r.t the third component
                Mass_pre(3,3,k,l) = Mass_pre(3,3,k,l)+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(-j/2)*P(-i/2)*Jacobian_c(k,l,n3_c)&
                   *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
                   +XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof&
                   /int_cof_f(2*k+i,2*l+j)
             end do
          end do
          !
          do j = -4,2,2
             ! first set equation w.r.t the first component
             Mass_pre(1,1,k,l) = Mass_pre(1,1,k,l)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
               +lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! first set equation w.r.t the second component
             Mass_pre(1,2,k,l) = Mass_pre(1,2,k,l)+ Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
               +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! first set equation w.r.t the third component
             Mass_pre(1,3,k,l) = Mass_pre(1,3,k,l)+ Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
               +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! second set equation w.r.t the first component
             Mass_pre(2,1,k,l) = Mass_pre(2,1,k,l)+ Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
               *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! second set equation w.r.t the second component
             Mass_pre(2,2,k,l) = Mass_pre(2,2,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
               +lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! second set equation w.r.t the third component
             Mass_pre(2,3,k,l) = Mass_pre(2,3,k,l) + Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
               *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! third set equation w.r.t the first component
             Mass_pre(3,1,k,l) = Mass_pre(3,1,k,l)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
               *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! third set equation w.r.t the second component
             Mass_pre(3,2,k,l) = Mass_pre(3,2,k,l)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
               * P(-j/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
               *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
             ! third set equation w.r.t the third component
             Mass_pre(3,3,k,l) = Mass_pre(3,3,k,l)+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
               * P(-j/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
               *XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
               +XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l+j)
          end do
          !
          do i = -4,2,2
             ! first set equation w.r.t the first component
             Mass_pre(1,1,k,l) = Mass_pre(1,1,k,l)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1) &
               * P(-i/2)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
               +lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
               *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l-1)
             ! first set equation w.r.t the second component
             Mass_pre(1,2,k,l) = Mass_pre(1,2,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1) &
               * P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
               +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l-1)
             ! first set equation w.r.t the third component
             Mass_pre(1,3,k,l) = Mass_pre(1,3,k,l) + Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1) &
               * P(-i/2)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
               +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l-1)
             ! second set equation w.r.t the first component
             Mass_pre(2,1,k,l) = Mass_pre(2,1,k,l)+ Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
               /int_cof_f(2*k+i,2*l-1)
             ! second set equation w.r.t the second component
             Mass_pre(2,2,k,l) = Mass_pre(2,2,k,l)+ Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2&
               +XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l-1)
             ! second set equation w.r.t the third component
             Mass_pre(2,3,k,l) = Mass_pre(2,3,k,l)+ Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
               /int_cof_f(2*k+i,2*l-1)
             ! third set equation w.r.t the first component
             Mass_pre(3,1,k,l) = Mass_pre(3,1,k,l)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
               /int_cof_f(2*k+i,2*l-1)
             ! third set equation w.r.t the second component
             Mass_pre(3,2,k,l) = Mass_pre(3,2,k,l)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof&
               /int_cof_f(2*k+i,2*l-1)
             ! third set equation w.r.t the third component
             Mass_pre(3,3,k,l) = Mass_pre(3,3,k,l)+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(-i/2)*Jacobian_c(k,l,n3_c)&
               *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2&
               +mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k+i,2*l-1)
          end do
          !
          ! first set equation w.r.t the first component
          Mass_pre(1,1,k,l) = Mass_pre(1,1,k,l)+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1) &
           * Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2&
           +mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! first set equation w.r.t the second component
          Mass_pre(1,2,k,l) = Mass_pre(1,2,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1) &
           * Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)&
           *XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! first set equation w.r.t the third component
          Mass_pre(1,3,k,l) = Mass_pre(1,3,k,l) + Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1) &
           * Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)&
           *XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! second set equation w.r.t the first component
          Mass_pre(2,1,k,l) = Mass_pre(2,1,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! second set equation w.r.t the second component
          Mass_pre(2,2,k,l) = Mass_pre(2,2,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! second set equation w.r.t the third component
          Mass_pre(2,3,k,l) = Mass_pre(2,3,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! third set equation w.r.t the first component
          Mass_pre(3,1,k,l) = Mass_pre(3,1,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! third set equation w.r.t the second component
          Mass_pre(3,2,k,l) = Mass_pre(3,2,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          ! third set equation w.r.t the third component
          Mass_pre(3,3,k,l) = Mass_pre(3,3,k,l)+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))&
           *XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*int_cof/int_cof_f(2*k-1,2*l-1)
          !
          ! first set equation w.r.t the first component
          Mass_pre(1,1,k,l) = Mass_pre(1,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
             *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
             *(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c/int_cof_c(k,l)
          ! first set equation w.r.t the second component
          Mass_pre(1,2,k,l) = Mass_pre(1,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! first set equation w.r.t the third component
          Mass_pre(1,3,k,l) = Mass_pre(1,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! second set equation w.r.t the first component
          Mass_pre(2,1,k,l) = Mass_pre(2,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! second set equation w.r.t the second component
          Mass_pre(2,2,k,l) = Mass_pre(2,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
           *(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c/int_cof_c(k,l)
          ! second set equation w.r.t the third component
          Mass_pre(2,3,k,l) = Mass_pre(2,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! third set equation w.r.t the first component
          Mass_pre(3,1,k,l) = Mass_pre(3,1,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! third set equation w.r.t the second component
          Mass_pre(3,2,k,l) = Mass_pre(3,2,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c/int_cof_c(k,l)
          ! third set equation w.r.t the third component
          Mass_pre(3,3,k,l) = Mass_pre(3,3,k,l)-Sb(0)*Jacobian_c(k,l,n3_c)&
          *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)&
          *(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/h3_c/int_cof_c(k,l)
       end do
    end do
  !
  end subroutine Interface_preconditioner
  !
  subroutine Interface_system(Mass)

    ! System of linear equations
    real(dp), dimension (1:n1_c*n2_c*3,1:n1_c*n2_c*3) :: Mass
    real(dp) :: int_cof
    real(dp), dimension (:,:), allocatable :: int_cof_c, int_cof_f
    real(dp), dimension (:,:,:,:), allocatable :: Mass_p_11,Mass_p_12,Mass_p_13,Mass_p_21,Mass_p_22
    real(dp), dimension (:,:,:,:), allocatable :: Mass_p_23,Mass_p_31,Mass_p_32,Mass_p_33
    real(dp), dimension (:,:,:,:), allocatable :: Mass_r_11,Mass_r_12,Mass_r_13,Mass_r_21,Mass_r_22
    real(dp), dimension (:,:,:,:), allocatable :: Mass_r_23,Mass_r_31,Mass_r_32,Mass_r_33
    !
    allocate(int_cof_c(1:n1_c,1:n2_c))
    allocate(int_cof_f(1-nrg:n1_f+nrg,1-nrg:n2_f+nrg))
    allocate(Mass_p_11(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_12(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_13(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_21(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_22(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_23(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_31(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_32(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_p_33(-2:n1_f+3,-2:n2_f+3,-2:n1_c+3,-2:n2_c+3))
    allocate(Mass_r_11(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_12(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_13(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_21(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_22(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_23(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_31(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_32(1:n1_c,1:n2_c,1:n1_c,1:n2_c))
    allocate(Mass_r_33(1:n1_c,1:n2_c,1:n1_c,1:n2_c))

    Mass =0.d0
    Mass_p_11 = 0.d0
    Mass_p_12 = 0.d0
    Mass_p_13 = 0.d0
    Mass_p_21 = 0.d0
    Mass_p_22 = 0.d0
    Mass_p_23 = 0.d0
    Mass_p_31 = 0.d0
    Mass_p_32 = 0.d0
    Mass_p_33 = 0.d0
    Mass_r_11 = 0.d0
    Mass_r_12 = 0.d0
    Mass_r_13 = 0.d0
    Mass_r_21 = 0.d0
    Mass_r_22 = 0.d0
    Mass_r_23 = 0.d0
    Mass_r_31 = 0.d0
    Mass_r_32 = 0.d0
    Mass_r_33 = 0.d0

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

    do k = 0,n2_c+1
      do i = 0,n1_c+1
        Mass_p_11(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*((2.d0*mu_c(i,k,n3_c)+lambda_c(i,k,n3_c)) &
            *XI13_c(i,k,n3_c)*XI13_c(i,k,n3_c)+mu_c(i,k,n3_c)*(XI23_c(i,k,n3_c)*XI23_c(i,k,n3_c) &
            +XI33_c(i,k,n3_c)*XI33_c(i,k,n3_c)))/rho_c(i,k,n3_c)
        Mass_p_12(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI13_c(i,k,n3_c)*XI23_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_13(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI13_c(i,k,n3_c)*XI33_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_21(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI13_c(i,k,n3_c)*XI23_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_22(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*((2.d0*mu_c(i,k,n3_c)+lambda_c(i,k,n3_c)) &
            *XI23_c(i,k,n3_c)*XI23_c(i,k,n3_c)+XI13_c(i,k,n3_c)*XI13_c(i,k,n3_c)*mu_c(i,k,n3_c) &
            +XI33_c(i,k,n3_c)*XI33_c(i,k,n3_c)*mu_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_23(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI23_c(i,k,n3_c)*XI33_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_31(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI13_c(i,k,n3_c)*XI33_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_32(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*XI23_c(i,k,n3_c)*XI33_c(i,k,n3_c) &
            *(mu_c(i,k,n3_c)+lambda_c(i,k,n3_c))/rho_c(i,k,n3_c)
        Mass_p_33(2*i-1,2*k-1,i,k) = Jacobian_c(i,k,n3_c)*((2.d0*mu_c(i,k,n3_c)+lambda_c(i,k,n3_c)) &
            *XI33_c(i,k,n3_c)*XI33_c(i,k,n3_c)+mu_c(i,k,n3_c)*(XI13_c(i,k,n3_c)*XI13_c(i,k,n3_c) &
            +XI23_c(i,k,n3_c)*XI23_c(i,k,n3_c)))/rho_c(i,k,n3_c)
      end do
    end do
    !
    do k = -1,n2_c+1
      do i = 0,n1_c+1
        do j = -1,2
          Mass_p_11(2*i-1,2*k,i,k+j) = Mass_p_11(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *((2.d0*mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))*XI13_c(i,k+j,n3_c)*XI13_c(i,k+j,n3_c) &
            +mu_c(i,k+j,n3_c)*(XI23_c(i,k+j,n3_c)*XI23_c(i,k+j,n3_c)+XI33_c(i,k+j,n3_c) &
            *XI33_c(i,k+j,n3_c)))/rho_c(i,k+j,n3_c)
          Mass_p_12(2*i-1,2*k,i,k+j) = Mass_p_12(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI13_c(i,k+j,n3_c)*XI23_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_13(2*i-1,2*k,i,k+j) = Mass_p_13(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI13_c(i,k+j,n3_c)*XI33_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_21(2*i-1,2*k,i,k+j) = Mass_p_21(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI13_c(i,k+j,n3_c)*XI23_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_22(2*i-1,2*k,i,k+j) = Mass_p_22(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *((2.d0*mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))*XI23_c(i,k+j,n3_c)*XI23_c(i,k+j,n3_c) &
            +XI13_c(i,k+j,n3_c)*XI13_c(i,k+j,n3_c)*mu_c(i,k+j,n3_c)+XI33_c(i,k+j,n3_c) &
            *XI33_c(i,k+j,n3_c)*mu_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_23(2*i-1,2*k,i,k+j) = Mass_p_23(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI23_c(i,k+j,n3_c)*XI33_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_31(2*i-1,2*k,i,k+j) = Mass_p_31(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI13_c(i,k+j,n3_c)*XI33_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_32(2*i-1,2*k,i,k+j) = Mass_p_32(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *XI23_c(i,k+j,n3_c)*XI33_c(i,k+j,n3_c)*(mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))/rho_c(i,k+j,n3_c)
          Mass_p_33(2*i-1,2*k,i,k+j) = Mass_p_33(2*i-1,2*k,i,k+j)+P(j)*Jacobian_c(i,k+j,n3_c) &
            *((2.d0*mu_c(i,k+j,n3_c)+lambda_c(i,k+j,n3_c))*XI33_c(i,k+j,n3_c)*XI33_c(i,k+j,n3_c) &
            +mu_c(i,k+j,n3_c)*(XI13_c(i,k+j,n3_c)*XI13_c(i,k+j,n3_c)+XI23_c(i,k+j,n3_c) &
            *XI23_c(i,k+j,n3_c)))/rho_c(i,k+j,n3_c)
        end do
      end do
    end do
    !
    do k = 0,n2_c+1
      do i = -1,n1_c+1
        do j = -1,2
          Mass_p_11(2*i,2*k-1,i+j,k) = Mass_p_11(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *((2.d0*mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))*XI13_c(i+j,k,n3_c)*XI13_c(i+j,k,n3_c) &
            +mu_c(i+j,k,n3_c)*(XI23_c(i+j,k,n3_c)*XI23_c(i+j,k,n3_c)+XI33_c(i+j,k,n3_c) &
            *XI33_c(i+j,k,n3_c)))/rho_c(i+j,k,n3_c)
          Mass_p_12(2*i,2*k-1,i+j,k) = Mass_p_12(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI13_c(i+j,k,n3_c)*XI23_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_13(2*i,2*k-1,i+j,k) = Mass_p_13(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI13_c(i+j,k,n3_c)*XI33_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_21(2*i,2*k-1,i+j,k) = Mass_p_21(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI13_c(i+j,k,n3_c)*XI23_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_22(2*i,2*k-1,i+j,k) = Mass_p_22(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *((2.d0*mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))*XI23_c(i+j,k,n3_c)*XI23_c(i+j,k,n3_c) &
            +XI13_c(i+j,k,n3_c)*XI13_c(i+j,k,n3_c)*mu_c(i+j,k,n3_c)+XI33_c(i+j,k,n3_c) &
            *XI33_c(i+j,k,n3_c)*mu_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_23(2*i,2*k-1,i+j,k) = Mass_p_23(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI23_c(i+j,k,n3_c)*XI33_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_31(2*i,2*k-1,i+j,k) = Mass_p_31(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI13_c(i+j,k,n3_c)*XI33_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_32(2*i,2*k-1,i+j,k) = Mass_p_32(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *XI23_c(i+j,k,n3_c)*XI33_c(i+j,k,n3_c)*(mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))/rho_c(i+j,k,n3_c)
          Mass_p_33(2*i,2*k-1,i+j,k) = Mass_p_33(2*i,2*k-1,i+j,k)+P(j)*Jacobian_c(i+j,k,n3_c) &
            *((2.d0*mu_c(i+j,k,n3_c)+lambda_c(i+j,k,n3_c))*XI33_c(i+j,k,n3_c)*XI33_c(i+j,k,n3_c) &
            +mu_c(i+j,k,n3_c)*(XI13_c(i+j,k,n3_c)*XI13_c(i+j,k,n3_c)+XI23_c(i+j,k,n3_c) &
            *XI23_c(i+j,k,n3_c)))/rho_c(i+j,k,n3_c)
        end do
      end do
    end do
    !
    do k = -1,n2_c+1
      do i = -1,n1_c+1
        do l = -1,2
          do j = -1,2
            Mass_p_11(2*i,2*k,i+j,k+l)=Mass_p_11(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *((2.d0*mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c))*XI13_c(i+j,k+l,n3_c) &
              *XI13_c(i+j,k+l,n3_c)+mu_c(i+j,k+l,n3_c)*(XI23_c(i+j,k+l,n3_c)*XI23_c(i+j,k+l,n3_c) &
              +XI33_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)))/rho_c(i+j,k+l,n3_c))
            Mass_p_12(2*i,2*k,i+j,k+l)=Mass_p_12(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI13_c(i+j,k+l,n3_c)*XI23_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_13(2*i,2*k,i+j,k+l)=Mass_p_13(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI13_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_21(2*i,2*k,i+j,k+l)=Mass_p_21(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI13_c(i+j,k+l,n3_c)*XI23_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_22(2*i,2*k,i+j,k+l)=Mass_p_22(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *((2.d0*mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c))*XI23_c(i+j,k+l,n3_c) &
              *XI23_c(i+j,k+l,n3_c)+XI13_c(i+j,k+l,n3_c)*XI13_c(i+j,k+l,n3_c)*mu_c(i+j,k+l,n3_c) &
              +XI33_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)*mu_c(i+j,k+l,n3_c))/rho_c(i+j,k+l,n3_c))
            Mass_p_23(2*i,2*k,i+j,k+l)=Mass_p_23(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI23_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_31(2*i,2*k,i+j,k+l)=Mass_p_31(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI13_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_32(2*i,2*k,i+j,k+l)=Mass_p_32(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *XI23_c(i+j,k+l,n3_c)*XI33_c(i+j,k+l,n3_c)*(mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c)) &
              /rho_c(i+j,k+l,n3_c))
            Mass_p_33(2*i,2*k,i+j,k+l)=Mass_p_33(2*i,2*k,i+j,k+l)+P(l)*(P(j)*Jacobian_c(i+j,k+l,n3_c) &
              *((2.d0*mu_c(i+j,k+l,n3_c)+lambda_c(i+j,k+l,n3_c))*XI33_c(i+j,k+l,n3_c) &
              *XI33_c(i+j,k+l,n3_c)+mu_c(i+j,k+l,n3_c)*(XI13_c(i+j,k+l,n3_c)*XI13_c(i+j,k+l,n3_c) &
              +XI23_c(i+j,k+l,n3_c)*XI23_c(i+j,k+l,n3_c)))/rho_c(i+j,k+l,n3_c))
          end do
        end do
      end do
    end do

    Mass_p_11(:,:,-2:0,:) = 0.d0
    Mass_p_11(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_11(:,:,:,-2:0) = 0.d0
    Mass_p_11(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_12(:,:,-2:0,:) = 0.d0
    Mass_p_12(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_12(:,:,:,-2:0) = 0.d0
    Mass_p_12(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_13(:,:,-2:0,:) = 0.d0
    Mass_p_13(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_13(:,:,:,-2:0) = 0.d0
    Mass_p_13(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_21(:,:,-2:0,:) = 0.d0
    Mass_p_21(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_21(:,:,:,-2:0) = 0.d0
    Mass_p_21(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_22(:,:,-2:0,:) = 0.d0
    Mass_p_22(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_22(:,:,:,-2:0) = 0.d0
    Mass_p_22(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_23(:,:,-2:0,:) = 0.d0
    Mass_p_23(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_23(:,:,:,-2:0) = 0.d0
    Mass_p_23(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_31(:,:,-2:0,:) = 0.d0
    Mass_p_31(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_31(:,:,:,-2:0) = 0.d0
    Mass_p_31(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_32(:,:,-2:0,:) = 0.d0
    Mass_p_32(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_32(:,:,:,-2:0) = 0.d0
    Mass_p_32(:,:,:,n2_c+1:n2_c+3) = 0.d0
    !
    Mass_p_33(:,:,-2:0,:) = 0.d0
    Mass_p_33(:,:,n1_c+1:n1_c+3,:) = 0.d0
    Mass_p_33(:,:,:,-2:0) = 0.d0
    Mass_p_33(:,:,:,n2_c+1:n2_c+3) = 0.d0

    do k = 1,n2_c
      do i = 1,n1_c
        do l = -4,2
          do j = -4,2
            Mass_r_11(i,k,1:n1_c,1:n2_c) = Mass_r_11(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_11(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_12(i,k,1:n1_c,1:n2_c) = Mass_r_12(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_12(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_13(i,k,1:n1_c,1:n2_c) = Mass_r_13(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_13(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_21(i,k,1:n1_c,1:n2_c) = Mass_r_21(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_21(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_22(i,k,1:n1_c,1:n2_c) = Mass_r_22(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_22(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_23(i,k,1:n1_c,1:n2_c) = Mass_r_23(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_23(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_31(i,k,1:n1_c,1:n2_c) = Mass_r_31(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_31(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_32(i,k,1:n1_c,1:n2_c) = Mass_r_32(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_32(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
            !
            Mass_r_33(i,k,1:n1_c,1:n2_c) = Mass_r_33(i,k,1:n1_c,1:n2_c) &
                + Rop(l)*(Rop(j)*rho_f(2*i+j,2*k+l,1)*Mass_p_33(2*i+j,2*k+l,1:n1_c,1:n2_c)&
                /int_cof_f(2*i+j,2*k+l))
          end do
        end do
      end do
    end do
    do l = 1,n2_c
      do k = 1,n1_c
        do j = 1,n2_c
          do i = 1,n1_c
            Mass((l-1)*3*n1_c+(k-1)*3+1,(j-1)*3*n1_c+(i-1)*3+1) = Mass_r_11(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+1,(j-1)*3*n1_c+(i-1)*3+2) = Mass_r_12(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+1,(j-1)*3*n1_c+(i-1)*3+3) = Mass_r_13(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+2,(j-1)*3*n1_c+(i-1)*3+1) = Mass_r_21(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+2,(j-1)*3*n1_c+(i-1)*3+2) = Mass_r_22(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+2,(j-1)*3*n1_c+(i-1)*3+3) = Mass_r_23(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+3,(j-1)*3*n1_c+(i-1)*3+1) = Mass_r_31(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+3,(j-1)*3*n1_c+(i-1)*3+2) = Mass_r_32(k,l,i,j)*int_cof
            Mass((l-1)*3*n1_c+(k-1)*3+3,(j-1)*3*n1_c+(i-1)*3+3) = Mass_r_33(k,l,i,j)*int_cof
          end do
        end do
      end do
    end do
    do j = 1, n2_c
      do i = 1, n1_c
        Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+1) = &
              Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+1) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c)) &
              *XI13_c(i,j,n3_c)*XI13_c(i,j,n3_c)+mu_c(i,j,n3_c)*(XI23_c(i,j,n3_c) &
              *XI23_c(i,j,n3_c)+XI33_c(i,j,n3_c)*XI33_c(i,j,n3_c)))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+2) = &
              Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+2) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+3) = &
              Mass((j-1)*3*n1_c+(i-1)*3+1,(j-1)*3*n1_c+(i-1)*3+3) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+1) = &
              Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+1) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI13_c(i,j,n3_c)*XI23_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+2) = &
              Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+2) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c)) &
              *XI23_c(i,j,n3_c)*XI23_c(i,j,n3_c)+XI13_c(i,j,n3_c)*XI13_c(i,j,n3_c) &
              *mu_c(i,j,n3_c)+XI33_c(i,j,n3_c)*XI33_c(i,j,n3_c)*mu_c(i,j,n3_c)) &
              /h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+3) = &
              Mass((j-1)*3*n1_c+(i-1)*3+2,(j-1)*3*n1_c+(i-1)*3+3) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+1) = &
              Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+1) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI13_c(i,j,n3_c)*XI33_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+2) = &
              Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+2) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*XI23_c(i,j,n3_c)*XI33_c(i,j,n3_c) &
              *(mu_c(i,j,n3_c)+lambda_c(i,j,n3_c))/h3_c/int_cof_c(i,j)
        !
        Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+3) = &
              Mass((j-1)*3*n1_c+(i-1)*3+3,(j-1)*3*n1_c+(i-1)*3+3) &
              -Sb(0)*Jacobian_c(i,j,n3_c)*((2.d0*mu_c(i,j,n3_c)+lambda_c(i,j,n3_c)) &
              *XI33_c(i,j,n3_c)*XI33_c(i,j,n3_c)+mu_c(i,j,n3_c)*(XI13_c(i,j,n3_c) &
              *XI13_c(i,j,n3_c)+XI23_c(i,j,n3_c)*XI23_c(i,j,n3_c)))/h3_c/int_cof_c(i,j)
      end do
    end do
    !
  end subroutine Interface_system
  !
  subroutine Interface_LHS_cg(u_temp_cg)
    real(dp) :: int_cof
    real(dp), dimension (1:3*n1_c*n2_c) :: u_temp_cg
    real(dp), dimension (:,:,:), allocatable :: u_temp
    real(dp), dimension (1:n1_c,1:n2_c) :: int_cof_c
    real(dp), dimension (1-nrg:n1_f+nrg,1-nrg:n2_f+nrg) :: int_cof_f

    ! allocate memory for temporary arrays
    allocate(u_temp(1-nrg:n1_c+nrg,1-nrg:n2_c+nrg,1:dim))
    !
    do j = 1,n2_c
       do i = 1,n1_c
          int_cof_c(i,j) = Jacobian_c(i,j,n3_c)*sqrt(XI13_c(i,j,n3_c)**2+XI23_c(i,j,n3_c)**2+XI33_c(i,j,n3_c)**2)
       end do
    end do
    !
    do j = 1-nrg,n2_f+nrg
       do i = 1-nrg,n1_f+nrg
          int_cof_f(i,j) = Jacobian_f(i,j,1)*sqrt(XI13_f(i,j,1)**2+XI23_f(i,j,1)**2+XI33_f(i,j,1)**2)
       end do
    end do
    !
    u_temp = 0.d0
    do j = 1,n2_c
       do i = 1,n1_c
          u_temp(i,j,1) = u_temp_cg((j-1)*3*n1_c+(i-1)*3+1)
          u_temp(i,j,2) = u_temp_cg((j-1)*3*n1_c+(i-1)*3+2)
          u_temp(i,j,3) = u_temp_cg((j-1)*3*n1_c+(i-1)*3+3)
       end do
    end do
    !
    int_cof = 17.d0/48.d0*h3_f*ghcof(1)/h3_c**2
    !
    LHS = 0.d0
    do l = 1,n2_c
       do k = 1,n1_c
          !
          do j = -4,2,2
             do i = -4,2,2
                do j1 = -1,2
                   do i1 = -1,2
                      ! first set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+(Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof &
                       +Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI13_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof+Rop(j)&
                       *Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI13_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof)&
                       /int_cof_f(2*k+i,2*l+j)
                      ! second set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+(Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)*XI23_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof &
                       + Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)*Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))*XI23_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof)&
                       /int_cof_f(2*k+i,2*l+j)
                      ! third set equation
                      LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+(Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI13_c(k+i/2+i1,l+j/2+j1,n3_c)*XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,1)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*(lambda_c(k+i/2+i1,l+j/2+j1,n3_c)+mu_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI23_c(k+i/2+i1,l+j/2+j1,n3_c)*XI33_c(k+i/2+i1,l+j/2+j1,n3_c)/rho_c(k+i/2+i1,l+j/2+j1,n3_c)&
                       *u_temp(k+i/2+i1,l+j/2+j1,2)*int_cof+Rop(j)*Rop(i)*rho_f(2*k+i,2*l+j,1)*P(j1)*P(i1)&
                       *Jacobian_c(k+i/2+i1,l+j/2+j1,n3_c)*((2.d0*mu_c(k+i/2+i1,l+j/2+j1,n3_c)+lambda_c(k+i/2+i1,l+j/2+j1,n3_c))&
                       *XI33_c(k+i/2+i1,l+j/2+j1,n3_c)**2+mu_c(k+i/2+i1,l+j/2+j1,n3_c)*(XI13_c(k+i/2+i1,l+j/2+j1,n3_c)**2&
                       +XI23_c(k+i/2+i1,l+j/2+j1,n3_c)**2))/rho_c(k+i/2+i1,l+j/2+j1,n3_c)*u_temp(k+i/2+i1,l+j/2+j1,3)*int_cof)&
                       /int_cof_f(2*k+i,2*l+j)
                   end do
                end do
             end do
          end do
          !
          do j = -4,2,2
             do j1 = -1,2
                ! first set equation
                LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+(Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI13_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI23_c(k,l+j/2+j1,n3_c)**2+XI33_c(k,l+j/2+j1,n3_c)**2))&
                 /rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI23_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,3)*int_cof)/int_cof_f(2*k-1,2*l+j)
                ! second set equation
                LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+(Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI23_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1) &
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI23_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI13_c(k,l+j/2+j1,n3_c)**2+XI33_c(k,l+j/2+j1,n3_c)**2))&
                 /rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI23_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)* u_temp(k,l+j/2+j1,3)*int_cof)/int_cof_f(2*k-1,2*l+j)
                ! third set equation
                LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+(Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)*P(j1)&
                 *Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI13_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,1)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 * P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*(lambda_c(k,l+j/2+j1,n3_c)+mu_c(k,l+j/2+j1,n3_c))*XI23_c(k,l+j/2+j1,n3_c)&
                 *XI33_c(k,l+j/2+j1,n3_c)/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,2)*int_cof+Rop(j)*Rop(-1)*rho_f(2*k-1,2*l+j,1)&
                 *P(j1)*Jacobian_c(k,l+j/2+j1,n3_c)*((2.d0*mu_c(k,l+j/2+j1,n3_c)+lambda_c(k,l+j/2+j1,n3_c))&
                 *XI33_c(k,l+j/2+j1,n3_c)**2+mu_c(k,l+j/2+j1,n3_c)*(XI13_c(k,l+j/2+j1,n3_c)**2&
                 +XI23_c(k,l+j/2+j1,n3_c)**2))/rho_c(k,l+j/2+j1,n3_c)*u_temp(k,l+j/2+j1,3)*int_cof)/int_cof_f(2*k-1,2*l+j)
             end do
          end do
          !
          do i = -4,2,2
             do i1 = -1,2
                ! first set equation
                LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+(Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)**2&
                 +mu_c(k+i/2+i1,l,n3_c)*(XI23_c(k+i/2+i1,l,n3_c)**2+XI33_c(k+i/2+i1,l,n3_c)**2))/rho_c(k+i/2+i1,l,n3_c)&
                 *u_temp(k+i/2+i1,l,1)*int_cof+ Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)&
                 *(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)*XI23_c(k+i/2+i1,l,n3_c)&
                 /rho_c(k+i/2+i1,l,n3_c)* u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)* u_temp(k+i/2+i1,l,3)*int_cof)/int_cof_f(2*k+i,2*l-1)
                ! second set equation
                LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+(Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI23_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,1)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))&
                 *XI23_c(k+i/2+i1,l,n3_c)**2+mu_c(k+i/2+i1,l,n3_c)*(XI13_c(k+i/2+i1,l,n3_c)**2+XI33_c(k+i/2+i1,l,n3_c)**2))&
                 /rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI23_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,3)*int_cof)/int_cof_f(2*k+i,2*l-1)
                ! third set equation
                LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+(Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)*P(i1)&
                 *Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI13_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,1)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*(lambda_c(k+i/2+i1,l,n3_c)+mu_c(k+i/2+i1,l,n3_c))*XI23_c(k+i/2+i1,l,n3_c)&
                 *XI33_c(k+i/2+i1,l,n3_c)/rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,2)*int_cof+Rop(-1)*Rop(i)*rho_f(2*k+i,2*l-1,1)&
                 *P(i1)*Jacobian_c(k+i/2+i1,l,n3_c)*((2.d0*mu_c(k+i/2+i1,l,n3_c)+lambda_c(k+i/2+i1,l,n3_c))&
                 *XI33_c(k+i/2+i1,l,n3_c)**2+mu_c(k+i/2+i1,l,n3_c)*(XI13_c(k+i/2+i1,l,n3_c)**2+XI23_c(k+i/2+i1,l,n3_c)**2))&
                 /rho_c(k+i/2+i1,l,n3_c)*u_temp(k+i/2+i1,l,3)*int_cof)/int_cof_f(2*k+i,2*l-1)
             end do
          end do
          !
          ! first set equation
          LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+(Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)&
           +mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)&
           *XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof)/int_cof_f(2*k-1,2*l-1)
          ! second set equation
          LHS((l-1)*3*n1_c+(k-1)*3+2) = LHS((l-1)*3*n1_c+(k-1)*3+2)+(Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof &
           +Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2&
           +mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)&
           *rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof)/int_cof_f(2*k-1,2*l-1)
          ! third set equation
          LHS((l-1)*3*n1_c+(k-1)*3+3) = LHS((l-1)*3*n1_c+(k-1)*3+3)+(Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,1)*int_cof &
           +Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)&
           *XI33_c(k,l,n3_c)/rho_c(k,l,n3_c)*u_temp(k,l,2)*int_cof+Rop(-1)*Rop(-1)*rho_f(2*k-1,2*l-1,1)*Jacobian_c(k,l,n3_c)&
           *((2.d0*mu_c(k,l,n3_c)+lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))&
           /rho_c(k,l,n3_c)*u_temp(k,l,3)*int_cof)/int_cof_f(2*k-1,2*l-1)
          !
          ! first set equation
          LHS((l-1)*3*n1_c+(k-1)*3+1) = LHS((l-1)*3*n1_c+(k-1)*3+1)+(-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI13_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI23_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,1) &
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c*u_temp(k,l,2) &
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,3))&
           /int_cof_c(k,l)
          ! second set equation
          LHS((l-1)*3*n1_c+(k-1)*3+2)=LHS((l-1)*3*n1_c+(k-1)*3+2)+(-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI23_c(k,l,n3_c)/h3_c*u_temp(k,l,1)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI23_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI33_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,2)&
           -Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))*XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,3))&
           /int_cof_c(k,l)
          ! third set equation
          LHS((l-1)*3*n1_c+(k-1)*3+3)=LHS((l-1)*3*n1_c+(k-1)*3+3)+(-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI13_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,1)-Sb(0)*Jacobian_c(k,l,n3_c)*(lambda_c(k,l,n3_c)+mu_c(k,l,n3_c))&
           *XI23_c(k,l,n3_c)*XI33_c(k,l,n3_c)/h3_c*u_temp(k,l,2)-Sb(0)*Jacobian_c(k,l,n3_c)*((2.d0*mu_c(k,l,n3_c)&
           +lambda_c(k,l,n3_c))*XI33_c(k,l,n3_c)**2+mu_c(k,l,n3_c)*(XI13_c(k,l,n3_c)**2+XI23_c(k,l,n3_c)**2))/h3_c*u_temp(k,l,3))&
           /int_cof_c(k,l)
       end do
    end do
  !
  end subroutine Interface_LHS_cg
  !
end program condition_number
