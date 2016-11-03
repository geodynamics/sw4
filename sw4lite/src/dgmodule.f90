module dgmodule

contains

  SUBROUTINE get_comm_sides(x_in_b,x_in_e,y_in_b,y_in_e,&
       v_in_all_faces,w_in_all_faces,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,w_in_all_faces
    REAL(dp), DIMENSION(0:nint,0:nint,3,jfirst:jlast,kfirst:klast,2) :: x_in_b,x_in_e
    REAL(dp), DIMENSION(0:nint,0:nint,3,ifirst:ilast,kfirst:klast,2) :: y_in_b,y_in_e
    integer :: i1,i2,i3,i4,i5

    !$OMP PARALLEL DO PRIVATE(i5,i4,i3,i2,i1)
    do i5 = kfirst,klast
     do i4 = jfirst,jlast
      do i3 = 1,3
       do i2 = 0,nint
        do i1 = 0,nint
         x_in_b(i1,i2,i3,i4,i5,1) = v_in_all_faces(i1,i2,1,i3,ifirst,i4,i5)
         x_in_b(i1,i2,i3,i4,i5,2) = w_in_all_faces(i1,i2,1,i3,ifirst,i4,i5)
         x_in_e(i1,i2,i3,i4,i5,1) = v_in_all_faces(i1,i2,2,i3,ilast,i4,i5)
         x_in_e(i1,i2,i3,i4,i5,2) = w_in_all_faces(i1,i2,2,i3,ilast,i4,i5)
        end do
       end do
      end do
     end do
     do i4 = ifirst,ilast
      do i3  = 1,3
       do i2 = 0,nint
        do i1 = 0,nint
         y_in_b(i1,i2,i3,i4,i5,1) = v_in_all_faces(i1,i2,3,i3,i4,jfirst,i5)
         y_in_b(i1,i2,i3,i4,i5,2) = w_in_all_faces(i1,i2,3,i3,i4,jfirst,i5)
         y_in_e(i1,i2,i3,i4,i5,1) = v_in_all_faces(i1,i2,4,i3,i4,jlast,i5)
         y_in_e(i1,i2,i3,i4,i5,2) = w_in_all_faces(i1,i2,4,i3,i4,jlast,i5)
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE get_comm_sides

  SUBROUTINE put_comm_sides(x_out_b,x_out_e,y_out_b,y_out_e,v_out_all_faces,w_out_all_faces,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_out_all_faces,w_out_all_faces
    REAL(dp), DIMENSION(0:nint,0:nint,3,jfirst:jlast,kfirst:klast,2) :: x_out_b,x_out_e
    REAL(dp), DIMENSION(0:nint,0:nint,3,ifirst:ilast,kfirst:klast,2) :: y_out_b,y_out_e
    integer :: i1,i2,i3,i4,i5

    !$OMP PARALLEL DO PRIVATE(i5,i4,i3,i2,i1)
    do i5 = kfirst,klast
     do i4 = jfirst,jlast
      do i3 = 1,3
       do i2 = 0,nint
        do i1 = 0,nint
         v_out_all_faces(i1,i2,1,i3,ifirst,i4,i5) = x_out_b(i1,i2,i3,i4,i5,1)
         w_out_all_faces(i1,i2,1,i3,ifirst,i4,i5) = x_out_b(i1,i2,i3,i4,i5,2)
         v_out_all_faces(i1,i2,2,i3,ilast,i4,i5) = x_out_e(i1,i2,i3,i4,i5,1)
         w_out_all_faces(i1,i2,2,i3,ilast,i4,i5) = x_out_e(i1,i2,i3,i4,i5,2)
        end do
       end do
      end do
     end do
     do i4 = ifirst,ilast
      do i3  = 1,3
       do i2 = 0,nint
        do i1 = 0,nint
         v_out_all_faces(i1,i2,3,i3,i4,jfirst,i5) = y_out_b(i1,i2,i3,i4,i5,1)
         w_out_all_faces(i1,i2,3,i3,i4,jfirst,i5) = y_out_b(i1,i2,i3,i4,i5,2)
         v_out_all_faces(i1,i2,4,i3,i4,jlast,i5) = y_out_e(i1,i2,i3,i4,i5,1)
         w_out_all_faces(i1,i2,4,i3,i4,jlast,i5) = y_out_e(i1,i2,i3,i4,i5,2)
        end do
       end do
      end do
     end do
    end do


  END SUBROUTINE put_comm_sides

  SUBROUTINE compute_time_derivatives(utdg,vtdg,udg,vdg,force_u,force_v,&
       MU,MV,SU,SV,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,q_v,q_u) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_v+1)**3) :: MV,SU
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_u+1)**3) :: MU,SV
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg,force_u,utdg
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg,force_v,vtdg
    REAL(dp) :: df
    INTEGER :: nrows,ncols,j1,j2,j3

    nrows = 3*(q_u+1)**3
    ncols = 3*(q_u+1)**3
    !$OMP PARALLEL DO PRIVATE(j3,j2,j1)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       ! Compute u_t
       ! Use the inverse of MU
       utdg(0:q_u,0:q_u,0:q_u,1:3,j1,j2,j3) = vdg(0:q_v,0:q_v,0:q_v,1:3,j1,j2,j3)
       call DGEMV('N',nrows,ncols,1.d0,MU,nrows,&
            force_u(0:q_u,0:q_u,0:q_u,1:3,j1,j2,j3),1,&
            1.d0,utdg(0:q_u,0:q_u,0:q_u,1:3,j1,j2,j3),1)
       ! Compute v_t
       ! first force_v <- SU*u + force_v
       call DGEMV('N',nrows,ncols,1.d0,SU,nrows,&
            udg(0:q_u,0:q_u,0:q_u,1:3,j1,j2,j3),1,&
            1.d0,force_v(0:q_v,0:q_v,0:q_v,1:3,j1,j2,j3),1)
       ! Use the inverse of MV
       vtdg(0:q_v,0:q_v,0:q_v,1:3,j1,j2,j3) = 0.0_dp
       call DGEMV('N',nrows,ncols,1.d0,MV,nrows,&
            force_v(0:q_u,0:q_u,0:q_u,1:3,j1,j2,j3),1,&
            1.d0,vtdg(0:q_v,0:q_v,0:q_v,1:3,j1,j2,j3),1)
      end do
     end do
    end do
  END SUBROUTINE compute_time_derivatives

  SUBROUTINE compute_numerical_fluxes(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,&
       v_star_all_faces,w_star_all_faces,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint,lambda_lame,mu_lame,rho) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp) :: lambda_lame,mu_lame,rho
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,v_out_all_faces,w_in_all_faces,w_out_all_faces, &
         v_star_all_faces,w_star_all_faces
    REAL(dp) :: alpha_flip, beta_up, tau_up

    alpha_flip = 0.5_dp
    beta_up = 0.5_dp
    tau_up  = 0.1_dp*rho/(2.0_dp*mu_lame+lambda_lame)

    !$OMP WORKSHARE
    ! Generalized averaging for the flux (allows for flip-flop).
    v_star_all_faces(:,:,1:5:2,:,:,:,:) = alpha_flip*v_in_all_faces(:,:,1:5:2,:,:,:,:)&
         +(1.0_dp-alpha_flip)*v_out_all_faces(:,:,1:5:2,:,:,:,:)
    v_star_all_faces(:,:,2:6:2,:,:,:,:) = (1.0_dp-alpha_flip)*v_in_all_faces(:,:,2:6:2,:,:,:,:)&
         +alpha_flip*v_out_all_faces(:,:,2:6:2,:,:,:,:)
    w_star_all_faces(:,:,1:5:2,:,:,:,:) = (1.0_dp-alpha_flip)*w_in_all_faces(:,:,1:5:2,:,:,:,:)&
         +alpha_flip*w_out_all_faces(:,:,1:5:2,:,:,:,:)
    w_star_all_faces(:,:,2:6:2,:,:,:,:) = alpha_flip*w_in_all_faces(:,:,2:6:2,:,:,:,:)&
         +(1.0_dp-alpha_flip)*w_out_all_faces(:,:,2:6:2,:,:,:,:)
    ! Upwind flux at 1,3,5 faces
    w_star_all_faces(:,:,1:5:2,:,:,:,:) = w_star_all_faces(:,:,1:5:2,:,:,:,:) &
         +beta_up*(v_in_all_faces(:,:,1:5:2,:,:,:,:)-v_out_all_faces(:,:,1:5:2,:,:,:,:))
    v_star_all_faces(:,:,1:5:2,:,:,:,:) = v_star_all_faces(:,:,1:5:2,:,:,:,:) &
         +tau_up*(w_in_all_faces(:,:,1:5:2,:,:,:,:)-w_out_all_faces(:,:,1:5:2,:,:,:,:))
    ! Upwind flux at 2,4,6 faces
    w_star_all_faces(:,:,2:6:2,:,:,:,:) = w_star_all_faces(:,:,2:6:2,:,:,:,:) &
         -beta_up*(v_in_all_faces(:,:,2:6:2,:,:,:,:)-v_out_all_faces(:,:,2:6:2,:,:,:,:))
    v_star_all_faces(:,:,2:6:2,:,:,:,:) = v_star_all_faces(:,:,2:6:2,:,:,:,:) &
         -tau_up*(w_in_all_faces(:,:,2:6:2,:,:,:,:)-w_out_all_faces(:,:,2:6:2,:,:,:,:))
    !$OMP END WORKSHARE

  END SUBROUTINE compute_numerical_fluxes


  SUBROUTINE pass_outside_fluxes(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,v_out_all_faces,w_in_all_faces,w_out_all_faces

    !$OMP WORKSHARE
    v_out_all_faces(:,:,1,:,ifirst+1:ilast,:,:) = v_in_all_faces(:,:,2,:,ifirst:ilast-1,:,:)
    v_out_all_faces(:,:,2,:,ifirst:ilast-1,:,:) = v_in_all_faces(:,:,1,:,ifirst+1:ilast,:,:)

    v_out_all_faces(:,:,3,:,:,jfirst+1:jlast,:) = v_in_all_faces(:,:,4,:,:,jfirst:jlast-1,:)
    v_out_all_faces(:,:,4,:,:,jfirst:jlast-1,:) = v_in_all_faces(:,:,3,:,:,jfirst+1:jlast,:)

    v_out_all_faces(:,:,5,:,:,:,kfirst+1:klast) = v_in_all_faces(:,:,6,:,:,:,kfirst:klast-1)
    v_out_all_faces(:,:,6,:,:,:,kfirst:klast-1) = v_in_all_faces(:,:,5,:,:,:,kfirst+1:klast)
    !$OMP END WORKSHARE

    !$OMP WORKSHARE
    w_out_all_faces(:,:,1,:,ifirst+1:ilast,:,:) = w_in_all_faces(:,:,2,:,ifirst:ilast-1,:,:)
    w_out_all_faces(:,:,2,:,ifirst:ilast-1,:,:) = w_in_all_faces(:,:,1,:,ifirst+1:ilast,:,:)

    w_out_all_faces(:,:,3,:,:,jfirst+1:jlast,:) = w_in_all_faces(:,:,4,:,:,jfirst:jlast-1,:)
    w_out_all_faces(:,:,4,:,:,jfirst:jlast-1,:) = w_in_all_faces(:,:,3,:,:,jfirst+1:jlast,:)

    w_out_all_faces(:,:,5,:,:,:,kfirst+1:klast) = w_in_all_faces(:,:,6,:,:,:,kfirst:klast-1)
    w_out_all_faces(:,:,6,:,:,:,kfirst:klast-1) = w_in_all_faces(:,:,5,:,:,:,kfirst+1:klast)
    !$OMP END WORKSHARE

  END SUBROUTINE pass_outside_fluxes

  SUBROUTINE set_boundary_conditions(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint,sbx_b,sbx_e,sby_b,sby_e,bc_type) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast,sbx_b,sbx_e,sby_b,sby_e,bc_type
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,v_out_all_faces, &
         w_in_all_faces,w_out_all_faces
    if(bc_type .eq. 1) then
     if (sbx_b.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,1,1,ifirst,:,:)     =  v_in_all_faces(:,:,1,1,ifirst,:,:)
      v_out_all_faces(:,:,1,2,ifirst,:,:)     = -v_in_all_faces(:,:,1,2,ifirst,:,:)
      v_out_all_faces(:,:,1,3,ifirst,:,:)     = -v_in_all_faces(:,:,1,3,ifirst,:,:)
      w_out_all_faces(:,:,1,1,ifirst,:,:)     = -w_in_all_faces(:,:,1,1,ifirst,:,:)
      w_out_all_faces(:,:,1,2,ifirst,:,:)     =  w_in_all_faces(:,:,1,2,ifirst,:,:)
      w_out_all_faces(:,:,1,3,ifirst,:,:)     =  w_in_all_faces(:,:,1,3,ifirst,:,:)
      !$OMP END WORKSHARE
     end if
     if (sbx_e.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,2,1,ilast,:,:)     =  v_in_all_faces(:,:,2,1,ilast,:,:)
      v_out_all_faces(:,:,2,2,ilast,:,:)     = -v_in_all_faces(:,:,2,2,ilast,:,:)
      v_out_all_faces(:,:,2,3,ilast,:,:)     = -v_in_all_faces(:,:,2,3,ilast,:,:)
      w_out_all_faces(:,:,2,1,ilast,:,:)     = -w_in_all_faces(:,:,2,1,ilast,:,:)
      w_out_all_faces(:,:,2,2,ilast,:,:)     =  w_in_all_faces(:,:,2,2,ilast,:,:)
      w_out_all_faces(:,:,2,3,ilast,:,:)     =  w_in_all_faces(:,:,2,3,ilast,:,:)
      !$OMP END WORKSHARE
     end if

     if (sby_b.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,3,1,:,jfirst,:)     = -v_in_all_faces(:,:,3,1,:,jfirst,:)
      v_out_all_faces(:,:,3,2,:,jfirst,:)     =  v_in_all_faces(:,:,3,2,:,jfirst,:)
      v_out_all_faces(:,:,3,3,:,jfirst,:)     = -v_in_all_faces(:,:,3,3,:,jfirst,:)
      w_out_all_faces(:,:,3,1,:,jfirst,:)     =  w_in_all_faces(:,:,3,1,:,jfirst,:)
      w_out_all_faces(:,:,3,2,:,jfirst,:)     = -w_in_all_faces(:,:,3,2,:,jfirst,:)
      w_out_all_faces(:,:,3,3,:,jfirst,:)     =  w_in_all_faces(:,:,3,3,:,jfirst,:)
      !$OMP END WORKSHARE
     end if
     if (sby_e.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,4,1,:,jlast,:)     = -v_in_all_faces(:,:,4,1,:,jlast,:)
      v_out_all_faces(:,:,4,2,:,jlast,:)     =  v_in_all_faces(:,:,4,2,:,jlast,:)
      v_out_all_faces(:,:,4,3,:,jlast,:)     = -v_in_all_faces(:,:,4,3,:,jlast,:)
      w_out_all_faces(:,:,4,1,:,jlast,:)     =  w_in_all_faces(:,:,4,1,:,jlast,:)
      w_out_all_faces(:,:,4,2,:,jlast,:)     = -w_in_all_faces(:,:,4,2,:,jlast,:)
      w_out_all_faces(:,:,4,3,:,jlast,:)     =  w_in_all_faces(:,:,4,3,:,jlast,:)
      !$OMP END WORKSHARE
     end if
     !$OMP WORKSHARE
     v_out_all_faces(:,:,5,1,:,:,kfirst)     = -v_in_all_faces(:,:,5,1,:,:,kfirst)
     v_out_all_faces(:,:,5,2,:,:,kfirst)     = -v_in_all_faces(:,:,5,2,:,:,kfirst)
     v_out_all_faces(:,:,5,3,:,:,kfirst)     =  v_in_all_faces(:,:,5,3,:,:,kfirst)
     v_out_all_faces(:,:,6,1,:,:,klast)     = -v_in_all_faces(:,:,6,1,:,:,klast)
     v_out_all_faces(:,:,6,2,:,:,klast)     = -v_in_all_faces(:,:,6,2,:,:,klast)
     v_out_all_faces(:,:,6,3,:,:,klast)     =  v_in_all_faces(:,:,6,3,:,:,klast)

     w_out_all_faces(:,:,5,1,:,:,kfirst)     =  w_in_all_faces(:,:,5,1,:,:,kfirst)
     w_out_all_faces(:,:,5,2,:,:,kfirst)     =  w_in_all_faces(:,:,5,2,:,:,kfirst)
     w_out_all_faces(:,:,5,3,:,:,kfirst)     = -w_in_all_faces(:,:,5,3,:,:,kfirst)
     w_out_all_faces(:,:,6,1,:,:,klast)     =  w_in_all_faces(:,:,6,1,:,:,klast)
     w_out_all_faces(:,:,6,2,:,:,klast)     =  w_in_all_faces(:,:,6,2,:,:,klast)
     w_out_all_faces(:,:,6,3,:,:,klast)     = -w_in_all_faces(:,:,6,3,:,:,klast)
     !$OMP END WORKSHARE
    else !! One free surface
     if (sbx_b.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,1,1,ifirst,:,:)     = -v_in_all_faces(:,:,1,1,ifirst,:,:)
      v_out_all_faces(:,:,1,2,ifirst,:,:)     = -v_in_all_faces(:,:,1,2,ifirst,:,:)
      v_out_all_faces(:,:,1,3,ifirst,:,:)     = -v_in_all_faces(:,:,1,3,ifirst,:,:)
      w_out_all_faces(:,:,1,1,ifirst,:,:)     =  w_in_all_faces(:,:,1,1,ifirst,:,:)
      w_out_all_faces(:,:,1,2,ifirst,:,:)     =  w_in_all_faces(:,:,1,2,ifirst,:,:)
      w_out_all_faces(:,:,1,3,ifirst,:,:)     =  w_in_all_faces(:,:,1,3,ifirst,:,:)
      !$OMP END WORKSHARE
     end if
     if (sbx_e.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,2,1,ilast,:,:)     = -v_in_all_faces(:,:,2,1,ilast,:,:)
      v_out_all_faces(:,:,2,2,ilast,:,:)     = -v_in_all_faces(:,:,2,2,ilast,:,:)
      v_out_all_faces(:,:,2,3,ilast,:,:)     = -v_in_all_faces(:,:,2,3,ilast,:,:)
      w_out_all_faces(:,:,2,1,ilast,:,:)     =  w_in_all_faces(:,:,2,1,ilast,:,:)
      w_out_all_faces(:,:,2,2,ilast,:,:)     =  w_in_all_faces(:,:,2,2,ilast,:,:)
      w_out_all_faces(:,:,2,3,ilast,:,:)     =  w_in_all_faces(:,:,2,3,ilast,:,:)
      !$OMP END WORKSHARE
     end if

     if (sby_b.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,3,1,:,jfirst,:)     = -v_in_all_faces(:,:,3,1,:,jfirst,:)
      v_out_all_faces(:,:,3,2,:,jfirst,:)     = -v_in_all_faces(:,:,3,2,:,jfirst,:)
      v_out_all_faces(:,:,3,3,:,jfirst,:)     = -v_in_all_faces(:,:,3,3,:,jfirst,:)
      w_out_all_faces(:,:,3,1,:,jfirst,:)     =  w_in_all_faces(:,:,3,1,:,jfirst,:)
      w_out_all_faces(:,:,3,2,:,jfirst,:)     =  w_in_all_faces(:,:,3,2,:,jfirst,:)
      w_out_all_faces(:,:,3,3,:,jfirst,:)     =  w_in_all_faces(:,:,3,3,:,jfirst,:)
      !$OMP END WORKSHARE
     end if
     if (sby_e.eq.1) then
      !$OMP WORKSHARE
      v_out_all_faces(:,:,4,1,:,jlast,:)     = -v_in_all_faces(:,:,4,1,:,jlast,:)
      v_out_all_faces(:,:,4,2,:,jlast,:)     = -v_in_all_faces(:,:,4,2,:,jlast,:)
      v_out_all_faces(:,:,4,3,:,jlast,:)     = -v_in_all_faces(:,:,4,3,:,jlast,:)
      w_out_all_faces(:,:,4,1,:,jlast,:)     =  w_in_all_faces(:,:,4,1,:,jlast,:)
      w_out_all_faces(:,:,4,2,:,jlast,:)     =  w_in_all_faces(:,:,4,2,:,jlast,:)
      w_out_all_faces(:,:,4,3,:,jlast,:)     =  w_in_all_faces(:,:,4,3,:,jlast,:)
      !$OMP END WORKSHARE
     end if
     !$OMP WORKSHARE
     v_out_all_faces(:,:,5,1,:,:,kfirst)     =  v_in_all_faces(:,:,5,1,:,:,kfirst)
     v_out_all_faces(:,:,5,2,:,:,kfirst)     =  v_in_all_faces(:,:,5,2,:,:,kfirst)
     v_out_all_faces(:,:,5,3,:,:,kfirst)     =  v_in_all_faces(:,:,5,3,:,:,kfirst)
     v_out_all_faces(:,:,6,1,:,:,klast)      = -v_in_all_faces(:,:,6,1,:,:,klast)
     v_out_all_faces(:,:,6,2,:,:,klast)      = -v_in_all_faces(:,:,6,2,:,:,klast)
     v_out_all_faces(:,:,6,3,:,:,klast)      = -v_in_all_faces(:,:,6,3,:,:,klast)

     w_out_all_faces(:,:,5,1,:,:,kfirst)     =  -w_in_all_faces(:,:,5,1,:,:,kfirst)
     w_out_all_faces(:,:,5,2,:,:,kfirst)     =  -w_in_all_faces(:,:,5,2,:,:,kfirst)
     w_out_all_faces(:,:,5,3,:,:,kfirst)     =  -w_in_all_faces(:,:,5,3,:,:,kfirst)
     w_out_all_faces(:,:,6,1,:,:,klast)      =  w_in_all_faces(:,:,6,1,:,:,klast)
     w_out_all_faces(:,:,6,2,:,:,klast)      =  w_in_all_faces(:,:,6,2,:,:,klast)
     w_out_all_faces(:,:,6,3,:,:,klast)      =  w_in_all_faces(:,:,6,3,:,:,klast)
     !$OMP END WORKSHARE
    end if
  END SUBROUTINE set_boundary_conditions

  SUBROUTINE taylor_swap(utdg,vtdg,updg,vpdg,udg,vdg,df,ifirst,ilast,jfirst,jlast,kfirst,klast,q_v,q_u) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg,updg,utdg
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg,vpdg,vtdg
    REAL(dp) :: df

    !$OMP WORKSHARE
    udg = utdg
    vdg = vtdg
    updg = updg + df*utdg
    vpdg = vpdg + df*vtdg
    !$OMP END WORKSHARE

  END SUBROUTINE taylor_swap

  SUBROUTINE start_next_step(updg,vpdg,udg,vdg,ifirst,ilast,jfirst,jlast,kfirst,klast,q_v,q_u) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg,updg
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg,vpdg

    !$OMP WORKSHARE
    updg = udg
    vpdg = vdg
    !$OMP END WORKSHARE

  END SUBROUTINE start_next_step

  SUBROUTINE get_initial_data(udg,vdg,ifirst,ilast,jfirst,jlast,kfirst,klast,&
       h,q_v,q_u,nint,id_type) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v,nint
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast,id_type
    REAL(dp), intent(in) :: h
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg
    integer :: i1,i2,i3
    integer :: l1,l2,l3
    integer :: k1,k2,k3
    integer :: ivar
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp), DIMENSION(0:nint,0:nint,0:nint) :: x_int,y_int,z_int,f_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp)  :: P_weights(0:q_u)

    !
    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    p_weights = 0.0_dp
    do l1 = 0,q_u
     do k1 = 0,nint
      P_weights(l1) = P_weights(l1) &
           + w_int(k1)*P_int(k1,l1)*P_int(k1,l1)
     end do
    end do

    ! do an l2 projection of the initial data
    udg = 0.0_dp
    vdg = 0.0_dp
    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,ivar,k3,k2,k1,x_int,y_int,z_int,f_int)
    do i3 = kfirst,klast
     do i2 = jfirst,jlast
      do i1 = ifirst,ilast
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          x_int(k1,k2,k3) = (real(i1,dp)-0.5_dp)*h + r_int(k1)/2.0_dp*h
          y_int(k1,k2,k3) = (real(i2,dp)-0.5_dp)*h + r_int(k2)/2.0_dp*h
          z_int(k1,k2,k3) = (real(i3,dp)-0.5_dp)*h + r_int(k3)/2.0_dp*h
         end do
        end do
       end do

       do ivar = 1,3
        if (id_type .eq. 0) then
         if (ivar.eq.1) f_int = 0.0_dp
         if (ivar.eq.2) f_int = 0.0_dp
         if (ivar.eq.3) f_int = 0.0_dp
        elseif(id_type .eq. 1) then
         if (ivar.eq.1) f_int = cos(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
         if (ivar.eq.2) f_int = sin(2.0_dp*pi*x_int)*cos(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
         if (ivar.eq.3) f_int = sin(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*cos(2.0_dp*pi*z_int)
        end if
        do l3 = 0,q_u
         do l2 = 0,q_u
          do l1 = 0,q_u
           do k3 = 0,nint
            do k2 = 0,nint
             do k1 = 0,nint
              udg(l1,l2,l3,ivar,i1,i2,i3) = udg(l1,l2,l3,ivar,i1,i2,i3)  &
                   + f_int(k1,k2,k3)&
                   *P_int(k1,l1)*P_int(k2,l2)*P_int(k3,l3) &
                   *w_int(k1)*w_int(k2)*w_int(k3)
             end do
            end do
           end do
           udg(l1,l2,l3,ivar,i1,i2,i3) = udg(l1,l2,l3,ivar,i1,i2,i3) &
                /(P_weights(l1)*P_weights(l2)*P_weights(l3))
          end do
         end do
        end do
       end do
      end do
     end do
    end do
    do i3 = kfirst,klast
     do i2 = jfirst,jlast
      do i1 = ifirst,ilast
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          x_int(k1,k2,k3) = (real(i1,dp)-0.5_dp)*h + r_int(k1)/2.0_dp*h
          y_int(k1,k2,k3) = (real(i2,dp)-0.5_dp)*h + r_int(k2)/2.0_dp*h
          z_int(k1,k2,k3) = (real(i3,dp)-0.5_dp)*h + r_int(k3)/2.0_dp*h
         end do
        end do
       end do
       do ivar = 1,3
        if (id_type .eq. 0) then
         if (ivar.eq.1) f_int = 0.0_dp
         if (ivar.eq.2) f_int = 0.0_dp
         if (ivar.eq.3) f_int = 0.0_dp
        elseif(id_type .eq. 1) then
         if (ivar.eq.1) f_int = 0.0_dp
         if (ivar.eq.2) f_int = 0.0_dp
         if (ivar.eq.3) f_int = 0.0_dp
        end if
        do l3 = 0,q_v
         do l2 = 0,q_v
          do l1 = 0,q_v
           do k3 = 0,nint
            do k2 = 0,nint
             do k1 = 0,nint
              vdg(l1,l2,l3,ivar,i1,i2,i3) = vdg(l1,l2,l3,ivar,i1,i2,i3)  &
                   + f_int(k1,k2,k3)&
                   *P_int(k1,l1)*P_int(k2,l2)*P_int(k3,l3) &
                   *w_int(k1)*w_int(k2)*w_int(k3)
             end do
            end do
           end do
           vdg(l1,l2,l3,ivar,i1,i2,i3) = vdg(l1,l2,l3,ivar,i1,i2,i3)&
                /(P_weights(l1)*P_weights(l2)*P_weights(l3))
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE get_initial_data

  SUBROUTINE ASSEMBLE(MU,MV,SU,SV,LU,LV,q_u,q_v,nint,hx) bind(c)
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in)  :: q_u,q_v,nint
    REAL(dp), intent(in) :: hx
    !
    ! MU*u_t = SV*v + "LU(v*-v)"
    ! MV*v_t = SU*u + "LV(n.u_x)"
    !
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_v+1)**3) :: MV
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_u+1)**3) :: MU
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_v+1)**3) :: SV
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_u+1)**3) :: SU
    REAL(dp), DIMENSION(0:nint,0:nint,(q_v+1)**3,6) :: LV
    REAL(dp), DIMENSION(0:nint,0:nint,(q_u+1)**3,3,3,6) :: LU
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp)  :: x_r,r_x,jac_loc,rho_loc
    INTEGER   :: i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,row,col,iface
    REAL(dp)  :: uxyz_loc(3,3),phixyz_loc(3,3),phi_x,phi_y,phi_z
    REAL(dp)  :: lam_loc,mu_loc,l2m ! Local Lame parameters

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    x_r = hx/2.0_dp
    r_x = 1.0_dp/x_r

    ! COMPUTE MU
    ! i represents the test function, i.e. the row of the matrix
    ! 4 is for the field component, 3,2,1 for dimensions z,y,x.
    lam_loc = 2.0_dp
    mu_loc = 1.0_dp
    l2m = lam_loc + 2.0_dp*mu_loc
    jac_loc = x_r**3

    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i3 = 0,q_u
     do i2 = 0,q_u
      do i1 = 0,q_u
       do i4 = 1,3
        row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
        ! j represents the expansion of the solution
        do j4 = 1,3
         do j3 = 0,q_u
          do j2 = 0,q_u
           do j1 = 0,q_u
            col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
            ! Loop over quadrature points
            MU(row,col) = 0.0_dp
            do k3 = 0,nint
             do k2 = 0,nint
              do k1 = 0,nint
               ! Compute various xyz derivatives...
               ! Note that this will be lengthier if we have a full metric
               ! but that the structure will be the same
               phixyz_loc = 0.0_dp
               phixyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
               phixyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
               phixyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
               uxyz_loc = 0.0_dp
               uxyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
               uxyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
               uxyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
               !
               MU(row,col) = MU(row,col) &
                    +jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                    +(l2m*phixyz_loc(1,1)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(3,3)))*uxyz_loc(1,1)&
                    +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(1,2)&
                    +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(1,3)&
                                !
                    +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(2,1)&
                    +(l2m*phixyz_loc(2,2)+lam_loc*(phixyz_loc(1,1)+phixyz_loc(3,3)))*uxyz_loc(2,2)&
                    +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(2,3)&
                                !
                    +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(3,1)&
                    +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(3,2)&
                    +(l2m*phixyz_loc(3,3)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(1,1)))*uxyz_loc(3,3)&
                    )
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to constants with independent equations.
    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i4 = 1,3
     i3 = 0
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     do j4 = 1,3
      do j3 = 0,q_u
       do j2 = 0,q_u
        do j1 = 0,q_u
         col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
         ! Loop over quadrature points
         MU(row,col) = 0.0_dp
         do k3 = 0,nint
          do k2 = 0,nint
           do k1 = 0,nint
            ! Compute the basis
            uxyz_loc = 0.0_dp
            uxyz_loc(j4,j4) = P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
            MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(i4,j4)
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to the remaining invariants
    i4 = 1
    i3 = 0
    i2 = 1
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 1 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          uxyz_loc(1,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,2)
         end do
        end do
       end do
       j4 = 2 ! Second component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(2,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,1)
         end do
        end do
       end do
      end do
     end do
    end do
    ! invariant 2
    i4 = 2
    i3 = 1
    i2 = 0
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 2 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(2,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,3)
         end do
        end do
       end do
       j4 = 3 ! second component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(3,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,2)
         end do
        end do
       end do
      end do
     end do
    end do
    ! invariant 3
    i4 = 3
    i3 = 0
    i2 = 0
    i1 = 1
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 1 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(1,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,3)
         end do
        end do
       end do
       j4 = 3 ! Last component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(3,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,1)
         end do
        end do
       end do
      end do
     end do
    end do


    ! NOW COMPUTE SV
    if(q_u .eq. q_v) then
     SV = MU
    else
     ! i represents the test function, i.e. the row of the matrix
     ! 4 is for the field component, 3,2,1 for dimensions z,y,x.
     do i4 = 1,3
      do i3 = 0,q_u
       do i2 = 0,q_u
        do i1 = 0,q_u
         row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
         ! j represents the expansion of the solution
         do j4 = 1,3
          do j3 = 0,q_v
           do j2 = 0,q_v
            do j1 = 0,q_v
             col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
             ! Loop over quadrature points
             SV(row,col) = 0.0_dp
             do k3 = 0,nint
              do k2 = 0,nint
               do k1 = 0,nint
                ! Compute various xyz derivatives...
                ! Note that this will be lengthier if we have a full metric
                ! but that the structure will be the same
                phixyz_loc = 0.0_dp
                phixyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
                phixyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
                phixyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
                uxyz_loc = 0.0_dp
                uxyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
                uxyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
                uxyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
                !
                lam_loc = 2.0_dp
                mu_loc = 1.0_dp
                l2m = lam_loc + 2.0_dp*mu_loc
                jac_loc = x_r**3
                SV(row,col) = SV(row,col) &
                     +jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                     +(l2m*phixyz_loc(1,1)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(3,3)))*uxyz_loc(1,1)&
                     +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(1,2)&
                     +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(1,3)&
                                !
                     +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(2,1)&
                     +(l2m*phixyz_loc(2,2)+lam_loc*(phixyz_loc(1,1)+phixyz_loc(3,3)))*uxyz_loc(2,2)&
                     +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(2,3)&
                                !
                     +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(3,1)&
                     +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(3,2)&
                     +(l2m*phixyz_loc(3,3)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(1,1)))*uxyz_loc(3,3)&
                     )
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do

     ! Replace the rows corresponding to constants with independent equations.
     do i4 = 1,3
      i3 = 0
      i2 = 0
      i1 = 0
      row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
      do j4 = 1,3
       do j3 = 0,q_v
        do j2 = 0,q_v
         do j1 = 0,q_v
          col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
          ! Loop over quadrature points
          SV(row,col) = 0.0_dp
          do k3 = 0,nint
           do k2 = 0,nint
            do k1 = 0,nint
             jac_loc = x_r**3
             ! Compute the basis
             uxyz_loc = 0.0_dp
             uxyz_loc(j4,j4) = P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
             SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(i4,j4)
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do

     ! Replace the rows corresponding to the remaining invariants

     i4 = 1
     i3 = 0
     i2 = 1
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 1 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           uxyz_loc(1,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,2)
          end do
         end do
        end do
        j4 = 2 ! Second component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(2,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,1)
          end do
         end do
        end do
       end do
      end do
     end do
     ! invariant 2
     i4 = 2
     i3 = 1
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 2 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(2,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,3)
          end do
         end do
        end do
        j4 = 3 ! second component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(3,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,2)
          end do
         end do
        end do
       end do
      end do
     end do
     ! invariant 3
     i4 = 3
     i3 = 0
     i2 = 0
     i1 = 1
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 1 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(1,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,3)
          end do
         end do
        end do
        j4 = 3 ! Last component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(3,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,1)
          end do
         end do
        end do
       end do
      end do
     end do
    end if

    ! NOW SU
    ! i represents the test function, i.e. the row of the matrix
    ! 4 is for the field component, 3,2,1 for dimensions z,y,x.

    lam_loc = 2.0_dp
    mu_loc = 1.0_dp
    l2m = lam_loc + 2.0_dp*mu_loc
    jac_loc = x_r**3

    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       do i4 = 1,3
        row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3 + (q_v+1)**3*(i4-1)
        ! j represents the expansion of the solution
        do j4 = 1,3
         do j3 = 0,q_u
          do j2 = 0,q_u
           do j1 = 0,q_u
            col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
            ! Loop over quadrature points
            SU(row,col) = 0.0_dp
            do k3 = 0,nint
             do k2 = 0,nint
              do k1 = 0,nint
               ! Compute various xyz derivatives...
               ! Note that this will be lengthier if we have a full metric
               ! but that the structure will be the same
               uxyz_loc = 0.0_dp
               uxyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
               uxyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
               uxyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
               phixyz_loc = 0.0_dp
               phixyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
               phixyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
               phixyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
               !
               SU(row,col) = SU(row,col) &
                    -jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                    +(l2m*uxyz_loc(1,1)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(3,3)))*phixyz_loc(1,1)&
                    +(mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2)))*phixyz_loc(1,2)&
                    +(mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3)))*phixyz_loc(1,3)&
                                !
                    +(mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2)))*phixyz_loc(2,1)&
                    +(l2m*uxyz_loc(2,2)+lam_loc*(uxyz_loc(1,1)+uxyz_loc(3,3)))*phixyz_loc(2,2)&
                    +(mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3)))*phixyz_loc(2,3)&
                                !
                    +(mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3)))*phixyz_loc(3,1)&
                    +(mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3)))*phixyz_loc(3,2)&
                    +(l2m*uxyz_loc(3,3)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(1,1)))*phixyz_loc(3,3)&
                    )
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    rho_loc = 1.0_dp
    jac_loc = x_r**3

    MV = 0.0_dp
    ! NOW MV

    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3
       do j3= 0,q_v
        do j2 = 0,q_v
         do j1 = 0,q_v
          col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3
          if(row.eq.col) then
           do k3= 0,nint
            do k2 = 0,nint
             do k1 = 0,nint
              MV(row,col) = MV(row,col) &
                   +rho_loc*jac_loc*w_int(k1)*w_int(k2)*w_int(k3)&
                   *P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)&
                   *P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
             end do
            end do
           end do
          else
           MV(row,col) = 0.0_dp
          end if
          ! copy the diagonal ...
          do k1 = 1,2
           MV(row+(q_v+1)**3*k1,col+(q_v+1)**3*k1) = MV(row,col)
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! NOW LV
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3
       ! x-sides
       do k3 = 0,nint
        do k2 = 0,nint
         k1 = 0
         LV(k2,k3,row,1) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k1 = nint
         LV(k2,k3,row,2) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
       ! y-sides
       do k3 = 0,nint
        do k1 = 0,nint
         k2 = 0
         LV(k1,k3,row,3) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k2 = nint
         LV(k1,k3,row,4) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
       ! z-sides
       do k2 = 0,nint
        do k1 = 0,nint
         k3 = 0
         LV(k1,k2,row,5) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k3 = nint
         LV(k1,k2,row,6) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
      end do
     end do
    end do

    ! NOW LU
    LU = 0.d0
    ! First set of test functions
    ! LU("quad points","row",1,"v_i^*-v_i",side)
    ! LU("quad points","row",1,1,x-sides) = (2*\mu+\lambda)*\phi_x
    ! LU("quad points","row",1,2,x-sides) = \mu*\phi_y
    ! LU("quad points","row",1,3,x-sides) = \mu*\phi_z
    !
    ! LU("quad points","row",1,1,y-sides) = \mu*\phi_y
    ! LU("quad points","row",1,2,y-sides) = \lambda*\phi_x
    ! LU("quad points","row",1,3,y-sides) = 0
    !
    ! LU("quad points","row",1,1,z-sides) = \mu*\phi_z
    ! LU("quad points","row",1,2,z-sides) = 0
    ! LU("quad points","row",1,3,z-sides) = \lambda*\phi_x

    ! Second set of test functions
    ! LU("quad points","row",2,1,x-sides) = \lambda*\phi_y
    ! LU("quad points","row",2,2,x-sides) = \mu*\phi_x
    ! LU("quad points","row",2,3,x-sides) = 0
    !
    ! LU("quad points","row",2,1,y-sides) = \mu*\phi_x
    ! LU("quad points","row",2,2,y-sides) = (2*\mu+\lambda)*\phi_y
    ! LU("quad points","row",2,3,y-sides) = \mu*\phi_z
    !
    ! LU("quad points","row",2,1,z-sides) = 0
    ! LU("quad points","row",2,2,z-sides) = \mu*\phi_z
    ! LU("quad points","row",2,3,z-sides) = \lambda*\phi_y

    ! Third set of test functions
    ! LU("quad points","row",3,1,x-sides) = \lambda*\phi_z
    ! LU("quad points","row",3,2,x-sides) = 0
    ! LU("quad points","row",3,3,x-sides) = \mu*\phi_x
    !
    ! LU("quad points","row",3,1,y-sides) = 0
    ! LU("quad points","row",3,2,y-sides) = \lambda*\phi_z
    ! LU("quad points","row",3,3,y-sides) = \mu*\phi_y
    !
    ! LU("quad points","row",3,1,z-sides) = \mu*\phi_x
    ! LU("quad points","row",3,2,z-sides) = \mu*\phi_y
    ! LU("quad points","row",3,3,z-sides) = (2*\mu+\lambda)*\phi_z

    do i3 = 0,q_u
     do i2 = 0,q_u
      do i1 = 0,q_u
       row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
       ! x-sides
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint,nint
          lam_loc = 2.0_dp
          mu_loc = 1.0_dp
          l2m = lam_loc + 2.0_dp*mu_loc
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 1 + k1/nint ! 1 or 2
          LU(k2,k3,row,1,1,iface) = l2m*phi_x
          LU(k2,k3,row,1,2,iface) = mu_loc*phi_y
          LU(k2,k3,row,1,3,iface) = mu_loc*phi_z
          LU(k2,k3,row,2,1,iface) = lam_loc*phi_y
          LU(k2,k3,row,2,2,iface) = mu_loc*phi_x
          LU(k2,k3,row,2,3,iface) = 0.0_dp
          LU(k2,k3,row,3,1,iface) = lam_loc*phi_z
          LU(k2,k3,row,3,2,iface) = 0.0_dp
          LU(k2,k3,row,3,3,iface) = mu_loc*phi_x
         end do
        end do
       end do
       ! y-sides
       do k3 = 0,nint
        do k2 = 0,nint,nint
         do k1 = 0,nint
          lam_loc = 2.0_dp
          mu_loc = 1.0_dp
          l2m = lam_loc + 2.0_dp*mu_loc
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 3 + k2/nint ! 3 or 4
          LU(k1,k3,row,1,1,iface) = mu_loc*phi_y
          LU(k1,k3,row,1,2,iface) = lam_loc*phi_x
          LU(k1,k3,row,1,3,iface) = 0.0_dp
          LU(k1,k3,row,2,1,iface) = mu_loc*phi_x
          LU(k1,k3,row,2,2,iface) = l2m*phi_y
          LU(k1,k3,row,2,3,iface) = mu_loc*phi_z
          LU(k1,k3,row,3,1,iface) = 0.0_dp
          LU(k1,k3,row,3,2,iface) = lam_loc*phi_z
          LU(k1,k3,row,3,3,iface) = mu_loc*phi_y
         end do
        end do
       end do
       ! z-sides
       do k3 = 0,nint,nint
        do k2 = 0,nint
         do k1 = 0,nint
          lam_loc = 2.0_dp
          mu_loc = 1.0_dp
          l2m = lam_loc + 2.0_dp*mu_loc
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 5 + k3/nint ! 5 or 6
          LU(k1,k2,row,1,1,iface) = mu_loc*phi_z
          LU(k1,k2,row,1,2,iface) = 0.0_dp
          LU(k1,k2,row,1,3,iface) = lam_loc*phi_x
          LU(k1,k2,row,2,1,iface) = 0.0_dp
          LU(k1,k2,row,2,2,iface) = mu_loc*phi_z
          LU(k1,k2,row,2,3,iface) = lam_loc*phi_y
          LU(k1,k2,row,3,1,iface) = mu_loc*phi_x
          LU(k1,k2,row,3,2,iface) = mu_loc*phi_y
          LU(k1,k2,row,3,3,iface) = l2m*phi_z
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to constants with zeros
    do i4 = 1,3
     i3 = 0
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
     LU(:,:,row,i4,:,:) = 0.0_dp
    end do
    i4 = 1
    i3 = 0
    i2 = 1
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp
    ! invariant 2
    i4 = 2
    i3 = 1
    i2 = 0
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp
    ! invariant 3
    i4 = 3
    i3 = 0
    i2 = 0
    i1 = 1
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp

  END SUBROUTINE ASSEMBLE

  SUBROUTINE ASSEMBLE_CONST_COEFF(MU,MV,SU,SV,LU,LV,q_u,q_v,nint,hx,lambda_lame,mu_lame,rho) bind(c)
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(inout)  :: q_u,q_v,nint
    REAL(dp), intent(inout) :: hx,lambda_lame,mu_lame,rho
    !
    ! MU*u_t = SV*v + "LU(v*-v)"
    ! MV*v_t = SU*u + "LV(n.u_x)"
    !
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_v+1)**3) :: MV
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_u+1)**3) :: MU
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_v+1)**3) :: SV
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_u+1)**3) :: SU
    REAL(dp), DIMENSION(0:nint,0:nint,(q_v+1)**3,6) :: LV
    REAL(dp), DIMENSION(0:nint,0:nint,(q_u+1)**3,3,3,6) :: LU
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp)  :: x_r,r_x,jac_loc,rho_loc
    INTEGER   :: i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,row,col,iface
    REAL(dp)  :: uxyz_loc(3,3),phixyz_loc(3,3),phi_x,phi_y,phi_z
    REAL(dp)  :: lam_loc,mu_loc,l2m ! Local Lame parameters

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    x_r = hx/2.0_dp
    r_x = 1.0_dp/x_r

    ! COMPUTE MU
    ! i represents the test function, i.e. the row of the matrix
    ! 4 is for the field component, 3,2,1 for dimensions z,y,x.
    lam_loc = lambda_lame
    mu_loc = mu_lame
    l2m = lam_loc + 2.0_dp*mu_loc
    rho_loc = rho
    jac_loc = x_r**3

    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i3 = 0,q_u
     do i2 = 0,q_u
      do i1 = 0,q_u
       do i4 = 1,3
        row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
        ! j represents the expansion of the solution
        do j4 = 1,3
         do j3 = 0,q_u
          do j2 = 0,q_u
           do j1 = 0,q_u
            col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
            ! Loop over quadrature points
            MU(row,col) = 0.0_dp
            do k3 = 0,nint
             do k2 = 0,nint
              do k1 = 0,nint
               ! Compute various xyz derivatives...
               ! Note that this will be lengthier if we have a full metric
               ! but that the structure will be the same
               phixyz_loc = 0.0_dp
               phixyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
               phixyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
               phixyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
               uxyz_loc = 0.0_dp
               uxyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
               uxyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
               uxyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
               !
               MU(row,col) = MU(row,col) &
                    +jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                    +(l2m*phixyz_loc(1,1)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(3,3)))*uxyz_loc(1,1)&
                    +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(1,2)&
                    +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(1,3)&
                                !
                    +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(2,1)&
                    +(l2m*phixyz_loc(2,2)+lam_loc*(phixyz_loc(1,1)+phixyz_loc(3,3)))*uxyz_loc(2,2)&
                    +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(2,3)&
                                !
                    +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(3,1)&
                    +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(3,2)&
                    +(l2m*phixyz_loc(3,3)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(1,1)))*uxyz_loc(3,3)&
                    )
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to constants with independent equations.
    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i4 = 1,3
     i3 = 0
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     do j4 = 1,3
      do j3 = 0,q_u
       do j2 = 0,q_u
        do j1 = 0,q_u
         col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
         ! Loop over quadrature points
         MU(row,col) = 0.0_dp
         do k3 = 0,nint
          do k2 = 0,nint
           do k1 = 0,nint
            ! Compute the basis
            uxyz_loc = 0.0_dp
            uxyz_loc(j4,j4) = P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
            MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(i4,j4)
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to the remaining invariants
    i4 = 1
    i3 = 0
    i2 = 1
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 1 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          uxyz_loc(1,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,2)
         end do
        end do
       end do
       j4 = 2 ! Second component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(2,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,1)
         end do
        end do
       end do
      end do
     end do
    end do
    ! invariant 2
    i4 = 2
    i3 = 1
    i2 = 0
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 2 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(2,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,3)
         end do
        end do
       end do
       j4 = 3 ! second component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(3,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,2)
         end do
        end do
       end do
      end do
     end do
    end do
    ! invariant 3
    i4 = 3
    i3 = 0
    i2 = 0
    i1 = 1
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
    MU(row,:) = 0.0_dp ! Zero out the full row!
    do j3 = 0,q_u
     do j2 = 0,q_u
      do j1 = 0,q_u
       j4 = 1 ! First component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(1,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
          MU(row,col) = MU(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,3)
         end do
        end do
       end do
       j4 = 3 ! Last component
       col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
       ! Loop over quadrature points
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint
          jac_loc = x_r**3
          ! Compute the basis
          uxyz_loc(3,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
          MU(row,col) = MU(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,1)
         end do
        end do
       end do
      end do
     end do
    end do


    ! NOW COMPUTE SV
    if(q_u .eq. q_v) then
     SV = MU
    else
     ! i represents the test function, i.e. the row of the matrix
     ! 4 is for the field component, 3,2,1 for dimensions z,y,x.
     do i4 = 1,3
      do i3 = 0,q_u
       do i2 = 0,q_u
        do i1 = 0,q_u
         row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
         ! j represents the expansion of the solution
         do j4 = 1,3
          do j3 = 0,q_v
           do j2 = 0,q_v
            do j1 = 0,q_v
             col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
             ! Loop over quadrature points
             SV(row,col) = 0.0_dp
             do k3 = 0,nint
              do k2 = 0,nint
               do k1 = 0,nint
                ! Compute various xyz derivatives...
                ! Note that this will be lengthier if we have a full metric
                ! but that the structure will be the same
                phixyz_loc = 0.0_dp
                phixyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
                phixyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
                phixyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
                uxyz_loc = 0.0_dp
                uxyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
                uxyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
                uxyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
                !
                SV(row,col) = SV(row,col) &
                     +jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                     +(l2m*phixyz_loc(1,1)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(3,3)))*uxyz_loc(1,1)&
                     +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(1,2)&
                     +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(1,3)&
                                !
                     +(mu_loc*(phixyz_loc(2,1)+phixyz_loc(1,2)))*uxyz_loc(2,1)&
                     +(l2m*phixyz_loc(2,2)+lam_loc*(phixyz_loc(1,1)+phixyz_loc(3,3)))*uxyz_loc(2,2)&
                     +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(2,3)&
                                !
                     +(mu_loc*(phixyz_loc(3,1)+phixyz_loc(1,3)))*uxyz_loc(3,1)&
                     +(mu_loc*(phixyz_loc(3,2)+phixyz_loc(2,3)))*uxyz_loc(3,2)&
                     +(l2m*phixyz_loc(3,3)+lam_loc*(phixyz_loc(2,2)+phixyz_loc(1,1)))*uxyz_loc(3,3)&
                     )
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do

     ! Replace the rows corresponding to constants with independent equations.
     do i4 = 1,3
      i3 = 0
      i2 = 0
      i1 = 0
      row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
      do j4 = 1,3
       do j3 = 0,q_v
        do j2 = 0,q_v
         do j1 = 0,q_v
          col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
          ! Loop over quadrature points
          SV(row,col) = 0.0_dp
          do k3 = 0,nint
           do k2 = 0,nint
            do k1 = 0,nint
             jac_loc = x_r**3
             ! Compute the basis
             uxyz_loc = 0.0_dp
             uxyz_loc(j4,j4) = P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
             SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(i4,j4)
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do

     ! Replace the rows corresponding to the remaining invariants

     i4 = 1
     i3 = 0
     i2 = 1
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 1 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           uxyz_loc(1,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,2)
          end do
         end do
        end do
        j4 = 2 ! Second component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(2,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,1)
          end do
         end do
        end do
       end do
      end do
     end do
     ! invariant 2
     i4 = 2
     i3 = 1
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 2 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(2,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(2,3)
          end do
         end do
        end do
        j4 = 3 ! second component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(3,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,2)
          end do
         end do
        end do
       end do
      end do
     end do
     ! invariant 3
     i4 = 3
     i3 = 0
     i2 = 0
     i1 = 1
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3 + (q_u+1)**3*(i4-1)
     SV(row,:) = 0.0_dp ! Zero out the full row!
     do j3 = 0,q_v
      do j2 = 0,q_v
       do j1 = 0,q_v
        j4 = 1 ! First component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(1,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
           SV(row,col) = SV(row,col) + jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(1,3)
          end do
         end do
        end do
        j4 = 3 ! Last component
        col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3 + (q_v+1)**3*(j4-1)
        ! Loop over quadrature points
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           jac_loc = x_r**3
           ! Compute the basis
           uxyz_loc(3,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
           SV(row,col) = SV(row,col) - jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*uxyz_loc(3,1)
          end do
         end do
        end do
       end do
      end do
     end do
    end if

    ! NOW SU
    ! i represents the test function, i.e. the row of the matrix
    ! 4 is for the field component, 3,2,1 for dimensions z,y,x.
    !$OMP PARALLEL DO PRIVATE(i4,i3,i2,i1,row,j4,j3,j2,j1,col,k1,k2,k3,phixyz_loc,uxyz_loc)
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       do i4 = 1,3
        row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3 + (q_v+1)**3*(i4-1)
        ! j represents the expansion of the solution
        do j4 = 1,3
         do j3 = 0,q_u
          do j2 = 0,q_u
           do j1 = 0,q_u
            col = 1 + j1 + (q_u+1)*j2 + (q_u+1)**2*j3 + (q_u+1)**3*(j4-1)
            ! Loop over quadrature points
            SU(row,col) = 0.0_dp
            do k3 = 0,nint
             do k2 = 0,nint
              do k1 = 0,nint
               ! Compute various xyz derivatives...
               ! Note that this will be lengthier if we have a full metric
               ! but that the structure will be the same
               uxyz_loc = 0.0_dp
               uxyz_loc(i4,1) = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
               uxyz_loc(i4,2) = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
               uxyz_loc(i4,3) = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
               phixyz_loc = 0.0_dp
               phixyz_loc(j4,1) = r_x*Pr_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
               phixyz_loc(j4,2) = P_int(k1,j1)*r_x*Pr_int(k2,j2)*P_int(k3,j3)
               phixyz_loc(j4,3) = P_int(k1,j1)*P_int(k2,j2)*r_x*Pr_int(k3,j3)
               !
               SU(row,col) = SU(row,col) &
                    -jac_loc*w_int(k1)*w_int(k2)*w_int(k3)*(&
                                !
                    +(l2m*uxyz_loc(1,1)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(3,3)))*phixyz_loc(1,1)&
                    +(mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2)))*phixyz_loc(1,2)&
                    +(mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3)))*phixyz_loc(1,3)&
                                !
                    +(mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2)))*phixyz_loc(2,1)&
                    +(l2m*uxyz_loc(2,2)+lam_loc*(uxyz_loc(1,1)+uxyz_loc(3,3)))*phixyz_loc(2,2)&
                    +(mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3)))*phixyz_loc(2,3)&
                                !
                    +(mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3)))*phixyz_loc(3,1)&
                    +(mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3)))*phixyz_loc(3,2)&
                    +(l2m*uxyz_loc(3,3)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(1,1)))*phixyz_loc(3,3)&
                    )
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do


    jac_loc = x_r**3

    ! NOW MV
    MV = 0.0_dp
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3
       do j3= 0,q_v
        do j2 = 0,q_v
         do j1 = 0,q_v
          col = 1 + j1 + (q_v+1)*j2 + (q_v+1)**2*j3
          if(row.eq.col) then
           do k3= 0,nint
            do k2 = 0,nint
             do k1 = 0,nint
              MV(row,col) = MV(row,col) &
                   +rho_loc*jac_loc*w_int(k1)*w_int(k2)*w_int(k3)&
                   *P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)&
                   *P_int(k1,j1)*P_int(k2,j2)*P_int(k3,j3)
             end do
            end do
           end do
          else
           MV(row,col) = 0.0_dp
          end if
          ! copy the diagonal ...
          do k1 = 1,2
           MV(row+(q_v+1)**3*k1,col+(q_v+1)**3*k1) = MV(row,col)
          end do
         end do
        end do
       end do
      end do
     end do
    end do

    ! NOW LV
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3
       ! x-sides
       do k3 = 0,nint
        do k2 = 0,nint
         k1 = 0
         LV(k2,k3,row,1) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k1 = nint
         LV(k2,k3,row,2) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
       ! y-sides
       do k3 = 0,nint
        do k1 = 0,nint
         k2 = 0
         LV(k1,k3,row,3) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k2 = nint
         LV(k1,k3,row,4) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
       ! z-sides
       do k2 = 0,nint
        do k1 = 0,nint
         k3 = 0
         LV(k1,k2,row,5) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
         k3 = nint
         LV(k1,k2,row,6) = P_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
        end do
       end do
      end do
     end do
    end do

    ! NOW LU
    LU = 0.d0
    ! First set of test functions
    ! LU("quad points","row",1,"v_i^*-v_i",side)
    ! LU("quad points","row",1,1,x-sides) = (2*\mu+\lambda)*\phi_x
    ! LU("quad points","row",1,2,x-sides) = \mu*\phi_y
    ! LU("quad points","row",1,3,x-sides) = \mu*\phi_z
    !
    ! LU("quad points","row",1,1,y-sides) = \mu*\phi_y
    ! LU("quad points","row",1,2,y-sides) = \lambda*\phi_x
    ! LU("quad points","row",1,3,y-sides) = 0
    !
    ! LU("quad points","row",1,1,z-sides) = \mu*\phi_z
    ! LU("quad points","row",1,2,z-sides) = 0
    ! LU("quad points","row",1,3,z-sides) = \lambda*\phi_x

    ! Second set of test functions
    ! LU("quad points","row",2,1,x-sides) = \lambda*\phi_y
    ! LU("quad points","row",2,2,x-sides) = \mu*\phi_x
    ! LU("quad points","row",2,3,x-sides) = 0
    !
    ! LU("quad points","row",2,1,y-sides) = \mu*\phi_x
    ! LU("quad points","row",2,2,y-sides) = (2*\mu+\lambda)*\phi_y
    ! LU("quad points","row",2,3,y-sides) = \mu*\phi_z
    !
    ! LU("quad points","row",2,1,z-sides) = 0
    ! LU("quad points","row",2,2,z-sides) = \mu*\phi_z
    ! LU("quad points","row",2,3,z-sides) = \lambda*\phi_y

    ! Third set of test functions
    ! LU("quad points","row",3,1,x-sides) = \lambda*\phi_z
    ! LU("quad points","row",3,2,x-sides) = 0
    ! LU("quad points","row",3,3,x-sides) = \mu*\phi_x
    !
    ! LU("quad points","row",3,1,y-sides) = 0
    ! LU("quad points","row",3,2,y-sides) = \lambda*\phi_z
    ! LU("quad points","row",3,3,y-sides) = \mu*\phi_y
    !
    ! LU("quad points","row",3,1,z-sides) = \mu*\phi_x
    ! LU("quad points","row",3,2,z-sides) = \mu*\phi_y
    ! LU("quad points","row",3,3,z-sides) = (2*\mu+\lambda)*\phi_z

    do i3 = 0,q_u
     do i2 = 0,q_u
      do i1 = 0,q_u
       row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
       ! x-sides
       do k3 = 0,nint
        do k2 = 0,nint
         do k1 = 0,nint,nint
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 1 + k1/nint ! 1 or 2
          LU(k2,k3,row,1,1,iface) = l2m*phi_x
          LU(k2,k3,row,1,2,iface) = mu_loc*phi_y
          LU(k2,k3,row,1,3,iface) = mu_loc*phi_z
          LU(k2,k3,row,2,1,iface) = lam_loc*phi_y
          LU(k2,k3,row,2,2,iface) = mu_loc*phi_x
          LU(k2,k3,row,2,3,iface) = 0.0_dp
          LU(k2,k3,row,3,1,iface) = lam_loc*phi_z
          LU(k2,k3,row,3,2,iface) = 0.0_dp
          LU(k2,k3,row,3,3,iface) = mu_loc*phi_x
         end do
        end do
       end do
       ! y-sides
       do k3 = 0,nint
        do k2 = 0,nint,nint
         do k1 = 0,nint
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 3 + k2/nint ! 3 or 4
          LU(k1,k3,row,1,1,iface) = mu_loc*phi_y
          LU(k1,k3,row,1,2,iface) = lam_loc*phi_x
          LU(k1,k3,row,1,3,iface) = 0.0_dp
          LU(k1,k3,row,2,1,iface) = mu_loc*phi_x
          LU(k1,k3,row,2,2,iface) = l2m*phi_y
          LU(k1,k3,row,2,3,iface) = mu_loc*phi_z
          LU(k1,k3,row,3,1,iface) = 0.0_dp
          LU(k1,k3,row,3,2,iface) = lam_loc*phi_z
          LU(k1,k3,row,3,3,iface) = mu_loc*phi_y
         end do
        end do
       end do
       ! z-sides
       do k3 = 0,nint,nint
        do k2 = 0,nint
         do k1 = 0,nint
          phi_x = r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
          phi_y = P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
          phi_z = P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
          iface = 5 + k3/nint ! 5 or 6
          LU(k1,k2,row,1,1,iface) = mu_loc*phi_z
          LU(k1,k2,row,1,2,iface) = 0.0_dp
          LU(k1,k2,row,1,3,iface) = lam_loc*phi_x
          LU(k1,k2,row,2,1,iface) = 0.0_dp
          LU(k1,k2,row,2,2,iface) = mu_loc*phi_z
          LU(k1,k2,row,2,3,iface) = lam_loc*phi_y
          LU(k1,k2,row,3,1,iface) = mu_loc*phi_x
          LU(k1,k2,row,3,2,iface) = mu_loc*phi_y
          LU(k1,k2,row,3,3,iface) = l2m*phi_z
         end do
        end do
       end do
      end do
     end do
    end do

    ! Replace the rows corresponding to constants with zeros
    do i4 = 1,3
     i3 = 0
     i2 = 0
     i1 = 0
     row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
     LU(:,:,row,i4,:,:) = 0.0_dp
    end do
    i4 = 1
    i3 = 0
    i2 = 1
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp
    ! invariant 2
    i4 = 2
    i3 = 1
    i2 = 0
    i1 = 0
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp
    ! invariant 3
    i4 = 3
    i3 = 0
    i2 = 0
    i1 = 1
    row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
    LU(:,:,row,i4,:,:) = 0.0_dp

  END SUBROUTINE ASSEMBLE_CONST_COEFF


  SUBROUTINE build_my_v(vdg,v_all_faces,ifirst,ilast,jfirst,jlast,kfirst,klast,&
       q_v,nint) bind(c)

    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER, intent(in) :: q_v,nint
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast)&
         :: v_all_faces

    ! (quad_points,quad_points,sides,fields,elements)
    integer :: i1,i2,i3
    integer :: j1,j2,j3
    integer :: k1,k2,k3
    integer :: iface
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_v),Pr_int(0:nint,0:q_v)
    REAL(dp)  :: p1,p2,p3,p_loc

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_v)
    !
    v_all_faces = 0.0_dp
    ! When done, w will hold the components of the stress tensor on all faces
    ! and on all elements

    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,j3,j2,j1,k1,k2,k3,p1,p2,p3,iface)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_v
        do i2 = 0,q_v
         do i1 = 0,q_v
          ! x - side
          do k1 = 0,nint,nint
           iface = 1 + k1/nint ! 1 or 2.
           p1 = P_int(k1,i1)
           do k3 = 0,nint
            p3 = P_int(k3,i3)*p1
            do k2 = 0,nint
             p_loc = P_int(k2,i2)*p3
             v_all_faces(k2,k3,iface,1:3,j1,j2,j3) = v_all_faces(k2,k3,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
          ! y-faces
          do k2 = 0,nint,nint
           iface = 3 + k2/nint ! 3 or 4.
           p2 = P_int(k2,i2)
           do k3 = 0,nint
            p3 = P_int(k3,i3)*p2
            do k1 = 0,nint
             p_loc = P_int(k1,i1)*p3
             v_all_faces(k1,k3,iface,1:3,j1,j2,j3) = v_all_faces(k1,k3,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
          ! z-faces
          do k3 = 0,nint,nint
           iface = 5 + k3/nint ! 5 or 6.
           p3 = P_int(k3,i3)
           do k2 = 0,nint
            p2 = P_int(k2,i2)*p3
            do k1 = 0,nint
             p_loc = P_int(k1,i1)*p2
             v_all_faces(k1,k2,iface,1:3,j1,j2,j3) = v_all_faces(k1,k2,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE build_my_v

  SUBROUTINE build_my_w(udg,w_all_faces,ifirst,ilast,jfirst,jlast,kfirst,klast,&
       h,q_u,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER,  intent(in) :: q_u,nint
    REAL(dp), intent(in) :: h
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast)&
         :: w_all_faces ! (quad_points,quad_points,sides,fields,elements)
    integer :: i1,i2,i3
    integer :: j1,j2,j3
    integer :: k1,k2,k3
    integer :: ivar,iface
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp) :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp) :: lam_loc,mu_loc,l2m ! Local Lame parameters
    REAL(dp) :: x_r,r_x
    REAL(dp) :: uxyz_loc(3,3),m_c
    !
    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    x_r = h/2.0_dp
    r_x = 1.0_dp/x_r
    !
    w_all_faces = 0.0_dp
    ! When done, w will hold the components of the stress tensor on all faces
    ! and on all elements
    lam_loc = 2.0_dp
    mu_loc = 1.0_dp
    l2m = lam_loc + 2.0_dp*mu_loc
    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,j3,j2,j1,k1,k2,k3,uxyz_loc,m_c,ivar,iface)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_u
        do i2 = 0,q_u
         do i1 = 0,q_u
          ! x-faces
          do k3 = 0,nint
           do k2 = 0,nint
            do k1 = 0,nint,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ! Note that this could be made much more efficient by moving
             ! the computations out in the loop hierarcy
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 1 + k1/nint ! 1 or 2.
             w_all_faces(k2,k3,iface,1,j1,j2,j3) = w_all_faces(k2,k3,iface,1,j1,j2,j3) &
                  +(l2m*uxyz_loc(1,1)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(3,3)))
             w_all_faces(k2,k3,iface,2,j1,j2,j3) = w_all_faces(k2,k3,iface,2,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2))
             w_all_faces(k2,k3,iface,3,j1,j2,j3) = w_all_faces(k2,k3,iface,3,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3))
            end do
           end do
          end do
          ! y-faces
          do k3 = 0,nint
           do k2 = 0,nint,nint
            do k1 = 0,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 3 + k2/nint ! 3 or 4.
             w_all_faces(k1,k3,iface,1,j1,j2,j3) = w_all_faces(k1,k3,iface,1,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2))
             w_all_faces(k1,k3,iface,2,j1,j2,j3) = w_all_faces(k1,k3,iface,2,j1,j2,j3) &
                  +(l2m*uxyz_loc(2,2)+lam_loc*(uxyz_loc(1,1)+uxyz_loc(3,3)))
             w_all_faces(k1,k3,iface,3,j1,j2,j3) = w_all_faces(k1,k3,iface,3,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3))
            end do
           end do
          end do
          ! z-faces
          do k3 = 0,nint,nint
           do k2 = 0,nint
            do k1 = 0,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 5 + k3/nint ! 5 or 6.
             w_all_faces(k1,k2,iface,1,j1,j2,j3) = w_all_faces(k1,k2,iface,1,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3))
             w_all_faces(k1,k2,iface,2,j1,j2,j3) = w_all_faces(k1,k2,iface,2,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3))
             w_all_faces(k1,k2,iface,3,j1,j2,j3) = w_all_faces(k1,k2,iface,3,j1,j2,j3) &
                  +(l2m*uxyz_loc(3,3)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(1,1)))
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE build_my_w

  SUBROUTINE build_my_v_const_coeff(vdg,v_all_faces,ifirst,ilast,jfirst,jlast,kfirst,klast,&
       q_v,nint) bind(c)

    USE type_defs
    IMPLICIT NONE
    INTEGER :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER :: q_v,nint
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast)&
         :: v_all_faces

    ! (quad_points,quad_points,sides,fields,elements)
    integer :: i1,i2,i3
    integer :: j1,j2,j3
    integer :: k1,k2,k3
    integer :: iface
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_v),Pr_int(0:nint,0:q_v)
    REAL(dp)  :: p1,p2,p3,p_loc

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_v)
    !
    v_all_faces = 0.0_dp
    ! When done, w will hold the components of the stress tensor on all faces
    ! and on all elements

    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,j3,j2,j1,k1,k2,k3,p1,p2,p3,iface)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_v
        do i2 = 0,q_v
         do i1 = 0,q_v
          ! x - side
          do k1 = 0,nint,nint
           iface = 1 + k1/nint ! 1 or 2.
           p1 = P_int(k1,i1)
           do k3 = 0,nint
            p3 = P_int(k3,i3)*p1
            do k2 = 0,nint
             p_loc = P_int(k2,i2)*p3
             v_all_faces(k2,k3,iface,1:3,j1,j2,j3) = v_all_faces(k2,k3,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
          ! y-faces
          do k2 = 0,nint,nint
           iface = 3 + k2/nint ! 3 or 4.
           p2 = P_int(k2,i2)
           do k3 = 0,nint
            p3 = P_int(k3,i3)*p2
            do k1 = 0,nint
             p_loc = P_int(k1,i1)*p3
             v_all_faces(k1,k3,iface,1:3,j1,j2,j3) = v_all_faces(k1,k3,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
          ! z-faces
          do k3 = 0,nint,nint
           iface = 5 + k3/nint ! 5 or 6.
           p3 = P_int(k3,i3)
           do k2 = 0,nint
            p2 = P_int(k2,i2)*p3
            do k1 = 0,nint
             p_loc = P_int(k1,i1)*p2
             v_all_faces(k1,k2,iface,1:3,j1,j2,j3) = v_all_faces(k1,k2,iface,1:3,j1,j2,j3) &
                  + vdg(i1,i2,i3,1:3,j1,j2,j3)*p_loc
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE build_my_v_const_coeff

  SUBROUTINE build_my_w_const_coeff(udg,w_all_faces,ifirst,ilast,jfirst,jlast,kfirst,klast,&
       h,q_u,nint,lambda_lame,mu_lame) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER,  intent(in) :: q_u,nint
    REAL(dp), intent(in) :: h,lambda_lame,mu_lame
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast)&
         :: w_all_faces ! (quad_points,quad_points,sides,fields,elements)
    integer :: i1,i2,i3
    integer :: j1,j2,j3
    integer :: k1,k2,k3
    integer :: ivar,iface
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp) :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp) :: lam_loc,mu_loc,l2m ! Local Lame parameters
    REAL(dp) :: x_r,r_x
    REAL(dp) :: uxyz_loc(3,3),m_c
    !
    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    x_r = h/2.0_dp
    r_x = 1.0_dp/x_r
    !
    w_all_faces = 0.0_dp
    ! When done, w will hold the components of the stress tensor on all faces
    ! and on all elements
    lam_loc = lambda_lame
    mu_loc = mu_lame
    l2m = lam_loc + 2.0_dp*mu_loc
    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,j3,j2,j1,k1,k2,k3,uxyz_loc,m_c,ivar,iface)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_u
        do i2 = 0,q_u
         do i1 = 0,q_u
          ! x-faces
          do k3 = 0,nint
           do k2 = 0,nint
            do k1 = 0,nint,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ! Note that this could be made much more efficient by moving
             ! the computations out in the loop hierarcy
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 1 + k1/nint ! 1 or 2.
             w_all_faces(k2,k3,iface,1,j1,j2,j3) = w_all_faces(k2,k3,iface,1,j1,j2,j3) &
                  +(l2m*uxyz_loc(1,1)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(3,3)))
             w_all_faces(k2,k3,iface,2,j1,j2,j3) = w_all_faces(k2,k3,iface,2,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2))
             w_all_faces(k2,k3,iface,3,j1,j2,j3) = w_all_faces(k2,k3,iface,3,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3))
            end do
           end do
          end do
          ! y-faces
          do k3 = 0,nint
           do k2 = 0,nint,nint
            do k1 = 0,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 3 + k2/nint ! 3 or 4.
             w_all_faces(k1,k3,iface,1,j1,j2,j3) = w_all_faces(k1,k3,iface,1,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(2,1)+uxyz_loc(1,2))
             w_all_faces(k1,k3,iface,2,j1,j2,j3) = w_all_faces(k1,k3,iface,2,j1,j2,j3) &
                  +(l2m*uxyz_loc(2,2)+lam_loc*(uxyz_loc(1,1)+uxyz_loc(3,3)))
             w_all_faces(k1,k3,iface,3,j1,j2,j3) = w_all_faces(k1,k3,iface,3,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3))
            end do
           end do
          end do
          ! z-faces
          do k3 = 0,nint,nint
           do k2 = 0,nint
            do k1 = 0,nint
             uxyz_loc = 0.0_dp
             ! A single mode of x,y,z derivatives.
             ivar = 1
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 2
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ivar = 3
             m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
             uxyz_loc(ivar,1) = m_c*r_x*Pr_int(k1,i1)*P_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,2) = m_c*P_int(k1,i1)*r_x*Pr_int(k2,i2)*P_int(k3,i3)
             uxyz_loc(ivar,3) = m_c*P_int(k1,i1)*P_int(k2,i2)*r_x*Pr_int(k3,i3)
             ! Local material parameters, can be fixed later
             !
             iface = 5 + k3/nint ! 5 or 6.
             w_all_faces(k1,k2,iface,1,j1,j2,j3) = w_all_faces(k1,k2,iface,1,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,1)+uxyz_loc(1,3))
             w_all_faces(k1,k2,iface,2,j1,j2,j3) = w_all_faces(k1,k2,iface,2,j1,j2,j3) &
                  +mu_loc*(uxyz_loc(3,2)+uxyz_loc(2,3))
             w_all_faces(k1,k2,iface,3,j1,j2,j3) = w_all_faces(k1,k2,iface,3,j1,j2,j3) &
                  +(l2m*uxyz_loc(3,3)+lam_loc*(uxyz_loc(2,2)+uxyz_loc(1,1)))
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE build_my_w_const_coeff


  SUBROUTINE compute_single_mode_error(l2_err,udg,ifirst,ilast,jfirst,jlast,&
       kfirst,klast,h,t,q_u,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER,  intent(in) :: q_u,nint
    REAL(dp), intent(in) :: h,t
    REAL(dp), intent(inout) :: l2_err(3)
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    integer :: j1,j2,j3,k1,k2,k3,ivar,i1,i2,i3
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp), DIMENSION(0:nint,0:nint,0:nint,3) :: u_elem
    REAL(dp), DIMENSION(0:nint,0:nint,0:nint) :: x_int,y_int,z_int,f_int
    REAL(dp) :: m_c,py,pz

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    w_int = w_int*h/2.0_dp
    l2_err = 0.0_dp

    ! Compute the error in u for each component
    do ivar = 1,3
     l2_err(ivar) = 0.0_dp
     do j3 = kfirst,klast
      do j2 = jfirst,jlast
       do j1 = ifirst,ilast
        ! For each element build u_elem
        u_elem = 0.0_dp
        do i3 = 0,q_u
         do i2 = 0,q_u
          do i1 = 0,q_u
           m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
           do k3 = 0,nint
            pz = m_c*P_int(k3,i3)
            do k2 = 0,nint
             py = pz*P_int(k2,i2)
             do k1 = 0,nint
              u_elem(k1,k2,k3,ivar) = u_elem(k1,k2,k3,ivar) + py*P_int(k1,i1)
             end do
            end do
           end do
          end do
         end do
        end do
        ! Build a local grid
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           x_int(k1,k2,k3) = (real(j1,dp)-0.5_dp)*h + r_int(k1)/2.0_dp*h
           y_int(k1,k2,k3) = (real(j2,dp)-0.5_dp)*h + r_int(k2)/2.0_dp*h
           z_int(k1,k2,k3) = (real(j3,dp)-0.5_dp)*h + r_int(k3)/2.0_dp*h
          end do
         end do
        end do
        ! Compute l2-err
        if (ivar.eq.1) f_int = cos(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        if (ivar.eq.2) f_int = sin(2.0_dp*pi*x_int)*cos(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        if (ivar.eq.3) f_int = sin(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*cos(2.0_dp*pi*z_int)
        f_int = (cos(t*sqrt(12.0_dp*(2.0_dp*pi)**2))*f_int - u_elem(:,:,:,ivar))**2
        !f_int = u_elem(:,:,:,ivar)
        do k3 = 0,nint
         pz = w_int(k3)
         do k2 = 0,nint
          py = pz*w_int(k2)
          do k1 = 0,nint
           l2_err(ivar) = l2_err(ivar) + f_int(k1,k2,k3)*w_int(k1)*py
          end do
         end do
        end do
       end do
      end do
     end do
    end do

  END SUBROUTINE compute_single_mode_error

  SUBROUTINE compute_point_dirac_error(l2_err,udg,ifirst,ilast,jfirst,jlast,&
       kfirst,klast,h,t,q_u,nint,parameters) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    INTEGER,  intent(in) :: q_u,nint
    REAL(dp), intent(in) :: h,t
    REAL(dp), intent(inout) :: l2_err(3),parameters(11)
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    integer :: j1,j2,j3,k1,k2,k3,ivar,i1,i2,i3
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp), DIMENSION(0:nint,0:nint,0:nint) :: u_elem,u_exact,f_int
    REAL(dp) :: m_c,py,pz,eps,rho,beta,alpha,x0,y0,z0,t0,fx,fy,fz,time,a,b,fr,x,y,z,r,&
         max_err,max_err_loc

    eps = 0.01_dp*h
    rho   = parameters(1)
    beta  = parameters(2)
    alpha = parameters(3)
    x0    = parameters(4)
    y0    = parameters(5)
    z0    = parameters(6)
    fr    = parameters(7)
    t0    = parameters(8)
    fx    = parameters(9)
    fy    = parameters(10)
    fz    = parameters(11)
    time  = (t-t0)*fr
    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    w_int = w_int*h/2.0_dp

    max_err = 0.0_dp

    ! Compute the error in u for each component
    l2_err = 0.0_dp
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do ivar = 1,3
        ! For each element build u_elem
        u_elem = 0.0_dp
        do i3 = 0,q_u
         do i2 = 0,q_u
          do i1 = 0,q_u
           m_c = udg(i1,i2,i3,ivar,j1,j2,j3)
           do k3 = 0,nint
            pz = m_c*P_int(k3,i3)
            do k2 = 0,nint
             py = pz*P_int(k2,i2)
             do k1 = 0,nint
              u_elem(k1,k2,k3) = u_elem(k1,k2,k3) + py*P_int(k1,i1)
             end do
            end do
           end do
          end do
         end do
        end do
        ! Build a local grid
        do k3 = 0,nint
         do k2 = 0,nint
          do k1 = 0,nint
           x = (real(j1,dp)-0.5_dp)*h + r_int(k1)/2.0_dp*h
           y = (real(j2,dp)-0.5_dp)*h + r_int(k2)/2.0_dp*h
           z = (real(j3,dp)-0.5_dp)*h + r_int(k3)/2.0_dp*h
           R = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
           if( R .lt. eps ) then
            u_exact(k1,k2,k3) = 0.0_dp
           else
            ! (tD == iC6SmoothBump)
            A = ((1.0_dp/alpha**2)*C6SmoothBump_f(time, fr*R, alpha)-1.0_dp/(beta**2)*C6SmoothBump_f(time,fr*R,beta) &
                 +3.0_dp/((fr*R)**2)*C6SmoothBump_x_T_Integral_f(time,fr*R,alpha,beta))/(4.0_dp*pi*rho*R*R*R)
            B = (1.0_dp/(beta**2)*C6SmoothBump_f(time,fr*R,beta)&
                 -1.0_dp/((fr*R)**2)*C6SmoothBump_x_T_Integral_f(time,fr*R,alpha,beta))/(4.0_dp*pi*rho*R)
            if (ivar.eq.1) u_exact(k1,k2,k3) = ( (x - x0)*(x - x0)*fx + (x - x0)*(y - y0)*fy + (x - x0)*(z - z0)*fz )*A + fx*B
            if (ivar.eq.2) u_exact(k1,k2,k3) = ( (y - y0)*(x - x0)*fx + (y - y0)*(y - y0)*fy + (y - y0)*(z - z0)*fz )*A + fy*B
            if (ivar.eq.3) u_exact(k1,k2,k3) = ( (z - z0)*(x - x0)*fx + (z - z0)*(y - y0)*fy + (z - z0)*(z - z0)*fz )*A + fz*B
           end if
          end do
         end do
        end do
        ! Compute l2-err
        f_int = (u_exact - u_elem)**2
        if( R .lt. 1.99_dp*h) then
         f_int = 0.0_dp
        end if
        if (ivar.eq.3) then
         max_err_loc = sqrt(maxval(f_int))
         if(max_err_loc .gt. max_err) max_err = max_err_loc
        end if
        do k3 = 0,nint
         pz = w_int(k3)
         do k2 = 0,nint
          py = pz*w_int(k2)
          do k1 = 0,nint
           l2_err(ivar) = l2_err(ivar) + f_int(k1,k2,k3)*w_int(k1)*py
          end do
         end do
        end do
       end do
      end do
     end do
    end do
    l2_err(1) = max_err
  END SUBROUTINE compute_point_dirac_error


  SUBROUTINE factor(MU,MV,q_u,q_v) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(inout)  :: q_u,q_v
    REAL(dp), DIMENSION(3*(q_v+1)**3,3*(q_v+1)**3) :: MV
    REAL(dp), DIMENSION(3*(q_u+1)**3,3*(q_u+1)**3) :: MU
    INTEGER :: nrows,ncols
    INTEGER :: INFO,LWORK
    INTEGER, DIMENSION(3*(q_v+1)**3) :: IPIVV
    INTEGER, DIMENSION(3*(q_u+1)**3) :: IPIVU
    REAL(dp), DIMENSION(5*3*(q_u+1)**3) :: UWORK

    nrows = 3*(q_v+1)**3
    ncols = 3*(q_v+1)**3
    CALL DGETRF(nrows,ncols,MV,nrows,IPIVV,INFO)
    if (info .ne. 0) write(*,*) "Warning LU factorization of MV did not work. INFO = ",INFO
    LWORK = 5*3*(q_u+1)**3
    call DGETRI(nrows,MV,nrows,IPIVV,UWORK,LWORK,INFO)
    if (info .ne. 0) write(*,*) "Inversion of MV did not work. INFO = ",INFO
    nrows = 3*(q_u+1)**3
    ncols = 3*(q_u+1)**3
    CALL DGETRF(nrows,ncols,MU,nrows,IPIVU,INFO)
    if (info .ne. 0) write(*,*) "Warning LU factorization of MU did not work. INFO = ",INFO
    LWORK = 5*3*(q_u+1)**3
    call DGETRI(nrows,MU,nrows,IPIVU,UWORK,LWORK,INFO)
    if (info .ne. 0) write(*,*) "Inversion of MU did not work. INFO = ",INFO
  END SUBROUTINE factor

  SUBROUTINE compute_surface_integrals(v_in_all_faces,v_star_all_faces,&
       w_star_all_faces,force_u,force_v,LU,LV,h,q_u,q_v,nint,&
       ifirst,ilast,jfirst,jlast,kfirst,klast) bind(c)
    !
    use type_defs
    IMPLICIT NONE
    !
    INTEGER :: q_u,q_v,nint
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
    real(dp) :: h
    REAL(dp), DIMENSION(0:nint,0:nint,(q_v+1)**3,6) :: LV
    REAL(dp), DIMENSION(0:nint,0:nint,(q_u+1)**3,3,3,6) :: LU
    !
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: force_u
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: force_v
    !
    ! (quad_points,quad_points,sides,fields,elements)
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,w_star_all_faces,v_star_all_faces
    REAL(dp), DIMENSION(0:nint) :: r_int,w_int
    REAL(dp)  :: P_int(0:nint,0:q_u),Pr_int(0:nint,0:q_u)
    REAL(dp)  :: x_r,r_x,jac_loc,int_loc
    integer :: i1,i2,i3
    integer :: j1,j2,j3
    integer :: k1,k2,k3
    integer :: ivar,row,i4

    force_u = 0.0_dp
    force_v = 0.0_dp

    call lglnodes(r_int,w_int,P_int,Pr_int,nint,q_u)
    x_r = h/2.0_dp
    r_x = 1.0_dp/x_r
    !
    ! Compute the forcing in "u_t = v".
    !

    !
    ! Note that the matrices LU and LV can absorb many of the
    ! computaitons below
    !
    jac_loc = x_r**2
    !$OMP PARALLEL DO PRIVATE(j3,j2,j1,i1,i2,i3,ivar,i4,row,int_loc,k1,k2,k3)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_u
        do i2 = 0,q_u
         do i1 = 0,q_u
          row = 1 + i1 + (q_u+1)*i2 + (q_u+1)**2*i3
          ! first set of test functions
          do ivar = 1,3
           int_loc = 0.d0
           do i4 = 1,3 ! Loop over v_i4^*-v_i4
            ! x-faces
            do k3 = 0,nint
             do k2 = 0,nint
              int_loc = int_loc + jac_loc*w_int(k2)*w_int(k3)*(&
                   +LU(k2,k3,row,ivar,i4,2)*(v_star_all_faces(k2,k3,2,i4,j1,j2,j3)-v_in_all_faces(k2,k3,2,i4,j1,j2,j3)) &
                   -LU(k2,k3,row,ivar,i4,1)*(v_star_all_faces(k2,k3,1,i4,j1,j2,j3)-v_in_all_faces(k2,k3,1,i4,j1,j2,j3)))
             end do
            end do
            ! y-faces
            do k3 = 0,nint
             do k1 = 0,nint
              int_loc = int_loc + jac_loc*w_int(k1)*w_int(k3)*(&
                   +LU(k1,k3,row,ivar,i4,4)*(v_star_all_faces(k1,k3,4,i4,j1,j2,j3)-v_in_all_faces(k1,k3,4,i4,j1,j2,j3)) &
                   -LU(k1,k3,row,ivar,i4,3)*(v_star_all_faces(k1,k3,3,i4,j1,j2,j3)-v_in_all_faces(k1,k3,3,i4,j1,j2,j3)))
             end do
            end do
            ! z-faces
            do k2 = 0,nint
             do k1 = 0,nint
              int_loc = int_loc + jac_loc*w_int(k1)*w_int(k2)*(&
                   +LU(k1,k2,row,ivar,i4,6)*(v_star_all_faces(k1,k2,6,i4,j1,j2,j3)-v_in_all_faces(k1,k2,6,i4,j1,j2,j3)) &
                   -LU(k1,k2,row,ivar,i4,5)*(v_star_all_faces(k1,k2,5,i4,j1,j2,j3)-v_in_all_faces(k1,k2,5,i4,j1,j2,j3)))
             end do
            end do
           end do
           force_u(i1,i2,i3,ivar,j1,j2,j3) = int_loc
          end do
         end do
        end do
       end do
      end do
     end do
    end do
    !
    ! Compute the forcing in "v_t = \div(T)".
    !
    !$OMP PARALLEL DO PRIVATE(j3,j2,j1,i1,i2,i3,i4,row,k1,k2,k3)
    do j3 = kfirst,klast
     do j2 = jfirst,jlast
      do j1 = ifirst,ilast
       do i3 = 0,q_v
        do i2 = 0,q_v
         do i1 = 0,q_v
          row = 1 + i1 + (q_v+1)*i2 + (q_v+1)**2*i3
          !
          do ivar = 1,3
           ! x-faces
           do k3 = 0,nint
            do k2 = 0,nint
             force_v(i1,i2,i3,ivar,j1,j2,j3) = force_v(i1,i2,i3,ivar,j1,j2,j3) &
                  +jac_loc*w_int(k2)*w_int(k3)*(&
                  (LV(k2,k3,row,2)*w_star_all_faces(k2,k3,2,ivar,j1,j2,j3) &
                  -LV(k2,k3,row,1)*w_star_all_faces(k2,k3,1,ivar,j1,j2,j3)))
            end do
           end do
           ! y-faces
           do k3 = 0,nint
            do k1 = 0,nint
             force_v(i1,i2,i3,ivar,j1,j2,j3) = force_v(i1,i2,i3,ivar,j1,j2,j3) &
                  +jac_loc*w_int(k1)*w_int(k3)*(&
                  (LV(k1,k3,row,4)*w_star_all_faces(k1,k3,4,ivar,j1,j2,j3) &
                  -LV(k1,k3,row,3)*w_star_all_faces(k1,k3,3,ivar,j1,j2,j3)))
            end do
           end do
           ! z-faces
           do k2 = 0,nint
            do k1 = 0,nint
             force_v(i1,i2,i3,ivar,j1,j2,j3) = force_v(i1,i2,i3,ivar,j1,j2,j3) &
                  +jac_loc*w_int(k1)*w_int(k2)*(&
                  (LV(k1,k2,row,6)*w_star_all_faces(k1,k2,6,ivar,j1,j2,j3) &
                  -LV(k1,k2,row,5)*w_star_all_faces(k1,k2,5,ivar,j1,j2,j3)))
            end do
           end do
          end do
         end do
        end do
       end do
      end do
     end do
    end do
    !

  END SUBROUTINE compute_surface_integrals

  SUBROUTINE EVAL_LEGENDRE(P,Pr,qu,x) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in)  :: qu
    REAL(dp) :: x,P(0:qu),Pr(0:qu)
    INTEGER :: k

    P(0) = 1.0_dp
    P(1) = x
    ! k*P_k(x)=(2*k-1)*x*P_{k-1}(x)-(k-1)*P_{k-2}(x)
    do k=2,qu
     P(k)=(real(2*k-1,dp)*x*P(k-1)-real(k-1,dp)*P(k-2))/real(k,dp)
    end do

    Pr(0) = 0.0_dp
    Pr(1) = 1.0_dp
    ! P_k'(x)=P_{k-2}'(x)+(2k-1)P_{k-1}(x)
    do k=2,qu
     Pr(k) = Pr(k-2)+real(2*k-1,dp)*P(k-1)
    end do
  END SUBROUTINE EVAL_LEGENDRE

  SUBROUTINE get_recorder(udg,vdg,urec,vrec,ifirst,ilast,jfirst,jlast,kfirst,klast,q_v,q_u,&
       Px,Prx,Py,Pry,Pz,Prz,ix,iy,iz,ncomponents) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast,ix,iy,iz,ncomponents
    REAL(dp), DIMENSION(0:q_u,0:q_u,0:q_u,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: udg
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: vdg
    REAL(dp) :: Px(0:q_u),Prx(0:q_u),Py(0:q_u),Pry(0:q_u),Pz(0:q_u),Prz(0:q_u),&
         urec(ncomponents),vrec(ncomponents)
    INTEGER :: k1,k2,k3,ivar
    real(dp) :: tz,ty

    vrec = 0.0_dp
    urec = 0.0_dp
    do ivar = 1,3
     do k3 = 0,q_u
      do k2 = 0,q_u
       do k1 = 0,q_u
        urec(ivar) = urec(ivar) + Px(k1)*Py(k2)*Pz(k3)*udg(k1,k2,k3,ivar,ix,iy,iz)
       end do
      end do
     end do
    end do
  END SUBROUTINE get_recorder

  SUBROUTINE get_dirac_source(point_src,px,py,pz,q_v) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_v
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v) :: point_src
    REAL(dp), DIMENSION(0:q_v) ::px,py,pz
    INTEGER :: i1,i2,i3
    do i3 = 0,q_v
     do i2 = 0,q_v
      do i1 = 0,q_v
       point_src(i1,i2,i3) = px(i1)*py(i2)*pz(i3)
      end do
     end do
    end do

  end SUBROUTINE get_dirac_source

  SUBROUTINE add_dirac_source(force_v,point_src,f_amp,source_tay_coeff,&
       ifirst,ilast,jfirst,jlast,kfirst,klast,q_v,&
       ix,iy,iz) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_v
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast,ix,iy,iz
    REAL(dp), DIMENSION(0:q_v,0:q_v,0:q_v,3,ifirst:ilast,jfirst:jlast,kfirst:klast) :: force_v
    REAL(dp) :: point_src(0:q_v,0:q_v,0:q_v),f_amp(3),source_tay_coeff
    INTEGER :: ivar

    do ivar = 1,3
     force_v(0:q_v,0:q_v,0:q_v,ivar,ix,iy,iz) = force_v(0:q_v,0:q_v,0:q_v,ivar,ix,iy,iz) &
          + point_src*source_tay_coeff*f_amp(ivar)
    end do
  END SUBROUTINE add_dirac_source

  SUBROUTINE set_tay_weights(tg,ct,nsrc,ntay) bind(c)
    use type_defs
    IMPLICIT NONE
    INTEGER :: ntay,nsrc,i
    real(KIND = dp) :: tg(0:nsrc),ct(0:nsrc,0:ntay)

    ! Cheb grid
    do i = 0,nsrc
     tg(i) = -cos(pi*real(i,dp)/real(nsrc,dp))
    end do
    tg = 0.5_dp*tg
    call weights_qp(0.d0,tg,nsrc,nsrc,ntay,ct)

  END SUBROUTINE set_tay_weights

  SUBROUTINE get_source_tay_coeff(source_tay_coeff,tg,ct,nsrc,ntay,t,dt) bind(c)
    use type_defs
    IMPLICIT NONE
    INTEGER :: ntay,nsrc,i
    real(KIND = dp) :: source_tay_coeff(0:ntay),fsrc(0:nsrc),tg(0:nsrc),ct(0:nsrc,0:ntay),dt,df,t

    ! Evaluate the source on a Chebyshev grid around the current time
    do i = 0,nsrc
     fsrc(i) = source_fun(t+dt*tg(i))
    end do
    df = 1.0_dp
    do i = 0,ntay
     source_tay_coeff(i) = sum(ct(:,i)*fsrc)*df
     df = df/dt
    end do

  END SUBROUTINE get_source_tay_coeff

  !// C6 smooth bump for time dependence for further testing of point force
  real(dp) function C6SmoothBump_f(t,R,c)
    use type_defs
    IMPLICIT NONE
    real(dp) :: t,R,c
    if(((t-R/c) .gt. 0.0_dp).and.( (t-R/c) .lt. 1.0_dp)) then
     C6SmoothBump_f = 51480.0_dp*((t-R/c)*(1.0_dp-t+R/c))**7
    else
     C6SmoothBump_f = 0.0_dp
    end if
  end function C6SmoothBump_f

  real(kind = dp) function C6SmoothBump_x_T_Integral_f(t,R,alpha,beta)
    ! // Integral of H(t-T)*H(1-t+T)*C6SmoothBump(t-T)*T from R/alpha to R/beta
    use type_defs
    IMPLICIT NONE
    real(dp) :: t,r,alpha,beta,hil,lowl,temp

    temp = R
    !//  lowL = where(R / alpha > t - 1, R/alpha, t - 1); hiL = where(R / beta < t, R / beta, t);
    if ( R/alpha .gt. (t-1.0_dp)) then
     lowL = R/alpha
    else
     lowL = t-1.0_dp
    end if
    if( R/beta .lt. t ) then
     hiL = R/beta
    else
     hiL = t
    end if
    !//  temp = where (lowL < t && hiL > t - 1, VSBTP(hiL, t) - VSBTP(lowL, t), 0.0);
    if( (lowL .lt. t).and.(hiL > (t-1.0_dp))) then
     temp = C6SBTP(hiL, t) - C6SBTP(lowL, t)
    else
     temp = 0.0_dp
    end if
    C6SmoothBump_x_T_Integral_f=temp
  end function C6SmoothBump_x_T_Integral_f

  real(kind = dp) function C6SBTP(Lim,t)
    ! // Primitive function (for T) of C6SmoothBump(t-T)*T
    use type_defs
    IMPLICIT NONE
    real(dp) :: t,lim,x
    x = t-Lim
    C6SBTP =  (x**8)*(-3217.5_dp*(x**8)+3432.0_dp*(7.0_dp+t)*(x**7)-25740.0_dp*(3.0_dp+t)*(x**6)&
         +27720.0_dp*(5.0_dp+3.0_dp*t)*(x**5)-150150.0_dp*(t+1.0_dp)*x*x*x*x &
         +32760.0_dp*(3.0_dp+5.0_dp*t)*x*x*x-36036.0_dp*(1.0_dp+3.0_dp*t)*x*x+5720.0_dp*(1.0_dp+7.0_dp*t)*x-6435.0_dp*t)
  end function C6SBTP

  real(kind = dp) function source_fun(t)

    use type_defs
    IMPLICIT NONE
    real(kind = dp) :: t
    real(kind = dp) :: sig2,t0
    real(kind = dp) :: omega
    t0 = 0.0_dp
    omega = 1.0_dp
    ! sig2 = 36.0_dp/(t0**2)
    ! source_fun = -2.0_dp*sig2*(t-t0)*exp(-sig2*(t-t0)**2)
    if (t.lt.t0) then
     source_fun = 0.0_dp
    elseif ((t0.le.t).and.(t.lt.(t0+1.0_dp/omega))) then
     source_fun = 51480.0_dp*omega**7*(t-t0)**7*(1.0_dp-omega*(t-t0))**7
    else
     source_fun = 0.0_dp
    end if
    source_fun = source_fun

  end function source_fun

  subroutine weights_qp(z,x,n,nd,m,c)
    !c---  This routine is from Bengt Fornbergs paper in SIAM Review
    !c---  Input Parameters
    !c---  z location where approximations are to be accurate,
    !c---  x(0:nd) grid point locations, found in x(0:n)
    !c---  n one less than total number of grid points; n must
    !c---  not exceed the parameter nd below,
    !c---  nd dimension of x- and c-arrays in calling program
    !c---  x(0:nd) and c(0:nd,0:m), respectively,
    !c---  m highest derivative for which weights are sought,

    !c---  Output Parameter
    !c---  c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
    !c---  of order 0:m, found in c(0:n,0:m)
    !c---
    use type_defs
    IMPLICIT NONE
    INTEGER :: nd,m,n,i,j,k,mn
    real(KIND = dp) :: x(0:nd),c(0:nd,0:m),c1,c2,c3,c4,c5,z

    c1 = 1.0_dp
    c4 = x(0)-z
    do k = 0,m
     do j = 0,n
      c(j,k) = 0.0_dp
     end do
    end do
    c(0,0) = 1.0_dp
    do i=1,n
     mn = min(i,m)
     c2 = 1.0_dp
     c5 = c4
     c4 = x(i)-z
     do j=0,i-1
      c3 = x(i)-x(j)
      c2 = c2*c3
      if (j.eq.i-1) then
       do  k=mn,1,-1
        c(i,k) = c1*(real(k,dp)*c(i-1,k-1)-c5*c(i-1,k))/c2
       end do
       c(i,0) = -c1*c5*c(i-1,0)/c2
      endif
      do k=mn,1,-1
       c(j,k) = (c4*c(j,k)-real(k,dp)*c(j,k-1))/c3
      end do
      c(j,0) = c4*c(j,0)/c3
     end do
     c1 = c2
    end do

  end subroutine weights_qp

end module dgmodule
