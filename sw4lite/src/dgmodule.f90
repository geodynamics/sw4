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
    end do

    do i5 = kfirst,klast
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
    end do

    do i5 = kfirst,klast
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
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,v_out_all_faces, &
         w_in_all_faces,w_out_all_faces, &
         v_star_all_faces,w_star_all_faces
    REAL(dp), parameter :: alpha_flip = 0.5_dp
    REAL(dp), parameter :: beta_up = 0.5_dp
    REAL(dp), parameter :: tau_up  = 0.1_dp

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
         v_in_all_faces,v_out_all_faces, &
         w_in_all_faces,w_out_all_faces

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
       ifirst,ilast,jfirst,jlast,kfirst,klast,nint,sbx_b,sbx_e,sby_b,sby_e) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: nint,ifirst,ilast,jfirst,jlast,kfirst,klast,sbx_b,sbx_e,sby_b,sby_e
    REAL(dp), DIMENSION(0:nint,0:nint,6,3,ifirst:ilast,jfirst:jlast,kfirst:klast) ::&
         v_in_all_faces,v_out_all_faces, &
         w_in_all_faces,w_out_all_faces

    if (sbx_b.eq.1) then
     v_out_all_faces(:,:,1,1,ifirst,:,:)     =  v_in_all_faces(:,:,1,1,ifirst,:,:)
     v_out_all_faces(:,:,1,2,ifirst,:,:)     = -v_in_all_faces(:,:,1,2,ifirst,:,:)
     v_out_all_faces(:,:,1,3,ifirst,:,:)     = -v_in_all_faces(:,:,1,3,ifirst,:,:)
     w_out_all_faces(:,:,1,1,ifirst,:,:)     = -w_in_all_faces(:,:,1,1,ifirst,:,:)
     w_out_all_faces(:,:,1,2,ifirst,:,:)     =  w_in_all_faces(:,:,1,2,ifirst,:,:)
     w_out_all_faces(:,:,1,3,ifirst,:,:)     =  w_in_all_faces(:,:,1,3,ifirst,:,:)
    end if
    if (sbx_e.eq.1) then
     v_out_all_faces(:,:,2,1,ilast,:,:)     =  v_in_all_faces(:,:,2,1,ilast,:,:)
     v_out_all_faces(:,:,2,2,ilast,:,:)     = -v_in_all_faces(:,:,2,2,ilast,:,:)
     v_out_all_faces(:,:,2,3,ilast,:,:)     = -v_in_all_faces(:,:,2,3,ilast,:,:)
     w_out_all_faces(:,:,2,1,ilast,:,:)     = -w_in_all_faces(:,:,2,1,ilast,:,:)
     w_out_all_faces(:,:,2,2,ilast,:,:)     =  w_in_all_faces(:,:,2,2,ilast,:,:)
     w_out_all_faces(:,:,2,3,ilast,:,:)     =  w_in_all_faces(:,:,2,3,ilast,:,:)
    end if

    if (sby_b.eq.1) then
     v_out_all_faces(:,:,3,1,:,jfirst,:)     = -v_in_all_faces(:,:,3,1,:,jfirst,:)
     v_out_all_faces(:,:,3,2,:,jfirst,:)     =  v_in_all_faces(:,:,3,2,:,jfirst,:)
     v_out_all_faces(:,:,3,3,:,jfirst,:)     = -v_in_all_faces(:,:,3,3,:,jfirst,:)
     w_out_all_faces(:,:,3,1,:,jfirst,:)     =  w_in_all_faces(:,:,3,1,:,jfirst,:)
     w_out_all_faces(:,:,3,2,:,jfirst,:)     = -w_in_all_faces(:,:,3,2,:,jfirst,:)
     w_out_all_faces(:,:,3,3,:,jfirst,:)     =  w_in_all_faces(:,:,3,3,:,jfirst,:)
    end if
    if (sby_e.eq.1) then
     v_out_all_faces(:,:,4,1,:,jlast,:)     = -v_in_all_faces(:,:,4,1,:,jlast,:)
     v_out_all_faces(:,:,4,2,:,jlast,:)     =  v_in_all_faces(:,:,4,2,:,jlast,:)
     v_out_all_faces(:,:,4,3,:,jlast,:)     = -v_in_all_faces(:,:,4,3,:,jlast,:)
     w_out_all_faces(:,:,4,1,:,jlast,:)     =  w_in_all_faces(:,:,4,1,:,jlast,:)
     w_out_all_faces(:,:,4,2,:,jlast,:)     = -w_in_all_faces(:,:,4,2,:,jlast,:)
     w_out_all_faces(:,:,4,3,:,jlast,:)     =  w_in_all_faces(:,:,4,3,:,jlast,:)
    end if
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
       h,q_v,q_u,nint) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in) :: q_u,q_v,nint
    INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast,kfirst,klast
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
        !if (ivar.eq.1) f_int = 1.0_dp
        !if (ivar.eq.2) f_int = 1.0_dp
        !if (ivar.eq.3) f_int = 1.0_dp
        if (ivar.eq.1) f_int = cos(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        if (ivar.eq.2) f_int = sin(2.0_dp*pi*x_int)*cos(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        if (ivar.eq.3) f_int = sin(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*cos(2.0_dp*pi*z_int)
        !f_int = x_int*y_int*z_int
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
        if (ivar.eq.1) f_int = 0.0_dp
        if (ivar.eq.2) f_int = 0.0_dp
        if (ivar.eq.3) f_int = 0.0_dp
        !if (ivar.eq.1) f_int = cos(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        !if (ivar.eq.2) f_int = sin(2.0_dp*pi*x_int)*cos(2.0_dp*pi*y_int)*sin(2.0_dp*pi*z_int)
        !if (ivar.eq.3) f_int = sin(2.0_dp*pi*x_int)*sin(2.0_dp*pi*y_int)*cos(2.0_dp*pi*z_int)
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

    ! NOW MV
    !$OMP PARALLEL DO PRIVATE(i3,i2,i1,row,j3,j2,j1,col,k1,k2,k3)
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
     !$OMP PARALLEL DO PRIVATE(j3,j2,j1,u_elem,i3,i2,i1,m_c,pz,py,x_int,y_int,z_int,f_int) reduction (+:l2_err)
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

  SUBROUTINE factor(MU,MV,q_u,q_v) bind(c)
    !
    USE type_defs
    IMPLICIT NONE
    INTEGER, intent(in)  :: q_u,q_v
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


end module dgmodule
