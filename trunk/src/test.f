      implicit none

      real*8 gettimec

      integer nx, ny, i, j, c, acc, seed, rn, bccnd(4), ini
      integer ifirst, ilast, jfirst, jlast, ni, nj
      integer maxsteps, st, nsteps, tacc, errfreq
      integer, dimension(:), allocatable :: seeds
      real*8 Lx, Ly, h, x, y, dt, err(2), li(2), l2(2), uex(2), t1
      real*8 l2size(2), lisize(2)
      real*8, allocatable, dimension(:,:,:) :: up, u, um, fo, um2, u2
      real*8, allocatable, dimension(:,:) :: mu, la, rho, fobnd
      real*8 acof(6,8,8), bop(4,6), ghcof, hnorm(4)
      real*8 dt2i, pi, dtpr
      real*8 locrad, maxrad, om, omm, ph, phm, th, la0, cfl
      real*8 energy, normfact, term, t, tfinal, cpcvrat, lamurat
      real*8 iop(5), iop2(5), bop2(4,6), gh2, s(0:4), RVEL
      real*8 cRel, Tperiod
      real*8 k1, k2, rp, rs, alpha, tmp1, tmp2, ifact, ang, muval, rpeff
      character*400 buf
      logical writesol, cpcsratset, cflset, lamuratset, onlyswave

      pi  = 4*ATAN(1d0)
c      L   = 5d0
      Ly   = 1
      Lx   = 1.5*Ly
      ny   = 51
      dt   = 1
      i    = 1
      acc  = 4
      seed = 298347
      bccnd(1) = 0
      bccnd(2) = 0
      bccnd(3) = 8
      bccnd(4) = 8
*** ini=1 --> Twilight test
***     2 --> Energy test
***     3 --> Surface wave
***     4 --> P-wave reflection
      ang = 45
      ini = 1
      cfl = 0.8d0
      cflset = .false.
      tfinal = 1.0
      maxsteps = 1000000
      cpcsratset = .false.
      lamuratset = .false.
      errfreq  = 10
      writesol = .true.
      onlyswave = .true.
      do while( i .le. IARGC() )
         call GETARG( i, buf )
         if( buf .eq. '-n' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') ny
            i = i+2
         elseif( buf.eq.'-acc' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') acc
            i = i+2
         elseif( buf.eq.'-maxsteps' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') maxsteps
            i = i+2
         elseif( buf.eq.'-tfinal' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') tfinal
            i = i+2
         elseif( buf.eq.'-cfl' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') cfl
            cflset = .true.
            i = i+2
         elseif( buf.eq.'-seed' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') seed
            i = i+2
         elseif( buf.eq.'-errfreq' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') errfreq
            i = i+2
         elseif( buf.eq.'-testno' )then
            call GETARG(i+1,buf)
            read(buf,'(i16)') ini
            i = i+2
         elseif( buf.eq.'-cpcsrat' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') cpcvrat
            cpcsratset = .true.
            i = i+2
         elseif( buf.eq.'-lamurat' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') lamurat
            lamuratset = .true.
            i = i+2
         elseif( buf.eq.'-ang' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') ang
            i = i+2
         elseif( buf.eq.'-reflpwave' )then
            onlyswave = .false.
            i = i+1
         elseif( buf.eq.'-h' )then
            write(*,*)'Example usage: ew -testno 2 -n 100 -acc 2 '//
     *       '-maxsteps 1000 -tfinal 2.5 -cfl 0.8 -seed 23984 '//
     *       ' -errfreq 10 -cpcsrat 30.4'
            write(*,*)
            write(*,*) '  -testno [1|2|3|4]   Problem number: '//
     *       ' 1-twilight test, 2-energy test, 3-surface wave,'//
     +       ' 4-reflected P-wave'
            write(*,*) '  -n 100   Number of grid points in y-dir'
            write(*,*) '  -acc [2|4]   Order of accuracy'
            write(*,*) '  -maxsteps 1000   Maximum no of time steps'
            write(*,*) '  -tfinal 2.5   Integrate to this time'
            write(*,*) '  -cfl 0.8   CFL-number'
            write(*,*) '  -seed 23984   Seed to random no generator'//
     *      ' for energy test'
            write(*,*) '  -errfreq 10  Print error every 10th step'//
     *      ' for surface waves'
            write(*,*) '  -cpcsrat 30.4  Set cp/cv ratio to 30.4'
            write(*,*) '  -lamurat 4.30  Set lambda/mu ratio to 4.30'
            write(*,*) '  -reflpwave Compute both the reflected P-'//
     +      ' and S-waves (only for -testno 4)'
            stop
         else
            write(*,*) 'WARNING: Option ',buf(1:20), ' not recognized'
            i = i+ 1
         endif
      enddo
      if( acc.ne.4 .and. acc.ne.2 )then
         write(*,*) 'acc must be two or four, not ',acc
         stop
      endif
      if( ini.eq.1 )then
         write(*,*)'Running twilight test'
      elseif( ini.eq.2 )then
         write(*,*) 'Running energy test NOT IMPLEMENTED'
         stop
      elseif( ini.eq.3 )then
         write(*,*) 'Running surface wave test NOT IMPLEMENTED'
         stop
      elseif( ini.eq.4 )then
         write(*,*) 'Running P-wave reflection test NOT IMPLEMENTED'
         stop
      else
         write(*,*) 'testno must be one, two, three, or four, not ',ini
         stop
      endif

      if( acc.eq.4 )then
         tacc = 4
      else
         tacc = 2
      endif

      if( ini.eq.1 )then
         Lx = 2*pi
         Ly = 2*pi
      elseif( ini.eq.4 )then
         ang = ang*pi/180
         Ly = 4*pi/sin(ang)
         Lx = 4*pi/cos(ang)
      endif

      if( .not.cpcsratset )then
        if( ini.eq.1 )then
          cpcvrat = 3
        elseif( ini.eq.2 )then
          cpcvrat = 100
        endif
      endif
      
      if(ini.eq.3 .or. ini.eq.4 )then
        if (.not. lamuratset) then
          lamurat = 100.0
        endif
        write(*,206) lamurat
 206    format(' Ratio lambda/mu =', f7.3)

      endif

      h  = Ly/(ny-1)
      nx = Lx/h + 1.5
      write(*,105) 'h = ', h
      write(*,105) 'Domain           = ', Lx, ' x ', Ly
 105  format(' ',a,g12.5,a,g12.5)
      write(*,104) 'Number of points = ', nx, ' x ', ny 
 104  format(' ',a,i5,a,i5)
      call random_seed(size=rn)
      allocate( seeds(rn) )
      do i=1,rn
         seeds(i) = seed/i+seed*10*sin(1d0*i*i)
      enddo
      call random_seed(put=seeds)

c starting and ending index (change to ilast=nx+2, etc)
      ifirst = -1
      ilast = nx+2
      jfirst = -1
      jlast = ny+2

c total number of grid points (ni=nx+4, nj=ny+4)
      ni = ilast-ifirst+1
      nj = jlast-jfirst+1

c allocate solution and material coefficient arrays
      allocate( up(2, ifirst:ilast, jfirst:jlast) )
      allocate( u(2,  ifirst:ilast, jfirst:jlast) )
      allocate( um(2, ifirst:ilast, jfirst:jlast) )
      allocate( rho(  ifirst:ilast, jfirst:jlast) )
      allocate( la(   ifirst:ilast, jfirst:jlast) )
      allocate( fo(2, ifirst:ilast, jfirst:jlast) )
      allocate( u2(2, ifirst:ilast, jfirst:jlast) )
      allocate( mu(   ifirst:ilast, jfirst:jlast) )
      allocate( um2(2, ifirst:ilast, jfirst:jlast) )
      allocate( fobnd(2, jfirst:jlast) )
c get coefficients for difference approximation of 2nd derivative with variable coefficients
      call VARCOEFFS4( acof, ghcof )
c get coefficients for difference approximation of 1st derivative
      call WAVEPROPBOP_4( iop, iop2, bop, bop2, gh2, hnorm, s )
      if( acc.eq.2 )then
         hnorm(1) = 0.5d0
         hnorm(2) = 1
         hnorm(3) = 1
         hnorm(4) = 1
      endif

c setup material and initial conditions
      if( ini.eq.1 )then
*** Twilight testing
*** bccnd=1 is homogeneous Dirichlet
*** bccnd=2 is Dirichlet with tw-forcing
*** bccnd=3 is Free surface (with tw-forcing)
         bccnd(1) = 3
c         bccnd(1) = 2
c         bccnd(2) = 3
         bccnd(2) = 2
         cpcvrat = 3
         ph = 0.3d0
         th = 0.2d0
         phm= 0.1d0
         om  = 1d0
         omm = 1d0
         la0 = 3*(cpcvrat**2-2)
         t   = 0
         st  = 0

         do j=jfirst,jlast
            y = (j-1)*h
            do i=ifirst,ilast
               x = (i-1)*h
               u(1,i,j) = cos(om*x+ph)*sin(om*y+th)*cos(t*t)
               u(2,i,j) = sin(om*x+th)*cos(om*y+ph)*sin(t)
               mu(i,j)  = sin(3*omm*x+phm)*sin(y*omm)+3
               la(i,j)  = cos(omm*x+phm)*sin(omm*3*y)**2+la0
               rho(i,j) = sin(om*x+ph)*sin(om*y-th)+2
            enddo
         enddo
      endif

      if (.not.cflset) then
        if (acc.eq.4) then
          cfl=1.3
        else
          cfl=0.9
        endif
      endif

      if( ini.eq.4 )then
         write(*,105) ' Points per wavelength = ',2*pi/h/sqrt(2+lamurat)
         write(*,105) ' Period in time = ',2*pi/om
         open(unit=11, file='errors.dat', status='unknown')
      endif

      write(*,202)  acc, cfl, errfreq
 202  format(' Order of accuracy =', i2, /, ' Cfl =', f7.3, /, 
     +     ' Err freq =', i4)
*** Determine time step
      maxrad = 0
      do j=1,ny
         do i=1,nx
            locrad = (3*mu(i,j)+la(i,j))/rho(i,j)
            if( locrad.gt.maxrad )then
               maxrad = locrad
            endif
         enddo
      enddo
      dt = cfl*h/SQRT(maxrad)
      nsteps = tfinal/dt
      dt = tfinal/nsteps

      dt2i = 1/(dt*dt)
*** Can not get twilight data at t=-dt until dt is computed.
      if( ini.eq.1 )then
         do j=jfirst,jlast
            y = (j-1)*h
            do i=ifirst,ilast
               x = (i-1)*h
               um(1,i,j) = cos(om*x+ph)*sin(om*y+th)*cos((t-dt)*(t-dt))
               um(2,i,j) = sin(om*x+th)*cos(om*y+ph)*sin(t-dt)
            enddo
         enddo
      endif

c test accuracy of spatial approximation

      call LC4TH( ni, nj, u2, u, mu, la, rho, dt, h )
      call exactRhs( ifirst, ilast, jfirst, jlast, t, fo, om, omm, 
     +     ph, phm, th, la0, h )

*** Compute errors at the final time
      if( ini.eq.1 )then
         do c=1,2
            li(c) = 0
            l2(c) = 0
         enddo
c this test includes all points
         do j=jfirst,jlast
            y = (j-1)*h
            do i=ifirst,ilast
               x = (i-1)*h
               uex(1)  = cos(om*x+ph)*sin(om*y+th)*cos(t*t)
               uex(2)  = sin(om*x+th)*cos(om*y+ph)*sin(t)
               do c=1,2
                  err(c) = ABS( u(c,i,j) - uex(c) )
               enddo
               if( li(1).lt.max(err(1),err(2)) )then
                 li(1) = max(err(1),err(2))
               endif
               l2(1) = l2(1) + h*h* (err(1)**2 + err(2)**2)
               um(1,i,j) = uex(1)-u(1,i,j)
               um(2,i,j) = uex(2)-u(2,i,j)
            enddo
         enddo
         write(*,*) 'Twilight errors max-norm and L2'
         write(*,101) li(1), SQRT(l2(1)/(Lx*Ly))
 101     format(' ', 2(g15.7,tr2))
c error in rhs
         do c=1,2
            li(c) = 0
            l2(c) = 0
         enddo
c this test only includes interior points
         do j=1,ny
            do i=1,nx
               do c=1,2
                  err(c) = ABS( fo(c,i,j) - u2(c,i,j) )
                  if( li(c).lt.err(c) )then
                    li(c) = err(c)
                  endif
               enddo
            end do
          end do
          write(*,*) 'Max errors in RHS: ', li(1), li(2)
      endif

*** Save errors and the solution (only interior points)
      if( ini.eq.1 .or. ini.eq.3 )then
         open(21,file='x.dat')
         open(22,file='y.dat')
         open(23,file='lue.dat')
         open(24,file='lve.dat')
         open(25,file='lu.dat')
         open(26,file='lv.dat')
         do j=1,ny
            write(21,*) ((i-1)*h,i=1,nx)
            write(22,*) ((j-1)*h,i=1,nx)
            write(23,*) (fo(1,i,j),i=1,nx)
            write(24,*) (fo(2,i,j),i=1,nx)
            write(25,*) (u2(1,i,j),i=1,nx)
            write(26,*) (u2(2,i,j),i=1,nx)
         enddo
         close(21)
         close(22)
         close(23)
         close(24)
         close(25)
         close(26)
      endif 

      if( writesol )then
         open(23,file='u.dat')
         open(24,file='v.dat')
         do j=1,ny
            write(23,*) (u(1,i,j),i=1,nx)
            write(24,*) (u(2,i,j),i=1,nx)
         enddo
         close(23)
         close(24)
      endif

      write(*,103) 'Taken ', st, ' steps to time ',t
 103  format(' ', a, i7, a, g12.5 )
      write(*,102) 'Loop execution time ', t1
 102  format(' ', a, g15.7)
      
      end

c-----------------------------------------------------------------------
      subroutine LC4TH( ni, nj, uacc, u, mu, la, rho, dt, h )

*** Only centered approximation of the right hand side of the elastic wave equation

      implicit none

      real*8 tf, i6, i144
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144 )

      integer ni, nj, i, j
      real*8 uacc(2,ni,nj), u(2,ni,nj)
      real*8 mu(ni,nj), la(ni,nj), rho(ni,nj)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 r1, r2, dt, h, dt2, cof

      dt2 = dt*dt
      cof = 1d0/(h*h)
c assume two ghost points
      do j=3,nj-2
         do i=3,ni-2
*** 29 + 46+54 +46+54 + 12 = 241 ops
            mux1 = mu(i-1,j)-tf*(mu(i,j)+mu(i-2,j))
            mux2 = mu(i-2,j)+mu(i+1,j)+3*(mu(i,j)+mu(i-1,j))
            mux3 = mu(i-1,j)+mu(i+2,j)+3*(mu(i+1,j)+mu(i,j))
            mux4 = mu(i+1,j)-tf*(mu(i,j)+mu(i+2,j))
            muy1 = mu(i,j-1)-tf*(mu(i,j)+mu(i,j-2))
            muy2 = mu(i,j-2)+mu(i,j+1)+3*(mu(i,j)+mu(i,j-1))
            muy3 = mu(i,j-1)+mu(i,j+2)+3*(mu(i,j+1)+mu(i,j))
            muy4 = mu(i,j+1)-tf*(mu(i,j)+mu(i,j+2))

*** xx and yy derivatives:
            r1 = i6*((2*mux1+la(i-1,j)-tf*(la(i,j)+la(i-2,j)))*
     *                         (u(1,i-2,j)-u(1,i,j))+
     *           (2*mux2+la(i-2,j)+la(i+1,j)+3*(la(i,j)+la(i-1,j)))*
     *                         (u(1,i-1,j)-u(1,i,j))+ 
     *           (2*mux3+la(i-1,j)+la(i+2,j)+3*(la(i+1,j)+la(i,j)))*
     *                         (u(1,i+1,j)-u(1,i,j))+
     *           (2*mux4+ la(i+1,j)-tf*(la(i,j)+la(i+2,j)))*
     *           (u(1,i+2,j)-u(1,i,j)) 
     *              + muy1*(u(1,i,j-2)-u(1,i,j)) + 
     *                muy2*(u(1,i,j-1)-u(1,i,j)) + 
     *                muy3*(u(1,i,j+1)-u(1,i,j)) +
     *                muy4*(u(1,i,j+2)-u(1,i,j)) )

            r2 = i6*(mux1*(u(2,i-2,j)-u(2,i,j)) + 
     *                 mux2*(u(2,i-1,j)-u(2,i,j)) + 
     *                 mux3*(u(2,i+1,j)-u(2,i,j)) +
     *                 mux4*(u(2,i+2,j)-u(2,i,j)) + 
     *             (2*muy1+la(i,j-1)-tf*(la(i,j)+la(i,j-2)))*
     *                     (u(2,i,j-2)-u(2,i,j))+
     *             (2*muy2+la(i,j-2)+la(i,j+1)+3*(la(i,j)+la(i,j-1)))*
     *                     (u(2,i,j-1)-u(2,i,j))+ 
     *             (2*muy3+la(i,j-1)+la(i,j+2)+3*(la(i,j+1)+la(i,j)))*
     *                     (u(2,i,j+1)-u(2,i,j))+
     *             (2*muy4+la(i,j+1)-tf*(la(i,j)+la(i,j+2)))*
     *                     (u(2,i,j+2)-u(2,i,j)))

*** Mixed derivatives:
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j)*(u(2,i-2,j-2)-u(2,i-2,j+2)+
     *                        8*(-u(2,i-2,j-1)+u(2,i-2,j+1))) - 8*(
     *                   la(i-1,j)*(u(2,i-1,j-2)-u(2,i-1,j+2)+
     *                        8*(-u(2,i-1,j-1)+u(2,i-1,j+1))) )+8*(
     *                   la(i+1,j)*(u(2,i+1,j-2)-u(2,i+1,j+2)+
     *                        8*(-u(2,i+1,j-1)+u(2,i+1,j+1))) ) - (
     *                   la(i+2,j)*(u(2,i+2,j-2)-u(2,i+2,j+2)+
     *                        8*(-u(2,i+2,j-1)+u(2,i+2,j+1))) ))
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j)*(u(1,i-2,j-2)-u(1,i-2,j+2)+
     *                        8*(-u(1,i-2,j-1)+u(1,i-2,j+1))) - 8*(
     *                   mu(i-1,j)*(u(1,i-1,j-2)-u(1,i-1,j+2)+
     *                        8*(-u(1,i-1,j-1)+u(1,i-1,j+1))) )+8*(
     *                   mu(i+1,j)*(u(1,i+1,j-2)-u(1,i+1,j+2)+
     *                        8*(-u(1,i+1,j-1)+u(1,i+1,j+1))) ) - (
     *                   mu(i+2,j)*(u(1,i+2,j-2)-u(1,i+2,j+2)+
     *                        8*(-u(1,i+2,j-1)+u(1,i+2,j+1))) )) 
*** (mu*v_x)_y
            r1 = r1 + i144*( mu(i,j-2)*(u(2,i-2,j-2)-u(2,i+2,j-2)+
     *                        8*(-u(2,i-1,j-2)+u(2,i+1,j-2))) - 8*(
     *                   mu(i,j-1)*(u(2,i-2,j-1)-u(2,i+2,j-1)+
     *                        8*(-u(2,i-1,j-1)+u(2,i+1,j-1))) )+8*(
     *                   mu(i,j+1)*(u(2,i-2,j+1)-u(2,i+2,j+1)+
     *                        8*(-u(2,i-1,j+1)+u(2,i+1,j+1))) ) - (
     *                   mu(i,j+2)*(u(2,i-2,j+2)-u(2,i+2,j+2)+
     *                        8*(-u(2,i-1,j+2)+u(2,i+1,j+2))) )) 
*** (la*u_x)_y
            r2 = r2 + i144*( la(i,j-2)*(u(1,i-2,j-2)-u(1,i+2,j-2)+
     *                        8*(-u(1,i-1,j-2)+u(1,i+1,j-2))) - 8*(
     *                   la(i,j-1)*(u(1,i-2,j-1)-u(1,i+2,j-1)+
     *                        8*(-u(1,i-1,j-1)+u(1,i+1,j-1))) )+8*(
     *                   la(i,j+1)*(u(1,i-2,j+1)-u(1,i+2,j+1)+
     *                        8*(-u(1,i-1,j+1)+u(1,i+1,j+1))) ) - (
     *                   la(i,j+2)*(u(1,i-2,j+2)-u(1,i+2,j+2)+
     *                        8*(-u(1,i-1,j+2)+u(1,i+1,j+2))) ))
            uacc(1,i,j) = cof*r1
            uacc(2,i,j) = cof*r2
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine INNER_LOOP4TH( ni, nj, up, u, um, mu, la, rho, dt, h, 
     *                          fo )

*** Only centered approximation 

      implicit none

      real*8 tf, i6, i144
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144 )

      integer ni, nj, i, j
      real*8 up(2,ni,nj), u(2,ni,nj), um(2,ni,nj)
      real*8 mu(ni,nj), la(ni,nj), fo(2,ni,nj), rho(ni,nj)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 r1, r2, dt, la2, h, dt2, cof

      dt2 = dt*dt
      la2 = dt2/(h*h)
      do j=3,nj-2
         do i=3,ni-2
*** 29 + 46+54 +46+54 + 12 = 241 ops
            cof = la2/rho(i,j)
            mux1 = mu(i-1,j)-tf*(mu(i,j)+mu(i-2,j))
            mux2 = mu(i-2,j)+mu(i+1,j)+3*(mu(i,j)+mu(i-1,j))
            mux3 = mu(i-1,j)+mu(i+2,j)+3*(mu(i+1,j)+mu(i,j))
            mux4 = mu(i+1,j)-tf*(mu(i,j)+mu(i+2,j))
            muy1 = mu(i,j-1)-tf*(mu(i,j)+mu(i,j-2))
            muy2 = mu(i,j-2)+mu(i,j+1)+3*(mu(i,j)+mu(i,j-1))
            muy3 = mu(i,j-1)+mu(i,j+2)+3*(mu(i,j+1)+mu(i,j))
            muy4 = mu(i,j+1)-tf*(mu(i,j)+mu(i,j+2))

*** xx and yy derivatives:
            r1 = i6*((2*mux1+la(i-1,j)-tf*(la(i,j)+la(i-2,j)))*
     *                         (u(1,i-2,j)-u(1,i,j))+
     *           (2*mux2+la(i-2,j)+la(i+1,j)+3*(la(i,j)+la(i-1,j)))*
     *                         (u(1,i-1,j)-u(1,i,j))+ 
     *           (2*mux3+la(i-1,j)+la(i+2,j)+3*(la(i+1,j)+la(i,j)))*
     *                         (u(1,i+1,j)-u(1,i,j))+
     *           (2*mux4+ la(i+1,j)-tf*(la(i,j)+la(i+2,j)))*
     *           (u(1,i+2,j)-u(1,i,j)) 
     *              + muy1*(u(1,i,j-2)-u(1,i,j)) + 
     *                muy2*(u(1,i,j-1)-u(1,i,j)) + 
     *                muy3*(u(1,i,j+1)-u(1,i,j)) +
     *                muy4*(u(1,i,j+2)-u(1,i,j)) )

            r2 = i6*(mux1*(u(2,i-2,j)-u(2,i,j)) + 
     *                 mux2*(u(2,i-1,j)-u(2,i,j)) + 
     *                 mux3*(u(2,i+1,j)-u(2,i,j)) +
     *                 mux4*(u(2,i+2,j)-u(2,i,j)) + 
     *             (2*muy1+la(i,j-1)-tf*(la(i,j)+la(i,j-2)))*
     *                     (u(2,i,j-2)-u(2,i,j))+
     *             (2*muy2+la(i,j-2)+la(i,j+1)+3*(la(i,j)+la(i,j-1)))*
     *                     (u(2,i,j-1)-u(2,i,j))+ 
     *             (2*muy3+la(i,j-1)+la(i,j+2)+3*(la(i,j+1)+la(i,j)))*
     *                     (u(2,i,j+1)-u(2,i,j))+
     *             (2*muy4+la(i,j+1)-tf*(la(i,j)+la(i,j+2)))*
     *                     (u(2,i,j+2)-u(2,i,j)))

*** Mixed derivatives:
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j)*(u(2,i-2,j-2)-u(2,i-2,j+2)+
     *                        8*(-u(2,i-2,j-1)+u(2,i-2,j+1))) - 8*(
     *                   la(i-1,j)*(u(2,i-1,j-2)-u(2,i-1,j+2)+
     *                        8*(-u(2,i-1,j-1)+u(2,i-1,j+1))) )+8*(
     *                   la(i+1,j)*(u(2,i+1,j-2)-u(2,i+1,j+2)+
     *                        8*(-u(2,i+1,j-1)+u(2,i+1,j+1))) ) - (
     *                   la(i+2,j)*(u(2,i+2,j-2)-u(2,i+2,j+2)+
     *                        8*(-u(2,i+2,j-1)+u(2,i+2,j+1))) ))
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j)*(u(1,i-2,j-2)-u(1,i-2,j+2)+
     *                        8*(-u(1,i-2,j-1)+u(1,i-2,j+1))) - 8*(
     *                   mu(i-1,j)*(u(1,i-1,j-2)-u(1,i-1,j+2)+
     *                        8*(-u(1,i-1,j-1)+u(1,i-1,j+1))) )+8*(
     *                   mu(i+1,j)*(u(1,i+1,j-2)-u(1,i+1,j+2)+
     *                        8*(-u(1,i+1,j-1)+u(1,i+1,j+1))) ) - (
     *                   mu(i+2,j)*(u(1,i+2,j-2)-u(1,i+2,j+2)+
     *                        8*(-u(1,i+2,j-1)+u(1,i+2,j+1))) )) 
*** (mu*v_x)_y
            r1 = r1 + i144*( mu(i,j-2)*(u(2,i-2,j-2)-u(2,i+2,j-2)+
     *                        8*(-u(2,i-1,j-2)+u(2,i+1,j-2))) - 8*(
     *                   mu(i,j-1)*(u(2,i-2,j-1)-u(2,i+2,j-1)+
     *                        8*(-u(2,i-1,j-1)+u(2,i+1,j-1))) )+8*(
     *                   mu(i,j+1)*(u(2,i-2,j+1)-u(2,i+2,j+1)+
     *                        8*(-u(2,i-1,j+1)+u(2,i+1,j+1))) ) - (
     *                   mu(i,j+2)*(u(2,i-2,j+2)-u(2,i+2,j+2)+
     *                        8*(-u(2,i-1,j+2)+u(2,i+1,j+2))) )) 
*** (la*u_x)_y
            r2 = r2 + i144*( la(i,j-2)*(u(1,i-2,j-2)-u(1,i+2,j-2)+
     *                        8*(-u(1,i-1,j-2)+u(1,i+1,j-2))) - 8*(
     *                   la(i,j-1)*(u(1,i-2,j-1)-u(1,i+2,j-1)+
     *                        8*(-u(1,i-1,j-1)+u(1,i+1,j-1))) )+8*(
     *                   la(i,j+1)*(u(1,i-2,j+1)-u(1,i+2,j+1)+
     *                        8*(-u(1,i-1,j+1)+u(1,i+1,j+1))) ) - (
     *                   la(i,j+2)*(u(1,i-2,j+2)-u(1,i+2,j+2)+
     *                        8*(-u(1,i-1,j+2)+u(1,i+1,j+2))) ))
            up(1,i,j) = 2*u(1,i,j) - um(1,i,j) + cof*r1 + dt2*fo(1,i,j)
            up(2,i,j) = 2*u(2,i,j) - um(2,i,j) + cof*r2 + dt2*fo(2,i,j)
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine INNER_LOOP2ND( ni, nj, up, u, um, mu, la, rho, dt, h, 
     *                          fo )

*** Only centered approximation 

      implicit none

      real*8 half, fourth
      parameter( half=1d0/2, fourth=1d0/4 )

      integer ni, nj, i, j
      real*8 up(2,ni,nj), u(2,ni,nj), um(2,ni,nj)
      real*8 mu(ni,nj), la(ni,nj), fo(2,ni,nj), rho(ni,nj)
      real*8 mupx, mumx, mupy, mumy
      real*8 dt, la2, h, dt2, cof

      dt2 = dt*dt
      la2 = dt2/(h*h)
      do j=2,nj-1
         do i=2,ni-1
*** 9+39+39 =85 ops
            cof = la2/rho(i,j)
            mupx = half*(mu(i+1,j)+mu(i,j))
            mumx = half*(mu(i-1,j)+mu(i,j))
            mupy = half*(mu(i,j+1)+mu(i,j))
            mumy = half*(mu(i,j-1)+mu(i,j))
            up(1,i,j) = 2*u(1,i,j)- um(1,i,j) + cof*(
     *         (2*mupx + half*(la(i,j)+la(i+1,j)))*
     *                                       (u(1,i+1,j)-u(1,i,j)) -
     *         (2*mumx + half*(la(i,j)+la(i-1,j)))*
     *                                       (u(1,i,j)-u(1,i-1,j)) +
     *          fourth*(la(i+1,j)*(u(2,i+1,j+1)-u(2,i+1,j-1))-
     *                  la(i-1,j)*(u(2,i-1,j+1)-u(2,i-1,j-1)) )+
     *          fourth*(mu(i,j+1)*(u(2,i+1,j+1)-u(2,i-1,j+1)) -
     *                  mu(i,j-1)*(u(2,i+1,j-1)-u(2,i-1,j-1))) +
     *        mupy*(u(1,i,j+1)-u(1,i,j))-mumy*(u(1,i,j)-u(1,i,j-1)))
     *        +  dt2*fo(1,i,j)
            up(2,i,j) = 2*u(2,i,j) - um(2,i,j) + cof*(
     *         (2*mupy + half*(la(i,j)+la(i,j+1)))*
     *                                       (u(2,i,j+1)-u(2,i,j)) -
     *         (2*mumy + half*(la(i,j)+la(i,j-1)))*
     *                                       (u(2,i,j)-u(2,i,j-1)) +
     *          fourth*(mu(i+1,j)*(u(1,i+1,j+1)-u(1,i+1,j-1))-
     *                  mu(i-1,j)*(u(1,i-1,j+1)-u(1,i-1,j-1)) )+
     *          fourth*(la(i,j+1)*(u(1,i+1,j+1)-u(1,i-1,j+1)) -
     *                  la(i,j-1)*(u(1,i+1,j-1)-u(1,i-1,j-1))) +
     *        mupx*(u(2,i+1,j)-u(2,i,j))-mumx*(u(2,i,j)-u(2,i-1,j)))
     *        +  dt2*fo(2,i,j)
         enddo
      enddo
*** One sided operators at left and right boundaries:
      i = 3
      do j=2,nj-1
         cof = la2/rho(i,j)
         mupx = half*(mu(i+1,j)+mu(i,j))
         mumx = half*(mu(i-1,j)+mu(i,j))
         mupy = half*(mu(i,j+1)+mu(i,j))
         mumy = half*(mu(i,j-1)+mu(i,j))
         up(1,i,j) = 2*u(1,i,j)- um(1,i,j) + cof*(
     *        (2*mupx + half*(la(i,j)+la(i+1,j)))*
     *        (u(1,i+1,j)-u(1,i,j)) -
     *        (2*mumx + half*(la(i,j)+la(i-1,j)))*
     *        (u(1,i,j)-u(1,i-1,j)) +
     *          half*(la(i+1,j)*(u(2,i+1,j+1)-u(2,i+1,j-1))-
     *                la(i,j)  *(u(2,i,  j+1)-u(2,i,  j-1)) )+
     *          half*(mu(i,j+1)*(u(2,i+1,j+1)-u(2,i,  j+1)) -
     *                mu(i,j-1)*(u(2,i+1,j-1)-u(2,i,  j-1)) )+
     *        mupy*(u(1,i,j+1)-u(1,i,j))-mumy*(u(1,i,j)-u(1,i,j-1)))
     *        +  dt2*fo(1,i,j)
         up(2,i,j) = 2*u(2,i,j) - um(2,i,j) + cof*(
     *         (2*mupy + half*(la(i,j)+la(i,j+1)))*
     *                                       (u(2,i,j+1)-u(2,i,j)) -
     *         (2*mumy + half*(la(i,j)+la(i,j-1)))*
     *                                       (u(2,i,j)-u(2,i,j-1)) +
     *          half*(mu(i+1,j)*(u(1,i+1,j+1)-u(1,i+1,j-1)) -
     *                mu(i,  j)*(u(1,i,  j+1)-u(1,i,  j-1)) )+
     *          half*(la(i,j+1)*(u(1,i+1,j+1)-u(1,i,j+1)) -
     *                la(i,j-1)*(u(1,i+1,j-1)-u(1,i,j-1)) ) +
     *        mupx*(u(2,i+1,j)-u(2,i,j))-mumx*(u(2,i,j)-u(2,i-1,j)) )
     *        +  dt2*fo(2,i,j)
      enddo
      i = ni-1
      do j=2,nj-1
         cof = la2/rho(i,j)
         mupx = half*(mu(i+1,j)+mu(i,j))
         mumx = half*(mu(i-1,j)+mu(i,j))
         mupy = half*(mu(i,j+1)+mu(i,j))
         mumy = half*(mu(i,j-1)+mu(i,j))
         up(1,i,j) = 2*u(1,i,j)- um(1,i,j) + cof*(
     *         (2*mupx + half*(la(i,j)+la(i+1,j)))*
     *                                       (u(1,i+1,j)-u(1,i,j)) -
     *         (2*mumx + half*(la(i,j)+la(i-1,j)))*
     *                                       (u(1,i,j)-u(1,i-1,j)) +
     *          half*(la(i,j)  *(u(2,i,j+1)  -u(2,i,j-1))-
     *                la(i-1,j)*(u(2,i-1,j+1)-u(2,i-1,j-1)) )+
     *          half*(mu(i,j+1)*(u(2,i,j+1)-u(2,i-1,j+1)) -
     *                mu(i,j-1)*(u(2,i,j-1)-u(2,i-1,j-1))) +
     *        mupy*(u(1,i,j+1)-u(1,i,j))-mumy*(u(1,i,j)-u(1,i,j-1)))
     *        +  dt2*fo(1,i,j)
         up(2,i,j) = 2*u(2,i,j) - um(2,i,j) + cof*(
     *         (2*mupy + half*(la(i,j)+la(i,j+1)))*
     *                                       (u(2,i,j+1)-u(2,i,j)) -
     *         (2*mumy + half*(la(i,j)+la(i,j-1)))*
     *                                       (u(2,i,j)-u(2,i,j-1)) +
     *          half*(mu(i,j)  *(u(1,i,j+1)  -u(1,i,j-1))-
     *                mu(i-1,j)*(u(1,i-1,j+1)-u(1,i-1,j-1)) )+
     *          half*(la(i,j+1)*(u(1,i,j+1)  -u(1,i-1,j+1)) -
     *                la(i,j-1)*(u(1,i,j-1)  -u(1,i-1,j-1))) +
     *        mupx*(u(2,i+1,j)-u(2,i,j))-mumx*(u(2,i,j)-u(2,i-1,j)))
     *        +  dt2*fo(2,i,j)
      enddo
      end

c-----------------------------------------------------------------------
      subroutine BC( nc, nx, ny, var, mu, la, s, bccnd, om, omm, ph, 
     *               phm, th, la0, t, h, tder, fobnd, ini, acc, k1, k2,
     *               rp, rs, alpha, onlyswave )
      implicit none
      integer nc, nx, ny, bccnd(4), c, i, j, tder, ini, acc
      real*8 var(nc,-1:nx+1,-1:ny+1), mu(-1:nx+1,-1:ny+1)
      real*8 la(-1:nx+1,-1:ny+1)
      real*8 s(0:4), ph, phm, th, omm, om, t, h, x, y
      real*8 fobnd(2,-1:ny+1), la0
      real*8 cRel, csw0, dxi, latilde, pr, xi, CSW, CSWP, amp1, amp2
      real*8 phase, xi2, cPhase, k1, k2, alpha, rp, rs, ifact, rpeff
      integer iter
      logical onlyswave

*** side j=1
      if( bccnd(3).eq.8 )then
*** Periodic in y
         do i=-1,nx+1
            do c=1,nc
               var(c,i,-1) = var(c,i,ny-2)
               var(c,i,0)  = var(c,i,ny-1)
               var(c,i,ny)  = var(c,i,1)
               var(c,i,ny+1)= var(c,i,2)
            enddo
         enddo
      elseif( bccnd(3).eq.1 )then
*** Dirichlet in y
         do i=-1,nx+1
            do c=1,nc
               var(c,i,1)  = 0
            enddo
         enddo
      endif

*** side j=nj
      if( bccnd(4).eq.1 )then
         do i=-1,nx+1
            do c=1,nc
               var(c,i,ny)  = 0
            enddo
         enddo
      endif

c all ok up to here

*** side i=1
      if( bccnd(1).eq.8 )then
*** Periodic in x
         do j=-1,ny+1
            do c=1,nc
               var(c,-1,j)  = var(c,nx-2,j)
               var(c,0,j)   = var(c,nx-1,j)
               var(c,nx,j)   = var(c,1,j)
               var(c,nx+1,j) = var(c,2,j)
            enddo
         enddo
      elseif( bccnd(1).eq.1 )then
         do j=-1,ny+1
            do c=1,nc
               var(c,1,j)  = 0
            enddo
         enddo
      elseif( bccnd(1).eq.2 )then
*** Twilight data
         x = 0
         if( tder.eq.0 )then
            do j=-1,ny+1
               y = (j-1)*h
               var(1,1,j)  = cos(om*x+ph)*sin(om*y+th)*cos(t*t)
               var(2,1,j)  = sin(om*x+th)*cos(om*y+ph)*sin(t)
            enddo
         elseif( tder.eq.2 )then
            do j=-1,ny+1
               y = (j-1)*h
               var(1,1,j)  = cos(om*x+ph)*sin(om*y+th)*
     *                          (-2*sin(t*t)-4*t*t*cos(t*t))
               var(2,1,j)  =-sin(om*x+th)*cos(om*y+ph)*sin(t)
            enddo
         endif            

      elseif( bccnd(1).eq.3 )then
         if( ini.eq.1 )then
            call BCTWFORCING( ny, h, t, 1, om, omm, ph, phm, 
     *               th, la0, fobnd )
         elseif( ini.eq.4 )then
            if( onlyswave )then
               rpeff = rp
            else
               rpeff = 0
            endif
            do j=-1,ny+1
               y = (j-1)*h
               fobnd(1,j) = -(2*mu(1,j)*k1*k1+la(1,j))*
     *                                   COS( om*t + k2*y )*(1+rpeff)
               fobnd(2,j) = -2*mu(1,j)*k1*k2*
     *                                   COS( om*t + k2*y )*(1-rpeff)
            enddo
         else
            do j=-1,ny+1
               fobnd(1,j) = 0
               fobnd(2,j) = 0
            enddo
         endif
         if( acc.eq.4 )then
            call BCFREE( nx, ny, var, mu, la, s, .true., .false., h, 
     *                   fobnd)
         else
            call BCFREE2ND( nx, ny, var, mu, la, .true., .false., h,
     *                      fobnd)
         endif
      endif

c all ok up to here

*** side i=n
      if( bccnd(2).eq.1 )then
         do j=-1,ny+1
            do c=1,nc
               var(c,nx,j)  = 0
            enddo
         enddo
      elseif( bccnd(2).eq.2 )then

c all ok up to here

c far right boundary
         x = (nx-1)*h
         if (ini.eq.1) then
*** Twilight data
           if( tder.eq.0 )then
             do j=-1,ny+1
               y = (j-1)*h
               var(1,nx,j)  = cos(om*x+ph)*sin(om*y+th)*cos(t*t)
               var(2,nx,j)  = sin(om*x+th)*cos(om*y+ph)*sin(t)
             enddo
           elseif( tder.eq.2 )then
             do j=-1,ny+1
               y = (j-1)*h
               var(1,nx,j)  = cos(om*x+ph)*sin(om*y+th)*
     *              (-2*sin(t*t)-4*t*t*cos(t*t))
               var(2,nx,j)  =-sin(om*x+th)*cos(om*y+ph)*sin(t)
             enddo
           endif
         else if (ini.eq.4) then
            ifact = 1/SQRT(alpha*alpha*k1*k1+k2*k2)
            if( onlyswave )then
               rpeff = 0
            else
               rpeff = rp
            endif
            do j=-1,ny+1
               y = (j-1)*h
               var(1,nx,j) = rpeff*(-k1)*sin(om*(t)-k1*x+k2*y) - 
     *                     rs*k2*ifact*sin(om*(t)-alpha*k1*x+k2*y)
               var(2,nx,j) = rpeff*k2*sin(om*(t)-k1*x+k2*y) -
     *                   rs*alpha*k1*ifact*sin(om*(t)-alpha*k1*x+k2*y)
            enddo
         endif
      elseif( bccnd(2).eq.3 )then
         if( ini.eq.1 )then
            call BCTWFORCING( ny, h, t, nx, om, omm, ph, phm,
     *                      th, la0, fobnd )
         else
            do j=-1,ny+1
               fobnd(1,j) = 0
               fobnd(2,j) = 0
            enddo
         endif
         if( acc.eq.4 )then
            call BCFREE( nx, ny, var, mu, la, s, .false., .true., h,
     *                   fobnd )
         else
            call BCFREE2ND( nx, ny, var, mu, la, .false., .true., h,
     *                      fobnd)
         endif
      endif

      if( bccnd(3).eq.8 )then
*** Periodic in y once again, to be on the safe side
         do i=-1,nx+1
            do c=1,nc
               var(c,i,-1) = var(c,i,ny-2)
               var(c,i,0)  = var(c,i,ny-1)
               var(c,i,ny)  = var(c,i,1)
               var(c,i,ny+1)= var(c,i,2)
            enddo
         enddo
      endif
      end

c-----------------------------------------------------------------------
      subroutine BCFREE( nx, ny, u, mu, la, s, lower, upper, h, fobnd )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=1d0/12 )
      integer nx, ny, j, k
      real*8 u(2,-1:nx+1,-1:ny+1), mu(-1:nx+1,-1:ny+1)
      real*8 la(-1:nx+1,-1:ny+1)
      real*8 s(0:4), d0yu, d0yv, h, fobnd(2,-1:ny+1)
      logical upper, lower

      if( lower )then
         do j=1,ny-1
            d0yu=d4b*(u(1,1,j-2)-u(1,1,j+2))+d4a*(u(1,1,j+1)-u(1,1,j-1))
     *                 -h*fobnd(2,j)/mu(1,j)
            d0yv=d4b*(u(2,1,j-2)-u(2,1,j+2))+d4a*(u(2,1,j+1)-u(2,1,j-1))
            d0yv = (d0yv*la(1,j)-h*fobnd(1,j))/(2*mu(1,j)+la(1,j))
            do k=1,4
               d0yu = d0yu + u(2,k,j)*s(k)
               d0yv = d0yv + u(1,k,j)*s(k)
            enddo
            u(1,0,j) = -d0yv/s(0)
            u(2,0,j) = -d0yu/s(0)
         enddo
      endif
      if( upper )then
         do j=1,ny-1
            d0yu=d4b*(u(1,nx,j-2)-u(1,nx,j+2))+
     *           d4a*(u(1,nx,j+1)-u(1,nx,j-1))
     *                 -h*fobnd(2,j)/mu(nx,j)
            d0yv=d4b*(u(2,nx,j-2)-u(2,nx,j+2))+
     *           d4a*(u(2,nx,j+1)-u(2,nx,j-1))
            d0yv = (d0yv*la(nx,j)-h*fobnd(1,j))/(2*mu(nx,j)+la(nx,j))
            do k=1,4
               d0yu = d0yu - u(2,nx-k+1,j)*s(k)
               d0yv = d0yv - u(1,nx-k+1,j)*s(k)
            enddo
            u(1,nx+1,j) = d0yv/s(0)
            u(2,nx+1,j) = d0yu/s(0)
         enddo
      endif
      end

c-----------------------------------------------------------------------
      subroutine BCFREE2ND( nx, ny, u, mu, la, lower, upper, h,fobnd)
      implicit none
      real*8 half
      parameter( half=1d0/2 )
      integer nx, ny, j
      real*8 u(2,-1:nx+1,-1:ny+1), mu(-1:nx+1,-1:ny+1)
      real*8 la(-1:nx+1,-1:ny+1)
      real*8 imy,  h, fobnd(2,-1:ny+1)
      logical upper, lower

      if( lower )then
         do j=1,ny
            imy = 1d0/(half*(la(1,j)+la(0,j))+mu(1,j)+mu(0,j))
            u(1,0,j) = u(1,1,j) + imy*( 
     *   (mu(1,j)+mu(2,j)+half*(la(1,j)+la(2,j)))*(u(1,2,j)-u(1,1,j))
     *        + la(1,j)*(u(2,1,j+1)-u(2,1,j-1)) - 2*h*fobnd(1,j) )
            imy = 1d0/(half*(mu(1,j)+mu(0,j)))
            u(2,0,j) = u(2,1,j) + imy*( 
     *           half*(mu(2,j)+mu(1,j))*(u(2,2,j)-u(2,1,j))
     *       + mu(1,j)*(u(1,1,j+1)-u(1,1,j-1))  -2*h*fobnd(2,j) )
         enddo
      endif
      if( upper )then
         do j=1,ny
            imy = 1d0/(half*(la(nx,j)+la(nx+1,j))+mu(nx,j)+mu(nx+1,j))
            u(1,nx+1,j) = u(1,nx,j) - imy*( 
     *   (mu(nx,j)+mu(nx-1,j)+half*(la(nx,j)+la(nx-1,j)))*
     *                                (u(1,nx,j)-u(1,nx-1,j))
     *        + la(nx,j)*(u(2,nx,j+1)-u(2,nx,j-1)) - 2*h*fobnd(1,j) )
            imy = 1d0/(half*(mu(nx,j)+mu(nx+1,j)))
            u(2,nx+1,j) = u(2,nx,j) - imy*( 
     *           half*(mu(nx-1,j)+mu(nx,j))*(u(2,nx,j)-u(2,nx-1,j))
     *       + mu(nx,j)*(u(1,nx,j+1)-u(1,nx,j-1))  -2*h*fobnd(2,j) )
         enddo
      endif
      end

c-----------------------------------------------------------------------
      subroutine exactRhs( ifirst, ilast, jfirst, jlast, t, fo, om, omm,
     +     ph, phm, th, la0, h )

      implicit none
      integer i, j, ifirst, ilast, jfirst, jlast
      real*8 fo(2,ifirst:ilast,jfirst:jlast), h

      doubleprecision x, y, t, om, ph, th, omm, phm, la0

      doubleprecision acc(2)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t11
      doubleprecision t13
      doubleprecision t14
      doubleprecision t15
      doubleprecision t16
      doubleprecision t18
      doubleprecision t19
      doubleprecision t20
      doubleprecision t22
      doubleprecision t23
      doubleprecision t24
      doubleprecision t26
      doubleprecision t27
      doubleprecision t3
      doubleprecision t30
      doubleprecision t31
      doubleprecision t33
      doubleprecision t34
      doubleprecision t35
      doubleprecision t36
      doubleprecision t38
      doubleprecision t4
      doubleprecision t40
      doubleprecision t42
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t47
      doubleprecision t51
      doubleprecision t52
      doubleprecision t55
      doubleprecision t57
      doubleprecision t59
      doubleprecision t6
      doubleprecision t61
      doubleprecision t63
      doubleprecision t65
      doubleprecision t68
      doubleprecision t7
      doubleprecision t71
      doubleprecision t79
      doubleprecision t8
      doubleprecision t85
      doubleprecision t88
      doubleprecision t90

      do j=jfirst,jlast
         y = (j-1)*h
         do i=ifirst,ilast
            x = (i-1)*h
        t1 = omm*x
        t3 = 3*t1+phm
        t4 = cos(t3)
        t6 = omm*y
        t7 = sin(t6)
        t8 = t4*omm*t7
        t10 = t1+phm
        t11 = sin(t10)
        t13 = 3*t6
        t14 = sin(t13)
        t15 = t14**2
        t16 = t11*omm*t15
        t18 = om*x
        t19 = t18+ph
        t20 = sin(t19)
        t22 = om*y
        t23 = t22+th
        t24 = sin(t23)
        t26 = t**2
        t27 = cos(t26)
        t30 = sin(t3)
        t31 = t30*t7
        t33 = cos(t10)
        t34 = t33*t15
        t35 = 2*t31+6+t34+la0
        t36 = cos(t19)
        t38 = om**2
        t40 = t38*t24*t27
        t42 = t18+th
        t43 = sin(t42)
        t44 = t22+ph
        t45 = sin(t44)
        t47 = sin(t)
        t51 = t34+la0
        t52 = cos(t42)
        t55 = t38*t45*t47
        t57 = cos(t6)
        t59 = t30*t57*omm
        t61 = cos(t44)
        t63 = t52*om*t61*t47
        t65 = t31+3
        t68 = cos(t23)
        t71 = t36*t68*om*t27
        acc(1) = -(6*t8-t16)*t20*om*t24*t27-t35*t36*t40+t16*t43*t45*om*t
     +       47-t51*t52*t55+t59*t63-t65*t52*t55+t59*t71-t65*t36*t40
        t79 = t38*t61*t47
        t85 = t38*t68*t27
        t88 = cos(t13)
        t90 = t33*t14*t88*omm
        acc(2) = 3*t8*t63-t65*t43*t79+3*t8*t71-t65*t20*t85-6*t90*t20*om*
     +      t24*t27-t51*t20*t85-(2*t59+6*t90)*t43*t45*om*t47-t35*t43*t79

        fo(1,i,j) = acc(1)
        fo(2,i,j) = acc(2)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine FORCING( nx, ny, t, fo, om, omm, ph, phm, th, la0, h )

      implicit none
      integer i, j, nx, ny
      real*8 t, fo(2,-1:nx+1,-1:ny+1), h
      doubleprecision x
      doubleprecision y
      doubleprecision om
      doubleprecision omm
      doubleprecision ph
      doubleprecision phm
      doubleprecision th
      doubleprecision la0


      doubleprecision eqs(2)
      doubleprecision t1
      doubleprecision t107
      doubleprecision t109
      doubleprecision t13
      doubleprecision t16
      doubleprecision t18
      doubleprecision t19
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t22
      doubleprecision t25
      doubleprecision t26
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t30
      doubleprecision t31
      doubleprecision t33
      doubleprecision t38
      doubleprecision t39
      doubleprecision t4
      doubleprecision t41
      doubleprecision t42
      doubleprecision t43
      doubleprecision t45
      doubleprecision t49
      doubleprecision t5
      doubleprecision t50
      doubleprecision t51
      doubleprecision t52
      doubleprecision t54
      doubleprecision t58
      doubleprecision t59
      doubleprecision t6
      doubleprecision t64
      doubleprecision t65
      doubleprecision t67
      doubleprecision t68
      doubleprecision t7
      doubleprecision t70
      doubleprecision t74
      doubleprecision t77
      doubleprecision t8
      doubleprecision t87
      doubleprecision t9
      doubleprecision t90
      do j=-1,ny+1
         y = (j-1)*h
         do i=-1,nx+1
            x = (i-1)*h
        t1 = om*x
        t2 = t1+ph
        t3 = cos(t2)
        t4 = om*y
        t5 = t4+th
        t6 = sin(t5)
        t7 = t3*t6
        t8 = t**2
        t9 = cos(t8)
        t13 = sin(t8)
        t16 = omm*x
        t18 = 3*t16+phm
        t19 = cos(t18)
        t20 = t19*omm
        t21 = omm*y
        t22 = sin(t21)
        t25 = t16+phm
        t26 = sin(t25)
        t28 = 3*t21
        t29 = sin(t28)
        t30 = t29**2
        t31 = t26*omm*t30
        t33 = sin(t2)
        t38 = sin(t18)
        t39 = t38*t22
        t41 = cos(t25)
        t42 = t41*t30
        t43 = 2*t39+6+t42+la0
        t45 = om**2
        t49 = t1+th
        t50 = sin(t49)
        t51 = t4+ph
        t52 = sin(t51)
        t54 = sin(t)
        t58 = t42+la0
        t59 = cos(t49)
        t64 = cos(t21)
        t65 = t38*t64
        t67 = cos(t51)
        t68 = t67*t54
        t70 = cos(t5)
        t74 = t59*om*t68+t3*t70*om*t9
        t77 = t39+3
        t87 = sin(t4-th)
        t90 = 1/(t33*t87+2)
        eqs(1) = -4*t7*t9*t8-2*t7*t13-(-(6*t20*t22-t31)*t33*om*t6*t9-t43
     #*t3*t45*t6*t9+t31*t50*t52*om*t54-t58*t59*t45*t52*t54+t65*omm*t74+t
     #77*(-t59*t45*t52*t54-t7*t45*t9))*t90
        t107 = cos(t28)
        t109 = t41*t29*t107*omm
        eqs(2) = -t50*t67*t54-(3*t20*t22*t74+t77*(-t50*t45*t68-t33*t45*t
     #70*t9)-(2*t65*omm+6*t109)*t50*t52*om*t54-t43*t50*t67*t45*t54-6*t10
     #9*t33*om*t6*t9-t58*t33*t45*t70*t9)*t90
        fo(1,i,j) = eqs(1)
        fo(2,i,j) = eqs(2)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine FORCINGTT( nx, ny, t, fo, om, omm, ph, phm, th, la0, h)
      implicit none

      integer i, j, nx, ny
      real*8 t, fo(2,-1:nx+1,-1:ny+1), h
      doubleprecision x
      doubleprecision y
      doubleprecision om
      doubleprecision omm
      doubleprecision ph
      doubleprecision phm
      doubleprecision th
      doubleprecision la0

      doubleprecision eqstt(2)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t111
      doubleprecision t114
      doubleprecision t123
      doubleprecision t125
      doubleprecision t136
      doubleprecision t138
      doubleprecision t14
      doubleprecision t149
      doubleprecision t157
      doubleprecision t2
      doubleprecision t20
      doubleprecision t22
      doubleprecision t23
      doubleprecision t24
      doubleprecision t25
      doubleprecision t26
      doubleprecision t29
      doubleprecision t3
      doubleprecision t30
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t35
      doubleprecision t37
      doubleprecision t38
      doubleprecision t4
      doubleprecision t41
      doubleprecision t48
      doubleprecision t49
      doubleprecision t5
      doubleprecision t51
      doubleprecision t52
      doubleprecision t53
      doubleprecision t54
      doubleprecision t55
      doubleprecision t6
      doubleprecision t63
      doubleprecision t64
      doubleprecision t65
      doubleprecision t66
      doubleprecision t68
      doubleprecision t7
      doubleprecision t72
      doubleprecision t73
      doubleprecision t78
      doubleprecision t79
      doubleprecision t8
      doubleprecision t81
      doubleprecision t82
      doubleprecision t84
      doubleprecision t85
      doubleprecision t9
      doubleprecision t93
      doubleprecision t96
      do j=-1,ny+1
         y = (j-1)*h
         do i=-1,nx+1
            x = (i-1)*h

        t1 = om*x
        t2 = t1+ph
        t3 = cos(t2)
        t4 = om*y
        t5 = t4+th
        t6 = sin(t5)
        t7 = t3*t6
        t8 = t**2
        t9 = cos(t8)
        t10 = t8**2
        t14 = sin(t8)
        t20 = omm*x
        t22 = 3*t20+phm
        t23 = cos(t22)
        t24 = t23*omm
        t25 = omm*y
        t26 = sin(t25)
        t29 = t20+phm
        t30 = sin(t29)
        t32 = 3*t25
        t33 = sin(t32)
        t34 = t33**2
        t35 = t30*omm*t34
        t37 = sin(t2)
        t38 = (6*t24*t26-t35)*t37
        t41 = t6*t9*t8
        t48 = sin(t22)
        t49 = t48*t26
        t51 = cos(t29)
        t52 = t51*t34
        t53 = 6+2*t49+t52+la0
        t54 = t53*t3
        t55 = om**2
        t63 = t1+th
        t64 = sin(t63)
        t65 = t4+ph
        t66 = sin(t65)
        t68 = sin(t)
        t72 = t52+la0
        t73 = cos(t63)
        t78 = cos(t25)
        t79 = t48*t78
        t81 = cos(t65)
        t82 = t81*t68
        t84 = cos(t5)
        t85 = t3*t84
        t93 = -t73*om*t82-4*t85*om*t9*t8-2*t85*om*t14
        t96 = 3+t49
        t111 = sin(t4-th)
        t114 = 1/(2+t37*t111)
        eqstt(1) = 16*t7*t9*t10+48*t7*t14*t8-12*t7*t9-(4*t38*om*t41+2*t3
     #8*om*t6*t14+4*t54*t55*t41+2*t54*t55*t6*t14-t35*t64*t66*om*t68+t72*
     #t73*t55*t66*t68+t79*omm*t93+t96*(t73*t55*t66*t68+4*t7*t55*t9*t8+2*
     #t7*t55*t14))*t114
        t123 = t37*t55
        t125 = t84*t9*t8
        t136 = cos(t32)
        t138 = t51*t33*t136*omm
        t149 = t37*om
        t157 = t72*t37
        eqstt(2) = t64*t81*t68-(3*t24*t26*t93+t96*(t64*t55*t82+4*t123*t1
     #25+2*t123*t84*t14)+(2*t79*omm+6*t138)*t64*t66*om*t68+t53*t64*t81*t
     #55*t68+24*t138*t149*t41+12*t138*t149*t6*t14+4*t157*t55*t125+2*t157
     #*t55*t84*t14)*t114
        fo(1,i,j) = eqstt(1)
        fo(2,i,j) = eqstt(2)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine BCTWFORCING( n, h, t, ind, om, omm, ph, phm, th, 
     *                        la0, fo )
      implicit none

      integer n, ind, j
      real*8  h, fo(2,-1:n+1)

      doubleprecision x
      doubleprecision y
      doubleprecision t
      doubleprecision om
      doubleprecision omm
      doubleprecision ph
      doubleprecision phm
      doubleprecision th
      doubleprecision la0

      doubleprecision bforce(2)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t12
      doubleprecision t13
      doubleprecision t14
      doubleprecision t16
      doubleprecision t17
      doubleprecision t18
      doubleprecision t20
      doubleprecision t21
      doubleprecision t22
      doubleprecision t24
      doubleprecision t25
      doubleprecision t29
      doubleprecision t30
      doubleprecision t32
      doubleprecision t33
      doubleprecision t35
      doubleprecision t38
      doubleprecision t39
      doubleprecision t4
      doubleprecision t41
      doubleprecision t45
      doubleprecision t47
      doubleprecision t5
      doubleprecision t6
      doubleprecision t7

      x = (ind-1)*h
      do j=-1,n+1
         y = (j-1)*h
        t1 = omm*x
        t4 = sin(3*t1+phm)
        t5 = omm*y
        t6 = sin(t5)
        t7 = t4*t6
        t10 = cos(t1+phm)
        t12 = sin(3*t5)
        t13 = t12**2
        t14 = t10*t13
        t16 = om*x
        t17 = t16+ph
        t18 = sin(t17)
        t20 = om*y
        t21 = t20+th
        t22 = sin(t21)
        t24 = t**2
        t25 = cos(t24)
        t29 = t16+th
        t30 = sin(t29)
        t32 = t20+ph
        t33 = sin(t32)
        t35 = sin(t)
        bforce(1) = -(2*t7+6+t14+la0)*t18*om*t22*t25-(t14+la0)*t30*t33*o
     #m*t35
        t38 = t7+3
        t39 = cos(t29)
        t41 = cos(t32)
        t45 = cos(t17)
        t47 = cos(t21)
        bforce(2) = t38*t39*om*t41*t35+t38*t45*t47*om*t25
        fo(1,j) = bforce(1)
        fo(2,j) = bforce(2)
      enddo
      end
