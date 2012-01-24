      implicit none

      real*8 gettimec

      integer nx, ny, nz, i, j, k, c, m, q, acc, seed, rn, bccnd(6), ini
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer ni, nj, nk, kb, mb, qb, kp
      integer maxsteps, st, nsteps, istep, tacc, errfreq
      integer onesided(6)
      integer, dimension(:), allocatable :: seeds
      real*8 Lx, Ly, Lz, h, x, y, z, dt, err(3), li(3), l2(3), uex(3)
      real*8 cpuTime
      real*8, allocatable, dimension(:,:,:,:) :: up,u,um,fo,lu,u2,um2
      real*8, allocatable, dimension(:,:,:) :: mu, la, rho, bforce
      real*8 acof(6,8,8), bop(4,6), bope(6,8), ghcof(6), hnorm(4)
      real*8 dt2, pi, dtpr, sbop(0:4)
      real*8 locrad, maxrad, om, omm, ph, phm, la0, cfl
      real*8 amplambda, ampmu, amprho
      real*8 energy, normfact, term, t, tfinal
      real*8 iop(5), iop2(5), bop2(4,6), gh2, s(0:4), RVEL
      real*8 cv, cRel, Tperiod
      real*8 k1, k2, rp, rs, alpha, tmp1, tmp2, ifact, ang, muval, rpeff
      real*8 i12, dt4, dt2i
      parameter(i12=1D0/12)

      character*400 buf
      logical writesol, cflset

      pi  = 4*ATAN(1d0)
c      L   = 5d0
      Lx   = 1.5
      Ly   = 1
      Lz   = 0.8
      ny   = 26
      dt   = 1
      i    = 1
      acc  = 4
      seed = 298347
c by default, all sides use centered formulas
      do c=1,6
        onesided(c) = 0
      enddo

c default boundary conditions
*** bccnd=1 is Dirichlet (homogeneous)
      do c=1,6
        bccnd(c) = 1
      enddo

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
      errfreq  = 10
      writesol = .true.
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
         elseif( buf.eq.'-ang' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') ang
            i = i+2
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
         Lx = 6.28
         Ly = 6.28
         Lz = 6.28
      endif

      h  = Ly/(ny-1)
      nx = Lx/h + 1.5
      Lx = (nx-1)*h
      nz = Lz/h + 1.5
      Lz = (nz-1)*h
      write(*,105) 'h = ', h
      write(*,105) 'Domain (Lx, Ly, Lz) = ', Lx, ' x ', Ly, ' x ', Lz
 105  format(' ',a,g12.5,a,g12.5,a,g12.5)
      write(*,104) 'Number of points = ', nx, ' x ', ny, ' x ', nz
 104  format(' ',a,i5,a,i5,a,i5)
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
c only 1 ghost point in z
c$$$      kfirst = 0
c$$$      klast = nz+1
c to be perfectly consistent with the C++ program, we use the 
c same number of ghost point in all directions
      kfirst = -1
      klast = nz+2

c total number of grid points (ni=nx+4, nj=ny+4)
      ni = ilast-ifirst+1
      nj = jlast-jfirst+1
      nk = klast-kfirst+1

c allocate solution and material coefficient arrays
      allocate( up(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( u(3,  ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( um(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( fo(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( u2(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( um2(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( lu(3, ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( rho(  ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( mu(   ifirst:ilast, jfirst:jlast, kfirst:klast) )
      allocate( la(   ifirst:ilast, jfirst:jlast, kfirst:klast) )

c temporary workspace
      allocate( bforce(3, ifirst:ilast, jfirst:jlast) )

c get coefficients for difference approximation of 2nd derivative with variable coefficients
      call VARCOEFFS4( acof, ghcof )
c get coefficients for difference approximation of 1st derivative
      call WAVEPROPBOP_4( iop, iop2, bop, bop2, gh2, hnorm, sbop )
c extend the definition of the 1st derivative tothe first 6 points
      call BOPEXT4TH( bop, bope )

c test 1-D operators
c$$$      call accur1d(h, nz, bope, acof, ghcof)

c setup material and initial conditions
      if( ini.eq.1 )then
*** Twilight testing
*** bccnd(1): low-i
*** bccnd(2): high-i
*** bccnd(3): low-j
*** bccnd(4): high-j
*** bccnd(5): low-k
*** bccnd(6): high-k
c
*** bccnd=1 is Dirichlet (homogeneous)
*** bccnd=2 is Dirichlet (with tw-forcing)
*** bccnd=3 is Free surface (homogeneous)
*** bccnd=4 is Free surface (with tw-forcing)
        do c=1,6
          bccnd(c) = 2
        enddo
c free surface on low-k and high-k boundary
        bccnd(5) = 4
c        bccnd(6) = 4

        ph = 0.3d0
        phm= 0.1d0
        om  = 1d0
        omm = 1d0
        cv = 1.3
        amprho = 1.0
        ampmu = 1.0
        amplambda = 1.0
        t   = 0
        st  = 0
c setup material model
        call exactMat( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       rho, mu, la, omm, phm, amprho, ampmu, amplambda, h )
      endif

c one-sided formulas are used together with free surface bc
      do c=1,6
        if (bccnd(c).eq.3 .or. bccnd(c).eq.4) onesided(c)=1
      enddo
      write(*,'( a, $)') ' Boundary conditions: '
      do c=1,6
        write(*,'(a, i1, a, i2, $)') ' bc(', c, ')= ', bccnd(c)
      end do
      write(*,*)
      write(*,'( a, $)') ' Onesided: '
      do c=1,6
        write(*,'(a, i1, a, i2, $)') ' onesided(', c, ')= ', onesided(c)
      end do
      write(*,*)

      if (.not.cflset) then
        if (tacc.eq.4) then
          cfl=1.3
        else
          cfl=0.9
        endif
      endif

*** Determine time step
      call setDT(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     mu, la, rho, h, cfl, dt)
      nsteps = tfinal/dt + 0.5
      write(*,202)  tacc, cfl, dt, errfreq
 202  format(' Order of accuracy in time =', i2, ' Cfl =', f7.3, 
     +     ' prel dt=', e14.3, ' Err freq =', i4)
c correct the time step to exactly end up at tfinal
      dt = tfinal/nsteps

      dt2 = dt*dt
      dt2i= 1d0/dt2
      dt4 = dt2*dt2

c assign initial data
      if( ini.eq.1 )then
        call exactSol( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       u, 0D0, om, cv, ph, h )
*** Can not assign initial data at t=-dt until dt is computed.
        call exactSol( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       um, -dt, om, cv, ph, h )
      endif

      t   = 0
c beginning of various tests
      write(*,*)
      write(*,*) 'Testing the accuracy of the spatial diff. approx.'
c test accuracy of spatial approximation
      call exactRhs3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, cv, ph, omm, phm, amprho, ampmu, amplambda, h)

c evaluate right hand side: L(u), result in array 'u2'
      call RHS4TH3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nz, onesided, acof, bope, ghcof,
     +     u2, u, mu, la, rho, h )

c error in rhs
      call rhsErr(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
     +     h, fo, u2, 0)

c test accuracy of forcing
      write(*,*) 'Testing accuracy of F = rho*utt - Lu'
      call forcingERR(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     nz, u, u2, um2, lu, fo, rho, mu, la,
     +     onesided, acof, bope, ghcof, t, dt, h, 
     +     om, cv, ph, omm, phm, amprho, ampmu, amplambda, nz/2, 1)

c impose bc on initial data (required before time stepping is started)
      call BC( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, 0d0,
     *     bforce, om, ph, cv ) 

      call BC( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nx, ny, nz,
     +     um, h, bccnd, sbop, mu, la, -dt,
     *     bforce, om, ph, cv ) 

c test error in ghost points. Put exact solution in array 'up'
      write(*,*)
      write(*,*) 'Testing error in ghost points after enforcing BC'
      call exactSol( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     up, 0D0, om, cv, ph, h )
      call solErr(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, up, u, -nz/2)
c end of tests
c tmp
c      stop

c time-stepping (time stepping) loop
      t = 0
      st = 0
      write(*,*)
      write(*,*) 'Starting time-stepping'
      cpuTime = gettimec()
      do istep=1,nsteps

c evaluate twilight forcing
        call FORCING( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, cv, ph, omm, phm, amprho, ampmu, amplambda, h)

c evaluate right hand side: L(u), result in array 'lu'
        call RHS4TH3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       nz, onesided, acof, bope, ghcof,
     +       lu, u, mu, la, rho, h )

c assign solution at t+dt
        call predictor(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, lu, fo, rho, dt2 )

c impose bc on up
        call BC( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       nx, ny, nz,
     +       up, h, bccnd, sbop, mu, la, t+dt,
     *       bforce, om, ph, cv ) 

*** Modified equation to get 4th order in time
        if( tacc.eq.4 )then
          call FORCINGTT( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +         fo, t, om, cv, ph, omm, phm, amprho, ampmu, amplambda, h)

c evaluate D+t D-t(u) at all grid points, result in array 'u2'
          call dpdmInTime(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, u2, dt2i)

c evaluate right hand side: L(u2), result in array 'lu'
          call RHS4TH3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +         nz, onesided, acof, bope, ghcof,
     +         lu, u2, mu, la, rho, h )

c correct solution at t+dt to 4th order accuracy
          call corrector(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +         up, lu, fo, rho, dt4 )

c impose bc on corrected solution
          call BC( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +         nx, ny, nz,
     +         up, h, bccnd, sbop, mu, la, t+dt,
     *         bforce, om, ph, cv ) 
        endif
        
        call cycleArrays(ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       up, u, um )

c increment time step and time
        st = st+1
        t = t + dt
      enddo
      cpuTime = gettimec()-cpuTime

      write(*,103) 'Taken ', st, ' steps to time ',t
 103  format(' ', a, i7, a, g12.5 )
      write(*,102) 'Loop execution time ', cpuTime, ' seconds'
 102  format(' ', a, g15.7, a)
      
*** Compute errors at the final time
      if( ini.eq.1 )then
c store exact solution in array 'up'
        call exactSol( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +       up, t, om, cv, ph, h )
        call solErr(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, up, u, nz/2)
      endif

c save the final solution
      if( writesol )then
        call saveSol(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +       nz/2, u)
      endif

      end

c-----------------------------------------------------------------------
      subroutine RHS4TH3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nz, onesided, acof, bope, ghcof,
     +     uacc, u, mu, la, rho, h )

*** Only centered approximation of the right hand side of the elastic wave equation

      implicit none

      real*8 tf, i6, i144, i12
      parameter( tf=3d0/4, i6=1d0/6, i144=1d0/144, i12=1d0/12 )

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      integer nz, onesided(6), m, q, kb, mb, qb, k1, k2
      real*8 acof(6,8,8), bope(6,8), ghcof(6)
      real*8 uacc(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mux1, mux2, mux3, mux4, muy1, muy2, muy3, muy4
      real*8 muz1, muz2, muz3, muz4
      real*8 mucof(8), mu1zz, mu2zz, lau2yz
      real*8 lap2mu(8), mu3zz, mu3xz, mu3yz, lau1xz
      real*8 lau3zx, u3zim2, u3zim1, u3zip1, u3zip2
      real*8 lau3zy, u3zjm2, u3zjm1, u3zjp1, u3zjp2
      real*8 mu1zx, u1zim2, u1zim1, u1zip1, u1zip2
      real*8 mu2zy, u2zjm2, u2zjm1, u2zjp1, u2zjp2
      real*8 r1, r2, r3, h, cof, d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )

      cof = 1d0/(h*h)

      k1 = kfirst+2
      if (onesided(5).eq.1) k1 = 7;
      k2 = klast-2
      if (onesided(6).eq.1) k2 = nz-6;
c the centered stencil can be evaluated 2 points away from the boundary
      do k=k1,k2
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))
c
            muz1 = mu(i,j,k-1)-tf*(mu(i,j,k)+mu(i,j,k-2))
            muz2 = mu(i,j,k-2)+mu(i,j,k+1)+3*(mu(i,j,k)+mu(i,j,k-1))
            muz3 = mu(i,j,k-1)+mu(i,j,k+2)+3*(mu(i,j,k+1)+mu(i,j,k))
            muz4 = mu(i,j,k+1)-tf*(mu(i,j,k)+mu(i,j,k+2))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) 
     *              + muz1*(u(1,i,j,k-2)-u(1,i,j,k)) + 
     *                muz2*(u(1,i,j,k-1)-u(1,i,j,k)) + 
     *                muz3*(u(1,i,j,k+1)-u(1,i,j,k)) +
     *                muz4*(u(1,i,j,k+2)-u(1,i,j,k)) )

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k))
     *              + muz1*(u(2,i,j,k-2)-u(2,i,j,k)) + 
     *                muz2*(u(2,i,j,k-1)-u(2,i,j,k)) + 
     *                muz3*(u(2,i,j,k+1)-u(2,i,j,k)) +
     *                muz4*(u(2,i,j,k+2)-u(2,i,j,k)) ) 

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) +
     *             (2*muz1+la(i,j,k-1)-tf*(la(i,j,k)+la(i,j,k-2)))*
     *                     (u(3,i,j,k-2)-u(3,i,j,k))+
     *      (2*muz2+la(i,j,k-2)+la(i,j,k+1)+3*(la(i,j,k)+la(i,j,k-1)))*
     *                     (u(3,i,j,k-1)-u(3,i,j,k))+ 
     *      (2*muz3+la(i,j,k-1)+la(i,j,k+2)+3*(la(i,j,k+1)+la(i,j,k)))*
     *                     (u(3,i,j,k+1)-u(3,i,j,k))+
     *             (2*muz4+la(i,j,k+1)-tf*(la(i,j,k)+la(i,j,k+2)))*
     *                     (u(3,i,j,k+2)-u(3,i,j,k)) )


*** Mixed derivatives:
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (la*w_z)_x
     *          + i144*( la(i-2,j,k)*(u(3,i-2,j,k-2)-u(3,i-2,j,k+2)+
     *                        8*(-u(3,i-2,j,k-1)+u(3,i-2,j,k+1))) - 8*(
     *                   la(i-1,j,k)*(u(3,i-1,j,k-2)-u(3,i-1,j,k+2)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i-1,j,k+1))) )+8*(
     *                   la(i+1,j,k)*(u(3,i+1,j,k-2)-u(3,i+1,j,k+2)+
     *                        8*(-u(3,i+1,j,k-1)+u(3,i+1,j,k+1))) ) - (
     *                   la(i+2,j,k)*(u(3,i+2,j,k-2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i+2,j,k-1)+u(3,i+2,j,k+1))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) )) 
***   (mu*w_x)_z
     *          + i144*( mu(i,j,k-2)*(u(3,i-2,j,k-2)-u(3,i+2,j,k-2)+
     *                        8*(-u(3,i-1,j,k-2)+u(3,i+1,j,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i-2,j,k-1)-u(3,i+2,j,k-1)+
     *                        8*(-u(3,i-1,j,k-1)+u(3,i+1,j,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i-2,j,k+1)-u(3,i+2,j,k+1)+
     *                        8*(-u(3,i-1,j,k+1)+u(3,i+1,j,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i-2,j,k+2)-u(3,i+2,j,k+2)+
     *                        8*(-u(3,i-1,j,k+2)+u(3,i+1,j,k+2))) )) 

***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) )) 
*** (la*w_z)_y
     *          + i144*( la(i,j-2,k)*(u(3,i,j-2,k-2)-u(3,i,j-2,k+2)+
     *                        8*(-u(3,i,j-2,k-1)+u(3,i,j-2,k+1))) - 8*(
     *                   la(i,j-1,k)*(u(3,i,j-1,k-2)-u(3,i,j-1,k+2)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j-1,k+1))) )+8*(
     *                   la(i,j+1,k)*(u(3,i,j+1,k-2)-u(3,i,j+1,k+2)+
     *                        8*(-u(3,i,j+1,k-1)+u(3,i,j+1,k+1))) ) - (
     *                   la(i,j+2,k)*(u(3,i,j+2,k-2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j+2,k-1)+u(3,i,j+2,k+1))) ))
*** (mu*w_y)_z
     *          + i144*( mu(i,j,k-2)*(u(3,i,j-2,k-2)-u(3,i,j+2,k-2)+
     *                        8*(-u(3,i,j-1,k-2)+u(3,i,j+1,k-2))) - 8*(
     *                   mu(i,j,k-1)*(u(3,i,j-2,k-1)-u(3,i,j+2,k-1)+
     *                        8*(-u(3,i,j-1,k-1)+u(3,i,j+1,k-1))) )+8*(
     *                   mu(i,j,k+1)*(u(3,i,j-2,k+1)-u(3,i,j+2,k+1)+
     *                        8*(-u(3,i,j-1,k+1)+u(3,i,j+1,k+1))) ) - (
     *                   mu(i,j,k+2)*(u(3,i,j-2,k+2)-u(3,i,j+2,k+2)+
     *                        8*(-u(3,i,j-1,k+2)+u(3,i,j+1,k+2))) )) 
***  (mu*u_z)_x
            r3 = r3 + i144*( mu(i-2,j,k)*(u(1,i-2,j,k-2)-u(1,i-2,j,k+2)+
     *                        8*(-u(1,i-2,j,k-1)+u(1,i-2,j,k+1))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j,k-2)-u(1,i-1,j,k+2)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i-1,j,k+1))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j,k-2)-u(1,i+1,j,k+2)+
     *                        8*(-u(1,i+1,j,k-1)+u(1,i+1,j,k+1))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j,k-2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i+2,j,k-1)+u(1,i+2,j,k+1))) )) 
*** (mu*v_z)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i,j-2,k-2)-u(2,i,j-2,k+2)+
     *                        8*(-u(2,i,j-2,k-1)+u(2,i,j-2,k+1))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i,j-1,k-2)-u(2,i,j-1,k+2)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j-1,k+1))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i,j+1,k-2)-u(2,i,j+1,k+2)+
     *                        8*(-u(2,i,j+1,k-1)+u(2,i,j+1,k+1))) ) - (
     *                   mu(i,j+2,k)*(u(2,i,j+2,k-2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j+2,k-1)+u(2,i,j+2,k+1))) ))
***   (la*u_x)_z
     *          + i144*( la(i,j,k-2)*(u(1,i-2,j,k-2)-u(1,i+2,j,k-2)+
     *                        8*(-u(1,i-1,j,k-2)+u(1,i+1,j,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(1,i-2,j,k-1)-u(1,i+2,j,k-1)+
     *                        8*(-u(1,i-1,j,k-1)+u(1,i+1,j,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(1,i-2,j,k+1)-u(1,i+2,j,k+1)+
     *                        8*(-u(1,i-1,j,k+1)+u(1,i+1,j,k+1))) ) - (
     *                   la(i,j,k+2)*(u(1,i-2,j,k+2)-u(1,i+2,j,k+2)+
     *                        8*(-u(1,i-1,j,k+2)+u(1,i+1,j,k+2))) )) 
*** (la*v_y)_z
     *          + i144*( la(i,j,k-2)*(u(2,i,j-2,k-2)-u(2,i,j+2,k-2)+
     *                        8*(-u(2,i,j-1,k-2)+u(2,i,j+1,k-2))) - 8*(
     *                   la(i,j,k-1)*(u(2,i,j-2,k-1)-u(2,i,j+2,k-1)+
     *                        8*(-u(2,i,j-1,k-1)+u(2,i,j+1,k-1))) )+8*(
     *                   la(i,j,k+1)*(u(2,i,j-2,k+1)-u(2,i,j+2,k+1)+
     *                        8*(-u(2,i,j-1,k+1)+u(2,i,j+1,k+1))) ) - (
     *                   la(i,j,k+2)*(u(2,i,j-2,k+2)-u(2,i,j+2,k+2)+
     *                        8*(-u(2,i,j-1,k+2)+u(2,i,j+1,k+2))) )) 

            uacc(1,i,j,k) = cof*r1
            uacc(2,i,j,k) = cof*r2
            uacc(3,i,j,k) = cof*r3

            enddo
         enddo
      enddo

c low-k boundary modified stencils
      if (onesided(5).eq.1) then
      do k=1,6
c the centered stencil can be used in the x- and y-directions
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) )
c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            do q=1,8
              mucof(q)=0
              do m=1,8
                mucof(q) = mucof(q)+acof(k,q,m)*mu(i,j,m)
              enddo
            end do
c computing the second derivative
            mu1zz = 0
            do q=1,8
              mu1zz = mu1zz + mucof(q)*u(1,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(k)*mu(i,j,1)*u(1,i,j,0))

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) )
c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
            mu2zz = 0
            do q=1,8
              mu2zz = mu2zz + mucof(q)*u(2,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(k)*mu(i,j,1)*u(2,i,j,0))

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
            do q=1,8
              lap2mu(q)=0
              do m=1,8
                lap2mu(q) = lap2mu(q)+acof(k,q,m)*
     +                                (la(i,j,m)+2*mu(i,j,m))
              enddo
            end do
c computing the second derivative
            mu3zz = 0
            do q=1,8
              mu3zz = mu3zz + lap2mu(q)*u(3,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(k)*(la(i,j,1)+2*mu(i,j,1))*
     +           u(3,i,j,0))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) ))
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do q=1,8
              u3zip2 = u3zip2 + bope(k,q)*u(3,i+2,j,q)
              u3zip1 = u3zip1 + bope(k,q)*u(3,i+1,j,q)
              u3zim1 = u3zim1 + bope(k,q)*u(3,i-1,j,q)
              u3zim2 = u3zim2 + bope(k,q)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do q=1,8
              mu3xz = mu3xz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) ))
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do q=1,8
              u3zjp2 = u3zjp2 + bope(k,q)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 + bope(k,q)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 + bope(k,q)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 + bope(k,q)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do q=1,8
              mu3yz = mu3yz + bope(k,q)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do q=1,8
              u1zip2 = u1zip2 + bope(k,q)*u(1,i+2,j,q)
              u1zip1 = u1zip1 + bope(k,q)*u(1,i+1,j,q)
              u1zim1 = u1zim1 + bope(k,q)*u(1,i-1,j,q)
              u1zim2 = u1zim2 + bope(k,q)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do q=1,8
              u2zjp2 = u2zjp2 + bope(k,q)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 + bope(k,q)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 + bope(k,q)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 + bope(k,q)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do q=1,8
              lau1xz = lau1xz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do q=1,8
              lau2yz = lau2yz + bope(k,q)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + lau2yz

            uacc(1,i,j,k) = cof*r1
            uacc(2,i,j,k) = cof*r2
            uacc(3,i,j,k) = cof*r3

            enddo
         enddo
      enddo
      endif

c high-k boundary
      if (onesided(6).eq.1) then
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
c from inner_loop_4a
            mux1 = mu(i-1,j,k)-tf*(mu(i,j,k)+mu(i-2,j,k))
            mux2 = mu(i-2,j,k)+mu(i+1,j,k)+3*(mu(i,j,k)+mu(i-1,j,k))
            mux3 = mu(i-1,j,k)+mu(i+2,j,k)+3*(mu(i+1,j,k)+mu(i,j,k))
            mux4 = mu(i+1,j,k)-tf*(mu(i,j,k)+mu(i+2,j,k))
c
            muy1 = mu(i,j-1,k)-tf*(mu(i,j,k)+mu(i,j-2,k))
            muy2 = mu(i,j-2,k)+mu(i,j+1,k)+3*(mu(i,j,k)+mu(i,j-1,k))
            muy3 = mu(i,j-1,k)+mu(i,j+2,k)+3*(mu(i,j+1,k)+mu(i,j,k))
            muy4 = mu(i,j+1,k)-tf*(mu(i,j,k)+mu(i,j+2,k))

*** xx, yy, and zz derivatives:
c note that we could have introduced intermediate variables for the average of lambda 
c in the same way as we did for mu
            r1 = i6*((2*mux1+la(i-1,j,k)-tf*(la(i,j,k)+la(i-2,j,k)))*
     *                         (u(1,i-2,j,k)-u(1,i,j,k))+
     *      (2*mux2+la(i-2,j,k)+la(i+1,j,k)+3*(la(i,j,k)+la(i-1,j,k)))*
     *                         (u(1,i-1,j,k)-u(1,i,j,k))+ 
     *      (2*mux3+la(i-1,j,k)+la(i+2,j,k)+3*(la(i+1,j,k)+la(i,j,k)))*
     *                         (u(1,i+1,j,k)-u(1,i,j,k))+
     *           (2*mux4+ la(i+1,j,k)-tf*(la(i,j,k)+la(i+2,j,k)))*
     *           (u(1,i+2,j,k)-u(1,i,j,k)) 
     *              + muy1*(u(1,i,j-2,k)-u(1,i,j,k)) + 
     *                muy2*(u(1,i,j-1,k)-u(1,i,j,k)) + 
     *                muy3*(u(1,i,j+1,k)-u(1,i,j,k)) +
     *                muy4*(u(1,i,j+2,k)-u(1,i,j,k)) )
c all indices ending with 'b' are indices relative to the boundary, going into the domain (1,2,3,...)
            kb = nz-k+1
c all coefficient arrays (acof, bope, ghcof) should be indexed with these indices
c all solution and material property arrays should be indexed with (i,j,k)

c (mu*uz)_z can not be centered
c second derivative (mu*u_z)_z at grid point z_k
c averaging the coefficient
            do qb=1,8
              mucof(qb)=0
              do mb=1,8
                m = nz-mb+1
                mucof(qb) = mucof(qb)+acof(kb,qb,mb)*mu(i,j,m)
              enddo
            end do
c computing the second derivative
            mu1zz = 0
            do qb=1,8
              q = nz-qb+1
              mu1zz = mu1zz + mucof(qb)*u(1,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r1 = r1 + (mu1zz + ghcof(kb)*mu(i,j,nz)*u(1,i,j,nz+1))

            r2 = i6*(mux1*(u(2,i-2,j,k)-u(2,i,j,k)) + 
     *                 mux2*(u(2,i-1,j,k)-u(2,i,j,k)) + 
     *                 mux3*(u(2,i+1,j,k)-u(2,i,j,k)) +
     *                 mux4*(u(2,i+2,j,k)-u(2,i,j,k)) + 
     *             (2*muy1+la(i,j-1,k)-tf*(la(i,j,k)+la(i,j-2,k)))*
     *                     (u(2,i,j-2,k)-u(2,i,j,k))+
     *      (2*muy2+la(i,j-2,k)+la(i,j+1,k)+3*(la(i,j,k)+la(i,j-1,k)))*
     *                     (u(2,i,j-1,k)-u(2,i,j,k))+ 
     *      (2*muy3+la(i,j-1,k)+la(i,j+2,k)+3*(la(i,j+1,k)+la(i,j,k)))*
     *                     (u(2,i,j+1,k)-u(2,i,j,k))+
     *             (2*muy4+la(i,j+1,k)-tf*(la(i,j,k)+la(i,j+2,k)))*
     *                     (u(2,i,j+2,k)-u(2,i,j,k)) )
c (mu*vz)_z can not be centered
c second derivative (mu*v_z)_z at grid point z_k
c averaging the coefficient: already done above
c computing the second derivative
            mu2zz = 0
            do qb=1,8
              q = nz-qb+1
              mu2zz = mu2zz + mucof(qb)*u(2,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r2 = r2 + (mu2zz + ghcof(kb)*mu(i,j,nz)*u(2,i,j,nz+1))

            r3 = i6*(mux1*(u(3,i-2,j,k)-u(3,i,j,k)) + 
     *                 mux2*(u(3,i-1,j,k)-u(3,i,j,k)) + 
     *                 mux3*(u(3,i+1,j,k)-u(3,i,j,k)) +
     *                 mux4*(u(3,i+2,j,k)-u(3,i,j,k))  
     *              + muy1*(u(3,i,j-2,k)-u(3,i,j,k)) + 
     *                muy2*(u(3,i,j-1,k)-u(3,i,j,k)) + 
     *                muy3*(u(3,i,j+1,k)-u(3,i,j,k)) +
     *                muy4*(u(3,i,j+2,k)-u(3,i,j,k)) )
c ((2*mu+lambda)*w_z)_z can not be centered
c averaging the coefficient
            do qb=1,8
              lap2mu(qb)=0
              do mb=1,8
                m = nz-mb+1
                lap2mu(qb) = lap2mu(qb)+acof(kb,qb,mb)*
     +                                (la(i,j,m)+2*mu(i,j,m))
              enddo
            end do
c computing the second derivative
            mu3zz = 0
            do qb=1,8
              q = nz-qb+1
              mu3zz = mu3zz + lap2mu(qb)*u(3,i,j,q)
            enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
            r3 = r3 + (mu3zz + ghcof(kb)*(la(i,j,nz)+2*mu(i,j,nz))*
     +           u(3,i,j,nz+1))

c General formula for the first derivative of u1 at z_k
c$$$        u1z=0
c$$$        do q=1,8
c$$$          u1z = u1z + bope(k,q)*u1(q)
c$$$        enddo

c cross-terms in first component of rhs
***   (la*v_y)_x
            r1 = r1 + i144*( la(i-2,j,k)*(u(2,i-2,j-2,k)-u(2,i-2,j+2,k)+
     *                        8*(-u(2,i-2,j-1,k)+u(2,i-2,j+1,k))) - 8*(
     *                   la(i-1,j,k)*(u(2,i-1,j-2,k)-u(2,i-1,j+2,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i-1,j+1,k))) )+8*(
     *                   la(i+1,j,k)*(u(2,i+1,j-2,k)-u(2,i+1,j+2,k)+
     *                        8*(-u(2,i+1,j-1,k)+u(2,i+1,j+1,k))) ) - (
     *                   la(i+2,j,k)*(u(2,i+2,j-2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i+2,j-1,k)+u(2,i+2,j+1,k))) )) 
***   (mu*v_x)_y
     *          + i144*( mu(i,j-2,k)*(u(2,i-2,j-2,k)-u(2,i+2,j-2,k)+
     *                        8*(-u(2,i-1,j-2,k)+u(2,i+1,j-2,k))) - 8*(
     *                   mu(i,j-1,k)*(u(2,i-2,j-1,k)-u(2,i+2,j-1,k)+
     *                        8*(-u(2,i-1,j-1,k)+u(2,i+1,j-1,k))) )+8*(
     *                   mu(i,j+1,k)*(u(2,i-2,j+1,k)-u(2,i+2,j+1,k)+
     *                        8*(-u(2,i-1,j+1,k)+u(2,i+1,j+1,k))) ) - (
     *                   mu(i,j+2,k)*(u(2,i-2,j+2,k)-u(2,i+2,j+2,k)+
     *                        8*(-u(2,i-1,j+2,k)+u(2,i+1,j+2,k))) ))
***   (la*w_z)_x: NOT CENTERED
            u3zip2=0
            u3zip1=0
            u3zim1=0
            u3zim2=0
            do qb=1,8
              q = nz-qb+1
              u3zip2 = u3zip2 - bope(kb,qb)*u(3,i+2,j,q)
              u3zip1 = u3zip1 - bope(kb,qb)*u(3,i+1,j,q)
              u3zim1 = u3zim1 - bope(kb,qb)*u(3,i-1,j,q)
              u3zim2 = u3zim2 - bope(kb,qb)*u(3,i-2,j,q)
            enddo
            lau3zx= i12*(-la(i+2,j,k)*u3zip2 + 8*la(i+1,j,k)*u3zip1
     +                   -8*la(i-1,j,k)*u3zim1 + la(i-2,j,k)*u3zim2)
            r1 = r1 + lau3zx

***   (mu*w_x)_z: NOT CENTERED
            mu3xz=0
            do qb=1,8
              q = nz-qb+1
              mu3xz = mu3xz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i+2,j,q) + 8*u(3,i+1,j,q)
     +              -8*u(3,i-1,j,q) + u(3,i-2,j,q)) )
            enddo

            r1 = r1 + mu3xz

c cross-terms in second component of rhs
***   (mu*u_y)_x
            r2 = r2 + i144*( mu(i-2,j,k)*(u(1,i-2,j-2,k)-u(1,i-2,j+2,k)+
     *                        8*(-u(1,i-2,j-1,k)+u(1,i-2,j+1,k))) - 8*(
     *                   mu(i-1,j,k)*(u(1,i-1,j-2,k)-u(1,i-1,j+2,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i-1,j+1,k))) )+8*(
     *                   mu(i+1,j,k)*(u(1,i+1,j-2,k)-u(1,i+1,j+2,k)+
     *                        8*(-u(1,i+1,j-1,k)+u(1,i+1,j+1,k))) ) - (
     *                   mu(i+2,j,k)*(u(1,i+2,j-2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i+2,j-1,k)+u(1,i+2,j+1,k))) )) 
*** (la*u_x)_y
     *          + i144*( la(i,j-2,k)*(u(1,i-2,j-2,k)-u(1,i+2,j-2,k)+
     *                        8*(-u(1,i-1,j-2,k)+u(1,i+1,j-2,k))) - 8*(
     *                   la(i,j-1,k)*(u(1,i-2,j-1,k)-u(1,i+2,j-1,k)+
     *                        8*(-u(1,i-1,j-1,k)+u(1,i+1,j-1,k))) )+8*(
     *                   la(i,j+1,k)*(u(1,i-2,j+1,k)-u(1,i+2,j+1,k)+
     *                        8*(-u(1,i-1,j+1,k)+u(1,i+1,j+1,k))) ) - (
     *                   la(i,j+2,k)*(u(1,i-2,j+2,k)-u(1,i+2,j+2,k)+
     *                        8*(-u(1,i-1,j+2,k)+u(1,i+1,j+2,k))) ))
*** (la*w_z)_y : NOT CENTERED
            u3zjp2=0
            u3zjp1=0
            u3zjm1=0
            u3zjm2=0
            do qb=1,8
              q = nz-qb+1
              u3zjp2 = u3zjp2 - bope(kb,qb)*u(3,i,j+2,q)
              u3zjp1 = u3zjp1 - bope(kb,qb)*u(3,i,j+1,q)
              u3zjm1 = u3zjm1 - bope(kb,qb)*u(3,i,j-1,q)
              u3zjm2 = u3zjm2 - bope(kb,qb)*u(3,i,j-2,q)
            enddo
            lau3zy= i12*(-la(i,j+2,k)*u3zjp2 + 8*la(i,j+1,k)*u3zjp1
     +                   -8*la(i,j-1,k)*u3zjm1 + la(i,j-2,k)*u3zjm2)

            r2 = r2 + lau3zy

*** (mu*w_y)_z: NOT CENTERED
            mu3yz=0
            do qb=1,8
              q = nz-qb+1
              mu3yz = mu3yz - bope(kb,qb)*( mu(i,j,q)*i12*
     +             (-u(3,i,j+2,q) + 8*u(3,i,j+1,q)
     +              -8*u(3,i,j-1,q) + u(3,i,j-2,q)) )
            enddo

            r2 = r2 + mu3yz

c No centered cross terms in r3
***  (mu*u_z)_x: NOT CENTERED
            u1zip2=0
            u1zip1=0
            u1zim1=0
            u1zim2=0
            do qb=1,8
              q = nz-qb+1
              u1zip2 = u1zip2 - bope(kb,qb)*u(1,i+2,j,q)
              u1zip1 = u1zip1 - bope(kb,qb)*u(1,i+1,j,q)
              u1zim1 = u1zim1 - bope(kb,qb)*u(1,i-1,j,q)
              u1zim2 = u1zim2 - bope(kb,qb)*u(1,i-2,j,q)
            enddo
            mu1zx= i12*(-mu(i+2,j,k)*u1zip2 + 8*mu(i+1,j,k)*u1zip1
     +                   -8*mu(i-1,j,k)*u1zim1 + mu(i-2,j,k)*u1zim2)
            r3 = r3 + mu1zx

*** (mu*v_z)_y: NOT CENTERED
            u2zjp2=0
            u2zjp1=0
            u2zjm1=0
            u2zjm2=0
            do qb=1,8
              q = nz-qb+1
              u2zjp2 = u2zjp2 - bope(kb,qb)*u(2,i,j+2,q)
              u2zjp1 = u2zjp1 - bope(kb,qb)*u(2,i,j+1,q)
              u2zjm1 = u2zjm1 - bope(kb,qb)*u(2,i,j-1,q)
              u2zjm2 = u2zjm2 - bope(kb,qb)*u(2,i,j-2,q)
            enddo
            mu2zy= i12*(-mu(i,j+2,k)*u2zjp2 + 8*mu(i,j+1,k)*u2zjp1
     +                   -8*mu(i,j-1,k)*u2zjm1 + mu(i,j-2,k)*u2zjm2)
            r3 = r3 + mu2zy

***   (la*u_x)_z: NOT CENTERED
            lau1xz=0
            do qb=1,8
              q = nz-qb+1
              lau1xz = lau1xz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(1,i+2,j,q) + 8*u(1,i+1,j,q)
     +              -8*u(1,i-1,j,q) + u(1,i-2,j,q)) )
            enddo

            r3 = r3 + lau1xz

*** (la*v_y)_z: NOT CENTERED
            lau2yz=0
            do qb=1,8
              q = nz-qb+1
              lau2yz = lau2yz - bope(kb,qb)*( la(i,j,q)*i12*
     +             (-u(2,i,j+2,q) + 8*u(2,i,j+1,q)
     +              -8*u(2,i,j-1,q) + u(2,i,j-2,q)) )
            enddo
            r3 = r3 + lau2yz

            uacc(1,i,j,k) = cof*r1
            uacc(2,i,j,k) = cof*r2
            uacc(3,i,j,k) = cof*r3

            enddo
         enddo
      enddo
      endif

      end

c-----------------------------------------------------------------------
      subroutine BC( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nx, ny, nz,
     +     u, h, bccnd, sbop, mu, la, t,
     *     bforce, om, ph, cv )
      implicit none
      real*8 d4a, d4b
      parameter( d4a=2d0/3, d4b=-1d0/12 )
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nx, ny, nz
      integer s, wind(6), i, j, k, bccnd(6), w, kl
      real*8 x, y, z, h, sbop(0:4)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 bforce(3,ifirst:ilast,jfirst:jlast)
      real*8 ux, vy, wx, wy, uz, vz, wz, t, om, ph, cv
      do s=1,6
        wind(1) = ifirst
        wind(2) = ilast
        wind(3) = jfirst
        wind(4) = jlast
        wind(5) = kfirst
        wind(6) = klast
c by default, always assign 2 ghost points
        if( s.eq.1 )then
          wind(2) = ifirst+1
        elseif( s.eq.2 )then
          wind(1) = ilast-1
        elseif( s.eq.3 )then
          wind(4) = jfirst+1
        elseif( s.eq.4 )then
          wind(3) = jlast-1
        elseif( s.eq.5 )then
          wind(6) = kfirst+1
        else
          wind(5) = klast-1
        endif

*** Homogeneous dirichlet condition
        if( bccnd(s).eq.1 )then
          do k=wind(5),wind(6)
            do j=wind(3),wind(4)
              do i=wind(1),wind(2)
                u(1,i,j,k) =0
                u(2,i,j,k) =0
                u(3,i,j,k) =0
              enddo
            enddo
          enddo
        elseif( bccnd(s).eq.2 )then
*** Twilight forced dirichlet condition
          do k=wind(5),wind(6)
            z = (k-1)*h
            do j=wind(3),wind(4)
              y = (j-1)*h
              do i=wind(1),wind(2)
                x = (i-1)*h
            u(1,i,j,k) = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph)
            u(2,i,j,k) = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph)
            u(3,i,j,k) = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t))
              enddo
            enddo
          enddo
*** Free surface condition, 3-homogeneous, 4-Twilight forced
        elseif( bccnd(s).eq.3 .or. bccnd(s).eq.4 )then
          if( s.ne.5 .and. s.ne.6 )then
            write(*,*) 
     *      'ERROR: Free surface condition not implemented for side ', s
            stop
          endif

c calculate boundary forcing
          if( bccnd(s).eq.4 )then
            if( s.eq.5 )then
              k = 1
            else
              k = nz
            endif
c update this call
            call TWFRSURFZ( ifirst, ilast, jfirst, jlast, kfirst, klast,
     +           nx, ny, nz, h, k, t, om, cv, ph, bforce, mu, la )
          else
            do j=jfirst+2,jlast-2
              do i=ifirst+2,ilast-2
                bforce(1,i,j) =0
                bforce(2,i,j) =0
                bforce(3,i,j) =0
              enddo
            enddo
          endif

**** Do the free surface condition (kl is the direction)
          if( s.eq.5 )then
            k = 1
            kl= 1
          else
            k = nz
            kl= -1
          endif

          do j=jfirst+2,jlast-2
            do i=ifirst+2,ilast-2
*** Compute wx, wy, ux, vy by centered differences along boundary
              wx = d4a*(u(3,i+1,j,k)-u(3,i-1,j,k))+
     *             d4b*(u(3,i+2,j,k)-u(3,i-2,j,k))
              ux = d4a*(u(1,i+1,j,k)-u(1,i-1,j,k))+
     *             d4b*(u(1,i+2,j,k)-u(1,i-2,j,k))

              wy = d4a*(u(3,i,j+1,k)-u(3,i,j-1,k))+
     *             d4b*(u(3,i,j+2,k)-u(3,i,j-2,k))
              vy = d4a*(u(2,i,j+1,k)-u(2,i,j-1,k))+
     *             d4b*(u(2,i,j+2,k)-u(2,i,j-2,k))

              uz = 0
              vz = 0
              wz = 0
c interior contribution to uz, vz, wz (kl is the direction)
              do w=1,4
                uz = uz + sbop(w)*u(1,i,j,k+kl*(w-1))
                vz = vz + sbop(w)*u(2,i,j,k+kl*(w-1))
                wz = wz + sbop(w)*u(3,i,j,k+kl*(w-1))
              enddo
              u(1,i,j,k-kl)=(-uz-kl*wx+kl*h*bforce(1,i,j)/mu(i,j,k))
     *             /sbop(0)
              u(2,i,j,k-kl)=(-vz-kl*wy+kl*h*bforce(2,i,j)/mu(i,j,k))
     *             /sbop(0)
              u(3,i,j,k-kl)=(-wz+(-kl*la(i,j,k)*(ux+vy)+
     *             kl*h*bforce(3,i,j))/(2*mu(i,j,k)+la(i,j,k)))
     *             /sbop(0)
            enddo
          enddo
        endif
      enddo
      end

c-----------------------------------------------------------------------
      subroutine predictor(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, lu, fo, rho, dt2 )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt2
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
c 2nd order accurate predictor of solution at t+dt
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              up(c,i,j,k) = 2*u(c,i,j,k) - um(c,i,j,k) + 
     +             dt2*(lu(c,i,j,k) + fo(c,i,j,k))/rho(i,j,k)
            enddo
          enddo
        enddo
      enddo
      return 
      end

c-----------------------------------------------------------------------
      subroutine corrector(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, lu, fo, rho, dt4 )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt4, i12
      parameter( i12=1d0/12 )
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
c correct solution to 4th order accuracy
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              up(c,i,j,k) = up(c,i,j,k) + i12*dt4*
     +             (lu(c,i,j,k) + fo(c,i,j,k))/rho(i,j,k)
            enddo
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine dpdmInTime(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     up, u, um, u2, dt2i)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 dt2i
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c evaluate 2nd divided time difference D+D-(u)
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              u2(c,i,j,k) = dt2i*
     +             (up(c,i,j,k) - 2*u(c,i,j,k) + um(c,i,j,k))
            enddo
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine cycleArrays(ifirst, ilast, jfirst, jlast, kfirst, klast
     +     , up, u, um )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k, c
      real*8 up(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
c cycle solution arrays (this could be done more efficiently with handles in C)
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            do c=1,3
              um(c,i,j,k) = u(c,i,j,k)
              u(c,i,j,k)  = up(c,i,j,k)
            enddo
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine exactSol( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     u, t, om, cv, ph, h )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
* v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
* w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, cv, ph
      doubleprecision ampmu, amplambda, h
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      do k=kfirst,klast
        z = (k-1)*h
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            u(1,i,j,k) = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph)
            u(2,i,j,k) = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph)
            u(3,i,j,k) = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t))
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine exactMat( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     rho, mu, la, omm, phm, amprho, ampmu, amplambda, h )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* rho    := amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
* mu     := ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*sin(omm*z+phm) );
* lambda := amplambda*(2 + sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) );
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda, h
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)

      do k=kfirst,klast
        z = (k-1)*h
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            rho(i,j,k) = amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*
     +           sin(omm*z+phm) );
            mu(i,j,k)  = ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*
     +           sin(omm*z+phm) )
            la(i,j,k)  = amplambda*(2 + 
     +           sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) )
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for right hand side forcing
c corresponding to twilight exact solution
c
      subroutine exactAcc( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     utt, t, om, c, ph, h )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
* v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
* w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph
      doubleprecision ampmu, amplambda, h
      real*8 utt(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision acc(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t14
      doubleprecision t19
      doubleprecision t22
      doubleprecision t30
      doubleprecision t4
      doubleprecision t5
      doubleprecision t7

      do k=kfirst,klast
        z = (k-1)*h
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
        t1 = c*t
        t4 = sin(om*(x-t1))
        t5 = om**2
        t7 = c**2
        t10 = sin(om*y+ph)
        t14 = sin(om*z+ph)
        acc(1) = -t4*t5*t7*t10*t14
        t19 = sin(om*x+ph)
        t22 = sin(om*(y-t1))
        acc(2) = -t19*t22*t5*t7*t14
        t30 = sin(om*(z-t1))
        acc(3) = -t19*t10*t30*t5*t7

        utt(1,i,j,k) = acc(1)
        utt(2,i,j,k) = acc(2)
        utt(3,i,j,k) = acc(3)
        enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for right hand side forcing
c corresponding to twilight exact solution
c
      subroutine exactRhs3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, h)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t101
      doubleprecision t102
      doubleprecision t103
      doubleprecision t105
      doubleprecision t106
      doubleprecision t108
      doubleprecision t11
      doubleprecision t114
      doubleprecision t115
      doubleprecision t118
      doubleprecision t119
      doubleprecision t12
      doubleprecision t120
      doubleprecision t122
      doubleprecision t125
      doubleprecision t129
      doubleprecision t135
      doubleprecision t143
      doubleprecision t15
      doubleprecision t150
      doubleprecision t152
      doubleprecision t159
      doubleprecision t16
      doubleprecision t166
      doubleprecision t168
      doubleprecision t17
      doubleprecision t173
      doubleprecision t175
      doubleprecision t18
      doubleprecision t2
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t38
      doubleprecision t4
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t47
      doubleprecision t48
      doubleprecision t51
      doubleprecision t53
      doubleprecision t54
      doubleprecision t55
      doubleprecision t57
      doubleprecision t58
      doubleprecision t6
      doubleprecision t60
      doubleprecision t63
      doubleprecision t64
      doubleprecision t69
      doubleprecision t7
      doubleprecision t70
      doubleprecision t71
      doubleprecision t72
      doubleprecision t77
      doubleprecision t78
      doubleprecision t79
      doubleprecision t8
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t83
      doubleprecision t85
      doubleprecision t87
      doubleprecision t89
      doubleprecision t92
      doubleprecision t95
      doubleprecision t97

      do k=kfirst,klast
         z = (k-1)*h
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t2 = omm*x+phm
        t3 = sin(t2)
        t4 = ampmu*t3
        t6 = omm*y+phm
        t7 = sin(t6)
        t8 = omm*t7
        t10 = omm*z+phm
        t11 = sin(t10)
        t12 = t8*t11
        t15 = cos(t2)
        t16 = amplambda*t15
        t17 = cos(t10)
        t18 = t8*t17
        t21 = c*t
        t23 = om*(x-t21)
        t24 = cos(t23)
        t27 = om*y+ph
        t28 = sin(t27)
        t31 = om*z+ph
        t32 = sin(t31)
        t38 = ampmu*(3+t15*t7*t11)
        t43 = amplambda*(2+t3*t7*t17)
        t44 = 2*t38+t43
        t45 = sin(t23)
        t47 = om**2
        t48 = t47*t28
        t51 = t16*t8
        t53 = om*x+ph
        t54 = sin(t53)
        t55 = t17*t54
        t57 = om*(y-t21)
        t58 = cos(t57)
        t60 = t58*om*t32
        t63 = cos(t53)
        t64 = t43*t63
        t69 = om*(z-t21)
        t70 = cos(t69)
        t71 = t28*t70
        t72 = t71*om
        t77 = ampmu*t15
        t78 = cos(t6)
        t79 = t77*t78
        t80 = omm*t11
        t81 = t63*om
        t82 = sin(t57)
        t83 = t82*t32
        t85 = cos(t27)
        t87 = om*t32
        t89 = t81*t83+t45*t85*t87
        t92 = t63*t47
        t95 = t45*t28
        t97 = t95*t47*t32
        t100 = t77*omm
        t101 = t7*t17
        t102 = sin(t69)
        t103 = t28*t102
        t105 = cos(t31)
        t106 = t105*om
        t108 = t81*t103+t95*t106
        forces(1) = (-2*t4*t12+t16*t18)*t24*om*t28*t32-t44*t45*t48*t32+t
     #51*t55*t60+t64*t47*t58*t32+t51*t55*t72+t64*t48*t70+t79*t80*t89+t38
     #*(t92*t58*t32-t97)+t100*t101*t108+t38*(t92*t71-t97)
        t114 = t4*omm
        t115 = t7*t11
        t118 = t54*t47
        t119 = t118*t83
        t120 = t24*t47
        t122 = t120*t85*t32
        t125 = t78*omm
        t129 = amplambda*t3
        t135 = t44*t54
        t143 = t24*om*t28*t32
        t150 = t54*t85
        t152 = t150*t47*t70
        t159 = t150*om*t102+t54*t82*t106
        forces(2) = -t114*t115*t89+t38*(-t119+t122)+(2*t77*t125*t11+t129
     #*t125*t17)*t54*t60-t135*t82*t47*t32+t129*t78*omm*t17*(t143+t54*t28
     #*t70*om)+t43*(t122+t152)+t100*t101*t159+t38*(t152-t119)
        t166 = t118*t103
        t168 = t120*t28*t105
        t173 = t54*t58
        t175 = t173*t47*t105
        forces(3) = -t114*t115*t108+t38*(-t166+t168)+t79*t80*t159+t38*(-
     #t166+t175)+(2*t77*t18-t129*t12)*t54*t72-t135*t103*t47-t129*omm*t11
     #5*(t143+t173*t87)+t43*(t168+t175)

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for twilight forcing functions for the elastic wave
c equation, corresponding to exactSol and exactMat
c
      subroutine FORCING( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, h)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t102
      doubleprecision t105
      doubleprecision t107
      doubleprecision t110
      doubleprecision t111
      doubleprecision t112
      doubleprecision t113
      doubleprecision t115
      doubleprecision t116
      doubleprecision t118
      doubleprecision t124
      doubleprecision t125
      doubleprecision t129
      doubleprecision t13
      doubleprecision t130
      doubleprecision t133
      doubleprecision t134
      doubleprecision t135
      doubleprecision t137
      doubleprecision t14
      doubleprecision t140
      doubleprecision t144
      doubleprecision t150
      doubleprecision t156
      doubleprecision t16
      doubleprecision t163
      doubleprecision t165
      doubleprecision t17
      doubleprecision t172
      doubleprecision t181
      doubleprecision t183
      doubleprecision t188
      doubleprecision t19
      doubleprecision t190
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t26
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t39
      doubleprecision t40
      doubleprecision t43
      doubleprecision t5
      doubleprecision t51
      doubleprecision t56
      doubleprecision t57
      doubleprecision t59
      doubleprecision t6
      doubleprecision t62
      doubleprecision t64
      doubleprecision t65
      doubleprecision t66
      doubleprecision t68
      doubleprecision t69
      doubleprecision t71
      doubleprecision t74
      doubleprecision t75
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t83
      doubleprecision t88
      doubleprecision t89
      doubleprecision t9
      doubleprecision t90
      doubleprecision t91
      doubleprecision t92
      doubleprecision t93
      doubleprecision t95
      doubleprecision t97
      doubleprecision t99

      do k=kfirst,klast
         z = (k-1)*h
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t2 = omm*x+phm
        t3 = sin(t2)
        t5 = omm*y+phm
        t6 = cos(t5)
        t9 = omm*z+phm
        t10 = sin(t9)
        t13 = amprho*(2+t3*t6*t10)
        t14 = c*t
        t16 = om*(x-t14)
        t17 = sin(t16)
        t19 = om**2
        t20 = c**2
        t21 = t19*t20
        t23 = om*y+ph
        t24 = sin(t23)
        t26 = om*z+ph
        t27 = sin(t26)
        t28 = t24*t27
        t31 = ampmu*t3
        t32 = sin(t5)
        t33 = omm*t32
        t34 = t33*t10
        t37 = cos(t2)
        t38 = amplambda*t37
        t39 = cos(t9)
        t40 = t33*t39
        t43 = cos(t16)
        t51 = ampmu*(3+t37*t32*t10)
        t56 = amplambda*(2+t3*t32*t39)
        t57 = 2*t51+t56
        t59 = t19*t24
        t62 = t38*t33
        t64 = om*x+ph
        t65 = sin(t64)
        t66 = t39*t65
        t68 = om*(y-t14)
        t69 = cos(t68)
        t71 = t69*om*t27
        t74 = cos(t64)
        t75 = t56*t74
        t80 = om*(z-t14)
        t81 = cos(t80)
        t82 = t24*t81
        t83 = t82*om
        t88 = ampmu*t37
        t89 = t88*t6
        t90 = omm*t10
        t91 = t74*om
        t92 = sin(t68)
        t93 = t92*t27
        t95 = cos(t23)
        t97 = om*t27
        t99 = t91*t93+t17*t95*t97
        t102 = t74*t19
        t105 = t17*t24
        t107 = t105*t19*t27
        t110 = t88*omm
        t111 = t32*t39
        t112 = sin(t80)
        t113 = t24*t112
        t115 = cos(t26)
        t116 = t115*om
        t118 = t91*t113+t105*t116
        forces(1) = -t13*t17*t21*t28-(-2*t31*t34+t38*t40)*t43*om*t24*t27
     #+t57*t17*t59*t27-t62*t66*t71-t75*t19*t69*t27-t62*t66*t83-t75*t59*t
     #81-t89*t90*t99-t51*(t102*t69*t27-t107)-t110*t111*t118-t51*(t102*t8
     #2-t107)
        t124 = t13*t65
        t125 = t92*t19
        t129 = t31*omm
        t130 = t32*t10
        t133 = t65*t19
        t134 = t133*t93
        t135 = t43*t19
        t137 = t135*t95*t27
        t140 = t6*omm
        t144 = amplambda*t3
        t150 = t57*t65
        t156 = t43*om*t28
        t163 = t65*t95
        t165 = t163*t19*t81
        t172 = t163*om*t112+t65*t92*t116
        forces(2) = -t124*t125*t20*t27+t129*t130*t99-t51*(-t134+t137)-(2
     #*t88*t140*t10+t144*t140*t39)*t65*t71+t150*t125*t27-t144*t6*omm*t39
     #*(t156+t65*t24*t81*om)-t56*(t137+t165)-t110*t111*t172-t51*(t165-t1
     #34)
        t181 = t133*t113
        t183 = t135*t24*t115
        t188 = t65*t69
        t190 = t188*t19*t115
        forces(3) = -t124*t113*t21+t129*t130*t118-t51*(-t181+t183)-t89*t
     #90*t172-t51*(-t181+t190)-(2*t88*t40-t144*t34)*t65*t83+t150*t113*t1
     #9+t144*omm*t130*(t156+t188*t97)-t56*(t183+t190)

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine FORCINGTT( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, h)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t102
      doubleprecision t103
      doubleprecision t105
      doubleprecision t107
      doubleprecision t110
      doubleprecision t115
      doubleprecision t118
      doubleprecision t119
      doubleprecision t120
      doubleprecision t121
      doubleprecision t122
      doubleprecision t124
      doubleprecision t125
      doubleprecision t127
      doubleprecision t13
      doubleprecision t135
      doubleprecision t14
      doubleprecision t140
      doubleprecision t141
      doubleprecision t144
      doubleprecision t145
      doubleprecision t146
      doubleprecision t147
      doubleprecision t150
      doubleprecision t154
      doubleprecision t16
      doubleprecision t161
      doubleprecision t163
      doubleprecision t169
      doubleprecision t17
      doubleprecision t173
      doubleprecision t176
      doubleprecision t185
      doubleprecision t19
      doubleprecision t194
      doubleprecision t195
      doubleprecision t2
      doubleprecision t20
      doubleprecision t201
      doubleprecision t21
      doubleprecision t22
      doubleprecision t23
      doubleprecision t25
      doubleprecision t26
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t33
      doubleprecision t34
      doubleprecision t35
      doubleprecision t36
      doubleprecision t39
      doubleprecision t41
      doubleprecision t42
      doubleprecision t43
      doubleprecision t45
      doubleprecision t47
      doubleprecision t49
      doubleprecision t5
      doubleprecision t50
      doubleprecision t55
      doubleprecision t6
      doubleprecision t60
      doubleprecision t61
      doubleprecision t66
      doubleprecision t67
      doubleprecision t69
      doubleprecision t70
      doubleprecision t71
      doubleprecision t72
      doubleprecision t73
      doubleprecision t74
      doubleprecision t76
      doubleprecision t77
      doubleprecision t84
      doubleprecision t85
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t94
      doubleprecision t95
      doubleprecision t96
      doubleprecision t97
      doubleprecision t98

      do k=kfirst,klast
         z = (k-1)*h
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t2 = omm*x+phm
        t3 = sin(t2)
        t5 = omm*y+phm
        t6 = cos(t5)
        t9 = omm*z+phm
        t10 = sin(t9)
        t13 = amprho*(2+t3*t6*t10)
        t14 = c*t
        t16 = om*(x-t14)
        t17 = sin(t16)
        t19 = om**2
        t20 = t19**2
        t21 = c**2
        t22 = t21**2
        t23 = t20*t22
        t25 = om*y+ph
        t26 = sin(t25)
        t28 = om*z+ph
        t29 = sin(t28)
        t33 = ampmu*t3
        t34 = sin(t5)
        t35 = omm*t34
        t36 = t35*t10
        t39 = cos(t2)
        t41 = cos(t9)
        t42 = t35*t41
        t43 = amplambda*t39*t42
        t45 = cos(t16)
        t47 = t19*om
        t49 = t21*t26
        t50 = t49*t29
        t55 = ampmu*(3+t39*t34*t10)
        t60 = amplambda*(2+t3*t34*t41)
        t61 = 2*t55+t60
        t66 = om*x+ph
        t67 = sin(t66)
        t69 = om*(y-t14)
        t70 = cos(t69)
        t71 = t67*t70
        t72 = t47*t21
        t73 = t72*t29
        t74 = t71*t73
        t76 = cos(t66)
        t77 = t60*t76
        t84 = om*(z-t14)
        t85 = cos(t84)
        t87 = t85*t47*t21
        t88 = t67*t26*t87
        t94 = ampmu*t39
        t95 = t94*t6
        t96 = omm*t10
        t97 = t76*t47
        t98 = sin(t69)
        t100 = t98*t21*t29
        t102 = t17*t47
        t103 = cos(t25)
        t105 = t21*t103*t29
        t107 = -t97*t100-t102*t105
        t110 = t76*t20
        t115 = t17*t20*t50
        t118 = t94*omm
        t119 = t34*t41
        t120 = sin(t84)
        t121 = t26*t120
        t122 = t121*t21
        t124 = cos(t28)
        t125 = t49*t124
        t127 = -t97*t122-t102*t125
        forces(1) = t13*t17*t23*t26*t29+(-2*t33*t36+t43)*t45*t47*t50-t61
     #*t17*t20*t50+t43*t74+t77*t20*t70*t21*t29+t43*t88+t77*t20*t26*t85*t
     #21-t95*t96*t107-t55*(-t110*t70*t21*t29+t115)-t118*t119*t127-t55*(-
     #t110*t26*t85*t21+t115)
        t135 = t13*t67
        t140 = t33*omm
        t141 = t34*t10
        t144 = t67*t20
        t145 = t144*t100
        t146 = t45*t20
        t147 = t146*t105
        t150 = t6*omm
        t154 = amplambda*t3
        t161 = t61*t67
        t163 = t20*t21
        t169 = t45*t47*t50
        t173 = t67*t103
        t176 = t173*t20*t85*t21
        t185 = -t173*t47*t120*t21-t67*t98*t72*t124
        forces(2) = t135*t98*t20*t22*t29+t140*t141*t107-t55*(t145-t147)+
     #(2*t94*t150*t10+t154*t150*t41)*t67*t70*t73-t161*t98*t163*t29-t154*
     #t6*omm*t41*(-t169-t88)-t60*(-t147-t176)-t118*t119*t185-t55*(-t176+
     #t145)
        t194 = t144*t122
        t195 = t146*t125
        t201 = t71*t163*t124
        forces(3) = t135*t121*t23+t140*t141*t127-t55*(t194-t195)-t95*t96
     #*t185-t55*(t194-t201)+(2*t94*t42-t154*t36)*t67*t26*t87-t161*t26*t1
     #20*t20*t21+t154*omm*t141*(-t169-t74)-t60*(-t195-t201)

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine TWFRSURFZ( ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     ni, nj, nk, h, kz, t, omega, c, phase, bforce,
     *                      mu, lambda )
c     *        momega,   mphase, ampmu, amplambda, fo )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 bforce(3,ifirst:ilast,jfirst:jlast), h
      integer ni, nj, nk, i, j, kz
      doubleprecision mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision lambda(ifirst:ilast,jfirst:jlast,kfirst:klast)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
c      doubleprecision momega
c      doubleprecision mphase
c      doubleprecision ampmu
c      doubleprecision amplambda

      doubleprecision forces(3)
c      doubleprecision t10
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t19
c      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t28
      doubleprecision t29
c      doubleprecision t3
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t43
      doubleprecision t44
      doubleprecision t49
c      doubleprecision t54
c      doubleprecision t56
c      doubleprecision t6
      doubleprecision t60
      doubleprecision t62
      doubleprecision t65
c      doubleprecision t9

      z = (kz-1)*h
      do j=1,nj
         y = (j-1)*h
         do i=1,ni
            x=(i-1)*h
c        t2 = momega*x+mphase
c        t3 = cos(t2)
c        t6 = sin(momega*y+mphase)
c        t9 = momega*z+mphase
c        t10 = sin(t9)
        t13 = mu(i,j,kz)
        t15 = omega*x+phase
        t16 = cos(t15)
        t19 = omega*y+phase
        t20 = sin(t19)
        t21 = c*t
        t23 = omega*(z-t21)
        t24 = sin(t23)
        t28 = omega*(x-t21)
        t29 = sin(t28)
        t32 = omega*z+phase
        t33 = cos(t32)
        t34 = t33*omega
        forces(1) = t13*(t16*omega*t20*t24+t29*t20*t34)
        t37 = sin(t15)
        t38 = cos(t19)
        t43 = omega*(y-t21)
        t44 = sin(t43)
        forces(2) = t13*(t37*t38*omega*t24+t37*t44*t34)
        t49 = cos(t23)
c        t54 = sin(t2)
c        t56 = cos(t9)
        t60 = cos(t28)
        t62 = sin(t32)
        t65 = cos(t43)
        forces(3) = 2*t13*t37*t20*t49*omega+lambda(i,j,kz)*(t6
     #0*omega*t20*t62+t37*t65*omega*t62+t37*t20*t49*omega)
        bforce(1,i,j) = forces(1)
        bforce(2,i,j) = forces(2)
        bforce(3,i,j) = forces(3)
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine accur1d(h, nz, bope, acof, ghcof)
      implicit none
      real*8 h, bope(6,8), acof(6,8,8), ghcof(6)
      integer nz
      integer k, q, m, kb, qb, mb
      real*8, allocatable,dimension(:)::u1,mu1, u1ze, mu1zze, u1z, mu1zz
      real*8 z, mucof(8), err1, err2, ph1, ph2, err(2)
c 1-d test arrays
      allocate( u1(0:nz+1) )
      allocate( mu1(0:nz+1) )
      allocate( u1ze(1:nz) )
      allocate( mu1zze(1:nz) )
      allocate( u1z(1:nz) )
      allocate( mu1zz(1:nz) )

      ph1=-0.2
      ph2=-0.5
      do k=0,nz+1
        z = (k-1)*h
        u1(k) = cos(2*z+ph1)
        mu1(k) = sin(z+ph2)
      enddo
      do k=1,nz
        z = (k-1)*h
        u1ze(k) = -2*sin(2*z+ph1)
        mu1zze(k) = -2*cos(z+ph2)*sin(2*z+ph1) 
     +       -4*sin(z+ph2)*cos(2*z+ph1)
      enddo
c numerical approximation, boundary modified in the first 6 points: low-k
      do k=1,6
        u1z(k)=0
        mu1zz(k)=0
c first derivative
        do q=1,8
          u1z(k) = u1z(k) + bope(k,q)*u1(q)
        enddo
        u1z(k) = u1z(k)/h
c second derivative (mu*u_z)_z
c averaging the coefficient
        do q=1,8
          mucof(q)=0
          do m=1,8
            mucof(q) = mucof(q)+acof(k,q,m)*mu1(m)
          enddo
        end do
c computing the second derivative
        do q=1,8
          mu1zz(k) = mu1zz(k) + mucof(q)*u1(q)
        enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
        mu1zz(k) = (mu1zz(k) + ghcof(k)*mu1(1)*u1(0))/(h*h)
      enddo

c numerical approximation, boundary modified in the first 6 points: high-k
      do k=nz,nz-5,-1
        kb = nz+1-k
        u1z(k)=0
        mu1zz(k)=0
c first derivative: note minus before bope!
        do qb=1,8
          q = nz+1-qb
          u1z(k) = u1z(k) - bope(kb,qb)*u1(q)
        enddo
        u1z(k) = u1z(k)/h
c second derivative (mu*u_z)_z
c averaging the coefficient
        do qb=1,8
          mucof(qb)=0
          do mb=1,8
            m = nz+1-mb
            mucof(qb) = mucof(qb)+acof(kb,qb,mb)*mu1(m)
          enddo
        end do
c computing the second derivative
        do qb=1,8
          q = nz+1-qb
          mu1zz(k) = mu1zz(k) + mucof(qb)*u1(q)
        enddo
c ghost point only influences the first point (k=1) because ghcof(k)=0 for k>=2
        mu1zz(k) = (mu1zz(k) + ghcof(kb)*mu1(nz)*u1(nz+1))/(h*h)
      enddo

c check accuracy and save data
      open(21,file='u1zz.dat')
      open(22,file='u1zze.dat')
      open(23,file='u1z.dat')
      open(24,file='u1ze.dat')
      write(21,*) (mu1zz(k), k=nz-5,nz)
      write(22,*) (mu1zze(k), k=nz-5,nz)
      write(23,*) (u1z(k), k=nz-5,nz)
      write(24,*) (u1ze(k), k=nz-5,nz)
      close(21)
      close(22)
      close(23)
      close(24)
c low-k
      err1=0
      err2=0
      do k=1,6
        err(1) = abs(mu1zz(k) - mu1zze(k))
        err(2) = abs(u1z(k)-u1ze(k))
        if (err(1) .gt. err1) err1 = err(1)
        if (err(2) .gt. err2) err2 = err(2)
      enddo
      write(*,*) 
      write(*,*) '1-D test of accuracy of boundary operators'
      write(*,*) 'Step size (h)=', h
      write(*,*) 'low-k: Max error in d/dz(mu*du/dz): ', err1
      write(*,*) 'low-k: Max error in du/dz: ', err2
c high-k
      err1=0
      err2=0
      do k=nz,nz-5,-1
        err(1) = abs(mu1zz(k) - mu1zze(k))
        err(2) = abs(u1z(k)-u1ze(k))
        if (err(1) .gt. err1) err1 = err(1)
        if (err(2) .gt. err2) err2 = err(2)
      enddo
      write(*,*) 'high-k: Max error in d/dz(mu*du/dz): ', err1
      write(*,*) 'high-k: Max error in du/dz: ', err2
      write(*,*) 
      return
      end

c----------------------------------------------------------
      subroutine setDT(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     mu, la, rho, h, cfl, dt)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 cfl, h, dt
      real*8 maxrad, locrad
      integer k, j, i
      maxrad = 0
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
c This formula works for 3D, assuming that la < 2*mu
            locrad = (4*mu(i,j,k)+la(i,j,k))/rho(i,j,k)
            if( locrad.gt.maxrad )then
              maxrad = locrad
            endif
          enddo
        enddo
      enddo
      dt = cfl*h/SQRT(maxrad)
      return
      end

c----------------------------------------------------------
      subroutine rhsErr(ifirst, ilast, jfirst, jlast, kfirst, klast, nz,
     +     h, fo, u2, saveData)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nz, saveData
      real*8 h
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li(3), l2(3), err(3)
      do c=1,3
        li(c) = 0
        l2(c) = 0
      enddo
c this test only includes low-k boundary points
      do k=1,6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( li(c).lt.err(c) )then
                li(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
      write(*,*) 'Max errors low-k boundary RHS:  ', 
     +     li(1), li(2), li(3)
c error in rhs
      do c=1,3
        li(c) = 0
        l2(c) = 0
      enddo
c this test only includes interior points
      do k=7,nz-6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( li(c).lt.err(c) )then
                li(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
c                'Max errors high-k boundary RHS: ', 
      write(*,*) 'Max errors interior RHS:        ', li(1), li(2), li(3)
c error in rhs
      do c=1,3
        li(c) = 0
        l2(c) = 0
      enddo
c this test only includes high-k boundary points
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              err(c) = ABS( fo(c,i,j,k) - u2(c,i,j,k) )
              if( li(c).lt.err(c) )then
                li(c) = err(c)
              endif
            enddo
          enddo
        enddo
      enddo
      write(*,*) 'Max errors high-k boundary RHS: ', 
     +     li(1), li(2), li(3)

      if (saveData.eq.1) then
*** Save errors (only interior points)
        open(21,file='x.dat')
        open(22,file='y.dat')
        open(23,file='lue.dat')
        open(24,file='lve.dat')
        open(25,file='lwe.dat')
        open(26,file='lu.dat')
        open(27,file='lv.dat')
        open(28,file='lw.dat')
c only output a slice k=const
        k=nz/2
        do j=jfirst+2,jlast-2
          write(21,*) ((i-1)*h,i=ifirst+2,ilast-2)
          write(22,*) ((j-1)*h,i=ifirst+2,ilast-2)
          write(23,*) (fo(1,i,j,k),i=ifirst+2,ilast-2)
          write(24,*) (fo(2,i,j,k),i=ifirst+2,ilast-2)
          write(25,*) (fo(3,i,j,k),i=ifirst+2,ilast-2)
          write(26,*) (u2(1,i,j,k),i=ifirst+2,ilast-2)
          write(27,*) (u2(2,i,j,k),i=ifirst+2,ilast-2)
          write(28,*) (u2(3,i,j,k),i=ifirst+2,ilast-2)
        enddo
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(27)
        close(28)
      endif

      write(*,*)

      return
      end

c------------------------------------------------------------
      subroutine solErr(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     h, uex, u, kp)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, kp
      real*8 h
      real*8 uex(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer c, k, j, i
      real*8 li, l2, err(3)

      li = 0
      l2 = 0
c this test includes all points in i,j, but only one ghost point in k
      do k=kfirst+1,klast-1
        do j=jfirst,jlast
          do i=ifirst,ilast
c exact solution in array 'uex'
            do c=1,3
              err(c) = ABS( u(c,i,j,k) - uex(c,i,j,k) )
            enddo
            if( li.lt.max(err(1),err(2),err(3)) )then
              li = max(err(1),err(2),err(3))
            endif
            l2 = l2 + 
     +           h*h*h* (err(1)**2 + err(2)**2 + err(3)**2)
          enddo
        enddo
      enddo
      write(*,101) 'Solution errors in max- and L2-norm: ', li, SQRT(l2)
 101  format(' ', a, 2(g15.7,tr2))

c save a slice of the error
      if (kp.ge. kfirst .and. kp.le.klast) then 
        open(23,file='uerr.dat')
        open(24,file='verr.dat')
        open(25,file='werr.dat')
        do j=jfirst, jlast
          write(23,*) (u(1,i,j,kp)-uex(1,i,j,kp),i=ifirst,ilast)
          write(24,*) (u(2,i,j,kp)-uex(2,i,j,kp),i=ifirst,ilast)
          write(25,*) (u(3,i,j,kp)-uex(3,i,j,kp),i=ifirst,ilast)
        enddo
        close(23)
        close(24)
        close(25)
        write(*,*)'Saved solution error to files uerr.dat, verr.dat,',
     +       ' werr.dat'
c      else
c        write(*,*)'Warning: saveSol: out of bounds: kplane= ', kp 
      endif
      return
      end

c--------------------------------------------------------------------
      subroutine saveSol(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     kp, u)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, kp
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
     
      integer j, i
c only output a slice
      if (kp.ge. kfirst .and. kp.le.klast) then 
        open(23,file='u.dat')
        open(24,file='v.dat')
        open(25,file='w.dat')
        do j=jfirst, jlast
          write(23,*) (u(1,i,j,kp),i=ifirst,ilast)
          write(24,*) (u(2,i,j,kp),i=ifirst,ilast)
          write(25,*) (u(3,i,j,kp),i=ifirst,ilast)
        enddo
        close(23)
        close(24)
        close(25)
        write(*,*)'Saved solution to files u.dat, v.dat, w.dat'
c      else
c        write(*,*)'Warning: saveSol: out of bounds: kplane= ', kp 
      endif
      return
      end


c-------------------------------------------------------------------
      subroutine forcingErr(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     nz, u, u2, um2, lu, fo, rho, mu, la,
     +     onesided, acof, bope, ghcof, t, dt, h, 
     +     om, cv, ph, omm, phm, amprho, ampmu, amplambda, kp, saveErr)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, kp, nz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 t, dt, om, cv, ph, omm, phm, amprho, ampmu, amplambda, h
      real*8 acof(6,8,8), bope(6,8), ghcof(6)
      real*8 li(3), err(3)
      integer c, i, j, k, saveErr, onesided(6)
c evaluate right hand side: L(u), result in array 'lu'
      call RHS4TH3( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     nz, onesided, acof, bope, ghcof,
     +     lu, u, mu, la, rho, h )
c evaluate twilight forcing
      call FORCING( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     fo, t, om, cv, ph, omm, phm, amprho, ampmu, amplambda, h)
c evaluate exact acceleration
      call exactAcc( ifirst, ilast, jfirst, jlast, kfirst, klast, 
     +     u2, t, om, cv, ph, h )
c compute the norm of rho*Utt - Lu - F
      call rhoUttLumF(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     nz, u, u2, um2, lu, fo, rho)

c save a slice of the truncation error
      if (saveErr.eq.1) then
        kp = nz/2
        if (kp.ge. kfirst+2 .and. kp.le.klast-2) then 
          open(23,file='tuerr.dat')
          open(24,file='tverr.dat')
          open(25,file='twerr.dat')
          do j=jfirst+2, jlast-2
            write(23,*) (um2(1,i,j,kp)-fo(1,i,j,kp),i=ifirst+2,ilast-2)
            write(24,*) (um2(2,i,j,kp)-fo(2,i,j,kp),i=ifirst+2,ilast-2)
            write(25,*) (um2(3,i,j,kp)-fo(3,i,j,kp),i=ifirst+2,ilast-2)
          enddo
          close(23)
          close(24)
          close(25)
          write(*,*)'Saved truncation errors to files tuerr.dat',
     +         ' tverr.dat, twerr.dat'
        else
          write(*,*)'Warning: saveSol: out of bounds: kplane= ', kp 
        endif
      endif
      return
      end


c-------------------------------------------------
      subroutine rhoUttLumF(ifirst, ilast, jfirst, jlast, kfirst, klast,
     +     nz, u, u2, um2, lu, fo, rho)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, kp, nz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 lu(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 u2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 um2(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 li(3), err(3)
      integer c, i, j, k

c  low-k boundary points
      do c=1,3
        li(c)=0
      enddo
      do k=1,6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              um2(c,i,j,k) = rho(i,j,k)*u2(c,i,j,k) - lu(c,i,j,k)
              err(c) = um2(c,i,j,k) - fo(c,i,j,k)
              if (abs(err(c)).gt.li(c)) li(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
      write(*,'(a, 3g14.7)')'Max spatial truncation error, low-k:    ', 
     +     li(1), li(2), li(3)

c evaluate error in the interior
      do c=1,3
        li(c)=0
      enddo
c interior points
      do k=7,nz-6
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              um2(c,i,j,k) = rho(i,j,k)*u2(c,i,j,k) - lu(c,i,j,k)
              err(c) = um2(c,i,j,k) - fo(c,i,j,k)
              if (abs(err(c)).gt.li(c)) li(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
      write(*,'(a, 3g14.7)')'Max spatial truncation error, interior: ', 
     +     li(1), li(2), li(3)
c this test only includes high-k boundary points
      do c=1,3
        li(c)=0
      enddo
      do k=nz-5,nz
        do j=jfirst+2,jlast-2
          do i=ifirst+2,ilast-2
            do c=1,3
              um2(c,i,j,k) = rho(i,j,k)*u2(c,i,j,k) - lu(c,i,j,k)
              err(c) = um2(c,i,j,k) - fo(c,i,j,k)
              if (abs(err(c)).gt.li(c)) li(c)=abs(err(c))
            enddo
          enddo
        enddo
      enddo
      write(*,'(a, 3g14.7)')'Max spatial truncation error, high-k:   ', 
     +     li(1), li(2), li(3)
      end
