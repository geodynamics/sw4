c$$$      implicit none
c$$$      real*8 gettimec
c$$$      integer ifirst, ilast, jfirst, jlast, kfirst, klast
c$$$      integer nig, njg, nkg, gh, iseed1, iseed2, iseed3
c$$$      integer randw(3)
c$$$      real*8 dist, distz, h, t1
c$$$      real*8, allocatable, dimension(:,:,:) :: w, wgh
c$$$
c$$$      ifirst = 10
c$$$      ilast  = 200
c$$$      jfirst = 80
c$$$      jlast  = 300
c$$$      kfirst = 1
c$$$      klast  = 110
c$$$      allocate( w(ifirst:ilast,jfirst:jlast,kfirst:klast) )
c$$$      allocate( wgh(ifirst:ilast,jfirst:jlast,kfirst:klast) )
c$$$      h = 0.1d0
c$$$      dist  = 0.5d0
c$$$      distz = 0.5d0
c$$$      gh = 2
c$$$      iseed1 = -1
c$$$      iseed2 = -1
c$$$      iseed3 = -1
c$$$      IF(ISEED1.LE.0.OR.ISEED1.GT.132362.OR.
c$$$     *     ISEED2.LE.0.OR.ISEED2.GT.131726.OR.
c$$$     *     ISEED3.LE.0.OR.ISEED3.GT.131656)THEN
c$$$         ISEED1=1234
c$$$         ISEED2=5678
c$$$         ISEED3=9876
c$$$      ENDIF
c$$$      randw(1) = iseed1
c$$$      randw(2) = iseed2
c$$$      randw(3) = iseed3
c$$$      nig = 1000
c$$$      njg = 1000
c$$$      nkg = 1000
c$$$      write(*,*) 'Give global dims, nx, ny, nz '
c$$$      read(*,*) nig, njg, nkg
c$$$      t1 = gettimec()
c$$$      call RANDOMFIELD3D( ifirst, ilast, jfirst, jlast, kfirst, klast,
c$$$     *       nig, njg, nkg, gh, w, wgh, dist, distz, h, randw ) 
c$$$      t1 = gettimec()-t1
c$$$      write(*,101) 'execution time ',t1, ' seconds '
c$$$ 101  format(' ',a,g15.5,a)
c$$$      end
      subroutine RANDOMFIELD3D( ifirst, ilast, jfirst, jlast, kfirst,
     *                          klast, nig, njg, nkg, gh, w, wgh, dist, 
     *                          distz, h, randw, savedrands, p, pz )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nig, njg, nkg
      integer i, j, k, randw(3), iseed1, iseed2, iseed3, krand, iz, gh
      integer p, pz, ii, jj, kk
c      integer n, nrand(20)
      real*8  w(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  wgh(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  savedrands(ifirst-p:ilast+p,jfirst-p:jlast+p,
     *                                        1-pz:1+pz)
      real*8  h, dist, distz
c      real*8  bnorm, bznorm
      real*8, allocatable, dimension(:) :: b, bz
      real    randno

      iseed1=randw(1)
      iseed2=randw(2)
      iseed3=randw(3)
c      p =int(dist/h)+1
c      pz=int(distz/h)+1
      allocate( b(-p:p), bz(-pz:pz) )
c      bnorm = 0
      do k=-p,p
         b(k) = exp(-k*k*h*h/(2*dist*dist))
c         bnorm = bnorm + b(k)*b(k)
      enddo
c      bznorm = 0
      do k=-pz,pz
         bz(k) = exp(-k*k*h*h/(2*distz*distz))
c         bznorm = bznorm + bz(k)*bz(k)
      enddo
c Normalize filter coefficients
c      bnorm  = SQRT(bnorm)
c      bznorm = SQRT(bznorm)
c      do k=-p,p
c         b(k) = b(k)/bnorm
c      enddo
c      do k=-pz,pz
c         bz(k) = bz(k)/bznorm
c      enddo

      w = 0
      wgh=0
c      nrand = 0
*** wgh is maintained because of boundary effects. The stencil is cut
*** at boundaries and will use fewer points there.
c      do k=1-gh,nkg+gh
      do k=1-pz,nkg+gh
         do j=1-gh,njg+gh
            do i=1-gh,nig+gh
*** Random number generator, expanded into loop
               krand=iseed1/206
               iseed1 = 157*(iseed1-krand*206)-krand*21
               if(iseed1.lt.0) iseed1=iseed1+32363
               krand=iseed2/217
               iseed2 = 146*(iseed2-krand*217)-krand*45
               if(iseed2.lt.0) iseed2=iseed2+31727
               krand=iseed3/222
               iseed3 = 142*(iseed3-krand*222)-krand*133
               if(iseed3.lt.0) iseed3=iseed3+31657
               iz = iseed1-iseed2
               if(iz.gt.706) iz = iz-32362
               iz = iz+ iseed3
               if(iz.lt.1) iz = iz+32362
*** First test if loop over stencil is necessary at all
               if( i+p.ge.ifirst .and. i-p.le.ilast .and. 
     *             j+p.ge.jfirst .and. j-p.le.jlast .and. 
     *             k+pz.ge.kfirst .and. k-pz.le.klast )then
*** Compute the random number in [-1,1], and loop over stencil
                  randno = 2*real(iz)*3.0899e-5 - 1
c                  randno = sqrt(3.0)*randno
c                  n = (randno+1)*10
c                  nrand(n+1) = nrand(n+1) + 1
                  if( k.le.1+pz )then
                     savedrands(i,j,k) =randno
c                     write(*,*) 'saved ',i,j,k,randno
                  endif
                  do kk=-pz,pz
                     if( k-kk.ge.kfirst .and. k-kk.le.klast )then 
                     do jj=-p,p
                        if( j-jj.ge.jfirst .and. j-jj.le.jlast )then
                        do ii=-p,p
                           if( i-ii.ge.ifirst .and. i-ii.le.ilast )then
                              w(i-ii,j-jj,k-kk) = w(i-ii,j-jj,k-kk) + 
     *                             b(ii)*b(jj)*bz(kk)*randno
                              wgh(i-ii,j-jj,k-kk) = wgh(i-ii,j-jj,k-kk)
     *                                   + b(ii)*b(jj)*bz(kk)
c                              wgh(i-ii,j-jj,k-kk) = wgh(i-ii,j-jj,k-kk)
c     *                           + b(ii)*b(jj)*bz(kk)*b(ii)*b(jj)*bz(kk)
                           endif
                        enddo
                     endif
                     enddo
                  endif
                  enddo
               endif
            enddo
         enddo
      enddo

c      write(*,*) 'number of numbers '
c      do n=1,20
c         write(*,*) n, ' ', nrand(n)
c      enddo
c      write(*,*) '-----------------------------'

*** Save state of random number generator
      randw(1)=iseed1    
      randw(2)=iseed2
      randw(3)=iseed3

*** Divide weights
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               w(i,j,k) = w(i,j,k)/wgh(i,j,k)
c               w(i,j,k) = w(i,j,k)/SQRT(wgh(i,j,k))
            enddo
         enddo
      enddo
      deallocate(b,bz)
      end

c-----------------------------------------------------------------------
      subroutine RANDOMFIELD3DC( ifirst, ilast, jfirst, jlast, kfirst,
     *                          klast, nig, njg, nkg, gh, w, wgh, dist, 
     *                          distz, h, z, randw, savedrands, p, pz )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, nig, njg, nkg
      integer i, j, k, randw(3), iseed1, iseed2, iseed3, krand, iz, gh
      integer p, pz, ii, jj, kk
      real*8  w(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  wgh(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  savedrands(ifirst-p:ilast+p,jfirst-p:jlast+p,
     *                                       nkg-pz:nkg+pz)
      real*8  h, dist, distz, i2dz, bz
      real*8, allocatable, dimension(:) :: b
      real    randno

      iseed1=randw(1)
      iseed2=randw(2)
      iseed3=randw(3)
c      p =int(dist/h)+1
c      pz=int(distz/h)+1
      allocate( b(-p:p) )
      do k=-p,p
         b(k) = exp(-k*k*h*h/(2*dist*dist))
      enddo
      i2dz = 1/(2*distz*distz)      
      w = 0
      wgh=0
      do k=1-gh,nkg+pz
         do j=1-gh,njg+gh
            do i=1-gh,nig+gh
*** Random number generator, expanded into loop
               krand=iseed1/206
               iseed1 = 157*(iseed1-krand*206)-krand*21
               if(iseed1.lt.0) iseed1=iseed1+32363
               krand=iseed2/217
               iseed2 = 146*(iseed2-krand*217)-krand*45
               if(iseed2.lt.0) iseed2=iseed2+31727
               krand=iseed3/222
               iseed3 = 142*(iseed3-krand*222)-krand*133
               if(iseed3.lt.0) iseed3=iseed3+31657
               iz = iseed1-iseed2
               if(iz.gt.706) iz = iz-32362
               iz = iz+ iseed3
               if(iz.lt.1) iz = iz+32362
*** First test if loop over stencil is necessary at all
               if( i+p.ge.ifirst  .and. i-p.le.ilast .and. 
     *             j+p.ge.jfirst  .and. j-p.le.jlast .and. 
     *             k+pz.ge.kfirst .and. k-pz.le.klast )then
*** Compute the random number in [-1,1], and loop over stencil
                  if( k.lt.nkg-pz )then
                     randno = 2*real(iz)*3.0899e-5-1
                  else
                     randno = savedrands(i,j,k)
                  endif
                  do kk=-pz,pz
                     if( k-kk.ge.kfirst .and. k-kk.le.klast )then 

                        do jj=-p,p
                           if( j-jj.ge.jfirst .and. j-jj.le.jlast )then
                           do ii=-p,p
                              if( i-ii.ge.ifirst .and. i-ii.le.ilast 
     *                                                          )then
                                 bz = EXP(-(z(i-ii,j-jj,k-kk)-
     *                                      z(i-ii,j-jj,k)    )**2*i2dz)
                                 w(i-ii,j-jj,k-kk) = w(i-ii,j-jj,k-kk) + 
     *                             b(ii)*b(jj)*bz*randno
                                 wgh(i-ii,j-jj,k-kk) = 
     *                         wgh(i-ii,j-jj,k-kk) + b(ii)*b(jj)*bz
                              endif
                           enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

*** Save and return state of random number generator
      randw(1)=iseed1    
      randw(2)=iseed2
      randw(3)=iseed3

*** Divide weights
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               w(i,j,k) = w(i,j,k)/wgh(i,j,k)
            enddo
         enddo
      enddo
      deallocate(b)
      end

c-----------------------------------------------------------------------
      subroutine PERTURBVELOCITY( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, vs, vp, per, amp, grad, zmin, h )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k
c      integer nrand(35), n
      real*8  per(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  vs(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  vp(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  amp, grad, a, z, zmin, h, sqrt3, maxper, minper

c      sqrt3 = SQRT(3d0)
c      nrand = 0
c      minper = 1d38
c      maxper = -1d38
c      write(*,*) 'amp = ',amp, ' grad= ',grad
      do k=kfirst,klast
         z = zmin + (k-1)*h
         A = amp + grad*z
         do j=jfirst,jlast
            do i=ifirst,ilast
               vs(i,j,k) = (1+A*per(i,j,k))*vs(i,j,k)
               vp(i,j,k) = (1+A*per(i,j,k))*vp(i,j,k)
c               n = (per(i,j,k)+sqrt3)*10
c               nrand(n+1) = nrand(n+1)+1
c               if( per(i,j,k).gt.maxper )then
c                  maxper = per(i,j,k)
c               endif
c               if( per(i,j,k).lt.minper )then
c                  minper = per(i,j,k)
c               endif
            enddo
         enddo
      enddo
c      write(*,*) 'number of numbers '
c      do n=1,35
c         write(*,*) n, ' ', nrand(n)
c      enddo
c      write(*,*) 'min and max perturbation ',minper,' ',maxper
c      write(*,*) '-----------------------------'

      end

c-----------------------------------------------------------------------
      subroutine PERTURBVELOCITYC( ifirst, ilast, jfirst, jlast, kfirst,
     *     klast, vs, vp, per, amp, grad, z )
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast
      integer i, j, k
      real*8  per(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  vs(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  vp(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  z(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8  amp, grad, a

      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               A = amp + grad*z(i,j,k)
               vs(i,j,k) = (1+A*per(i,j,k))*vs(i,j,k)
               vp(i,j,k) = (1+A*per(i,j,k))*vp(i,j,k)
            enddo
         enddo
      enddo
      end
