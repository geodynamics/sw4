!  SW4 LICENSE
! # ----------------------------------------------------------------------
! # SW4 - Seismic Waves, 4th order
! # ----------------------------------------------------------------------
! # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
! # Produced at the Lawrence Livermore National Laboratory. 
! # 
! # Written by:
! # N. Anders Petersson (petersson1@llnl.gov)
! # Bjorn Sjogreen      (sjogreen2@llnl.gov)
! # 
! # LLNL-CODE-643337 
! # 
! # All rights reserved. 
! # 
! # This file is part of SW4, Version: 1.0
! # 
! # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
! # 
! # This program is free software; you can redistribute it and/or modify
! # it under the terms of the GNU General Public License (as published by
! # the Free Software Foundation) version 2, dated June 1991. 
! # 
! # This program is distributed in the hope that it will be useful, but
! # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
! # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
! # conditions of the GNU General Public License for more details. 
! # 
! # You should have received a copy of the GNU General Public License
! # along with this program; if not, write to the Free Software
! # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
      program lambmain
      implicit none
      real*8 t, r, uex3, mu, cs, fz, rho, tmax, dt, rocs
      integer tfun, k, nt, i
      character(360) buf
      common /funpars/ t, rocs

c VerySmoothBump with tfun=1
c C6SmoothBump with tfun=2
      tfun = 1
      cs = 1
      rho = 1.5
      mu = cs*cs*rho
      fz = 1d0
c only seems to work for fz=1.0
c      fz = 0.5d0
c max time and number of time step
      tmax=20.0
c      nt = 1032
c      nt = 1677
      nt = 4695

c source to receiver distance
      r=10d0

c read  command line arguments
      i      = 1
      do while( i .lt. IARGC() )
         call GETARG( i, buf )
         if( buf .eq. '-tmax' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') tmax
            i = i+2
         elseif( buf.eq.'-nsteps' )then
            call GETARG(i+1,buf)
            read(buf,'(i10)') nt
            i = i+2
         elseif( buf.eq.'-dist' )then
            call GETARG(i+1,buf)
            read(buf,'(f16.0)') r
            i = i+2
         elseif( buf.eq.'-tfunc' )then
            call GETARG(i+1,buf)
            read(buf,'(i10)') tfun
            i = i+2
         else
            i = i+ 1
         endif
      enddo

      write(*,100)'Fz: ', fz
 100  format(' ', a, es12.5)
      write(*,"(' ', 'Cp/Cs=sqrt(3) (hardcoded)')")
      write(*,101)'Shear speed: ',cs,' Density: ',rho,
     +     ' Time function: ', tfun,' (1=VerySmoothBump, 2=C6SmoothBmp)'
      write(*,101)'Source-rec dist: ', r,' Max time: ',tmax,
     +     ' # time steps: ', nt
 101  format(' ', a, es12.5, a, es12.5, a, i7, a)


c basic checks
      if (.not.(tfun .eq. 1 .or. tfun .eq. 2)) then
        write(*,*)'Unknown time function, tfun=', tfun
        stop
      endif

      if (.not.(nt .gt. 0)) then
        write(*,*)'# time steps must be positive, not:', nt
        stop
      endif

      if (.not.(tmax .gt. 0)) then
        write(*,*)'End time must be positive, not:', tmax
        stop
      endif

      if (.not.(r .gt. 0)) then
        write(*,*)'src-rec distance must be positive, not:', r
        stop
      endif
      

c name of file
      open(10,file='uzex.dat',status='unknown')

c should not need to change anything beyond this point
      dt = tmax/(nt)
      rocs = r/cs

      do k=0,nt
        t = dt*k
        call LAMBONEPOINT(r,uex3, mu, cs, fz, tfun)
        write(10,102) t, uex3
 102    format(' ',es23.16,tr2,es23.16)
c testing
c        write(*,*) t, uex3
      enddo
      close(10)
      end

      subroutine LAMBONEPOINT(r, uex3, mu, cs, fz, tfun )
      implicit none

      external G1FUN, G2FUN, G2FUNNW
      external G1FUNS, G2FUNS, G2FUNNWS

      real*8 cpocs, gamma, gammac, pi, c21, c22, c12, c11, uex3
      parameter( cpocs=1.73205080756888d0 )
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
      parameter( pi= 3.14159265358979d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 G1FUN, G2FUN, G2FUNNW, INTFCN3
      real*8 G1FUNS, G2FUNS, G2FUNNWS, INTFCN3S

      integer tfun
      real*8 t, x0, y0, cs, uzex, r, cp
      real*8 a1, b1, a2, b2, a3, b3, al1, al2, al3
      real*8 fz, mu, csor, rocs, tim

      integer neval, ier, limit, lenw, last, integr
      integer, allocatable, dimension(:) :: iwork
      real*8 a, b, epsabs, epsrel, result
      real*8 abserr, exact, alfa, beta
      real*8, allocatable, dimension(:) :: work

      common /funpars/ t, rocs
c      save /funpars/

      if( tfun.ne.1 .and. tfun.ne.2 )then
         write(*,*) 'ERROR LAMBEXACTNUMQUAD_ only implemented '//
     * ' for VerySmoothBump and C6SmoothBump '
         return
      endif
c now passing t through the common block
c      tim = t
      limit = 100
      lenw = 4*limit
      epsabs = 1d-12
      epsrel = 1d-12
      allocate( iwork(limit+10), work(lenw+10))      

c$$$      do j=jfirst,jlast
c$$$         do i=ifirst,ilast
c$$$            x = (i-1)*h
c$$$            y = (j-1)*h
c$$$            r = SQRT( (x-x0)*(x-x0)+(y-y0)*(y-y0) )
            cp = cpocs*cs
            if( r.ge.cp*t )then
               uex3 = 0
            elseif( r.ge.1d-10 )then
               csor = cs/r
c               rocs = r/cs
               a1 = MAX(csor*(t-1),1/cpocs)
               b1 = MIN(csor*t,1d0)
               a2 = MAX(csor*(t-1),1d0)
               b2 = MIN(csor*t,gamma)
               a3 = MAX(csor*(t-1),gamma)
               b3 = csor*t
               uzex = 0
               if( b1.gt.a1+1e-12 )then
c     write(*,*) 'Doing 1'
                  if( tfun.eq.1 )then
                  call DQAGS( G1FUN, a1, b1, epsabs, epsrel, result, 
     *          abserr, neval, ier, limit, lenw, last, iwork, work )
                  else
                  call DQAGS( G1FUNS, a1, b1, epsabs, epsrel, result, 
     *          abserr, neval, ier, limit, lenw, last, iwork, work )
                  endif
                  if( ier.ne.0 )then
                     write(*,*)'ERROR: call to DQAGS returns ier = ',ier
                  endif
                  if( abserr.gt.10*epsabs )then
                     write(*,*) 'WARNING: DQAGS returns error ',abserr,
     *                    ' input tolerance ', epsabs           
                  endif
                  if( tfun.eq.1 )then
                  uzex = rocs*result + 
     *                 0.5d0*c21*INTFCN3(t-a1*rocs,t-b1*rocs) 
                  else
                  uzex = rocs*result + 
     *                 0.5d0*c21*INTFCN3S(t-a1*rocs,t-b1*rocs) 
                  endif
               endif
               if( b2.gt.a2+1e-12 )then
c            write(*,*) 'Doing 2'
                  alfa = 0
                  beta = -0.5d0
                  integr = 1
                  if( abs(b2-gamma).lt.1d-12 )then
*** Singular integral
                     if( tfun.eq.1 )then
                     call DQAWS( G2FUN, a2, b2, alfa, beta, integr, 
     *            epsabs, epsrel, result, abserr, neval, ier, limit, 
     *            lenw, last, iwork, work )
                     else
                     call DQAWS( G2FUNS, a2, b2, alfa, beta, integr, 
     *            epsabs, epsrel, result, abserr, neval, ier, limit, 
     *            lenw, last, iwork, work )
                     endif
                  else
*** Away from singularity
                     if( tfun.eq.1 )then
                     call DQAGS( G2FUNNW, a2, b2, epsabs, epsrel,result, 
     *                    abserr, neval, ier, limit, lenw, last, iwork, 
     *                    work )
                     else
                     call DQAGS( G2FUNNWS,a2, b2, epsabs, epsrel,result, 
     *                    abserr, neval, ier, limit, lenw, last, iwork, 
     *                    work )
                     endif
                  endif               
                  if( ier.ne.0 )then
                     write(*,*)'ERROR: call to DQAWS returns ier = ',ier
                  endif
                  if( abserr.gt.10*epsabs )then
                     write(*,*) 'WARNING: DQAWS returns error ',abserr,
     *                    ' input tolerance ', epsabs           
                  endif
                  uzex = uzex + rocs*result
                  if( tfun.eq.1 )then
                     uzex = uzex - c21*INTFCN3(t-b2*rocs,t-a2*rocs)
                  else
                     uzex = uzex - c21*INTFCN3S(t-b2*rocs,t-a2*rocs)
                  endif
               endif
               if( b3.gt.a3 )then
c            write(*,*) 'Doing 3'
                  if( tfun.eq.1 )then
                     uzex = uzex - c21*INTFCN3(0d0,t-a3*rocs)
                  else
                     uzex = uzex - c21*INTFCN3S(0d0,t-a3*rocs)
                  endif
               endif
               if( r.ge.1e-10 )then
                  uex3 = -uzex*fz/(pi*pi*mu*r)*cpocs*cpocs
                else
                  uex3=0
               endif
            endif
c$$$         enddo
c$$$      enddo

      deallocate( iwork, work )

      return
      end

c-----------------------------------------------------------------------
      real*8 function INTFCN3( a, b )
      implicit none
      real*8 cof
      parameter( cof = 1024d0 )
      real*8 a, b
      INTFCN3 = cof*(b**5*(1-b)**5-a**5*(1-a)**5)
      end

c-----------------------------------------------------------------------
      real*8 function INTFCN3S( a, b )
      implicit none
      real*8 cof
c      parameter( cof = 16384d0 )
      parameter( cof = 51480d0 )
      real*8 a, b
      INTFCN3S = cof*(b**7*(1-b)**7-a**7*(1-a)**7)
      end

c-----------------------------------------------------------------------
      real*8 function G1FUN( x )
      implicit none
      real*8 x
      real*8 cof
      parameter( cof = 1024d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/

      G1FUN = 5*cof*(t-rocs*x)**4*(1-t+rocs*x)**4*(1-2*(t-rocs*x))*
     * ( 0.5d0*c22/SQRT(gamma**2-x*x)-c11/SQRT(x*x-gammac**2)
     * + c12/sqrt(x*x-0.25d0) )

      end

c-----------------------------------------------------------------------
      real*8 function G1FUNS( x )
      implicit none
      real*8 x
      real*8 cof
c      parameter( cof = 16384d0 )
      parameter( cof = 51480d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/

      G1FUNS = 7*cof*(t-rocs*x)**6*(1-t+rocs*x)**6*(1-2*(t-rocs*x))*
     * ( 0.5d0*c22/SQRT(gamma**2-x*x)-c11/SQRT(x*x-gammac**2)
     * + c12/sqrt(x*x-0.25d0) )

      end

c-----------------------------------------------------------------------
      real*8 function G2FUN( x )
      implicit none
      real*8 x
      real*8 cof
      parameter( cof = 1024d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/

      G2FUN = 5*cof*(t-rocs*x)**4*(1-t+rocs*x)**4*(1-2*(t-rocs*x))*
     *  c22/SQRT(gamma+x)
      end

c-----------------------------------------------------------------------
      real*8 function G2FUNS( x )
      implicit none
      real*8 x
      real*8 cof
c      parameter( cof = 16384d0 )
      parameter( cof = 51480d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/

      G2FUNS = 7*cof*(t-rocs*x)**6*(1-t+rocs*x)**6*(1-2*(t-rocs*x))*
     *  c22/SQRT(gamma+x)
      end

c-----------------------------------------------------------------------
      real*8 function G2FUNNW( x )
      implicit none
      real*8 x
      real*8 cof
      parameter( cof = 1024d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/
      G2FUNNW = 5*cof*(t-rocs*x)**4*(1-t+rocs*x)**4*(1-2*(t-rocs*x))*
     *  c22/SQRT((gamma-x)*(gamma+x))
      end

c-----------------------------------------------------------------------
      real*8 function G2FUNNWS( x )
      implicit none
      real*8 x
      real*8 cof
c      parameter( cof = 16384d0 )
      parameter( cof = 51480d0 )
      real*8 c11, c12, c21, c22, gamma, gammac
      parameter( gamma=1.08766387358054d0 )
      parameter( gammac=0.563016250305247d0 )
*** pi/8
      parameter( c21=0.392699081698724d0 ) 
*** pi/48*sqrt(3*sqrt(3)+5)
      parameter( c22=0.208990620247103d0 )
*** pi/96*sqrt(3)
      parameter( c12=0.0566812301323193d0)
*** pi/96*sqrt(3*sqrt(3)-5)
      parameter( c11=0.0144935735221071d0)

      real*8 t, rocs
      common /funpars/ t, rocs
c      save /funpars/
      G2FUNNWS = 7*cof*(t-rocs*x)**6*(1-t+rocs*x)**6*(1-2*(t-rocs*x))*
     *  c22/SQRT((gamma-x)*(gamma+x))
      end
