      program lambmain
      implicit none
      real*8 t, r, uex3, mu, cs, fz, rho, tmax, dt
      integer tfun, k, nt

c VerySmoothBump with tfun=1
c C6SmoothBump with tfun=2
      tfun = 1
      cs = 1
      rho = 1.5
      mu = cs*cs*rho

c max time and number of time step
      tmax=30
      nt = 1290

c source to receiver distance
      r=1.0

      write(*,101)'Shear speed: ',cs,' Density: ',rho,
     +     ' Time function: ', tfun,' (1=VerySmoothBump, 2=C6SmoothBmp)'
      write(*,101)'Source-rec dist: ', r,' Max time: ',tmax,
     +     ' # time steps: ', nt
 101  format(' ', a, es12.5, a, es12.5, a, i7, a)

c name of file
      open(10,file='uz1e.dat',status='unknown')

c should not need to change anything beyond this point
      dt = tmax/(nt)

      do k=0,nt
        t = dt*k
        call LAMBONEPOINT(t,r,uex3, mu, cs, fz, tfun)
c testing
        write(*,*) t, uex3
        write(10,*) t, uex3
      enddo
      close(10)
      end

      subroutine LAMBONEPOINT(t, r, uex3, mu, cs, fz, tfun )
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

      common /funpars/ tim, rocs
      save /funpars/

      if( tfun.ne.1 .and. tfun.ne.2 )then
         write(*,*) 'ERROR LAMBEXACTNUMQUAD_ only implemented '//
     * ' for VerySmoothBump and C6SmoothBump '
         return
      endif
      tim = t
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
               uzex = 0
            elseif( r.ge.1d-10 )then
               csor = cs/r
               rocs = r/cs
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
      save /funpars/

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
      save /funpars/

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
      save /funpars/

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
      save /funpars/

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
      save /funpars/
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
      save /funpars/
      G2FUNNWS = 7*cof*(t-rocs*x)**6*(1-t+rocs*x)**6*(1-2*(t-rocs*x))*
     *  c22/SQRT((gamma-x)*(gamma+x))
      end
