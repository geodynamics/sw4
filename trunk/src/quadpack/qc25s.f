      subroutine qc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result,
     *   abserr,resasc,integr,nev)
c***begin prologue  qc25s
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (bl,br), with error
c            estimate, where the weight function w has a singular
c            behaviour of algebraico-logarithmic type at the points
c            a and/or b. (bl,br) is a part of (a,b).
c***description
c
c        integration rules for integrands having algebraico-logarithmic
c        end point singularities
c        standard fortran subroutine
c        real version
c
c        parameters
c           f      - real
c                    function subprogram defining the integrand
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - real
c                    left end point of the original interval
c
c           b      - real
c                    right end point of the original interval, b.gt.a
c
c           bl     - real
c                    lower limit of integration, bl.ge.a
c
c           br     - real
c                    upper limit of integration, br.le.b
c
c           alfa   - real
c                    parameter in the weight function
c
c           beta   - real
c                    parameter in the weight function
c
c           ri,rj,rg,rh - real
c                    modified chebyshev moments for the application
c                    of the generalized clenshaw-curtis
c                    method (computed in subroutine dqmomo)
c
c           result - real
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if b1 = a or br = b.
c                    in all other cases the 15-point kronrod
c                    rule is applied, obtained by optimal addition of
c                    abscissae to the 7-point gauss rule.
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           resasc - real
c                    approximation to the integral of abs(f*w-i/(b-a))
c
c           integr - integer
c                    which determines the weight function
c                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
c                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
c                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
c                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
c                                 log(b-x)
c
c           nev    - integer
c                    number of integrand evaluations
c
c***references  (none)
c***routines called  qcheb,qk15w
c***end prologue  qc25s
c
      real a,abserr,alfa,b,beta,bl,br,centr,cheb12,cheb24,
     *  dc,f,factor,fix,fval,hlgth,resabs,resasc,
     *  result,res12,res24,rg,rh,ri,rj,u,qwgts,x
      integer i,integr,isym,nev
c
      dimension cheb12(13),cheb24(25),fval(25),rg(25),rh(25),ri(25),
     *  rj(25),x(11)
c
      external f,qwgts
c
c           the vector x contains the values cos(k*pi/24)
c           k = 1, ..., 11, to be used for the computation of the
c           chebyshev series expansion of f.
c
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),
     *  x(11)/
     *     0.9914448613738104e+00,     0.9659258262890683e+00,
     *     0.9238795325112868e+00,     0.8660254037844386e+00,
     *     0.7933533402912352e+00,     0.7071067811865475e+00,
     *     0.6087614290087206e+00,     0.5000000000000000e+00,
     *     0.3826834323650898e+00,     0.2588190451025208e+00,
     *     0.1305261922200516e+00/
c
c           list of major variables
c           -----------------------
c
c           fval   - value of the function f at the points
c                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
c                    k = 0, ..., 24
c           cheb12 - coefficients of the chebyshev series expansion
c                    of degree 12, for the function f, in the
c                    interval (bl,br)
c           cheb24 - coefficients of the chebyshev series expansion
c                    of degree 24, for the function f, in the
c                    interval (bl,br)
c           res12  - approximation to the integral obtained from cheb12
c           res24  - approximation to the integral obtained from cheb24
c           qwgts - external function subprogram defining
c                    the four possible weight functions
c           hlgth  - half-length of the interval (bl,br)
c           centr  - mid point of the interval (bl,br)
c
c***first executable statement  qc25s
      nev = 25
      if(bl.eq.a.and.(alfa.ne.0.0e+00.or.integr.eq.2.or.integr.eq.4))
     * go to 10
      if(br.eq.b.and.(beta.ne.0.0e+00.or.integr.eq.3.or.integr.eq.4))
     * go to 140
c
c           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod
c           scheme.
c
c
      call qk15w(f,qwgts,a,b,alfa,beta,integr,bl,br,
     *    result,abserr,resabs,resasc)
      nev = 15
      go to 270
c
c           this part of the program is executed only if a = bl.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
c                  *f(0.5*(br-a)*x+0.5*(br+a))
c
   10 hlgth = 0.5e+00*(br-bl)
      centr = 0.5e+00*(br+bl)
      fix = b-centr
      fval(1) = 0.5e+00*f(hlgth+centr)*(fix-hlgth)**beta
      fval(13) = f(centr)*(fix**beta)
      fval(25) = 0.5e+00*f(centr-hlgth)*(fix+hlgth)**beta
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix-u)**beta
        fval(isym) = f(centr-u)*(fix+u)**beta
   20 continue
      factor = hlgth**(alfa+0.1e+01)
      result = 0.0e+00
      abserr = 0.0e+00
      res12 = 0.0e+00
      res24 = 0.0e+00
      if(integr.gt.2) go to 70
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 1  (or 2)
c
      do 30 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   30 continue
      do 40 i=14,25
        res24 = res24+cheb24(i)*ri(i)
   40 continue
      if(integr.eq.1) go to 130
c
c           integr = 2
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 50 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res12+cheb24(i)*rg(i)
   50 continue
      do 60 i=14,25
        res24 = res24+cheb24(i)*rg(i)
   60 continue
      go to 130
c
c           compute the chebyshev series expansion of the
c           following function
c           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
c
   70 fval(1) = fval(1)*alog(fix-hlgth)
      fval(13) = fval(13)*alog(fix)
      fval(25) = fval(25)*alog(fix+hlgth)
      do 80 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*alog(fix-u)
        fval(isym) = fval(isym)*alog(fix+u)
   80 continue
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 3  (or 4)
c
      do 90 i=1,13
        res12 = res12+cheb12(i)*ri(i)
        res24 = res24+cheb24(i)*ri(i)
   90 continue
      do 100 i=14,25
        res24 = res24+cheb24(i)*ri(i)
  100 continue
      if(integr.eq.3) go to 130
c
c           integr = 4
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 110 i=1,13
        res12 = res12+cheb12(i)*rg(i)
        res24 = res24+cheb24(i)*rg(i)
  110 continue
      do 120 i=14,25
        res24 = res24+cheb24(i)*rg(i)
  120 continue
  130 result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor
      go to 270
c
c           this part of the program is executed only if b = br.
c           ----------------------------------------------------
c
c           compute the chebyshev series expansion of the
c           following function
c           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
c                *f(0.5*(b-bl)*x+0.5*(b+bl))
c
  140 hlgth = 0.5e+00*(br-bl)
      centr = 0.5e+00*(br+bl)
      fix = centr-a
      fval(1) = 0.5e+00*f(hlgth+centr)*(fix+hlgth)**alfa
      fval(13) = f(centr)*(fix**alfa)
      fval(25) = 0.5e+00*f(centr-hlgth)*(fix-hlgth)**alfa
      do 150 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = f(u+centr)*(fix+u)**alfa
        fval(isym) = f(centr-u)*(fix-u)**alfa
  150 continue
      factor = hlgth**(beta+0.1e+01)
      result = 0.0e+00
      abserr = 0.0e+00
      res12 = 0.0e+00
      res24 = 0.0e+00
      if(integr.eq.2.or.integr.eq.4) go to 200
c
c           integr = 1  (or 3)
c
      call qcheb(x,fval,cheb12,cheb24)
      do 160 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  160 continue
      do 170 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  170 continue
      if(integr.eq.1) go to 260
c
c           integr = 3
c
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
      do 180 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  180 continue
      do 190 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  190 continue
      go to 260
c
c           compute the chebyshev series expansion of the
c           following function
c           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
c
  200 fval(1) = fval(1)*alog(hlgth+fix)
      fval(13) = fval(13)*alog(fix)
      fval(25) = fval(25)*alog(fix-hlgth)
      do 210 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
        fval(i) = fval(i)*alog(u+fix)
        fval(isym) = fval(isym)*alog(fix-u)
  210 continue
      call qcheb(x,fval,cheb12,cheb24)
c
c           integr = 2  (or 4)
c
      do 220 i=1,13
        res12 = res12+cheb12(i)*rj(i)
        res24 = res24+cheb24(i)*rj(i)
  220 continue
      do 230 i=14,25
        res24 = res24+cheb24(i)*rj(i)
  230 continue
      if(integr.eq.2) go to 260
      dc = alog(br-bl)
      result = res24*dc
      abserr = abs((res24-res12)*dc)
      res12 = 0.0e+00
      res24 = 0.0e+00
c
c           integr = 4
c
      do 240 i=1,13
        res12 = res12+cheb12(i)*rh(i)
        res24 = res24+cheb24(i)*rh(i)
  240 continue
      do 250 i=14,25
        res24 = res24+cheb24(i)*rh(i)
  250 continue
  260 result = (result+res24)*factor
      abserr = (abserr+abs(res24-res12))*factor
  270 return
      end
