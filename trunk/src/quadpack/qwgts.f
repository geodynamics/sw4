      real function qwgts(x,a,b,alfa,beta,integr)
c***begin prologue  qwgts
c***refer to qk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, algebraico-logarithmic
c             end-point singularities
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine qaws and defines the weight function.
c***end prologue  qwgts
c
      real a,alfa,b,beta,bmx,x,xma
      integer integr
c***first executable statement
      xma = x-a
      bmx = b-x
      qwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 qwgts = qwgts*alog(xma)
      go to 40
   20 qwgts = qwgts*alog(bmx)
      go to 40
   30 qwgts = qwgts*alog(xma)*alog(bmx)
   40 return
      end
