      subroutine facc(x,y,z,t,om,c,ph,omm,phm,amprho,ampmu,amplambda,cre
     #a_par)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision om
      doubleprecision c
      doubleprecision ph
      doubleprecision omm
      doubleprecision phm
      doubleprecision amprho
      doubleprecision ampmu
      doubleprecision amplambda
      doubleprecision crea_par(3)

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
        crea_par(1) = acc(1)
        crea_par(2) = acc(2)
        crea_par(3) = acc(3)
        return
        return
      end
