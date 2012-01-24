      subroutine ffsurfz(x,y,z,t,om,c,ph,omm,phm,amprho,ampmu,amplambda,
     #crea_par)
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

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t19
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t43
      doubleprecision t44
      doubleprecision t49
      doubleprecision t54
      doubleprecision t56
      doubleprecision t6
      doubleprecision t60
      doubleprecision t62
      doubleprecision t65
      doubleprecision t9

        t2 = omm*x+phm
        t3 = cos(t2)
        t6 = sin(omm*y+phm)
        t9 = omm*z+phm
        t10 = sin(t9)
        t13 = ampmu*(3+t3*t6*t10)
        t15 = om*x+ph
        t16 = cos(t15)
        t19 = om*y+ph
        t20 = sin(t19)
        t21 = c*t
        t23 = om*(z-t21)
        t24 = sin(t23)
        t28 = om*(x-t21)
        t29 = sin(t28)
        t32 = om*z+ph
        t33 = cos(t32)
        t34 = t33*om
        forces(1) = t13*(t16*om*t20*t24+t29*t20*t34)
        t37 = sin(t15)
        t38 = cos(t19)
        t43 = om*(y-t21)
        t44 = sin(t43)
        forces(2) = t13*(t37*t38*om*t24+t37*t44*t34)
        t49 = cos(t23)
        t54 = sin(t2)
        t56 = cos(t9)
        t60 = cos(t28)
        t62 = sin(t32)
        t65 = cos(t43)
        forces(3) = 2*t13*t37*t20*t49*om+amplambda*(2+t54*t6*t56)*(t60*o
     #m*t20*t62+t37*t65*om*t62+t37*t20*t49*om)
        crea_par(1) = forces(1)
        crea_par(2) = forces(2)
        crea_par(3) = forces(3)
        return
        return
      end
