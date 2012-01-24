      subroutine exactRhs(x,y,t,om,ph,th,omm,phm,la0,crea_par)
      doubleprecision x
      doubleprecision y
      doubleprecision t
      doubleprecision om
      doubleprecision ph
      doubleprecision th
      doubleprecision omm
      doubleprecision phm
      doubleprecision la0
      doubleprecision crea_par(2)

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
        crea_par(1) = acc(1)
        crea_par(2) = acc(2)
        return
        return
      end

