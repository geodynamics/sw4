      subroutine fgtt(x,y,z,t,om,c,ph,omm,phm,amprho,ampmu,amplambda,cre
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
        crea_par(1) = forces(1)
        crea_par(2) = forces(2)
        crea_par(3) = forces(3)
        return
        return
      end
