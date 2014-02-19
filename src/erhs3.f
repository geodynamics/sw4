      subroutine fgt(x,y,z,t,omega,c,phase,momega,mphase,amprho,ampmu,am
     #plambda,crea_par)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
      doubleprecision momega
      doubleprecision mphase
      doubleprecision amprho
      doubleprecision ampmu
      doubleprecision amplambda
      doubleprecision crea_par(3)

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t105
      doubleprecision t108
      doubleprecision t109
      doubleprecision t11
      doubleprecision t110
      doubleprecision t112
      doubleprecision t114
      doubleprecision t115
      doubleprecision t117
      doubleprecision t12
      doubleprecision t125
      doubleprecision t126
      doubleprecision t129
      doubleprecision t130
      doubleprecision t131
      doubleprecision t132
      doubleprecision t135
      doubleprecision t139
      doubleprecision t146
      doubleprecision t148
      doubleprecision t15
      doubleprecision t154
      doubleprecision t158
      doubleprecision t161
      doubleprecision t17
      doubleprecision t170
      doubleprecision t177
      doubleprecision t178
      doubleprecision t18
      doubleprecision t184
      doubleprecision t19
      doubleprecision t2
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t26
      doubleprecision t27
      doubleprecision t29
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t33
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t4
      doubleprecision t42
      doubleprecision t47
      doubleprecision t48
      doubleprecision t49
      doubleprecision t51
      doubleprecision t55
      doubleprecision t56
      doubleprecision t58
      doubleprecision t59
      doubleprecision t6
      doubleprecision t60
      doubleprecision t61
      doubleprecision t62
      doubleprecision t63
      doubleprecision t65
      doubleprecision t66
      doubleprecision t7
      doubleprecision t73
      doubleprecision t74
      doubleprecision t76
      doubleprecision t77
      doubleprecision t8
      doubleprecision t83
      doubleprecision t84
      doubleprecision t85
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t90
      doubleprecision t92
      doubleprecision t93
      doubleprecision t95
      doubleprecision t97

        t2 = momega*x+mphase
        t3 = sin(t2)
        t4 = ampmu*t3
        t6 = momega*y+mphase
        t7 = sin(t6)
        t8 = momega*t7
        t10 = momega*z+mphase
        t11 = sin(t10)
        t12 = t8*t11
        t15 = cos(t2)
        t17 = cos(t10)
        t18 = t8*t17
        t19 = amplambda*t15*t18
        t21 = c*t
        t23 = omega*(x-t21)
        t24 = cos(t23)
        t26 = omega**2
        t27 = t26*omega
        t29 = c**2
        t31 = omega*y+phase
        t32 = sin(t31)
        t33 = t29*t32
        t35 = omega*z+phase
        t36 = sin(t35)
        t37 = t33*t36
        t42 = ampmu*(3+t15*t7*t11)
        t47 = amplambda*(2+t3*t7*t17)
        t48 = 2*t42+t47
        t49 = sin(t23)
        t51 = t26**2
        t55 = omega*x+phase
        t56 = sin(t55)
        t58 = omega*(y-t21)
        t59 = cos(t58)
        t60 = t56*t59
        t61 = t27*t29
        t62 = t61*t36
        t63 = t60*t62
        t65 = cos(t55)
        t66 = t47*t65
        t73 = omega*(z-t21)
        t74 = cos(t73)
        t76 = t74*t27*t29
        t77 = t56*t32*t76
        t83 = ampmu*t15
        t84 = cos(t6)
        t85 = t83*t84
        t86 = momega*t11
        t87 = t65*t27
        t88 = sin(t58)
        t90 = t88*t29*t36
        t92 = t49*t27
        t93 = cos(t31)
        t95 = t29*t93*t36
        t97 = -t87*t90-t92*t95
        t100 = t65*t51
        t105 = t49*t51*t37
        t108 = t83*momega
        t109 = t7*t17
        t110 = sin(t73)
        t112 = t32*t110*t29
        t114 = cos(t35)
        t115 = t33*t114
        t117 = -t87*t112-t92*t115
        forces(1) = -(-2*t4*t12+t19)*t24*t27*t37+t48*t49*t51*t37-t19*t63
     #-t66*t51*t59*t29*t36-t19*t77-t66*t51*t32*t74*t29+t85*t86*t97+t42*(
     #-t100*t59*t29*t36+t105)+t108*t109*t117+t42*(-t100*t32*t74*t29+t105
     #)
        t125 = t4*momega
        t126 = t7*t11
        t129 = t56*t51
        t130 = t129*t90
        t131 = t24*t51
        t132 = t131*t95
        t135 = t84*momega
        t139 = amplambda*t3
        t146 = t48*t56
        t148 = t51*t29
        t154 = t24*t27*t37
        t158 = t56*t93
        t161 = t158*t51*t74*t29
        t170 = -t158*t27*t110*t29-t56*t88*t61*t114
        forces(2) = -t125*t126*t97+t42*(t130-t132)-(2*t83*t135*t11+t139*
     #t135*t17)*t56*t59*t62+t146*t88*t148*t36+t139*t84*momega*t17*(-t154
     #-t77)+t47*(-t132-t161)+t108*t109*t170+t42*(-t161+t130)
        t177 = t129*t112
        t178 = t131*t115
        t184 = t60*t148*t114
        forces(3) = -t125*t126*t117+t42*(t177-t178)+t85*t86*t170+t42*(t1
     #77-t184)-(2*t83*t18-t139*t12)*t56*t32*t76+t146*t32*t110*t51*t29-t1
     #39*momega*t126*(-t154-t63)+t47*(-t178-t184)
        crea_par(1) = forces(1)
        crea_par(2) = forces(2)
        crea_par(3) = forces(3)
        return
        return
      end
