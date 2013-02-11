      subroutine forcingfortsg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     h, zmin, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho, zmin
      doubleprecision ampmu, amplambda, omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t102
      doubleprecision t103
      doubleprecision t104
      doubleprecision t105
      doubleprecision t106
      doubleprecision t113
      doubleprecision t114
      doubleprecision t115
      doubleprecision t116
      doubleprecision t117
      doubleprecision t119
      doubleprecision t121
      doubleprecision t122
      doubleprecision t123
      doubleprecision t124
      doubleprecision t126
      doubleprecision t129
      doubleprecision t13
      doubleprecision t132
      doubleprecision t133
      doubleprecision t14
      doubleprecision t142
      doubleprecision t143
      doubleprecision t144
      doubleprecision t145
      doubleprecision t147
      doubleprecision t148
      doubleprecision t150
      doubleprecision t152
      doubleprecision t157
      doubleprecision t158
      doubleprecision t16
      doubleprecision t167
      doubleprecision t168
      doubleprecision t17
      doubleprecision t172
      doubleprecision t173
      doubleprecision t177
      doubleprecision t180
      doubleprecision t181
      doubleprecision t184
      doubleprecision t185
      doubleprecision t19
      doubleprecision t191
      doubleprecision t195
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t214
      doubleprecision t215
      doubleprecision t220
      doubleprecision t226
      doubleprecision t229
      doubleprecision t23
      doubleprecision t231
      doubleprecision t24
      doubleprecision t250
      doubleprecision t253
      doubleprecision t26
      doubleprecision t265
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t34
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t41
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t48
      doubleprecision t5
      doubleprecision t50
      doubleprecision t51
      doubleprecision t56
      doubleprecision t6
      doubleprecision t61
      doubleprecision t62
      doubleprecision t63
      doubleprecision t72
      doubleprecision t73
      doubleprecision t75
      doubleprecision t76
      doubleprecision t78
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t84
      doubleprecision t85
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t91
      doubleprecision t92
      doubleprecision t96
      doubleprecision t97
      doubleprecision t99

      do k=kfirst,klast
         z = (k-1)*h + zmin
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

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
        t20 = c**2
        t21 = t19*t20
        t23 = om*y+ph
        t24 = sin(t23)
        t26 = om*z+ph
        t27 = sin(t26)
        t28 = t24*t27
        t31 = omstrx*x
        t32 = sin(t31)
        t34 = 1+t32/2
        t35 = ampmu*t3
        t36 = sin(t5)
        t37 = omm*t36
        t38 = t37*t10
        t41 = cos(t2)
        t43 = cos(t9)
        t44 = t37*t43
        t45 = amplambda*t41*t44
        t48 = cos(t16)
        t50 = om*t24
        t51 = t50*t27
        t56 = ampmu*(3+t41*t36*t10)
        t61 = amplambda*(2+t3*t36*t43)
        t62 = 2*t56+t61
        t63 = cos(t31)
        t72 = t19*t24
        t73 = t72*t27
        t75 = omstry*y
        t76 = sin(t75)
        t78 = 1+t76/2
        t80 = om*x+ph
        t81 = sin(t80)
        t82 = t78*t81
        t84 = om*(y-t14)
        t85 = cos(t84)
        t87 = t85*om*t27
        t88 = t82*t87
        t91 = cos(t80)
        t92 = t91*t19
        t96 = omstrz*z
        t97 = sin(t96)
        t99 = 1+t97/2
        t100 = t99*t81
        t102 = om*(z-t14)
        t103 = cos(t102)
        t104 = t24*t103
        t105 = t104*om
        t106 = t100*t105
        t113 = ampmu*t41
        t114 = t113*t6
        t115 = omm*t10
        t116 = t34*t91
        t117 = sin(t84)
        t119 = om*t117*t27
        t121 = t78*t17
        t122 = cos(t23)
        t123 = t122*om
        t124 = t123*t27
        t126 = t116*t119+t121*t124
        t129 = t19*t85
        t132 = cos(t75)
        t133 = t132*omstry
        t142 = t113*omm
        t143 = t36*t43
        t144 = sin(t102)
        t145 = t50*t144
        t147 = t99*t17
        t148 = cos(t26)
        t150 = t24*t148*om
        t152 = t116*t145+t147*t150
        t157 = cos(t96)
        t158 = t157*omstrz
        forces(1) = -t13*t17*t21*t28-t34*((-2*t35*t38+t45)*t34*t48*t51+t
     #62*t63*omstrx*t48*om*t28/2-t62*t34*t17*t73+t45*t88+t61*t78*t92*t85
     #*t27+t45*t106+t61*t99*t92*t104)-t78*(t114*t115*t126+t56*(t116*t129
     #*t27+t133*t17*t124/2-t121*t73))-t99*(t142*t143*t152+t56*(t116*t72*
     #t103+t158*t17*t150/2-t147*t73))
        t167 = t13*t81
        t168 = t117*t19
        t172 = t35*omm
        t173 = t36*t10
        t177 = t63*omstrx*t91
        t180 = t34*t81
        t181 = t168*t27
        t184 = t19*t122
        t185 = t184*t27
        t191 = t6*omm
        t195 = amplambda*t3
        t214 = t34*t48
        t215 = t214*t51
        t220 = t184*t103
        t226 = t123*t144
        t229 = t117*t148*om
        t231 = t82*t226+t100*t229
        forces(2) = -t167*t168*t20*t27-t34*(-t172*t173*t126+t56*(t177*t1
     #19/2-t180*t181+t78*t48*t185))-t78*((2*t113*t191*t10+t195*t191*t43)
     #*t78*t81*t87+t62*t132*omstry*t81*t85*om*t27/2-t62*t78*t81*t181+t19
     #5*t6*omm*t43*(t215+t106)+t61*(t214*t185+t100*t220))-t99*(t142*t143
     #*t231+t56*(t82*t220+t158*t81*t229/2-t100*t181))
        t250 = t72*t144
        t253 = t72*t148
        t265 = t129*t148
        forces(3) = -t167*t24*t144*t21-t34*(-t172*t173*t152+t56*(t177*t1
     #45/2-t180*t250+t99*t48*t253))-t78*(t114*t115*t231+t56*(t133*t81*t2
     #26/2-t82*t250+t100*t265))-t99*((2*t113*t44-t195*t38)*t99*t81*t105+
     #t62*t157*omstrz*t81*t24*t103*om/2-t62*t99*t81*t250-t195*omm*t173*(
     #t215+t88)+t61*(t214*t253+t82*t265))
        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)
       enddo
       enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine forcingttfortsg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     h, zmin, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho, zmin
      doubleprecision ampmu, amplambda, omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t102
      doubleprecision t103
      doubleprecision t106
      doubleprecision t107
      doubleprecision t110
      doubleprecision t120
      doubleprecision t121
      doubleprecision t122
      doubleprecision t123
      doubleprecision t124
      doubleprecision t125
      doubleprecision t127
      doubleprecision t129
      doubleprecision t13
      doubleprecision t131
      doubleprecision t133
      doubleprecision t135
      doubleprecision t138
      doubleprecision t14
      doubleprecision t142
      doubleprecision t143
      doubleprecision t150
      doubleprecision t151
      doubleprecision t157
      doubleprecision t158
      doubleprecision t159
      doubleprecision t16
      doubleprecision t160
      doubleprecision t161
      doubleprecision t163
      doubleprecision t165
      doubleprecision t166
      doubleprecision t168
      doubleprecision t17
      doubleprecision t171
      doubleprecision t174
      doubleprecision t175
      doubleprecision t187
      doubleprecision t188
      doubleprecision t19
      doubleprecision t192
      doubleprecision t193
      doubleprecision t197
      doubleprecision t198
      doubleprecision t2
      doubleprecision t20
      doubleprecision t203
      doubleprecision t21
      doubleprecision t212
      doubleprecision t216
      doubleprecision t22
      doubleprecision t223
      doubleprecision t23
      doubleprecision t236
      doubleprecision t238
      doubleprecision t242
      doubleprecision t246
      doubleprecision t25
      doubleprecision t252
      doubleprecision t256
      doubleprecision t259
      doubleprecision t26
      doubleprecision t279
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t30
      doubleprecision t303
      doubleprecision t315
      doubleprecision t33
      doubleprecision t34
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t39
      doubleprecision t40
      doubleprecision t43
      doubleprecision t45
      doubleprecision t46
      doubleprecision t47
      doubleprecision t5
      doubleprecision t50
      doubleprecision t52
      doubleprecision t53
      doubleprecision t54
      doubleprecision t59
      doubleprecision t6
      doubleprecision t64
      doubleprecision t65
      doubleprecision t66
      doubleprecision t74
      doubleprecision t77
      doubleprecision t78
      doubleprecision t80
      doubleprecision t82
      doubleprecision t83
      doubleprecision t84
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t90
      doubleprecision t92
      doubleprecision t96
      doubleprecision t99

      do k=kfirst,klast
         z = (k-1)*h + zmin
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h
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
        t30 = t26*t29
        t33 = omstrx*x
        t34 = sin(t33)
        t36 = 1+t34/2
        t37 = ampmu*t3
        t38 = sin(t5)
        t39 = omm*t38
        t40 = t39*t10
        t43 = cos(t2)
        t45 = cos(t9)
        t46 = t39*t45
        t47 = amplambda*t43*t46
        t50 = cos(t16)
        t52 = t19*om
        t53 = t52*t21
        t54 = t53*t30
        t59 = ampmu*(3+t43*t38*t10)
        t64 = amplambda*(2+t3*t38*t45)
        t65 = 2*t59+t64
        t66 = cos(t33)
        t74 = t20*t21
        t77 = omstry*y
        t78 = sin(t77)
        t80 = 1+t78/2
        t82 = om*x+ph
        t83 = sin(t82)
        t84 = t80*t83
        t86 = om*(y-t14)
        t87 = cos(t86)
        t88 = t84*t87
        t90 = t88*t53*t29
        t92 = cos(t82)
        t96 = t21*t29
        t99 = omstrz*z
        t100 = sin(t99)
        t102 = 1+t100/2
        t103 = t102*t83
        t106 = om*(z-t14)
        t107 = cos(t106)
        t110 = t103*t26*t107*t52*t21
        t120 = ampmu*t43
        t121 = t120*t6
        t122 = omm*t10
        t123 = t36*t92
        t124 = t123*t52
        t125 = sin(t86)
        t127 = t125*t21*t29
        t129 = t80*t17
        t131 = cos(t25)
        t133 = t21*t131*t29
        t135 = -t124*t127-t129*t52*t133
        t138 = t123*t20
        t142 = cos(t77)
        t143 = t142*omstry
        t150 = t21*t26
        t151 = t150*t29
        t157 = t120*omm
        t158 = t38*t45
        t159 = sin(t106)
        t160 = t26*t159
        t161 = t160*t21
        t163 = t102*t17
        t165 = cos(t28)
        t166 = t150*t165
        t168 = -t124*t161-t163*t52*t166
        t171 = t26*t107
        t174 = cos(t99)
        t175 = t174*omstrz
        forces(1) = t13*t17*t23*t30-t36*(-(-2*t37*t40+t47)*t36*t50*t54-t
     #65*t66*omstrx*t50*t54/2+t65*t36*t17*t74*t30-t47*t90-t64*t80*t92*t2
     #0*t87*t96-t47*t110-t64*t102*t92*t20*t26*t107*t21)-t80*(t121*t122*t
     #135+t59*(-t138*t87*t21*t29-t143*t17*t53*t131*t29/2+t129*t20*t151))
     #-t102*(t157*t158*t168+t59*(-t138*t171*t21-t175*t17*t53*t26*t165/2+
     #t163*t20*t151))
        t187 = t13*t83
        t188 = t125*t20
        t192 = t37*omm
        t193 = t38*t10
        t197 = t66*omstrx*t92
        t198 = t52*t125
        t203 = t36*t83*t20
        t212 = t6*omm
        t216 = amplambda*t3
        t223 = t87*t52*t96
        t236 = t36*t50
        t238 = t236*t52*t151
        t242 = t236*t20
        t246 = t20*t107*t21
        t252 = t84*t131
        t256 = t103*t125
        t259 = -t252*t52*t159*t21-t256*t53*t165
        forces(2) = t187*t188*t22*t29-t36*(-t192*t193*t135+t59*(-t197*t1
     #98*t96/2+t203*t127-t80*t50*t20*t133))-t80*(-(2*t120*t212*t10+t216*
     #t212*t45)*t80*t83*t223-t65*t142*omstry*t83*t223/2+t65*t80*t83*t188
     #*t96+t216*t6*omm*t45*(-t238-t110)+t64*(-t242*t133-t103*t131*t246))
     #-t102*(t157*t158*t259+t59*(-t252*t246-t175*t83*t198*t21*t165/2+t25
     #6*t74*t29))
        t279 = t159*t21
        t303 = t74*t165
        t315 = t171*t53
        forces(3) = t187*t160*t23-t36*(-t192*t193*t168+t59*(-t197*t52*t2
     #6*t279/2+t203*t161-t102*t50*t20*t166))-t80*(t121*t122*t259+t59*(-t
     #143*t83*t131*t52*t279/2+t84*t26*t20*t159*t21-t103*t87*t303))-t102*
     #(-(2*t120*t46-t216*t40)*t102*t83*t315-t65*t174*omstrz*t83*t315/2+t
     #65*t102*t83*t160*t74-t216*omm*t193*(-t238-t90)+t64*(-t242*t166-t88
     #*t303))
        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)
      enddo
      enddo
      enddo
        return
      end

c-----------------------------------------------------------------------
      subroutine forcingfortcsg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     xx, yy, zz, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda, omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t102
      doubleprecision t103
      doubleprecision t104
      doubleprecision t105
      doubleprecision t106
      doubleprecision t113
      doubleprecision t114
      doubleprecision t115
      doubleprecision t116
      doubleprecision t117
      doubleprecision t119
      doubleprecision t121
      doubleprecision t122
      doubleprecision t123
      doubleprecision t124
      doubleprecision t126
      doubleprecision t129
      doubleprecision t13
      doubleprecision t132
      doubleprecision t133
      doubleprecision t14
      doubleprecision t142
      doubleprecision t143
      doubleprecision t144
      doubleprecision t145
      doubleprecision t147
      doubleprecision t148
      doubleprecision t150
      doubleprecision t152
      doubleprecision t157
      doubleprecision t158
      doubleprecision t16
      doubleprecision t167
      doubleprecision t168
      doubleprecision t17
      doubleprecision t172
      doubleprecision t173
      doubleprecision t177
      doubleprecision t180
      doubleprecision t181
      doubleprecision t184
      doubleprecision t185
      doubleprecision t19
      doubleprecision t191
      doubleprecision t195
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t214
      doubleprecision t215
      doubleprecision t220
      doubleprecision t226
      doubleprecision t229
      doubleprecision t23
      doubleprecision t231
      doubleprecision t24
      doubleprecision t250
      doubleprecision t253
      doubleprecision t26
      doubleprecision t265
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t34
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t41
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t48
      doubleprecision t5
      doubleprecision t50
      doubleprecision t51
      doubleprecision t56
      doubleprecision t6
      doubleprecision t61
      doubleprecision t62
      doubleprecision t63
      doubleprecision t72
      doubleprecision t73
      doubleprecision t75
      doubleprecision t76
      doubleprecision t78
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t84
      doubleprecision t85
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t91
      doubleprecision t92
      doubleprecision t96
      doubleprecision t97
      doubleprecision t99

      do k=kfirst,klast
         do j=jfirst,jlast
           do i=ifirst,ilast
             x = xx(i,j,k)
             y = yy(i,j,k)
             z = zz(i,j,k)

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
        t20 = c**2
        t21 = t19*t20
        t23 = om*y+ph
        t24 = sin(t23)
        t26 = om*z+ph
        t27 = sin(t26)
        t28 = t24*t27
        t31 = omstrx*x
        t32 = sin(t31)
        t34 = 1+t32/2
        t35 = ampmu*t3
        t36 = sin(t5)
        t37 = omm*t36
        t38 = t37*t10
        t41 = cos(t2)
        t43 = cos(t9)
        t44 = t37*t43
        t45 = amplambda*t41*t44
        t48 = cos(t16)
        t50 = om*t24
        t51 = t50*t27
        t56 = ampmu*(3+t41*t36*t10)
        t61 = amplambda*(2+t3*t36*t43)
        t62 = 2*t56+t61
        t63 = cos(t31)
        t72 = t19*t24
        t73 = t72*t27
        t75 = omstry*y
        t76 = sin(t75)
        t78 = 1+t76/2
        t80 = om*x+ph
        t81 = sin(t80)
        t82 = t78*t81
        t84 = om*(y-t14)
        t85 = cos(t84)
        t87 = t85*om*t27
        t88 = t82*t87
        t91 = cos(t80)
        t92 = t91*t19
        t96 = omstrz*z
        t97 = sin(t96)
        t99 = 1+t97/2
        t100 = t99*t81
        t102 = om*(z-t14)
        t103 = cos(t102)
        t104 = t24*t103
        t105 = t104*om
        t106 = t100*t105
        t113 = ampmu*t41
        t114 = t113*t6
        t115 = omm*t10
        t116 = t34*t91
        t117 = sin(t84)
        t119 = om*t117*t27
        t121 = t78*t17
        t122 = cos(t23)
        t123 = t122*om
        t124 = t123*t27
        t126 = t116*t119+t121*t124
        t129 = t19*t85
        t132 = cos(t75)
        t133 = t132*omstry
        t142 = t113*omm
        t143 = t36*t43
        t144 = sin(t102)
        t145 = t50*t144
        t147 = t99*t17
        t148 = cos(t26)
        t150 = t24*t148*om
        t152 = t116*t145+t147*t150
        t157 = cos(t96)
        t158 = t157*omstrz
        forces(1) = -t13*t17*t21*t28-t34*((-2*t35*t38+t45)*t34*t48*t51+t
     #62*t63*omstrx*t48*om*t28/2-t62*t34*t17*t73+t45*t88+t61*t78*t92*t85
     #*t27+t45*t106+t61*t99*t92*t104)-t78*(t114*t115*t126+t56*(t116*t129
     #*t27+t133*t17*t124/2-t121*t73))-t99*(t142*t143*t152+t56*(t116*t72*
     #t103+t158*t17*t150/2-t147*t73))
        t167 = t13*t81
        t168 = t117*t19
        t172 = t35*omm
        t173 = t36*t10
        t177 = t63*omstrx*t91
        t180 = t34*t81
        t181 = t168*t27
        t184 = t19*t122
        t185 = t184*t27
        t191 = t6*omm
        t195 = amplambda*t3
        t214 = t34*t48
        t215 = t214*t51
        t220 = t184*t103
        t226 = t123*t144
        t229 = t117*t148*om
        t231 = t82*t226+t100*t229
        forces(2) = -t167*t168*t20*t27-t34*(-t172*t173*t126+t56*(t177*t1
     #19/2-t180*t181+t78*t48*t185))-t78*((2*t113*t191*t10+t195*t191*t43)
     #*t78*t81*t87+t62*t132*omstry*t81*t85*om*t27/2-t62*t78*t81*t181+t19
     #5*t6*omm*t43*(t215+t106)+t61*(t214*t185+t100*t220))-t99*(t142*t143
     #*t231+t56*(t82*t220+t158*t81*t229/2-t100*t181))
        t250 = t72*t144
        t253 = t72*t148
        t265 = t129*t148
        forces(3) = -t167*t24*t144*t21-t34*(-t172*t173*t152+t56*(t177*t1
     #45/2-t180*t250+t99*t48*t253))-t78*(t114*t115*t231+t56*(t133*t81*t2
     #26/2-t82*t250+t100*t265))-t99*((2*t113*t44-t195*t38)*t99*t81*t105+
     #t62*t157*omstrz*t81*t24*t103*om/2-t62*t99*t81*t250-t195*omm*t173*(
     #t215+t88)+t61*(t214*t253+t82*t265))
        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)
       enddo
       enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine forcingttfortcsg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     xx, yy, zz, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda, omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t102
      doubleprecision t103
      doubleprecision t106
      doubleprecision t107
      doubleprecision t110
      doubleprecision t120
      doubleprecision t121
      doubleprecision t122
      doubleprecision t123
      doubleprecision t124
      doubleprecision t125
      doubleprecision t127
      doubleprecision t129
      doubleprecision t13
      doubleprecision t131
      doubleprecision t133
      doubleprecision t135
      doubleprecision t138
      doubleprecision t14
      doubleprecision t142
      doubleprecision t143
      doubleprecision t150
      doubleprecision t151
      doubleprecision t157
      doubleprecision t158
      doubleprecision t159
      doubleprecision t16
      doubleprecision t160
      doubleprecision t161
      doubleprecision t163
      doubleprecision t165
      doubleprecision t166
      doubleprecision t168
      doubleprecision t17
      doubleprecision t171
      doubleprecision t174
      doubleprecision t175
      doubleprecision t187
      doubleprecision t188
      doubleprecision t19
      doubleprecision t192
      doubleprecision t193
      doubleprecision t197
      doubleprecision t198
      doubleprecision t2
      doubleprecision t20
      doubleprecision t203
      doubleprecision t21
      doubleprecision t212
      doubleprecision t216
      doubleprecision t22
      doubleprecision t223
      doubleprecision t23
      doubleprecision t236
      doubleprecision t238
      doubleprecision t242
      doubleprecision t246
      doubleprecision t25
      doubleprecision t252
      doubleprecision t256
      doubleprecision t259
      doubleprecision t26
      doubleprecision t279
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t30
      doubleprecision t303
      doubleprecision t315
      doubleprecision t33
      doubleprecision t34
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t39
      doubleprecision t40
      doubleprecision t43
      doubleprecision t45
      doubleprecision t46
      doubleprecision t47
      doubleprecision t5
      doubleprecision t50
      doubleprecision t52
      doubleprecision t53
      doubleprecision t54
      doubleprecision t59
      doubleprecision t6
      doubleprecision t64
      doubleprecision t65
      doubleprecision t66
      doubleprecision t74
      doubleprecision t77
      doubleprecision t78
      doubleprecision t80
      doubleprecision t82
      doubleprecision t83
      doubleprecision t84
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t90
      doubleprecision t92
      doubleprecision t96
      doubleprecision t99

      do k=kfirst,klast
         do j=jfirst,jlast
           do i=ifirst,ilast
             x = xx(i,j,k)
             y = yy(i,j,k)
             z = zz(i,j,k)
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
        t30 = t26*t29
        t33 = omstrx*x
        t34 = sin(t33)
        t36 = 1+t34/2
        t37 = ampmu*t3
        t38 = sin(t5)
        t39 = omm*t38
        t40 = t39*t10
        t43 = cos(t2)
        t45 = cos(t9)
        t46 = t39*t45
        t47 = amplambda*t43*t46
        t50 = cos(t16)
        t52 = t19*om
        t53 = t52*t21
        t54 = t53*t30
        t59 = ampmu*(3+t43*t38*t10)
        t64 = amplambda*(2+t3*t38*t45)
        t65 = 2*t59+t64
        t66 = cos(t33)
        t74 = t20*t21
        t77 = omstry*y
        t78 = sin(t77)
        t80 = 1+t78/2
        t82 = om*x+ph
        t83 = sin(t82)
        t84 = t80*t83
        t86 = om*(y-t14)
        t87 = cos(t86)
        t88 = t84*t87
        t90 = t88*t53*t29
        t92 = cos(t82)
        t96 = t21*t29
        t99 = omstrz*z
        t100 = sin(t99)
        t102 = 1+t100/2
        t103 = t102*t83
        t106 = om*(z-t14)
        t107 = cos(t106)
        t110 = t103*t26*t107*t52*t21
        t120 = ampmu*t43
        t121 = t120*t6
        t122 = omm*t10
        t123 = t36*t92
        t124 = t123*t52
        t125 = sin(t86)
        t127 = t125*t21*t29
        t129 = t80*t17
        t131 = cos(t25)
        t133 = t21*t131*t29
        t135 = -t124*t127-t129*t52*t133
        t138 = t123*t20
        t142 = cos(t77)
        t143 = t142*omstry
        t150 = t21*t26
        t151 = t150*t29
        t157 = t120*omm
        t158 = t38*t45
        t159 = sin(t106)
        t160 = t26*t159
        t161 = t160*t21
        t163 = t102*t17
        t165 = cos(t28)
        t166 = t150*t165
        t168 = -t124*t161-t163*t52*t166
        t171 = t26*t107
        t174 = cos(t99)
        t175 = t174*omstrz
        forces(1) = t13*t17*t23*t30-t36*(-(-2*t37*t40+t47)*t36*t50*t54-t
     #65*t66*omstrx*t50*t54/2+t65*t36*t17*t74*t30-t47*t90-t64*t80*t92*t2
     #0*t87*t96-t47*t110-t64*t102*t92*t20*t26*t107*t21)-t80*(t121*t122*t
     #135+t59*(-t138*t87*t21*t29-t143*t17*t53*t131*t29/2+t129*t20*t151))
     #-t102*(t157*t158*t168+t59*(-t138*t171*t21-t175*t17*t53*t26*t165/2+
     #t163*t20*t151))
        t187 = t13*t83
        t188 = t125*t20
        t192 = t37*omm
        t193 = t38*t10
        t197 = t66*omstrx*t92
        t198 = t52*t125
        t203 = t36*t83*t20
        t212 = t6*omm
        t216 = amplambda*t3
        t223 = t87*t52*t96
        t236 = t36*t50
        t238 = t236*t52*t151
        t242 = t236*t20
        t246 = t20*t107*t21
        t252 = t84*t131
        t256 = t103*t125
        t259 = -t252*t52*t159*t21-t256*t53*t165
        forces(2) = t187*t188*t22*t29-t36*(-t192*t193*t135+t59*(-t197*t1
     #98*t96/2+t203*t127-t80*t50*t20*t133))-t80*(-(2*t120*t212*t10+t216*
     #t212*t45)*t80*t83*t223-t65*t142*omstry*t83*t223/2+t65*t80*t83*t188
     #*t96+t216*t6*omm*t45*(-t238-t110)+t64*(-t242*t133-t103*t131*t246))
     #-t102*(t157*t158*t259+t59*(-t252*t246-t175*t83*t198*t21*t165/2+t25
     #6*t74*t29))
        t279 = t159*t21
        t303 = t74*t165
        t315 = t171*t53
        forces(3) = t187*t160*t23-t36*(-t192*t193*t168+t59*(-t197*t52*t2
     #6*t279/2+t203*t161-t102*t50*t20*t166))-t80*(t121*t122*t259+t59*(-t
     #143*t83*t131*t52*t279/2+t84*t26*t20*t159*t21-t103*t87*t303))-t102*
     #(-(2*t120*t46-t216*t40)*t102*t83*t315-t65*t174*omstrz*t83*t315/2+t
     #65*t102*t83*t160*t74-t216*omm*t193*(-t238-t90)+t64*(-t242*t166-t88
     #*t303))
        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)
      enddo
      enddo
      enddo
        return
      end
