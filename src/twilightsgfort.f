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
      subroutine forcingfortsg( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     h, zmin, omstrx, omstry, omstrz ) bind(c)

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

c-----------------------------------------------------------------------
      subroutine forcingfortsgatt( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, omega, c, phase, momega, mphase, amprho, ampmu,
     +      amplambda, h, zmin, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, omega, c, phase, momega, mphase
      doubleprecision amprho, zmin, ampmu, amplambda, omstrx, omstry
      doubleprecision omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t100
      doubleprecision t103
      doubleprecision t104
      doubleprecision t105
      doubleprecision t106
      doubleprecision t107
      doubleprecision t11
      doubleprecision t110
      doubleprecision t111
      doubleprecision t116
      doubleprecision t117
      doubleprecision t119
      doubleprecision t120
      doubleprecision t122
      doubleprecision t125
      doubleprecision t128
      doubleprecision t132
      doubleprecision t133
      doubleprecision t134
      doubleprecision t135
      doubleprecision t137
      doubleprecision t138
      doubleprecision t14
      doubleprecision t140
      doubleprecision t141
      doubleprecision t142
      doubleprecision t144
      doubleprecision t148
      doubleprecision t15
      doubleprecision t151
      doubleprecision t152
      doubleprecision t16
      doubleprecision t164
      doubleprecision t167
      doubleprecision t17
      doubleprecision t176
      doubleprecision t18
      doubleprecision t183
      doubleprecision t194
      doubleprecision t197
      doubleprecision t198
      doubleprecision t2
      doubleprecision t20
      doubleprecision t203
      doubleprecision t210
      doubleprecision t212
      doubleprecision t214
      doubleprecision t215
      doubleprecision t217
      doubleprecision t219
      doubleprecision t24
      doubleprecision t243
      doubleprecision t26
      doubleprecision t265
      doubleprecision t27
      doubleprecision t272
      doubleprecision t28
      doubleprecision t31
      doubleprecision t32
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t41
      doubleprecision t42
      doubleprecision t44
      doubleprecision t50
      doubleprecision t56
      doubleprecision t57
      doubleprecision t58
      doubleprecision t6
      doubleprecision t64
      doubleprecision t67
      doubleprecision t7
      doubleprecision t73
      doubleprecision t74
      doubleprecision t76
      doubleprecision t77
      doubleprecision t78
      doubleprecision t8
      doubleprecision t81
      doubleprecision t82
      doubleprecision t85
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t92
      doubleprecision t95
      doubleprecision t97
      doubleprecision t98

      do k=kfirst,klast
         z = (k-1)*h + zmin
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t1 = omstrx*x
        t2 = sin(t1)
        t4 = 1+t2/2
        t6 = momega*x+mphase
        t7 = sin(t6)
        t8 = ampmu*t7
        t10 = momega*y+mphase
        t11 = cos(t10)
        t14 = momega*z+mphase
        t15 = sin(t14)
        t16 = momega*t11*t15
        t17 = t8*t16
        t18 = cos(t6)
        t20 = amplambda*t18*t16
        t24 = c*t
        t26 = omega*(x-t24)
        t27 = -t26-phase
        t28 = sin(t27)
        t31 = omega*x+phase
        t32 = sin(t31)
        t35 = -omega*(z-t24)-phase
        t36 = cos(t35)
        t37 = t32*t36
        t39 = cos(t27)
        t40 = cos(t31)
        t41 = t39*t40
        t42 = omega*t36
        t44 = t28*omega*t37+t41*t42
        t50 = ampmu*(3.D0/2.D0+t18*t11*t15/2)
        t56 = amplambda*(1.D0/2.D0+t7*t11*t15/4)
        t57 = 2*t50+t56
        t58 = cos(t1)
        t64 = omega**2
        t67 = t28*t64
        t73 = omstry*y
        t74 = sin(t73)
        t76 = 1+t74/2
        t77 = sin(t26)
        t78 = t76*t77
        t81 = -omega*(y-t24)-phase
        t82 = sin(t81)
        t85 = omega*z+phase
        t86 = cos(t85)
        t87 = t82*omega*t86
        t88 = t78*t87
        t92 = cos(t26)
        t95 = t92*t64*t82*t86
        t97 = omstrz*z
        t98 = sin(t97)
        t100 = 1+t98/2
        t103 = omega*y+phase
        t104 = cos(t103)
        t105 = t104*t36
        t106 = t105*omega
        t107 = t100*t40*t106
        t110 = t56*t100
        t111 = t32*t64
        t116 = ampmu*t18
        t117 = sin(t10)
        t119 = t117*momega*t15
        t120 = t116*t119
        t122 = cos(t81)
        t125 = t4*t92*omega*t122*t86
        t128 = t50*t4
        t132 = t116*t11
        t133 = cos(t14)
        t134 = t133*momega
        t135 = t4*t32
        t137 = sin(t35)
        t138 = omega*t104*t137
        t140 = t100*t39
        t141 = t32*t137
        t142 = t141*omega
        t144 = t135*t138+t140*t142
        t148 = t64*t104
        t151 = cos(t97)
        t152 = t151*omstrz
        forces(1) = t4*((-t17+t20/4)*t4*t44+t57*t58*omstrx*t44/2+t57*t4*
     #(-2*t39*t64*t37+2*t67*t40*t36)+t20*t88/4+t56*t76*t95+t20*t107/4-t1
     #10*t111*t105)+t76*(-t120*t125/2+t128*t95)+t100*(t132*t134*t144/2+t
     #50*(-t135*t148*t36+t152*t39*t142/2-t140*t37*t64))
        t164 = t58*omstrx
        t167 = t122*t86
        t176 = amplambda*t7
        t183 = cos(t73)
        t194 = t122*t64*t86
        t197 = momega*t15
        t198 = t4*t44
        t203 = sin(t103)
        t210 = t76*t40
        t212 = t203*omega*t137
        t214 = t100*t77
        t215 = sin(t85)
        t217 = t122*t215*omega
        t219 = t210*t212-t214*t217
        forces(2) = t4*(-t17*t125/2+t50*t164*t92*omega*t167/2-t128*t77*t
     #64*t167)+t76*((-t120-t176*t119/4)*t76*t77*t87+t57*t183*omstry*t77*
     #t82*omega*t86/2-t57*t76*t77*t194-t176*t117*t197*(t198+t107)/4-t110
     #*t40*t203*t64*t36)+t100*(t132*t134*t219/2+t50*(-t210*t203*t64*t36-
     #t152*t77*t217/2-t214*t194))
        t243 = t148*t137
        t265 = t82*t64*t215
        t272 = t11*t133*momega
        forces(3) = t4*(-t8*momega*t11*t15*t144/2+t50*(t164*t32*t138/2+t
     #4*t40*t243+t100*t28*t111*t137+t140*t40*t64*t137))+t76*(-t116*t117*
     #t197*t219/2+t50*(t183*omstry*t40*t212/2+t210*t243-t214*t265))+t100
     #*((t116*t272+t176*t272/4)*t100*t40*t106+t57*t151*omstrz*t40*t104*t
     #42/2+t57*t100*t40*t243+t176*t11*t134*(t198+t88)/4+t56*(t4*(t67*t14
     #1+t41*t64*t137)-t78*t265))

        fo(1,i,j,k) = fo(1,i,j,k)+forces(1)
        fo(2,i,j,k) = fo(2,i,j,k)+forces(2)
        fo(3,i,j,k) = fo(3,i,j,k)+forces(3)
      enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine forcingttattfortsg( ifirst, ilast, jfirst, jlast, 
     +    kfirst, klast, fo, t, omega, c, phase, momega, mphase, amprho,
     +     ampmu, amplambda, h, zmin, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, omega, c, phase, momega, mphase 
      doubleprecision amprho, zmin, ampmu, amplambda
      doubleprecision omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t100
      doubleprecision t101
      doubleprecision t104
      doubleprecision t105
      doubleprecision t107
      doubleprecision t108
      doubleprecision t109
      doubleprecision t11
      doubleprecision t110
      doubleprecision t114
      doubleprecision t115
      doubleprecision t116
      doubleprecision t120
      doubleprecision t121
      doubleprecision t124
      doubleprecision t125
      doubleprecision t127
      doubleprecision t128
      doubleprecision t130
      doubleprecision t131
      doubleprecision t135
      doubleprecision t14
      doubleprecision t141
      doubleprecision t146
      doubleprecision t147
      doubleprecision t149
      doubleprecision t15
      doubleprecision t150
      doubleprecision t151
      doubleprecision t153
      doubleprecision t155
      doubleprecision t157
      doubleprecision t16
      doubleprecision t161
      doubleprecision t164
      doubleprecision t167
      doubleprecision t168
      doubleprecision t169
      doubleprecision t17
      doubleprecision t170
      doubleprecision t172
      doubleprecision t173
      doubleprecision t175
      doubleprecision t179
      doubleprecision t18
      doubleprecision t183
      doubleprecision t188
      doubleprecision t191
      doubleprecision t192
      doubleprecision t194
      doubleprecision t2
      doubleprecision t20
      doubleprecision t202
      doubleprecision t204
      doubleprecision t205
      doubleprecision t207
      doubleprecision t214
      doubleprecision t215
      doubleprecision t224
      doubleprecision t228
      doubleprecision t230
      doubleprecision t234
      doubleprecision t237
      doubleprecision t238
      doubleprecision t24
      doubleprecision t245
      doubleprecision t253
      doubleprecision t254
      doubleprecision t26
      doubleprecision t260
      doubleprecision t266
      doubleprecision t267
      doubleprecision t27
      doubleprecision t271
      doubleprecision t273
      doubleprecision t274
      doubleprecision t277
      doubleprecision t279
      doubleprecision t28
      doubleprecision t282
      doubleprecision t29
      doubleprecision t297
      doubleprecision t30
      doubleprecision t300
      doubleprecision t31
      doubleprecision t314
      doubleprecision t32
      doubleprecision t34
      doubleprecision t35
      doubleprecision t352
      doubleprecision t359
      doubleprecision t36
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t41
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t46
      doubleprecision t48
      doubleprecision t49
      doubleprecision t50
      doubleprecision t52
      doubleprecision t55
      doubleprecision t6
      doubleprecision t61
      doubleprecision t67
      doubleprecision t68
      doubleprecision t69
      doubleprecision t7
      doubleprecision t75
      doubleprecision t76
      doubleprecision t78
      doubleprecision t8
      doubleprecision t82
      doubleprecision t85
      doubleprecision t86
      doubleprecision t88
      doubleprecision t89
      doubleprecision t90
      doubleprecision t94
      doubleprecision t95
      doubleprecision t96
      doubleprecision t98
      doubleprecision t99

      do k=kfirst,klast
         z = (k-1)*h + zmin
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t1 = omstrx*x
        t2 = sin(t1)
        t4 = 1+t2/2
        t6 = momega*x+mphase
        t7 = sin(t6)
        t8 = ampmu*t7
        t10 = momega*y+mphase
        t11 = cos(t10)
        t14 = momega*z+mphase
        t15 = sin(t14)
        t16 = momega*t11*t15
        t17 = t8*t16
        t18 = cos(t6)
        t20 = amplambda*t18*t16
        t24 = c*t
        t26 = omega*(x-t24)
        t27 = -t26-phase
        t28 = sin(t27)
        t29 = omega**2
        t30 = t29*omega
        t31 = t28*t30
        t32 = c**2
        t34 = omega*x+phase
        t35 = sin(t34)
        t36 = t32*t35
        t39 = -omega*(z-t24)-phase
        t40 = cos(t39)
        t41 = t36*t40
        t43 = cos(t27)
        t44 = t43*t30
        t45 = sin(t39)
        t46 = t36*t45
        t48 = cos(t34)
        t49 = t32*t48
        t50 = t49*t40
        t52 = t49*t45
        t55 = -2*t31*t41-2*t44*t46-2*t44*t50+2*t31*t52
        t61 = ampmu*(3.D0/2.D0+t18*t11*t15/2)
        t67 = amplambda*(1.D0/2.D0+t7*t11*t15/4)
        t68 = 2*t61+t67
        t69 = cos(t1)
        t75 = t29**2
        t76 = t43*t75
        t78 = t28*t75
        t82 = t76*t41-t78*t46-t78*t50-t76*t52
        t85 = omstry*y
        t86 = sin(t85)
        t88 = 1+t86/2
        t89 = sin(t26)
        t90 = t88*t89
        t94 = -omega*(y-t24)-phase
        t95 = sin(t94)
        t96 = t32*t95
        t98 = omega*z+phase
        t99 = cos(t98)
        t100 = t96*t99
        t101 = t90*t30*t100
        t104 = cos(t26)
        t105 = t88*t104
        t107 = cos(t94)
        t108 = t32*t107
        t109 = t108*t99
        t110 = t105*t30*t109
        t114 = t75*t32
        t115 = t95*t99
        t116 = t114*t115
        t120 = t107*t99
        t121 = t114*t120
        t124 = omstrz*z
        t125 = sin(t124)
        t127 = 1+t125/2
        t128 = t127*t48
        t130 = omega*y+phase
        t131 = cos(t130)
        t135 = t128*t131*t40*t30*t32
        t141 = t40*t32
        t146 = ampmu*t18
        t147 = sin(t10)
        t149 = t147*momega*t15
        t150 = t146*t149
        t151 = t4*t104
        t153 = t151*t30*t109
        t155 = t4*t89
        t157 = t155*t30*t100
        t161 = 2*t61*t151*t116
        t164 = 2*t61*t155*t121
        t167 = t146*t11
        t168 = cos(t14)
        t169 = t168*momega
        t170 = t4*t35
        t172 = t131*t45
        t173 = t172*t32
        t175 = t127*t43
        t179 = t127*t28
        t183 = -t170*t30*t173-2*t175*t30*t46-2*t179*t30*t41
        t188 = t131*t40
        t191 = cos(t124)
        t192 = t191*omstrz
        t194 = t30*t32
        t202 = t175*t75
        t204 = 2*t202*t41
        t205 = t179*t75
        t207 = 2*t205*t46
        forces(1) = t4*((-t17+t20/4)*t4*t55+t68*t69*omstrx*t55/2+4*t68*t
     #4*t82-t20*t101/2-t20*t110/2-2*t67*t105*t116+2*t67*t90*t121-t20*t13
     #5/4+t67*t127*t35*t75*t131*t141)+t88*(t150*t153+t150*t157-t161+t164
     #)+t127*(t167*t169*t183/2+t61*(t170*t75*t188*t32-t192*t43*t194*t35*
     #t45-t192*t28*t194*t35*t40+t204-t207))
        t214 = t69*omstrx
        t215 = t61*t214
        t224 = amplambda*t7
        t228 = (-t150-t224*t149/4)*t88
        t230 = t194*t115
        t234 = t194*t120
        t237 = cos(t85)
        t238 = t68*t237
        t245 = t68*t88
        t253 = momega*t15
        t254 = t4*t55
        t260 = sin(t130)
        t266 = t88*t48
        t267 = t266*t260
        t271 = t127*t89
        t273 = sin(t98)
        t274 = t108*t273
        t277 = t127*t104
        t279 = t96*t273
        t282 = -t267*t30*t45*t32+2*t271*t30*t274-2*t277*t30*t279
        t297 = t271*t75
        t300 = t277*t75
        forces(2) = t4*(t17*t153+t17*t157-t215*t104*t30*t109-t215*t89*t3
     #0*t100-t161+t164)+t88*(-2*t228*t89*t230-2*t228*t104*t234-t238*omst
     #ry*t89*t230-t238*omstry*t104*t234+2*t245*t89*t121-2*t245*t104*t116
     #-t224*t147*t253*(t254-t135)/4+t67*t128*t260*t75*t141)+t127*(t167*t
     #169*t282/2+t61*(t267*t75*t40*t32+t192*t89*t194*t107*t273-t192*t104
     #*t194*t95*t273+2*t297*t109-2*t300*t100))
        t314 = t45*t32
        t352 = t11*t168*momega
        t359 = t188*t194
        forces(3) = t4*(-t8*momega*t11*t15*t183/2+t61*(-t214*t35*t30*t13
     #1*t314/2-t4*t48*t75*t173+t204-t207-2*t202*t52-2*t205*t50))+t88*(-t
     #146*t147*t253*t282/2+t61*(-t237*omstry*t48*t260*t30*t314/2-t266*t1
     #31*t75*t45*t32+2*t297*t279+2*t300*t274))+t127*(-(t146*t352+t224*t3
     #52/4)*t127*t48*t359-t68*t191*omstrz*t48*t359/2-t68*t127*t48*t172*t
     #114+t224*t11*t169*(t254-2*t101-2*t110)/4+t67*(2*t4*t82+2*t90*t75*t
     #279+2*t105*t75*t274))
        fo(1,i,j,k) = fo(1,i,j,k)+forces(1)
        fo(2,i,j,k) = fo(2,i,j,k)+forces(2)
        fo(3,i,j,k) = fo(3,i,j,k)+forces(3)
      enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine forcingfortsgattc(ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, omega, c, phase, momega, mphase, amprho, ampmu,
     +      amplambda, xx, yy, zz, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, omega, c, phase, momega, mphase
      doubleprecision amprho, ampmu, amplambda, omstrx, omstry
      doubleprecision omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t100
      doubleprecision t103
      doubleprecision t104
      doubleprecision t105
      doubleprecision t106
      doubleprecision t107
      doubleprecision t11
      doubleprecision t110
      doubleprecision t111
      doubleprecision t116
      doubleprecision t117
      doubleprecision t119
      doubleprecision t120
      doubleprecision t122
      doubleprecision t125
      doubleprecision t128
      doubleprecision t132
      doubleprecision t133
      doubleprecision t134
      doubleprecision t135
      doubleprecision t137
      doubleprecision t138
      doubleprecision t14
      doubleprecision t140
      doubleprecision t141
      doubleprecision t142
      doubleprecision t144
      doubleprecision t148
      doubleprecision t15
      doubleprecision t151
      doubleprecision t152
      doubleprecision t16
      doubleprecision t164
      doubleprecision t167
      doubleprecision t17
      doubleprecision t176
      doubleprecision t18
      doubleprecision t183
      doubleprecision t194
      doubleprecision t197
      doubleprecision t198
      doubleprecision t2
      doubleprecision t20
      doubleprecision t203
      doubleprecision t210
      doubleprecision t212
      doubleprecision t214
      doubleprecision t215
      doubleprecision t217
      doubleprecision t219
      doubleprecision t24
      doubleprecision t243
      doubleprecision t26
      doubleprecision t265
      doubleprecision t27
      doubleprecision t272
      doubleprecision t28
      doubleprecision t31
      doubleprecision t32
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t41
      doubleprecision t42
      doubleprecision t44
      doubleprecision t50
      doubleprecision t56
      doubleprecision t57
      doubleprecision t58
      doubleprecision t6
      doubleprecision t64
      doubleprecision t67
      doubleprecision t7
      doubleprecision t73
      doubleprecision t74
      doubleprecision t76
      doubleprecision t77
      doubleprecision t78
      doubleprecision t8
      doubleprecision t81
      doubleprecision t82
      doubleprecision t85
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t92
      doubleprecision t95
      doubleprecision t97
      doubleprecision t98

      do k=kfirst,klast
         do j=jfirst,jlast
           do i=ifirst,ilast
             x = xx(i,j,k)
             y = yy(i,j,k)
             z = zz(i,j,k)

        t1 = omstrx*x
        t2 = sin(t1)
        t4 = 1+t2/2
        t6 = momega*x+mphase
        t7 = sin(t6)
        t8 = ampmu*t7
        t10 = momega*y+mphase
        t11 = cos(t10)
        t14 = momega*z+mphase
        t15 = sin(t14)
        t16 = momega*t11*t15
        t17 = t8*t16
        t18 = cos(t6)
        t20 = amplambda*t18*t16
        t24 = c*t
        t26 = omega*(x-t24)
        t27 = -t26-phase
        t28 = sin(t27)
        t31 = omega*x+phase
        t32 = sin(t31)
        t35 = -omega*(z-t24)-phase
        t36 = cos(t35)
        t37 = t32*t36
        t39 = cos(t27)
        t40 = cos(t31)
        t41 = t39*t40
        t42 = omega*t36
        t44 = t28*omega*t37+t41*t42
        t50 = ampmu*(3.D0/2.D0+t18*t11*t15/2)
        t56 = amplambda*(1.D0/2.D0+t7*t11*t15/4)
        t57 = 2*t50+t56
        t58 = cos(t1)
        t64 = omega**2
        t67 = t28*t64
        t73 = omstry*y
        t74 = sin(t73)
        t76 = 1+t74/2
        t77 = sin(t26)
        t78 = t76*t77
        t81 = -omega*(y-t24)-phase
        t82 = sin(t81)
        t85 = omega*z+phase
        t86 = cos(t85)
        t87 = t82*omega*t86
        t88 = t78*t87
        t92 = cos(t26)
        t95 = t92*t64*t82*t86
        t97 = omstrz*z
        t98 = sin(t97)
        t100 = 1+t98/2
        t103 = omega*y+phase
        t104 = cos(t103)
        t105 = t104*t36
        t106 = t105*omega
        t107 = t100*t40*t106
        t110 = t56*t100
        t111 = t32*t64
        t116 = ampmu*t18
        t117 = sin(t10)
        t119 = t117*momega*t15
        t120 = t116*t119
        t122 = cos(t81)
        t125 = t4*t92*omega*t122*t86
        t128 = t50*t4
        t132 = t116*t11
        t133 = cos(t14)
        t134 = t133*momega
        t135 = t4*t32
        t137 = sin(t35)
        t138 = omega*t104*t137
        t140 = t100*t39
        t141 = t32*t137
        t142 = t141*omega
        t144 = t135*t138+t140*t142
        t148 = t64*t104
        t151 = cos(t97)
        t152 = t151*omstrz
        forces(1) = t4*((-t17+t20/4)*t4*t44+t57*t58*omstrx*t44/2+t57*t4*
     #(-2*t39*t64*t37+2*t67*t40*t36)+t20*t88/4+t56*t76*t95+t20*t107/4-t1
     #10*t111*t105)+t76*(-t120*t125/2+t128*t95)+t100*(t132*t134*t144/2+t
     #50*(-t135*t148*t36+t152*t39*t142/2-t140*t37*t64))
        t164 = t58*omstrx
        t167 = t122*t86
        t176 = amplambda*t7
        t183 = cos(t73)
        t194 = t122*t64*t86
        t197 = momega*t15
        t198 = t4*t44
        t203 = sin(t103)
        t210 = t76*t40
        t212 = t203*omega*t137
        t214 = t100*t77
        t215 = sin(t85)
        t217 = t122*t215*omega
        t219 = t210*t212-t214*t217
        forces(2) = t4*(-t17*t125/2+t50*t164*t92*omega*t167/2-t128*t77*t
     #64*t167)+t76*((-t120-t176*t119/4)*t76*t77*t87+t57*t183*omstry*t77*
     #t82*omega*t86/2-t57*t76*t77*t194-t176*t117*t197*(t198+t107)/4-t110
     #*t40*t203*t64*t36)+t100*(t132*t134*t219/2+t50*(-t210*t203*t64*t36-
     #t152*t77*t217/2-t214*t194))
        t243 = t148*t137
        t265 = t82*t64*t215
        t272 = t11*t133*momega
        forces(3) = t4*(-t8*momega*t11*t15*t144/2+t50*(t164*t32*t138/2+t
     #4*t40*t243+t100*t28*t111*t137+t140*t40*t64*t137))+t76*(-t116*t117*
     #t197*t219/2+t50*(t183*omstry*t40*t212/2+t210*t243-t214*t265))+t100
     #*((t116*t272+t176*t272/4)*t100*t40*t106+t57*t151*omstrz*t40*t104*t
     #42/2+t57*t100*t40*t243+t176*t11*t134*(t198+t88)/4+t56*(t4*(t67*t14
     #1+t41*t64*t137)-t78*t265))

        fo(1,i,j,k) = fo(1,i,j,k)+forces(1)
        fo(2,i,j,k) = fo(2,i,j,k)+forces(2)
        fo(3,i,j,k) = fo(3,i,j,k)+forces(3)
      enddo
      enddo
      enddo
      end

c-----------------------------------------------------------------------
      subroutine forcingttattfortsgc( ifirst, ilast, jfirst, jlast, 
     +    kfirst, klast, fo, t, omega, c, phase, momega, mphase, amprho,
     +     ampmu, amplambda, xx, yy, zz, omstrx, omstry, omstrz )

      implicit none

      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, omega, c, phase, momega, mphase 
      doubleprecision amprho, ampmu, amplambda
      doubleprecision omstrx, omstry, omstrz
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 zz(ifirst:ilast,jfirst:jlast,kfirst:klast)

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t100
      doubleprecision t101
      doubleprecision t104
      doubleprecision t105
      doubleprecision t107
      doubleprecision t108
      doubleprecision t109
      doubleprecision t11
      doubleprecision t110
      doubleprecision t114
      doubleprecision t115
      doubleprecision t116
      doubleprecision t120
      doubleprecision t121
      doubleprecision t124
      doubleprecision t125
      doubleprecision t127
      doubleprecision t128
      doubleprecision t130
      doubleprecision t131
      doubleprecision t135
      doubleprecision t14
      doubleprecision t141
      doubleprecision t146
      doubleprecision t147
      doubleprecision t149
      doubleprecision t15
      doubleprecision t150
      doubleprecision t151
      doubleprecision t153
      doubleprecision t155
      doubleprecision t157
      doubleprecision t16
      doubleprecision t161
      doubleprecision t164
      doubleprecision t167
      doubleprecision t168
      doubleprecision t169
      doubleprecision t17
      doubleprecision t170
      doubleprecision t172
      doubleprecision t173
      doubleprecision t175
      doubleprecision t179
      doubleprecision t18
      doubleprecision t183
      doubleprecision t188
      doubleprecision t191
      doubleprecision t192
      doubleprecision t194
      doubleprecision t2
      doubleprecision t20
      doubleprecision t202
      doubleprecision t204
      doubleprecision t205
      doubleprecision t207
      doubleprecision t214
      doubleprecision t215
      doubleprecision t224
      doubleprecision t228
      doubleprecision t230
      doubleprecision t234
      doubleprecision t237
      doubleprecision t238
      doubleprecision t24
      doubleprecision t245
      doubleprecision t253
      doubleprecision t254
      doubleprecision t26
      doubleprecision t260
      doubleprecision t266
      doubleprecision t267
      doubleprecision t27
      doubleprecision t271
      doubleprecision t273
      doubleprecision t274
      doubleprecision t277
      doubleprecision t279
      doubleprecision t28
      doubleprecision t282
      doubleprecision t29
      doubleprecision t297
      doubleprecision t30
      doubleprecision t300
      doubleprecision t31
      doubleprecision t314
      doubleprecision t32
      doubleprecision t34
      doubleprecision t35
      doubleprecision t352
      doubleprecision t359
      doubleprecision t36
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t41
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t46
      doubleprecision t48
      doubleprecision t49
      doubleprecision t50
      doubleprecision t52
      doubleprecision t55
      doubleprecision t6
      doubleprecision t61
      doubleprecision t67
      doubleprecision t68
      doubleprecision t69
      doubleprecision t7
      doubleprecision t75
      doubleprecision t76
      doubleprecision t78
      doubleprecision t8
      doubleprecision t82
      doubleprecision t85
      doubleprecision t86
      doubleprecision t88
      doubleprecision t89
      doubleprecision t90
      doubleprecision t94
      doubleprecision t95
      doubleprecision t96
      doubleprecision t98
      doubleprecision t99

      do k=kfirst,klast
         do j=jfirst,jlast
           do i=ifirst,ilast
             x = xx(i,j,k)
             y = yy(i,j,k)
             z = zz(i,j,k)

        t1 = omstrx*x
        t2 = sin(t1)
        t4 = 1+t2/2
        t6 = momega*x+mphase
        t7 = sin(t6)
        t8 = ampmu*t7
        t10 = momega*y+mphase
        t11 = cos(t10)
        t14 = momega*z+mphase
        t15 = sin(t14)
        t16 = momega*t11*t15
        t17 = t8*t16
        t18 = cos(t6)
        t20 = amplambda*t18*t16
        t24 = c*t
        t26 = omega*(x-t24)
        t27 = -t26-phase
        t28 = sin(t27)
        t29 = omega**2
        t30 = t29*omega
        t31 = t28*t30
        t32 = c**2
        t34 = omega*x+phase
        t35 = sin(t34)
        t36 = t32*t35
        t39 = -omega*(z-t24)-phase
        t40 = cos(t39)
        t41 = t36*t40
        t43 = cos(t27)
        t44 = t43*t30
        t45 = sin(t39)
        t46 = t36*t45
        t48 = cos(t34)
        t49 = t32*t48
        t50 = t49*t40
        t52 = t49*t45
        t55 = -2*t31*t41-2*t44*t46-2*t44*t50+2*t31*t52
        t61 = ampmu*(3.D0/2.D0+t18*t11*t15/2)
        t67 = amplambda*(1.D0/2.D0+t7*t11*t15/4)
        t68 = 2*t61+t67
        t69 = cos(t1)
        t75 = t29**2
        t76 = t43*t75
        t78 = t28*t75
        t82 = t76*t41-t78*t46-t78*t50-t76*t52
        t85 = omstry*y
        t86 = sin(t85)
        t88 = 1+t86/2
        t89 = sin(t26)
        t90 = t88*t89
        t94 = -omega*(y-t24)-phase
        t95 = sin(t94)
        t96 = t32*t95
        t98 = omega*z+phase
        t99 = cos(t98)
        t100 = t96*t99
        t101 = t90*t30*t100
        t104 = cos(t26)
        t105 = t88*t104
        t107 = cos(t94)
        t108 = t32*t107
        t109 = t108*t99
        t110 = t105*t30*t109
        t114 = t75*t32
        t115 = t95*t99
        t116 = t114*t115
        t120 = t107*t99
        t121 = t114*t120
        t124 = omstrz*z
        t125 = sin(t124)
        t127 = 1+t125/2
        t128 = t127*t48
        t130 = omega*y+phase
        t131 = cos(t130)
        t135 = t128*t131*t40*t30*t32
        t141 = t40*t32
        t146 = ampmu*t18
        t147 = sin(t10)
        t149 = t147*momega*t15
        t150 = t146*t149
        t151 = t4*t104
        t153 = t151*t30*t109
        t155 = t4*t89
        t157 = t155*t30*t100
        t161 = 2*t61*t151*t116
        t164 = 2*t61*t155*t121
        t167 = t146*t11
        t168 = cos(t14)
        t169 = t168*momega
        t170 = t4*t35
        t172 = t131*t45
        t173 = t172*t32
        t175 = t127*t43
        t179 = t127*t28
        t183 = -t170*t30*t173-2*t175*t30*t46-2*t179*t30*t41
        t188 = t131*t40
        t191 = cos(t124)
        t192 = t191*omstrz
        t194 = t30*t32
        t202 = t175*t75
        t204 = 2*t202*t41
        t205 = t179*t75
        t207 = 2*t205*t46
        forces(1) = t4*((-t17+t20/4)*t4*t55+t68*t69*omstrx*t55/2+4*t68*t
     #4*t82-t20*t101/2-t20*t110/2-2*t67*t105*t116+2*t67*t90*t121-t20*t13
     #5/4+t67*t127*t35*t75*t131*t141)+t88*(t150*t153+t150*t157-t161+t164
     #)+t127*(t167*t169*t183/2+t61*(t170*t75*t188*t32-t192*t43*t194*t35*
     #t45-t192*t28*t194*t35*t40+t204-t207))
        t214 = t69*omstrx
        t215 = t61*t214
        t224 = amplambda*t7
        t228 = (-t150-t224*t149/4)*t88
        t230 = t194*t115
        t234 = t194*t120
        t237 = cos(t85)
        t238 = t68*t237
        t245 = t68*t88
        t253 = momega*t15
        t254 = t4*t55
        t260 = sin(t130)
        t266 = t88*t48
        t267 = t266*t260
        t271 = t127*t89
        t273 = sin(t98)
        t274 = t108*t273
        t277 = t127*t104
        t279 = t96*t273
        t282 = -t267*t30*t45*t32+2*t271*t30*t274-2*t277*t30*t279
        t297 = t271*t75
        t300 = t277*t75
        forces(2) = t4*(t17*t153+t17*t157-t215*t104*t30*t109-t215*t89*t3
     #0*t100-t161+t164)+t88*(-2*t228*t89*t230-2*t228*t104*t234-t238*omst
     #ry*t89*t230-t238*omstry*t104*t234+2*t245*t89*t121-2*t245*t104*t116
     #-t224*t147*t253*(t254-t135)/4+t67*t128*t260*t75*t141)+t127*(t167*t
     #169*t282/2+t61*(t267*t75*t40*t32+t192*t89*t194*t107*t273-t192*t104
     #*t194*t95*t273+2*t297*t109-2*t300*t100))
        t314 = t45*t32
        t352 = t11*t168*momega
        t359 = t188*t194
        forces(3) = t4*(-t8*momega*t11*t15*t183/2+t61*(-t214*t35*t30*t13
     #1*t314/2-t4*t48*t75*t173+t204-t207-2*t202*t52-2*t205*t50))+t88*(-t
     #146*t147*t253*t282/2+t61*(-t237*omstry*t48*t260*t30*t314/2-t266*t1
     #31*t75*t45*t32+2*t297*t279+2*t300*t274))+t127*(-(t146*t352+t224*t3
     #52/4)*t127*t48*t359-t68*t191*omstrz*t48*t359/2-t68*t127*t48*t172*t
     #114+t224*t11*t169*(t254-2*t101-2*t110)/4+t67*(2*t4*t82+2*t90*t75*t
     #279+2*t105*t75*t274))
        fo(1,i,j,k) = fo(1,i,j,k)+forces(1)
        fo(2,i,j,k) = fo(2,i,j,k)+forces(2)
        fo(3,i,j,k) = fo(3,i,j,k)+forces(3)
      enddo
      enddo
      enddo
      end
