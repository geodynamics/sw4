c-----------------------------------------------------------------------
      subroutine twilightfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, u, t, om, cv, ph, h, zmin )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
* v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
* w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, cv, ph, zmin
      doubleprecision ampmu, amplambda, h
      real*8 u(3,ifirst:ilast,jfirst:jlast,kfirst:klast)
      do k=kfirst,klast
        z = (k-1)*h + zmin
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            u(1,i,j,k) = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph)
            u(2,i,j,k) = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph)
            u(3,i,j,k) = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t))
          enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for right hand side forcing
c corresponding to twilight exact solution
c
      subroutine exactrhsfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda,
     +     h, zmin)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho, zmin
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t100
      doubleprecision t101
      doubleprecision t102
      doubleprecision t103
      doubleprecision t105
      doubleprecision t106
      doubleprecision t108
      doubleprecision t11
      doubleprecision t114
      doubleprecision t115
      doubleprecision t118
      doubleprecision t119
      doubleprecision t12
      doubleprecision t120
      doubleprecision t122
      doubleprecision t125
      doubleprecision t129
      doubleprecision t135
      doubleprecision t143
      doubleprecision t15
      doubleprecision t150
      doubleprecision t152
      doubleprecision t159
      doubleprecision t16
      doubleprecision t166
      doubleprecision t168
      doubleprecision t17
      doubleprecision t173
      doubleprecision t175
      doubleprecision t18
      doubleprecision t2
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t38
      doubleprecision t4
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t47
      doubleprecision t48
      doubleprecision t51
      doubleprecision t53
      doubleprecision t54
      doubleprecision t55
      doubleprecision t57
      doubleprecision t58
      doubleprecision t6
      doubleprecision t60
      doubleprecision t63
      doubleprecision t64
      doubleprecision t69
      doubleprecision t7
      doubleprecision t70
      doubleprecision t71
      doubleprecision t72
      doubleprecision t77
      doubleprecision t78
      doubleprecision t79
      doubleprecision t8
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t83
      doubleprecision t85
      doubleprecision t87
      doubleprecision t89
      doubleprecision t92
      doubleprecision t95
      doubleprecision t97

      do k=kfirst,klast
         z = (k-1)*h + zmin
         do j=jfirst,jlast
           y = (j-1)*h
           do i=ifirst,ilast
             x = (i-1)*h

        t2 = omm*x+phm
        t3 = sin(t2)
        t4 = ampmu*t3
        t6 = omm*y+phm
        t7 = sin(t6)
        t8 = omm*t7
        t10 = omm*z+phm
        t11 = sin(t10)
        t12 = t8*t11
        t15 = cos(t2)
        t16 = amplambda*t15
        t17 = cos(t10)
        t18 = t8*t17
        t21 = c*t
        t23 = om*(x-t21)
        t24 = cos(t23)
        t27 = om*y+ph
        t28 = sin(t27)
        t31 = om*z+ph
        t32 = sin(t31)
        t38 = ampmu*(3+t15*t7*t11)
        t43 = amplambda*(2+t3*t7*t17)
        t44 = 2*t38+t43
        t45 = sin(t23)
        t47 = om**2
        t48 = t47*t28
        t51 = t16*t8
        t53 = om*x+ph
        t54 = sin(t53)
        t55 = t17*t54
        t57 = om*(y-t21)
        t58 = cos(t57)
        t60 = t58*om*t32
        t63 = cos(t53)
        t64 = t43*t63
        t69 = om*(z-t21)
        t70 = cos(t69)
        t71 = t28*t70
        t72 = t71*om
        t77 = ampmu*t15
        t78 = cos(t6)
        t79 = t77*t78
        t80 = omm*t11
        t81 = t63*om
        t82 = sin(t57)
        t83 = t82*t32
        t85 = cos(t27)
        t87 = om*t32
        t89 = t81*t83+t45*t85*t87
        t92 = t63*t47
        t95 = t45*t28
        t97 = t95*t47*t32
        t100 = t77*omm
        t101 = t7*t17
        t102 = sin(t69)
        t103 = t28*t102
        t105 = cos(t31)
        t106 = t105*om
        t108 = t81*t103+t95*t106
        forces(1) = (-2*t4*t12+t16*t18)*t24*om*t28*t32-t44*t45*t48*t32+t
     #51*t55*t60+t64*t47*t58*t32+t51*t55*t72+t64*t48*t70+t79*t80*t89+t38
     #*(t92*t58*t32-t97)+t100*t101*t108+t38*(t92*t71-t97)
        t114 = t4*omm
        t115 = t7*t11
        t118 = t54*t47
        t119 = t118*t83
        t120 = t24*t47
        t122 = t120*t85*t32
        t125 = t78*omm
        t129 = amplambda*t3
        t135 = t44*t54
        t143 = t24*om*t28*t32
        t150 = t54*t85
        t152 = t150*t47*t70
        t159 = t150*om*t102+t54*t82*t106
        forces(2) = -t114*t115*t89+t38*(-t119+t122)+(2*t77*t125*t11+t129
     #*t125*t17)*t54*t60-t135*t82*t47*t32+t129*t78*omm*t17*(t143+t54*t28
     #*t70*om)+t43*(t122+t152)+t100*t101*t159+t38*(t152-t119)
        t166 = t118*t103
        t168 = t120*t28*t105
        t173 = t54*t58
        t175 = t173*t47*t105
        forces(3) = -t114*t115*t108+t38*(-t166+t168)+t79*t80*t159+t38*(-
     #t166+t175)+(2*t77*t18-t129*t12)*t54*t72-t135*t103*t47-t129*omm*t11
     #5*(t143+t173*t87)+t43*(t168+t175)

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for right hand side forcing
c corresponding to twilight exact solution
c
      subroutine exactaccfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, utt, t, om, c, ph, h, zmin )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
* v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
* w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, zmin
      doubleprecision ampmu, amplambda, h
      real*8 utt(3,ifirst:ilast,jfirst:jlast,kfirst:klast)

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

      do k=kfirst,klast
        z = (k-1)*h + zmin
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
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

        utt(1,i,j,k) = acc(1)
        utt(2,i,j,k) = acc(2)
        utt(3,i,j,k) = acc(3)
        enddo
        enddo
      enddo
      return
      end

c-----------------------------------------------------------------------
c
c Maple auto-generated code for twilight forcing functions for the elastic wave
c equation, corresponding to exactSol and exactMat
c
      subroutine forcingfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     h, zmin)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho, zmin
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t102
      doubleprecision t105
      doubleprecision t107
      doubleprecision t110
      doubleprecision t111
      doubleprecision t112
      doubleprecision t113
      doubleprecision t115
      doubleprecision t116
      doubleprecision t118
      doubleprecision t124
      doubleprecision t125
      doubleprecision t129
      doubleprecision t13
      doubleprecision t130
      doubleprecision t133
      doubleprecision t134
      doubleprecision t135
      doubleprecision t137
      doubleprecision t14
      doubleprecision t140
      doubleprecision t144
      doubleprecision t150
      doubleprecision t156
      doubleprecision t16
      doubleprecision t163
      doubleprecision t165
      doubleprecision t17
      doubleprecision t172
      doubleprecision t181
      doubleprecision t183
      doubleprecision t188
      doubleprecision t19
      doubleprecision t190
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t26
      doubleprecision t27
      doubleprecision t28
      doubleprecision t3
      doubleprecision t31
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t39
      doubleprecision t40
      doubleprecision t43
      doubleprecision t5
      doubleprecision t51
      doubleprecision t56
      doubleprecision t57
      doubleprecision t59
      doubleprecision t6
      doubleprecision t62
      doubleprecision t64
      doubleprecision t65
      doubleprecision t66
      doubleprecision t68
      doubleprecision t69
      doubleprecision t71
      doubleprecision t74
      doubleprecision t75
      doubleprecision t80
      doubleprecision t81
      doubleprecision t82
      doubleprecision t83
      doubleprecision t88
      doubleprecision t89
      doubleprecision t9
      doubleprecision t90
      doubleprecision t91
      doubleprecision t92
      doubleprecision t93
      doubleprecision t95
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
        t31 = ampmu*t3
        t32 = sin(t5)
        t33 = omm*t32
        t34 = t33*t10
        t37 = cos(t2)
        t38 = amplambda*t37
        t39 = cos(t9)
        t40 = t33*t39
        t43 = cos(t16)
        t51 = ampmu*(3+t37*t32*t10)
        t56 = amplambda*(2+t3*t32*t39)
        t57 = 2*t51+t56
        t59 = t19*t24
        t62 = t38*t33
        t64 = om*x+ph
        t65 = sin(t64)
        t66 = t39*t65
        t68 = om*(y-t14)
        t69 = cos(t68)
        t71 = t69*om*t27
        t74 = cos(t64)
        t75 = t56*t74
        t80 = om*(z-t14)
        t81 = cos(t80)
        t82 = t24*t81
        t83 = t82*om
        t88 = ampmu*t37
        t89 = t88*t6
        t90 = omm*t10
        t91 = t74*om
        t92 = sin(t68)
        t93 = t92*t27
        t95 = cos(t23)
        t97 = om*t27
        t99 = t91*t93+t17*t95*t97
        t102 = t74*t19
        t105 = t17*t24
        t107 = t105*t19*t27
        t110 = t88*omm
        t111 = t32*t39
        t112 = sin(t80)
        t113 = t24*t112
        t115 = cos(t26)
        t116 = t115*om
        t118 = t91*t113+t105*t116
        forces(1) = -t13*t17*t21*t28-(-2*t31*t34+t38*t40)*t43*om*t24*t27
     #+t57*t17*t59*t27-t62*t66*t71-t75*t19*t69*t27-t62*t66*t83-t75*t59*t
     #81-t89*t90*t99-t51*(t102*t69*t27-t107)-t110*t111*t118-t51*(t102*t8
     #2-t107)
        t124 = t13*t65
        t125 = t92*t19
        t129 = t31*omm
        t130 = t32*t10
        t133 = t65*t19
        t134 = t133*t93
        t135 = t43*t19
        t137 = t135*t95*t27
        t140 = t6*omm
        t144 = amplambda*t3
        t150 = t57*t65
        t156 = t43*om*t28
        t163 = t65*t95
        t165 = t163*t19*t81
        t172 = t163*om*t112+t65*t92*t116
        forces(2) = -t124*t125*t20*t27+t129*t130*t99-t51*(-t134+t137)-(2
     #*t88*t140*t10+t144*t140*t39)*t65*t71+t150*t125*t27-t144*t6*omm*t39
     #*(t156+t65*t24*t81*om)-t56*(t137+t165)-t110*t111*t172-t51*(t165-t1
     #34)
        t181 = t133*t113
        t183 = t135*t24*t115
        t188 = t65*t69
        t190 = t188*t19*t115
        forces(3) = -t124*t113*t21+t129*t130*t118-t51*(-t181+t183)-t89*t
     #90*t172-t51*(-t181+t190)-(2*t88*t40-t144*t34)*t65*t83+t150*t113*t1
     #9+t144*omm*t130*(t156+t188*t97)-t56*(t183+t190)

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine forcingttfort( ifirst, ilast, jfirst, jlast, kfirst,
     +     klast, fo, t, om, c, ph, omm, phm, amprho, ampmu, amplambda, 
     +     h, zmin)
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho, zmin
      doubleprecision ampmu, amplambda
      real*8 fo(3,ifirst:ilast,jfirst:jlast,kfirst:klast), h

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

        fo(1,i,j,k) = forces(1)
        fo(2,i,j,k) = forces(2)
        fo(3,i,j,k) = forces(3)

      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine exactmatfort( ifirst, ilast, jfirst, jlast, kfirst, 
     +     klast, rho, mu, la, omm, phm, amprho, ampmu, amplambda, h, 
     +     zmin )
* new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
* rho    := amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
* mu     := ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*sin(omm*z+phm) );
* lambda := amplambda*(2 + sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) );
      implicit none
      integer ifirst, ilast, jfirst, jlast, kfirst, klast, i, j, k
      doubleprecision x, y, z, t, om, c, ph, omm, phm, amprho
      doubleprecision ampmu, amplambda, h, zmin
      real*8 rho(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 mu(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real*8 la(ifirst:ilast,jfirst:jlast,kfirst:klast)

      do k=kfirst,klast
        z = (k-1)*h + zmin
        do j=jfirst,jlast
          y = (j-1)*h
          do i=ifirst,ilast
            x = (i-1)*h
            rho(i,j,k) = amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*
     +           sin(omm*z+phm) );
            mu(i,j,k)  = ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*
     +           sin(omm*z+phm) )
            la(i,j,k)  = amplambda*(2 + 
     +           sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) )
          enddo
        enddo
      enddo
      return
      end
