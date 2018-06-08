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
      subroutine fg(x,y,z,t,om,c,ph,omm,phm,amprho,ampmu,amplambda,crea_
     #par)
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
        crea_par(1) = forces(1)
        crea_par(2) = forces(2)
        crea_par(3) = forces(3)
        return
        return
      end
