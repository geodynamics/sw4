      subroutine tw_aniso_force_tt(ifirst, ilast, jfirst, jlast, kfirst, 
     *     klast, fo, t,om,cv,ph,omm,phm,amprho,phc,h, zmin) bind(c)
      use iso_c_binding
      implicit none
      integer, parameter:: dp = c_double
c
      integer, intent(in), value:: ifirst, ilast, jfirst, jlast, kfirst,
     *     klast
      real(dp), intent(in), value:: t
      real(dp), intent(in), value:: om
      real(dp), intent(in), value:: cv
      real(dp), intent(in), value:: ph
      real(dp), intent(in), value:: omm
      real(dp), intent(in), value:: phm
      real(dp), intent(in), value:: amprho
      real(dp), intent(in):: phc(21)
      real(dp), intent(out):: fo(3,ifirst:ilast,jfirst:jlast,
     *     kfirst:klast)
      real(dp), intent(in), value:: h, zmin

c local variables      
      doubleprecision x
      doubleprecision y
      doubleprecision z

      doubleprecision forces(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t100
      doubleprecision t101
      doubleprecision t103
      doubleprecision t107
      doubleprecision t108
      doubleprecision t109
      doubleprecision t110
      doubleprecision t111
      doubleprecision t112
      doubleprecision t114
      doubleprecision t115
      doubleprecision t116
      doubleprecision t117
      doubleprecision t118
      doubleprecision t119
      doubleprecision t12
      doubleprecision t120
      doubleprecision t121
      doubleprecision t122
      doubleprecision t123
      doubleprecision t125
      doubleprecision t127
      doubleprecision t130
      doubleprecision t131
      doubleprecision t134
      doubleprecision t137
      doubleprecision t139
      doubleprecision t14
      doubleprecision t140
      doubleprecision t141
      doubleprecision t142
      doubleprecision t143
      doubleprecision t144
      doubleprecision t145
      doubleprecision t146
      doubleprecision t147
      doubleprecision t149
      doubleprecision t15
      doubleprecision t152
      doubleprecision t154
      doubleprecision t155
      doubleprecision t156
      doubleprecision t157
      doubleprecision t158
      doubleprecision t159
      doubleprecision t160
      doubleprecision t162
      doubleprecision t163
      doubleprecision t166
      doubleprecision t167
      doubleprecision t169
      doubleprecision t17
      doubleprecision t170
      doubleprecision t171
      doubleprecision t172
      doubleprecision t174
      doubleprecision t176
      doubleprecision t181
      doubleprecision t183
      doubleprecision t185
      doubleprecision t186
      doubleprecision t187
      doubleprecision t188
      doubleprecision t189
      doubleprecision t19
      doubleprecision t191
      doubleprecision t193
      doubleprecision t194
      doubleprecision t197
      doubleprecision t199
      doubleprecision t2
      doubleprecision t20
      doubleprecision t200
      doubleprecision t203
      doubleprecision t205
      doubleprecision t207
      doubleprecision t208
      doubleprecision t213
      doubleprecision t214
      doubleprecision t215
      doubleprecision t216
      doubleprecision t217
      doubleprecision t22
      doubleprecision t222
      doubleprecision t224
      doubleprecision t226
      doubleprecision t227
      doubleprecision t23
      doubleprecision t230
      doubleprecision t233
      doubleprecision t234
      doubleprecision t236
      doubleprecision t239
      doubleprecision t24
      doubleprecision t244
      doubleprecision t246
      doubleprecision t248
      doubleprecision t25
      doubleprecision t250
      doubleprecision t254
      doubleprecision t255
      doubleprecision t262
      doubleprecision t263
      doubleprecision t265
      doubleprecision t266
      doubleprecision t267
      doubleprecision t269
      doubleprecision t27
      doubleprecision t273
      doubleprecision t28
      doubleprecision t283
      doubleprecision t284
      doubleprecision t286
      doubleprecision t29
      doubleprecision t291
      doubleprecision t299
      doubleprecision t3
      doubleprecision t301
      doubleprecision t303
      doubleprecision t306
      doubleprecision t307
      doubleprecision t308
      doubleprecision t31
      doubleprecision t310
      doubleprecision t314
      doubleprecision t317
      doubleprecision t32
      doubleprecision t320
      doubleprecision t323
      doubleprecision t327
      doubleprecision t33
      doubleprecision t334
      doubleprecision t34
      doubleprecision t341
      doubleprecision t347
      doubleprecision t349
      doubleprecision t35
      doubleprecision t352
      doubleprecision t355
      doubleprecision t357
      doubleprecision t358
      doubleprecision t36
      doubleprecision t366
      doubleprecision t369
      doubleprecision t37
      doubleprecision t372
      doubleprecision t373
      doubleprecision t374
      doubleprecision t38
      doubleprecision t386
      doubleprecision t388
      doubleprecision t389
      doubleprecision t392
      doubleprecision t393
      doubleprecision t394
      doubleprecision t395
      doubleprecision t396
      doubleprecision t397
      doubleprecision t4
      doubleprecision t400
      doubleprecision t404
      doubleprecision t407
      doubleprecision t408
      doubleprecision t409
      doubleprecision t41
      doubleprecision t412
      doubleprecision t415
      doubleprecision t42
      doubleprecision t420
      doubleprecision t423
      doubleprecision t425
      doubleprecision t43
      doubleprecision t431
      doubleprecision t433
      doubleprecision t438
      doubleprecision t441
      doubleprecision t444
      doubleprecision t447
      doubleprecision t450
      doubleprecision t452
      doubleprecision t46
      doubleprecision t461
      doubleprecision t464
      doubleprecision t468
      doubleprecision t47
      doubleprecision t470
      doubleprecision t472
      doubleprecision t474
      doubleprecision t478
      doubleprecision t48
      doubleprecision t481
      doubleprecision t484
      doubleprecision t486
      doubleprecision t489
      doubleprecision t49
      doubleprecision t494
      doubleprecision t5
      doubleprecision t50
      doubleprecision t502
      doubleprecision t510
      doubleprecision t512
      doubleprecision t515
      doubleprecision t516
      doubleprecision t517
      doubleprecision t518
      doubleprecision t52
      doubleprecision t521
      doubleprecision t525
      doubleprecision t53
      doubleprecision t536
      doubleprecision t538
      doubleprecision t54
      doubleprecision t543
      doubleprecision t55
      doubleprecision t554
      doubleprecision t556
      doubleprecision t562
      doubleprecision t573
      doubleprecision t58
      doubleprecision t581
      doubleprecision t587
      doubleprecision t588
      doubleprecision t592
      doubleprecision t598
      doubleprecision t6
      doubleprecision t60
      doubleprecision t600
      doubleprecision t603
      doubleprecision t608
      doubleprecision t61
      doubleprecision t619
      doubleprecision t625
      doubleprecision t628
      doubleprecision t637
      doubleprecision t638
      doubleprecision t639
      doubleprecision t640
      doubleprecision t65
      doubleprecision t652
      doubleprecision t66
      doubleprecision t663
      doubleprecision t669
      doubleprecision t67
      doubleprecision t672
      doubleprecision t684
      doubleprecision t69
      doubleprecision t70
      doubleprecision t72
      doubleprecision t73
      doubleprecision t76
      doubleprecision t77
      doubleprecision t8
      doubleprecision t80
      doubleprecision t82
      doubleprecision t83
      doubleprecision t85
      doubleprecision t86
      doubleprecision t87
      doubleprecision t88
      doubleprecision t9
      doubleprecision t90
      doubleprecision t91
      doubleprecision t92
      doubleprecision t93
      doubleprecision t94
      doubleprecision t96
      doubleprecision t97
      doubleprecision t98
      doubleprecision t99

      real(dp) ph1
      real(dp) ph2
      real(dp) ph3
      real(dp) ph4
      real(dp) ph5
      real(dp) ph6
      real(dp) ph7
      real(dp) ph8
      real(dp) ph9
      real(dp) ph10
      real(dp) ph11
      real(dp) ph12
      real(dp) ph13
      real(dp) ph14
      real(dp) ph15
      real(dp) ph16
      real(dp) ph17
      real(dp) ph18
      real(dp) ph19
      real(dp) ph20
      real(dp) ph21

      integer i, j, k
c extract all the phase angles for the stress matrix
      ph1 = phc(1)
      ph2 = phc(2)
      ph3 = phc(3)
      ph4 = phc(4)
      ph5 = phc(5)
      ph6 = phc(6)
      ph7 = phc(7)
      ph8 = phc(8)
      ph9 = phc(9)
      ph10 = phc(10)
      ph11 = phc(11)
      ph12 = phc(12)
      ph13 = phc(13)
      ph14 = phc(14)
      ph15 = phc(15)
      ph16 = phc(16)
      ph17 = phc(17)
      ph18 = phc(18)
      ph19 = phc(19)
      ph20 = phc(20)
      ph21 = phc(21)
      
      do k=kfirst,klast
        do j=jfirst,jlast
          do i=ifirst,ilast
            x=(i-1)*h
            y=(j-1)*h
            z=zmin+(k-1)*h

        t1 = omm*x
        t2 = t1+ph9
        t3 = sin(t2)
        t4 = omm*y
        t5 = t4+ph5
        t6 = cos(t5)
        t8 = omm*z
        t9 = t8+ph9
        t10 = sin(t9)
        t12 = 2+t3*t6*t10
        t14 = om*x+ph
        t15 = sin(t14)
        t17 = cv*t
        t19 = om*(y-t17)
        t20 = sin(t19)
        t22 = om**2
        t23 = t22**2
        t24 = cv**2
        t25 = t23*t24
        t27 = om*z+ph
        t28 = sin(t27)
        t29 = t25*t28
        t31 = t1+ph3
        t32 = sin(t31)
        t33 = t4+ph2
        t34 = cos(t33)
        t35 = t32*t34
        t36 = t8+ph3
        t37 = cos(t36)
        t38 = t37*omm
        t41 = om*(x-t17)
        t42 = cos(t41)
        t43 = t22*om
        t46 = om*y+ph
        t47 = sin(t46)
        t48 = t24*t47
        t49 = t48*t28
        t50 = t42*t43*t49
        t52 = t1+ph7
        t53 = sin(t52)
        t54 = t4+ph4
        t55 = cos(t54)
        t58 = sin(t8+ph7)
        t60 = 2+t53*t55*t58
        t61 = sin(t41)
        t65 = t1+ph10
        t66 = sin(t65)
        t67 = sin(t5)
        t69 = t8+ph1
        t70 = sin(t69)
        t72 = t66*t67*omm*t70
        t73 = cos(t46)
        t76 = om*(z-t17)
        t77 = sin(t76)
        t80 = t15*t73*t43*t77*t24
        t82 = t1+ph2
        t83 = cos(t82)
        t85 = t4+ph1
        t86 = cos(t85)
        t87 = t8+ph2
        t88 = sin(t87)
        t90 = t83*omm*t86*t88
        t91 = t61*t43
        t92 = t24*t73
        t93 = t92*t28
        t94 = t91*t93
        t96 = t1+ph12
        t97 = sin(t96)
        t98 = t4+ph6
        t99 = cos(t98)
        t100 = t97*t99
        t101 = sin(t36)
        t103 = 2+t100*t101
        t107 = t100*t38
        t108 = cos(t14)
        t109 = t108*t43
        t110 = t47*t77
        t111 = t110*t24
        t112 = t109*t111
        t114 = t1+ph8
        t115 = sin(t114)
        t116 = t115*t55
        t117 = t8+ph8
        t118 = cos(t117)
        t119 = t118*omm
        t120 = t116*t119
        t121 = t20*t24
        t122 = t121*t28
        t123 = t109*t122
        t125 = sin(t117)
        t127 = 2+t116*t125
        t130 = cos(t27)
        t131 = t92*t130
        t134 = cos(t31)
        t137 = t134*omm*t34*t101
        t139 = t1+ph15
        t140 = sin(t139)
        t141 = t4+ph8
        t142 = cos(t141)
        t143 = t140*t142
        t144 = t8+ph6
        t145 = cos(t144)
        t146 = t145*omm
        t147 = t143*t146
        t149 = cos(t76)
        t152 = t15*t47*t149*t43*t24
        t154 = t1+ph14
        t155 = sin(t154)
        t156 = t4+ph7
        t157 = cos(t156)
        t158 = t155*t157
        t159 = t8+ph5
        t160 = sin(t159)
        t162 = 2+t158*t160
        t163 = t162*t15
        t166 = t48*t130
        t167 = t91*t166
        t169 = -t12*t15*t20*t29+t35*t38*t50-t60*t61*t23*t49-t72*t80+t90*
     #t94-t103*t61*t23*t49+t107*t112+t120*t123+2*t127*t61*t23*t131+t137*
     #t112+t147*t152-t163*t20*t29+t137*t167
        t170 = cos(t159)
        t171 = t170*omm
        t172 = t158*t171
        t174 = t43*t24
        t176 = t15*t20*t174*t130
        t181 = cos(t19)
        t183 = t181*t24*t28
        t185 = t1+ph5
        t186 = sin(t185)
        t187 = t4+ph3
        t188 = cos(t187)
        t189 = t186*t188
        t191 = 2+t189*t160
        t193 = t191*t108*t23
        t194 = t121*t130
        t197 = sin(t54)
        t199 = omm*t125
        t200 = t115*t197*t199
        t203 = t127*t108*t23
        t205 = t73*t77*t24
        t207 = t1+ph1
        t208 = cos(t207)
        t213 = t1+ph13
        t214 = sin(t213)
        t215 = t214*t157
        t216 = t8+ph4
        t217 = cos(t216)
        t222 = t15*t181*t174*t28
        t224 = t66*t6
        t226 = 2+t224*t70
        t227 = t226*t15
        t230 = t23*t77*t24
        t233 = omm*t58
        t234 = t53*t197*t233
        t236 = sin(t82)
        t239 = 2+t236*t86*t88
        t244 = t25*t130
        t246 = t172*t176-t72*t176+t60*t108*t23*t183+t193*t194+t107*t167-
     #t200*t112+t203*t205+t208*omm*t86*t70*t50+t215*t217*omm*t222-t227*t
     #47*t230-t234*t123-t239*t15*t23*t122+t227*t181*t244
        t248 = sin(t144)
        t250 = 2+t143*t248
        t254 = t1+ph6
        t255 = cos(t254)
        t262 = omm*t10
        t263 = t3*t67*t262
        t265 = t1+ph11
        t266 = sin(t265)
        t267 = t266*t99
        t269 = 2+t267*t88
        t273 = t23*t149*t24
        t283 = t1+ph4
        t284 = cos(t283)
        t286 = sin(t216)
        t291 = 2+t35*t101
        t299 = t47*t149*t24
        t301 = sin(t85)
        t303 = omm*t88
        t306 = -t250*t15*t47*t230+t255*omm*t188*t248*t152+t90*t123-t263*
     #t222+t269*t15*t73*t273+t172*t80+t163*t73*t273+t203*t194+2*t239*t42
     #*t23*t93+t284*omm*t34*t286*t222+2*t291*t42*t23*t166+t103*t108*t23*
     #t299-t236*t301*t303*t50
        t307 = sin(t254)
        t308 = t307*t188
        t310 = 2+t308*t248
        t314 = cos(t185)
        t317 = t314*omm*t188*t160
        t320 = sin(t283)
        t323 = 2+t320*t34*t286
        t327 = sin(t207)
        t334 = sin(t98)
        t341 = 2+t215*t286
        t347 = sin(t1+phm)
        t349 = cos(t4+phm)
        t352 = sin(t8+phm)
        t355 = amprho*(2+t347*t349*t352)
        t357 = t24**2
        t358 = t23*t357
        t366 = t310*t108*t23*t299+t317*t176+t317*t80+t323*t108*t23*t183-
     #(2+t327*t86*t70)*t61*t23*t49-t266*t334*t303*t152+t193*t205-t234*t9
     #4+t341*t15*t181*t244-t200*t167+t355*t61*t358*t47*t28+t120*t94-t291
     #*t15*t23*t111
        forces(1) = t169+t246+t306+t366
        t369 = t162*t108*t23
        t372 = sin(t1+ph18)
        t373 = t4+ph9
        t374 = sin(t373)
        t386 = sin(t156)
        t388 = omm*t286
        t389 = t214*t386*t388
        t392 = sin(t1+ph19)
        t393 = t4+ph10
        t394 = cos(t393)
        t395 = t392*t394
        t396 = t8+ph10
        t397 = sin(t396)
        t400 = (2+t395*t397)*t15
        t404 = t226*t61*t23
        t407 = sin(t1+ph17)
        t408 = cos(t373)
        t409 = t407*t408
        t412 = (2+t409*t125)*t15
        t415 = cos(t265)
        t420 = cos(t52)
        t423 = t420*omm*t55*t58
        t425 = t355*t15
        t431 = t226*t108*t23
        t433 = t369*t299-t372*t374*t262*t152+2*t12*t108*t23*t183+t172*t1
     #67-t239*t61*t23*t49-t389*t167+t400*t73*t273+t404*t131-t412*t47*t23
     #0+t415*omm*t99*t88*t152+t423*t94+t425*t20*t23*t357*t28+t431*t205
        t438 = t191*t42*t23
        t441 = t341*t108*t23
        t444 = cos(t114)
        t447 = t444*omm*t55*t125
        t450 = cos(t69)
        t452 = t224*t450*omm
        t461 = t127*t42*t23
        t464 = t407*t374*t199
        t468 = t269*t108*t23
        t470 = -t60*t15*t23*t122+t438*t166+t441*t205-t389*t112+t447*t167
     #+t447*t112+t452*t123+t323*t42*t23*t93+2*t412*t181*t244+t461*t166-t
     #464*t80-t263*t94+t468*t299
        t472 = cos(t396)
        t474 = t395*t472*omm
        t478 = cos(t65)
        t481 = t478*omm*t6*t70
        t484 = t162*t61*t23
        t486 = t372*t408
        t489 = (2+t486*t10)*t15
        t494 = cos(t2)
        t502 = sin(t33)
        t510 = t474*t80-t400*t20*t29+t481*t80-t484*t49+t489*t73*t273-t26
     #3*t123+t474*t176+t494*omm*t6*t10*t222+t172*t112+t189*t171*t50-t320
     #*t502*t388*t50+t341*t61*t23*t131+t423*t123
        t512 = t127*t15*t23
        t515 = sin(t1+ph20)
        t516 = t515*t394
        t517 = t8+ph11
        t518 = sin(t517)
        t521 = (2+t516*t518)*t15
        t525 = sin(t1+ph16)
        t536 = cos(t517)
        t538 = t516*t536*omm
        t543 = sin(t141)
        t554 = -t512*t111-t521*t47*t230-(2+t525*t142*t58)*t15*t20*t29+t4
     #81*t176-t12*t61*t23*t49+t538*t152+2*t431*t194+t452*t94-t525*t543*t
     #233*t222+t90*t50+t60*t42*t23*t93-t464*t176+t409*t119*t222
        forces(2) = t433+t470+t510+t554
        t556 = cos(t213)
        t562 = sin(t393)
        t573 = t392*t562*omm*t397
        t581 = cos(t139)
        t587 = omm*t160
        t588 = t155*t386*t587
        t592 = t556*omm*t157*t286*t222+t484*t131-t515*t562*omm*t518*t152
     #+t461*t93+t425*t110*t358-t512*t122-t573*t80-t72*t123+2*t250*t108*t
     #23*t299+t468*t194+t581*omm*t142*t248*t152-t588*t167+t308*t146*t50
        t598 = cos(t87)
        t600 = t267*t598*omm
        t603 = sin(t187)
        t608 = cos(t9)
        t619 = t447*t94+t538*t176-t291*t61*t23*t49+t600*t94+t438*t93-t18
     #6*t603*t587*t50+t431*t183+t486*t608*omm*t222+t369*t194-t521*t20*t2
     #9+t538*t80+t600*t123+2*t369*t205
        t625 = cos(t96)
        t628 = t625*omm*t99*t101
        t637 = sin(t1+ph21)
        t638 = t637*t394
        t639 = t8+ph12
        t640 = cos(t639)
        t652 = -t573*t176+t441*t183+t489*t181*t244+t628*t112-t250*t61*t2
     #3*t49-t404*t49+t400*t181*t244+t638*t640*omm*t152-t72*t94+t447*t123
     #-t400*t47*t230-t464*t222+t269*t61*t23*t131
        t663 = sin(t639)
        t669 = cos(t154)
        t672 = t669*omm*t157*t160
        t684 = t103*t42*t23*t166+t147*t112+2*t521*t73*t273+t310*t42*t23*
     #t166-(2+t638*t663)*t15*t47*t230+t672*t176+t672*t80+t628*t167-t103*
     #t15*t23*t111+t147*t167-t412*t20*t29-t588*t112+t137*t50
        forces(3) = t592+t619+t652+t684

      fo(1,i,j,k) = forces(1)
      fo(2,i,j,k) = forces(2)
      fo(3,i,j,k) = forces(3)
      enddo
      enddo
      enddo

      return
      end
