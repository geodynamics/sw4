      subroutine tw_aniso_force(ifirst, ilast, jfirst, jlast, kfirst, 
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
c      real(dp) crea_par(3)

c local variables
      real(dp) x
      real(dp) y
      real(dp) z

      real(dp) t1
      real(dp) t10
      real(dp) t102
      real(dp) t104
      real(dp) t106
      real(dp) t109
      real(dp) t111
      real(dp) t112
      real(dp) t114
      real(dp) t115
      real(dp) t117
      real(dp) t118
      real(dp) t119
      real(dp) t120
      real(dp) t122
      real(dp) t124
      real(dp) t126
      real(dp) t127
      real(dp) t129
      real(dp) t13
      real(dp) t130
      real(dp) t131
      real(dp) t134
      real(dp) t135
      real(dp) t136
      real(dp) t137
      real(dp) t139
      real(dp) t141
      real(dp) t145
      real(dp) t147
      real(dp) t149
      real(dp) t15
      real(dp) t150
      real(dp) t152
      real(dp) t153
      real(dp) t157
      real(dp) t159
      real(dp) t16
      real(dp) t161
      real(dp) t162
      real(dp) t164
      real(dp) t166
      real(dp) t168
      real(dp) t170
      real(dp) t173
      real(dp) t174
      real(dp) t175
      real(dp) t176
      real(dp) t179
      real(dp) t180
      real(dp) t181
      real(dp) t183
      real(dp) t186
      real(dp) t189
      real(dp) t19
      real(dp) t190
      real(dp) t192
      real(dp) t193
      real(dp) t194
      real(dp) t195
      real(dp) t197
      real(dp) t199
      real(dp) t2
      real(dp) t20
      real(dp) t201
      real(dp) t207
      real(dp) t208
      real(dp) t209
      real(dp) t210
      real(dp) t212
      real(dp) t213
      real(dp) t214
      real(dp) t215
      real(dp) t217
      real(dp) t219
      real(dp) t22
      real(dp) t222
      real(dp) t225
      real(dp) t226
      real(dp) t229
      real(dp) t23
      real(dp) t231
      real(dp) t233
      real(dp) t234
      real(dp) t235
      real(dp) t238
      real(dp) t24
      real(dp) t242
      real(dp) t243
      real(dp) t244
      real(dp) t245
      real(dp) t249
      real(dp) t25
      real(dp) t250
      real(dp) t252
      real(dp) t254
      real(dp) t256
      real(dp) t257
      real(dp) t260
      real(dp) t261
      real(dp) t262
      real(dp) t264
      real(dp) t265
      real(dp) t266
      real(dp) t267
      real(dp) t27
      real(dp) t270
      real(dp) t273
      real(dp) t276
      real(dp) t277
      real(dp) t278
      real(dp) t279
      real(dp) t280
      real(dp) t281
      real(dp) t286
      real(dp) t289
      real(dp) t290
      real(dp) t291
      real(dp) t292
      real(dp) t293
      real(dp) t294
      real(dp) t297
      real(dp) t298
      real(dp) t3
      real(dp) t305
      real(dp) t306
      real(dp) t307
      real(dp) t308
      real(dp) t309
      real(dp) t31
      real(dp) t310
      real(dp) t311
      real(dp) t312
      real(dp) t315
      real(dp) t319
      real(dp) t321
      real(dp) t324
      real(dp) t327
      real(dp) t329
      real(dp) t33
      real(dp) t330
      real(dp) t333
      real(dp) t337
      real(dp) t34
      real(dp) t340
      real(dp) t344
      real(dp) t349
      real(dp) t35
      real(dp) t352
      real(dp) t355
      real(dp) t363
      real(dp) t366
      real(dp) t368
      real(dp) t37
      real(dp) t372
      real(dp) t373
      real(dp) t377
      real(dp) t38
      real(dp) t382
      real(dp) t384
      real(dp) t386
      real(dp) t394
      real(dp) t396
      real(dp) t40
      real(dp) t402
      real(dp) t403
      real(dp) t407
      real(dp) t41
      real(dp) t415
      real(dp) t416
      real(dp) t417
      real(dp) t419
      real(dp) t421
      real(dp) t422
      real(dp) t425
      real(dp) t429
      real(dp) t43
      real(dp) t434
      real(dp) t438
      real(dp) t44
      real(dp) t441
      real(dp) t445
      real(dp) t447
      real(dp) t449
      real(dp) t451
      real(dp) t453
      real(dp) t455
      real(dp) t46
      real(dp) t461
      real(dp) t462
      real(dp) t463
      real(dp) t464
      real(dp) t465
      real(dp) t466
      real(dp) t468
      real(dp) t47
      real(dp) t470
      real(dp) t473
      real(dp) t476
      real(dp) t481
      real(dp) t482
      real(dp) t483
      real(dp) t484
      real(dp) t486
      real(dp) t488
      real(dp) t49
      real(dp) t491
      real(dp) t493
      real(dp) t497
      real(dp) t499
      real(dp) t5
      real(dp) t512
      real(dp) t515
      real(dp) t519
      real(dp) t52
      real(dp) t522
      real(dp) t526
      real(dp) t532
      real(dp) t54
      real(dp) t540
      real(dp) t542
      real(dp) t549
      real(dp) t55
      real(dp) t554
      real(dp) t558
      real(dp) t561
      real(dp) t577
      real(dp) t579
      real(dp) t58
      real(dp) t583
      real(dp) t589
      real(dp) t59
      real(dp) t599
      real(dp) t6
      real(dp) t600
      real(dp) t601
      real(dp) t602
      real(dp) t606
      real(dp) t61
      real(dp) t613
      real(dp) t62
      real(dp) t63
      real(dp) t64
      real(dp) t66
      real(dp) t68
      real(dp) t69
      real(dp) t7
      real(dp) t70
      real(dp) t72
      real(dp) t73
      real(dp) t75
      real(dp) t77
      real(dp) t8
      real(dp) t81
      real(dp) t82
      real(dp) t83
      real(dp) t85
      real(dp) t86
      real(dp) t88
      real(dp) t9
      real(dp) t90
      real(dp) t92
      real(dp) t93
      real(dp) t95
      real(dp) t96
      real(dp) t98
      real(dp) t99
      real(dp) threecomp(3)

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
        t2 = t1+ph1
        t3 = cos(t2)
        t5 = omm*y
        t6 = t5+ph1
        t7 = cos(t6)
        t8 = omm*z
        t9 = t8+ph1
        t10 = sin(t9)
        t13 = cv*t
        t15 = om*(x-t13)
        t16 = cos(t15)
        t19 = om*y+ph
        t20 = sin(t19)
        t22 = om*z+ph
        t23 = sin(t22)
        t24 = t20*t23
        t25 = t16*om*t24
        t27 = sin(t2)
        t31 = sin(t15)
        t33 = om**2
        t34 = t33*t20
        t35 = t34*t23
        t37 = t1+ph2
        t38 = cos(t37)
        t40 = t8+ph2
        t41 = sin(t40)
        t43 = t38*omm*t7*t41
        t44 = cos(t19)
        t46 = om*t23
        t47 = t31*t44*t46
        t49 = sin(t37)
        t52 = 2+t49*t7*t41
        t54 = t33*t44
        t55 = t54*t23
        t58 = t1+ph3
        t59 = cos(t58)
        t61 = t5+ph2
        t62 = cos(t61)
        t63 = t8+ph3
        t64 = sin(t63)
        t66 = t59*omm*t62*t64
        t68 = cos(t22)
        t69 = t68*om
        t70 = t31*t20*t69
        t72 = sin(t58)
        t73 = t72*t62
        t75 = 2+t73*t64
        t77 = t34*t68
        t81 = om*x+ph
        t82 = cos(t81)
        t83 = t82*om
        t85 = om*(y-t13)
        t86 = sin(t85)
        t88 = t83*t86*t23
        t90 = sin(t81)
        t92 = t33*t86
        t93 = t92*t23
        t95 = t1+ph4
        t96 = cos(t95)
        t98 = t8+ph4
        t99 = sin(t98)
        t102 = cos(t85)
        t104 = t90*t102*t46
        t106 = sin(t95)
        t109 = 2+t106*t62*t99
        t111 = t33*t102
        t112 = t111*t23
        t114 = t1+ph5
        t115 = cos(t114)
        t117 = t5+ph3
        t118 = cos(t117)
        t119 = t8+ph5
        t120 = sin(t119)
        t122 = t115*omm*t118*t120
        t124 = t90*t86*t69
        t126 = sin(t114)
        t127 = t126*t118
        t129 = 2+t127*t120
        t130 = t129*t82
        t131 = t92*t68
        t134 = om*(z-t13)
        t135 = sin(t134)
        t136 = t20*t135
        t137 = t83*t136
        t139 = -t3*omm*t7*t10*t25+(10+t27*t7*t10)*t31*t35-t43*t47-2*t52*
     #t16*t55-t66*t70-2*t75*t16*t77-t43*t88+t52*t90*t93-t96*omm*t62*t99*
     #t104-t109*t82*t112-t122*t124-t130*t131-t66*t137
        t141 = t34*t135
        t145 = t90*t44*om*t135
        t147 = t54*t135
        t149 = t1+ph6
        t150 = cos(t149)
        t152 = t8+ph6
        t153 = sin(t152)
        t157 = cos(t134)
        t159 = t90*t20*t157*om
        t161 = sin(t149)
        t162 = t161*t118
        t164 = 2+t162*t153
        t166 = t34*t157
        t168 = sin(t6)
        t170 = omm*t41
        t173 = t1+ph7
        t174 = sin(t173)
        t175 = t5+ph4
        t176 = sin(t175)
        t179 = sin(t8+ph7)
        t180 = omm*t179
        t181 = t174*t176*t180
        t183 = cos(t175)
        t186 = 10+t174*t183*t179
        t189 = t1+ph8
        t190 = sin(t189)
        t192 = t8+ph8
        t193 = sin(t192)
        t194 = omm*t193
        t195 = t190*t176*t194
        t197 = t190*t183
        t199 = 2+t197*t193
        t201 = t54*t68
        t207 = t1+ph9
        t208 = sin(t207)
        t209 = t5+ph5
        t210 = sin(t209)
        t212 = t8+ph9
        t213 = sin(t212)
        t214 = omm*t213
        t215 = t208*t210*t214
        t217 = t75*t90*t141-t122*t145-t130*t147-t150*omm*t118*t153*t159-
     #t164*t82*t166+t49*t168*t170*t25+t181*t47+t186*t31*t35+t195*t70-2*t
     #199*t31*t201+t181*t88-t186*t82*t112+t215*t104
        t219 = cos(t209)
        t222 = 2+t208*t219*t213
        t225 = t1+ph10
        t226 = sin(t225)
        t229 = t226*t210*omm*t10
        t231 = t226*t219
        t233 = 2+t231*t10
        t234 = t233*t90
        t235 = t111*t68
        t238 = t199*t82
        t242 = t1+ph11
        t243 = sin(t242)
        t244 = t5+ph6
        t245 = sin(t244)
        t249 = cos(t244)
        t250 = t243*t249
        t252 = 2+t250*t41
        t254 = t54*t157
        t256 = cos(t63)
        t257 = t256*omm
        t260 = cos(t192)
        t261 = t260*omm
        t262 = t197*t261
        t264 = t1+ph12
        t265 = sin(t264)
        t266 = t265*t249
        t267 = t266*t257
        t270 = 10+t266*t64
        t273 = t222*t90*t93+t229*t124-t234*t235+t195*t137-t238*t147+t229
     #*t145+t234*t141+t243*t245*t170*t159-t252*t90*t254-t73*t257*t25-t26
     #2*t47-t267*t70+t270*t31*t35
        t276 = t1+ph13
        t277 = sin(t276)
        t278 = t5+ph7
        t279 = cos(t278)
        t280 = t277*t279
        t281 = cos(t98)
        t286 = 2+t280*t99
        t289 = t1+ph14
        t290 = sin(t289)
        t291 = t290*t279
        t292 = cos(t119)
        t293 = t292*omm
        t294 = t291*t293
        t297 = 2+t291*t120
        t298 = t297*t90
        t305 = t1+ph15
        t306 = sin(t305)
        t307 = t5+ph8
        t308 = cos(t307)
        t309 = t306*t308
        t310 = cos(t152)
        t311 = t310*omm
        t312 = t309*t311
        t315 = 2+t309*t153
        t319 = sin(t1+phm)
        t321 = cos(t5+phm)
        t324 = sin(t8+phm)
        t327 = amprho*(2+t319*t321*t324)
        t329 = cv**2
        t330 = t33*t329
        t333 = -t262*t88-t238*t131-t280*t281*omm*t104-t286*t90*t235-t294
     #*t124+t298*t93-t267*t137-t270*t82*t166-t294*t145-t298*t254-t312*t1
     #59+t315*t90*t141-t327*t31*t330*t24
        threecomp(1) = t139+t217+t273+t333
        t337 = cos(t173)
        t340 = t337*omm*t183*t179
        t344 = t199*t16
        t349 = cos(t189)
        t352 = t349*omm*t183*t193
        t355 = cos(t207)
        t363 = cos(t225)
        t366 = t363*omm*t219*t10
        t368 = t233*t82
        t372 = t186*t90*t93-t340*t88+t52*t31*t35-t344*t77-t43*t25-t186*t
     #16*t55-t352*t70-t340*t47-t355*omm*t219*t213*t104-2*t222*t82*t112-t
     #366*t124-2*t368*t131-t352*t137
        t373 = t199*t90
        t377 = cos(t242)
        t382 = t252*t82
        t384 = sin(t61)
        t386 = omm*t99
        t394 = sin(t278)
        t396 = t277*t394*t386
        t402 = sin(t1+ph16)
        t403 = sin(t307)
        t407 = t373*t141-t366*t145-t368*t147-t377*omm*t249*t41*t159-t382
     #*t166+t106*t384*t386*t25-t109*t16*t55+t215*t47+t222*t31*t35+t396*t
     #70-t286*t31*t201+t215*t88+t402*t403*t180*t104
        t415 = sin(t1+ph17)
        t416 = t5+ph9
        t417 = sin(t416)
        t419 = t415*t417*t194
        t421 = cos(t416)
        t422 = t415*t421
        t425 = (2+t422*t193)*t90
        t429 = t286*t82
        t434 = sin(t1+ph18)
        t438 = t434*t421
        t441 = (2+t438*t213)*t90
        t445 = t129*t16
        t447 = cos(t9)
        t449 = t231*t447*omm
        t451 = t233*t31
        t453 = (10+t402*t308*t179)*t90*t93+t419*t124-2*t425*t235+t396*t1
     #37-t429*t147+t419*t145+t425*t141+t434*t417*t214*t159-t441*t254-t12
     #7*t293*t25-t445*t77-t449*t47-t451*t201
        t455 = t297*t31
        t461 = sin(t1+ph19)
        t462 = t5+ph10
        t463 = cos(t462)
        t464 = t461*t463
        t465 = t8+ph10
        t466 = cos(t465)
        t468 = t464*t466*omm
        t470 = sin(t465)
        t473 = (10+t464*t470)*t90
        t476 = t297*t82
        t481 = sin(t1+ph20)
        t482 = t481*t463
        t483 = t8+ph11
        t484 = cos(t483)
        t486 = t482*t484*omm
        t488 = sin(t483)
        t491 = (2+t482*t488)*t90
        t493 = t327*t90
        t497 = -t294*t70+t455*t35-t449*t88-t422*t261*t104-t468*t124+t473
     #*t93-t294*t137-t476*t166-t468*t145-t473*t254-t486*t159+t491*t141-t
     #493*t92*t329*t23
        threecomp(2) = t372+t407+t453+t497
        t499 = cos(t276)
        t512 = cos(t264)
        t515 = t512*omm*t249*t64
        t519 = cos(t289)
        t522 = t519*omm*t279*t120
        t526 = -t499*omm*t279*t99*t104-t429*t112+t75*t31*t35-t352*t88-t3
     #52*t47-t66*t25-t270*t16*t77-t515*t70-t344*t55+t373*t93-t522*t124-t
     #476*t131-t515*t137
        t532 = cos(t305)
        t540 = sin(t117)
        t542 = omm*t120
        t549 = t290*t394*t542
        t554 = t270*t90*t141-t522*t145-2*t476*t147-t532*omm*t308*t153*t1
     #59-2*t315*t82*t166+t126*t540*t542*t25-t445*t55+t229*t47+t451*t35+t
     #549*t70-t455*t201+t229*t88-t368*t112
        t558 = sin(t462)
        t561 = t461*t558*omm*t470
        t577 = cos(t40)
        t579 = t250*t577*omm
        t583 = t419*t104+t425*t93+t561*t124-t473*t235+t549*t137+t561*t14
     #5+t473*t141+t481*t558*omm*t488*t159-2*t491*t254-t162*t311*t25-t164
     #*t16*t77-t579*t47-t252*t31*t201
        t589 = cos(t212)
        t599 = sin(t1+ph21)
        t600 = t599*t463
        t601 = t8+ph12
        t602 = cos(t601)
        t606 = sin(t601)
        t613 = -t312*t70+t315*t31*t35-t579*t88-t382*t131-t438*t589*omm*t
     #104-t441*t235-t486*t124+t491*t93-t312*t137-t486*t145-t600*t602*omm
     #*t159+(10+t600*t606)*t90*t141-t493*t136*t330
        threecomp(3) = t526+t554+t583+t613

      fo(1,i,j,k) = threecomp(1)
      fo(2,i,j,k) = threecomp(2)
      fo(3,i,j,k) = threecomp(3)
      enddo
      enddo
      enddo

      return
      end


c curvilinear case
      subroutine tw_aniso_curvi_force(ifirst, ilast, jfirst, jlast, 
     *     kfirst, klast, fo, t,om,cv,ph,omm,phm,amprho,phc,
     *     xx, yy, zz) bind(c)
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
      real(dp), intent(in):: xx(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real(dp), intent(in):: yy(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real(dp), intent(in):: zz(ifirst:ilast,jfirst:jlast,kfirst:klast)

c local variables
      real(dp) x
      real(dp) y
      real(dp) z

      real(dp) t1
      real(dp) t10
      real(dp) t102
      real(dp) t104
      real(dp) t106
      real(dp) t109
      real(dp) t111
      real(dp) t112
      real(dp) t114
      real(dp) t115
      real(dp) t117
      real(dp) t118
      real(dp) t119
      real(dp) t120
      real(dp) t122
      real(dp) t124
      real(dp) t126
      real(dp) t127
      real(dp) t129
      real(dp) t13
      real(dp) t130
      real(dp) t131
      real(dp) t134
      real(dp) t135
      real(dp) t136
      real(dp) t137
      real(dp) t139
      real(dp) t141
      real(dp) t145
      real(dp) t147
      real(dp) t149
      real(dp) t15
      real(dp) t150
      real(dp) t152
      real(dp) t153
      real(dp) t157
      real(dp) t159
      real(dp) t16
      real(dp) t161
      real(dp) t162
      real(dp) t164
      real(dp) t166
      real(dp) t168
      real(dp) t170
      real(dp) t173
      real(dp) t174
      real(dp) t175
      real(dp) t176
      real(dp) t179
      real(dp) t180
      real(dp) t181
      real(dp) t183
      real(dp) t186
      real(dp) t189
      real(dp) t19
      real(dp) t190
      real(dp) t192
      real(dp) t193
      real(dp) t194
      real(dp) t195
      real(dp) t197
      real(dp) t199
      real(dp) t2
      real(dp) t20
      real(dp) t201
      real(dp) t207
      real(dp) t208
      real(dp) t209
      real(dp) t210
      real(dp) t212
      real(dp) t213
      real(dp) t214
      real(dp) t215
      real(dp) t217
      real(dp) t219
      real(dp) t22
      real(dp) t222
      real(dp) t225
      real(dp) t226
      real(dp) t229
      real(dp) t23
      real(dp) t231
      real(dp) t233
      real(dp) t234
      real(dp) t235
      real(dp) t238
      real(dp) t24
      real(dp) t242
      real(dp) t243
      real(dp) t244
      real(dp) t245
      real(dp) t249
      real(dp) t25
      real(dp) t250
      real(dp) t252
      real(dp) t254
      real(dp) t256
      real(dp) t257
      real(dp) t260
      real(dp) t261
      real(dp) t262
      real(dp) t264
      real(dp) t265
      real(dp) t266
      real(dp) t267
      real(dp) t27
      real(dp) t270
      real(dp) t273
      real(dp) t276
      real(dp) t277
      real(dp) t278
      real(dp) t279
      real(dp) t280
      real(dp) t281
      real(dp) t286
      real(dp) t289
      real(dp) t290
      real(dp) t291
      real(dp) t292
      real(dp) t293
      real(dp) t294
      real(dp) t297
      real(dp) t298
      real(dp) t3
      real(dp) t305
      real(dp) t306
      real(dp) t307
      real(dp) t308
      real(dp) t309
      real(dp) t31
      real(dp) t310
      real(dp) t311
      real(dp) t312
      real(dp) t315
      real(dp) t319
      real(dp) t321
      real(dp) t324
      real(dp) t327
      real(dp) t329
      real(dp) t33
      real(dp) t330
      real(dp) t333
      real(dp) t337
      real(dp) t34
      real(dp) t340
      real(dp) t344
      real(dp) t349
      real(dp) t35
      real(dp) t352
      real(dp) t355
      real(dp) t363
      real(dp) t366
      real(dp) t368
      real(dp) t37
      real(dp) t372
      real(dp) t373
      real(dp) t377
      real(dp) t38
      real(dp) t382
      real(dp) t384
      real(dp) t386
      real(dp) t394
      real(dp) t396
      real(dp) t40
      real(dp) t402
      real(dp) t403
      real(dp) t407
      real(dp) t41
      real(dp) t415
      real(dp) t416
      real(dp) t417
      real(dp) t419
      real(dp) t421
      real(dp) t422
      real(dp) t425
      real(dp) t429
      real(dp) t43
      real(dp) t434
      real(dp) t438
      real(dp) t44
      real(dp) t441
      real(dp) t445
      real(dp) t447
      real(dp) t449
      real(dp) t451
      real(dp) t453
      real(dp) t455
      real(dp) t46
      real(dp) t461
      real(dp) t462
      real(dp) t463
      real(dp) t464
      real(dp) t465
      real(dp) t466
      real(dp) t468
      real(dp) t47
      real(dp) t470
      real(dp) t473
      real(dp) t476
      real(dp) t481
      real(dp) t482
      real(dp) t483
      real(dp) t484
      real(dp) t486
      real(dp) t488
      real(dp) t49
      real(dp) t491
      real(dp) t493
      real(dp) t497
      real(dp) t499
      real(dp) t5
      real(dp) t512
      real(dp) t515
      real(dp) t519
      real(dp) t52
      real(dp) t522
      real(dp) t526
      real(dp) t532
      real(dp) t54
      real(dp) t540
      real(dp) t542
      real(dp) t549
      real(dp) t55
      real(dp) t554
      real(dp) t558
      real(dp) t561
      real(dp) t577
      real(dp) t579
      real(dp) t58
      real(dp) t583
      real(dp) t589
      real(dp) t59
      real(dp) t599
      real(dp) t6
      real(dp) t600
      real(dp) t601
      real(dp) t602
      real(dp) t606
      real(dp) t61
      real(dp) t613
      real(dp) t62
      real(dp) t63
      real(dp) t64
      real(dp) t66
      real(dp) t68
      real(dp) t69
      real(dp) t7
      real(dp) t70
      real(dp) t72
      real(dp) t73
      real(dp) t75
      real(dp) t77
      real(dp) t8
      real(dp) t81
      real(dp) t82
      real(dp) t83
      real(dp) t85
      real(dp) t86
      real(dp) t88
      real(dp) t9
      real(dp) t90
      real(dp) t92
      real(dp) t93
      real(dp) t95
      real(dp) t96
      real(dp) t98
      real(dp) t99
      real(dp) threecomp(3)

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
             x = xx(i,j,k)
             y = yy(i,j,k)
             z = zz(i,j,k)

        t1 = omm*x
        t2 = t1+ph1
        t3 = cos(t2)
        t5 = omm*y
        t6 = t5+ph1
        t7 = cos(t6)
        t8 = omm*z
        t9 = t8+ph1
        t10 = sin(t9)
        t13 = cv*t
        t15 = om*(x-t13)
        t16 = cos(t15)
        t19 = om*y+ph
        t20 = sin(t19)
        t22 = om*z+ph
        t23 = sin(t22)
        t24 = t20*t23
        t25 = t16*om*t24
        t27 = sin(t2)
        t31 = sin(t15)
        t33 = om**2
        t34 = t33*t20
        t35 = t34*t23
        t37 = t1+ph2
        t38 = cos(t37)
        t40 = t8+ph2
        t41 = sin(t40)
        t43 = t38*omm*t7*t41
        t44 = cos(t19)
        t46 = om*t23
        t47 = t31*t44*t46
        t49 = sin(t37)
        t52 = 2+t49*t7*t41
        t54 = t33*t44
        t55 = t54*t23
        t58 = t1+ph3
        t59 = cos(t58)
        t61 = t5+ph2
        t62 = cos(t61)
        t63 = t8+ph3
        t64 = sin(t63)
        t66 = t59*omm*t62*t64
        t68 = cos(t22)
        t69 = t68*om
        t70 = t31*t20*t69
        t72 = sin(t58)
        t73 = t72*t62
        t75 = 2+t73*t64
        t77 = t34*t68
        t81 = om*x+ph
        t82 = cos(t81)
        t83 = t82*om
        t85 = om*(y-t13)
        t86 = sin(t85)
        t88 = t83*t86*t23
        t90 = sin(t81)
        t92 = t33*t86
        t93 = t92*t23
        t95 = t1+ph4
        t96 = cos(t95)
        t98 = t8+ph4
        t99 = sin(t98)
        t102 = cos(t85)
        t104 = t90*t102*t46
        t106 = sin(t95)
        t109 = 2+t106*t62*t99
        t111 = t33*t102
        t112 = t111*t23
        t114 = t1+ph5
        t115 = cos(t114)
        t117 = t5+ph3
        t118 = cos(t117)
        t119 = t8+ph5
        t120 = sin(t119)
        t122 = t115*omm*t118*t120
        t124 = t90*t86*t69
        t126 = sin(t114)
        t127 = t126*t118
        t129 = 2+t127*t120
        t130 = t129*t82
        t131 = t92*t68
        t134 = om*(z-t13)
        t135 = sin(t134)
        t136 = t20*t135
        t137 = t83*t136
        t139 = -t3*omm*t7*t10*t25+(10+t27*t7*t10)*t31*t35-t43*t47-2*t52*
     #t16*t55-t66*t70-2*t75*t16*t77-t43*t88+t52*t90*t93-t96*omm*t62*t99*
     #t104-t109*t82*t112-t122*t124-t130*t131-t66*t137
        t141 = t34*t135
        t145 = t90*t44*om*t135
        t147 = t54*t135
        t149 = t1+ph6
        t150 = cos(t149)
        t152 = t8+ph6
        t153 = sin(t152)
        t157 = cos(t134)
        t159 = t90*t20*t157*om
        t161 = sin(t149)
        t162 = t161*t118
        t164 = 2+t162*t153
        t166 = t34*t157
        t168 = sin(t6)
        t170 = omm*t41
        t173 = t1+ph7
        t174 = sin(t173)
        t175 = t5+ph4
        t176 = sin(t175)
        t179 = sin(t8+ph7)
        t180 = omm*t179
        t181 = t174*t176*t180
        t183 = cos(t175)
        t186 = 10+t174*t183*t179
        t189 = t1+ph8
        t190 = sin(t189)
        t192 = t8+ph8
        t193 = sin(t192)
        t194 = omm*t193
        t195 = t190*t176*t194
        t197 = t190*t183
        t199 = 2+t197*t193
        t201 = t54*t68
        t207 = t1+ph9
        t208 = sin(t207)
        t209 = t5+ph5
        t210 = sin(t209)
        t212 = t8+ph9
        t213 = sin(t212)
        t214 = omm*t213
        t215 = t208*t210*t214
        t217 = t75*t90*t141-t122*t145-t130*t147-t150*omm*t118*t153*t159-
     #t164*t82*t166+t49*t168*t170*t25+t181*t47+t186*t31*t35+t195*t70-2*t
     #199*t31*t201+t181*t88-t186*t82*t112+t215*t104
        t219 = cos(t209)
        t222 = 2+t208*t219*t213
        t225 = t1+ph10
        t226 = sin(t225)
        t229 = t226*t210*omm*t10
        t231 = t226*t219
        t233 = 2+t231*t10
        t234 = t233*t90
        t235 = t111*t68
        t238 = t199*t82
        t242 = t1+ph11
        t243 = sin(t242)
        t244 = t5+ph6
        t245 = sin(t244)
        t249 = cos(t244)
        t250 = t243*t249
        t252 = 2+t250*t41
        t254 = t54*t157
        t256 = cos(t63)
        t257 = t256*omm
        t260 = cos(t192)
        t261 = t260*omm
        t262 = t197*t261
        t264 = t1+ph12
        t265 = sin(t264)
        t266 = t265*t249
        t267 = t266*t257
        t270 = 10+t266*t64
        t273 = t222*t90*t93+t229*t124-t234*t235+t195*t137-t238*t147+t229
     #*t145+t234*t141+t243*t245*t170*t159-t252*t90*t254-t73*t257*t25-t26
     #2*t47-t267*t70+t270*t31*t35
        t276 = t1+ph13
        t277 = sin(t276)
        t278 = t5+ph7
        t279 = cos(t278)
        t280 = t277*t279
        t281 = cos(t98)
        t286 = 2+t280*t99
        t289 = t1+ph14
        t290 = sin(t289)
        t291 = t290*t279
        t292 = cos(t119)
        t293 = t292*omm
        t294 = t291*t293
        t297 = 2+t291*t120
        t298 = t297*t90
        t305 = t1+ph15
        t306 = sin(t305)
        t307 = t5+ph8
        t308 = cos(t307)
        t309 = t306*t308
        t310 = cos(t152)
        t311 = t310*omm
        t312 = t309*t311
        t315 = 2+t309*t153
        t319 = sin(t1+phm)
        t321 = cos(t5+phm)
        t324 = sin(t8+phm)
        t327 = amprho*(2+t319*t321*t324)
        t329 = cv**2
        t330 = t33*t329
        t333 = -t262*t88-t238*t131-t280*t281*omm*t104-t286*t90*t235-t294
     #*t124+t298*t93-t267*t137-t270*t82*t166-t294*t145-t298*t254-t312*t1
     #59+t315*t90*t141-t327*t31*t330*t24
        threecomp(1) = t139+t217+t273+t333
        t337 = cos(t173)
        t340 = t337*omm*t183*t179
        t344 = t199*t16
        t349 = cos(t189)
        t352 = t349*omm*t183*t193
        t355 = cos(t207)
        t363 = cos(t225)
        t366 = t363*omm*t219*t10
        t368 = t233*t82
        t372 = t186*t90*t93-t340*t88+t52*t31*t35-t344*t77-t43*t25-t186*t
     #16*t55-t352*t70-t340*t47-t355*omm*t219*t213*t104-2*t222*t82*t112-t
     #366*t124-2*t368*t131-t352*t137
        t373 = t199*t90
        t377 = cos(t242)
        t382 = t252*t82
        t384 = sin(t61)
        t386 = omm*t99
        t394 = sin(t278)
        t396 = t277*t394*t386
        t402 = sin(t1+ph16)
        t403 = sin(t307)
        t407 = t373*t141-t366*t145-t368*t147-t377*omm*t249*t41*t159-t382
     #*t166+t106*t384*t386*t25-t109*t16*t55+t215*t47+t222*t31*t35+t396*t
     #70-t286*t31*t201+t215*t88+t402*t403*t180*t104
        t415 = sin(t1+ph17)
        t416 = t5+ph9
        t417 = sin(t416)
        t419 = t415*t417*t194
        t421 = cos(t416)
        t422 = t415*t421
        t425 = (2+t422*t193)*t90
        t429 = t286*t82
        t434 = sin(t1+ph18)
        t438 = t434*t421
        t441 = (2+t438*t213)*t90
        t445 = t129*t16
        t447 = cos(t9)
        t449 = t231*t447*omm
        t451 = t233*t31
        t453 = (10+t402*t308*t179)*t90*t93+t419*t124-2*t425*t235+t396*t1
     #37-t429*t147+t419*t145+t425*t141+t434*t417*t214*t159-t441*t254-t12
     #7*t293*t25-t445*t77-t449*t47-t451*t201
        t455 = t297*t31
        t461 = sin(t1+ph19)
        t462 = t5+ph10
        t463 = cos(t462)
        t464 = t461*t463
        t465 = t8+ph10
        t466 = cos(t465)
        t468 = t464*t466*omm
        t470 = sin(t465)
        t473 = (10+t464*t470)*t90
        t476 = t297*t82
        t481 = sin(t1+ph20)
        t482 = t481*t463
        t483 = t8+ph11
        t484 = cos(t483)
        t486 = t482*t484*omm
        t488 = sin(t483)
        t491 = (2+t482*t488)*t90
        t493 = t327*t90
        t497 = -t294*t70+t455*t35-t449*t88-t422*t261*t104-t468*t124+t473
     #*t93-t294*t137-t476*t166-t468*t145-t473*t254-t486*t159+t491*t141-t
     #493*t92*t329*t23
        threecomp(2) = t372+t407+t453+t497
        t499 = cos(t276)
        t512 = cos(t264)
        t515 = t512*omm*t249*t64
        t519 = cos(t289)
        t522 = t519*omm*t279*t120
        t526 = -t499*omm*t279*t99*t104-t429*t112+t75*t31*t35-t352*t88-t3
     #52*t47-t66*t25-t270*t16*t77-t515*t70-t344*t55+t373*t93-t522*t124-t
     #476*t131-t515*t137
        t532 = cos(t305)
        t540 = sin(t117)
        t542 = omm*t120
        t549 = t290*t394*t542
        t554 = t270*t90*t141-t522*t145-2*t476*t147-t532*omm*t308*t153*t1
     #59-2*t315*t82*t166+t126*t540*t542*t25-t445*t55+t229*t47+t451*t35+t
     #549*t70-t455*t201+t229*t88-t368*t112
        t558 = sin(t462)
        t561 = t461*t558*omm*t470
        t577 = cos(t40)
        t579 = t250*t577*omm
        t583 = t419*t104+t425*t93+t561*t124-t473*t235+t549*t137+t561*t14
     #5+t473*t141+t481*t558*omm*t488*t159-2*t491*t254-t162*t311*t25-t164
     #*t16*t77-t579*t47-t252*t31*t201
        t589 = cos(t212)
        t599 = sin(t1+ph21)
        t600 = t599*t463
        t601 = t8+ph12
        t602 = cos(t601)
        t606 = sin(t601)
        t613 = -t312*t70+t315*t31*t35-t579*t88-t382*t131-t438*t589*omm*t
     #104-t441*t235-t486*t124+t491*t93-t312*t137-t486*t145-t600*t602*omm
     #*t159+(10+t600*t606)*t90*t141-t493*t136*t330
        threecomp(3) = t526+t554+t583+t613

      fo(1,i,j,k) = threecomp(1)
      fo(2,i,j,k) = threecomp(2)
      fo(3,i,j,k) = threecomp(3)
      enddo
      enddo
      enddo

      return
      end
      
