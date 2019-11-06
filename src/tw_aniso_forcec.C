#include "EW.h"
#include "caliper.h"
#include "sw4.h"
void EW::tw_aniso_force_ci(int ifirst, int ilast, int jfirst, int jlast,
                           int kfirst, int klast, float_sw4* fo, float_sw4 t,
                           float_sw4 om, float_sw4 cv, float_sw4 ph,
                           float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                           float_sw4 phc[21], float_sw4 h, float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);
#pragma omp parallel
  {
    float_sw4 t1, t10, t102, t104, t106, t109, t111, t112, t114, t115, t117,
        t118, t119, t120, t122, t124, t126, t127, t129, t13, t130, t131, t134,
        t135, t136, t137, t139, t141, t145, t147, t149, t15, t150, t152, t153,
        t157, t159, t16, t161, t162, t164, t166, t168, t170, t173, t174, t175,
        t176, t179, t180, t181, t183, t186, t189, t19, t190, t192, t193, t194,
        t195, t197, t199, t2, t20, t201, t207, t208, t209, t210, t212, t213,
        t214, t215, t217, t219, t22, t222, t225, t226, t229, t23, t231, t233,
        t234, t235, t238, t24, t242, t243, t244, t245, t249, t25, t250, t252,
        t254, t256, t257, t260, t261, t262, t264, t265, t266, t267, t27, t270,
        t273, t276, t277, t278, t279, t280, t281, t286, t289, t290, t291, t292,
        t293, t294, t297, t298, t3, t305, t306, t307, t308, t309, t31, t310,
        t311, t312, t315, t319, t321, t324, t327, t329, t33, t330, t333, t337,
        t34, t340, t344, t349, t35, t352, t355, t363, t366, t368, t37, t372,
        t373, t377, t38, t382, t384, t386, t394, t396, t40, t402, t403, t407,
        t41, t415, t416, t417, t419, t421, t422, t425, t429, t43, t434, t438,
        t44, t441, t445, t447, t449, t451, t453, t455, t46, t461, t462, t463,
        t464, t465, t466, t468, t47, t470, t473, t476, t481, t482, t483, t484,
        t486, t488, t49, t491, t493, t497, t499, t5, t512, t515, t519, t52,
        t522, t526, t532, t54, t540, t542, t549, t55, t554, t558, t561, t577,
        t579, t58, t583, t589, t59, t599, t6, t600, t601, t602, t606, t61, t613,
        t62, t63, t64, t66, t68, t69, t7, t70, t72, t73, t75, t77, t8, t81, t82,
        t83, t85, t86, t88, t9, t90, t92, t93, t95, t96, t98, t99, threecomp[3],
        ph1, ph2, ph3, ph4, ph5, ph6, ph7, ph8, ph9, ph10, ph11, ph12, ph13,
        ph14, ph15, ph16, ph17, ph18, ph19, ph20, ph21;
    // extract all the phase angles for the stress matrix
    ph1 = phc[0];
    ph2 = phc[1];
    ph3 = phc[2];
    ph4 = phc[3];
    ph5 = phc[4];
    ph6 = phc[5];
    ph7 = phc[6];
    ph8 = phc[7];
    ph9 = phc[8];
    ph10 = phc[9];
    ph11 = phc[10];
    ph12 = phc[11];
    ph13 = phc[12];
    ph14 = phc[13];
    ph15 = phc[14];
    ph16 = phc[15];
    ph17 = phc[16];
    ph18 = phc[17];
    ph19 = phc[18];
    ph20 = phc[19];
    ph21 = phc[20];

#pragma omp for
    for (int k = kfirst; k <= klast; k++)
      for (int j = jfirst; j <= jlast; j++)
        for (int i = ifirst; i <= ilast; i++) {
          float_sw4 x = (i - 1) * h;
          float_sw4 y = (j - 1) * h;
          float_sw4 z = zmin + (k - 1) * h;

          t1 = omm * x;
          t2 = t1 + ph1;
          t3 = cos(t2);
          t5 = omm * y;
          t6 = t5 + ph1;
          t7 = cos(t6);
          t8 = omm * z;
          t9 = t8 + ph1;
          t10 = sin(t9);
          t13 = cv * t;
          t15 = om * (x - t13);
          t16 = cos(t15);
          t19 = om * y + ph;
          t20 = sin(t19);
          t22 = om * z + ph;
          t23 = sin(t22);
          t24 = t20 * t23;
          t25 = t16 * om * t24;
          t27 = sin(t2);
          t31 = sin(t15);
          t33 = om * om;
          t34 = t33 * t20;
          t35 = t34 * t23;
          t37 = t1 + ph2;
          t38 = cos(t37);
          t40 = t8 + ph2;
          t41 = sin(t40);
          t43 = t38 * omm * t7 * t41;
          t44 = cos(t19);
          t46 = om * t23;
          t47 = t31 * t44 * t46;
          t49 = sin(t37);
          t52 = 2 + t49 * t7 * t41;
          t54 = t33 * t44;
          t55 = t54 * t23;
          t58 = t1 + ph3;
          t59 = cos(t58);
          t61 = t5 + ph2;
          t62 = cos(t61);
          t63 = t8 + ph3;
          t64 = sin(t63);
          t66 = t59 * omm * t62 * t64;
          t68 = cos(t22);
          t69 = t68 * om;
          t70 = t31 * t20 * t69;
          t72 = sin(t58);
          t73 = t72 * t62;
          t75 = 2 + t73 * t64;
          t77 = t34 * t68;
          t81 = om * x + ph;
          t82 = cos(t81);
          t83 = t82 * om;
          t85 = om * (y - t13);
          t86 = sin(t85);
          t88 = t83 * t86 * t23;
          t90 = sin(t81);
          t92 = t33 * t86;
          t93 = t92 * t23;
          t95 = t1 + ph4;
          t96 = cos(t95);
          t98 = t8 + ph4;
          t99 = sin(t98);
          t102 = cos(t85);
          t104 = t90 * t102 * t46;
          t106 = sin(t95);
          t109 = 2 + t106 * t62 * t99;
          t111 = t33 * t102;
          t112 = t111 * t23;
          t114 = t1 + ph5;
          t115 = cos(t114);
          t117 = t5 + ph3;
          t118 = cos(t117);
          t119 = t8 + ph5;
          t120 = sin(t119);
          t122 = t115 * omm * t118 * t120;
          t124 = t90 * t86 * t69;
          t126 = sin(t114);
          t127 = t126 * t118;
          t129 = 2 + t127 * t120;
          t130 = t129 * t82;
          t131 = t92 * t68;
          t134 = om * (z - t13);
          t135 = sin(t134);
          t136 = t20 * t135;
          t137 = t83 * t136;
          t139 = -t3 * omm * t7 * t10 * t25 +
                 (10 + t27 * t7 * t10) * t31 * t35 - t43 * t47 -
                 2 * t52 * t16 * t55 - t66 * t70 - 2 * t75 * t16 * t77 -
                 t43 * t88 + t52 * t90 * t93 - t96 * omm * t62 * t99 * t104 -
                 t109 * t82 * t112 - t122 * t124 - t130 * t131 - t66 * t137;
          t141 = t34 * t135;
          t145 = t90 * t44 * om * t135;
          t147 = t54 * t135;
          t149 = t1 + ph6;
          t150 = cos(t149);
          t152 = t8 + ph6;
          t153 = sin(t152);
          t157 = cos(t134);
          t159 = t90 * t20 * t157 * om;
          t161 = sin(t149);
          t162 = t161 * t118;
          t164 = 2 + t162 * t153;
          t166 = t34 * t157;
          t168 = sin(t6);
          t170 = omm * t41;
          t173 = t1 + ph7;
          t174 = sin(t173);
          t175 = t5 + ph4;
          t176 = sin(t175);
          t179 = sin(t8 + ph7);
          t180 = omm * t179;
          t181 = t174 * t176 * t180;
          t183 = cos(t175);
          t186 = 10 + t174 * t183 * t179;
          t189 = t1 + ph8;
          t190 = sin(t189);
          t192 = t8 + ph8;
          t193 = sin(t192);
          t194 = omm * t193;
          t195 = t190 * t176 * t194;
          t197 = t190 * t183;
          t199 = 2 + t197 * t193;
          t201 = t54 * t68;
          t207 = t1 + ph9;
          t208 = sin(t207);
          t209 = t5 + ph5;
          t210 = sin(t209);
          t212 = t8 + ph9;
          t213 = sin(t212);
          t214 = omm * t213;
          t215 = t208 * t210 * t214;
          t217 = t75 * t90 * t141 - t122 * t145 - t130 * t147 -
                 t150 * omm * t118 * t153 * t159 - t164 * t82 * t166 +
                 t49 * t168 * t170 * t25 + t181 * t47 + t186 * t31 * t35 +
                 t195 * t70 - 2 * t199 * t31 * t201 + t181 * t88 -
                 t186 * t82 * t112 + t215 * t104;
          t219 = cos(t209);
          t222 = 2 + t208 * t219 * t213;
          t225 = t1 + ph10;
          t226 = sin(t225);
          t229 = t226 * t210 * omm * t10;
          t231 = t226 * t219;
          t233 = 2 + t231 * t10;
          t234 = t233 * t90;
          t235 = t111 * t68;
          t238 = t199 * t82;
          t242 = t1 + ph11;
          t243 = sin(t242);
          t244 = t5 + ph6;
          t245 = sin(t244);
          t249 = cos(t244);
          t250 = t243 * t249;
          t252 = 2 + t250 * t41;
          t254 = t54 * t157;
          t256 = cos(t63);
          t257 = t256 * omm;
          t260 = cos(t192);
          t261 = t260 * omm;
          t262 = t197 * t261;
          t264 = t1 + ph12;
          t265 = sin(t264);
          t266 = t265 * t249;
          t267 = t266 * t257;
          t270 = 10 + t266 * t64;
          t273 = t222 * t90 * t93 + t229 * t124 - t234 * t235 + t195 * t137 -
                 t238 * t147 + t229 * t145 + t234 * t141 +
                 t243 * t245 * t170 * t159 - t252 * t90 * t254 -
                 t73 * t257 * t25 - t262 * t47 - t267 * t70 + t270 * t31 * t35;
          t276 = t1 + ph13;
          t277 = sin(t276);
          t278 = t5 + ph7;
          t279 = cos(t278);
          t280 = t277 * t279;
          t281 = cos(t98);
          t286 = 2 + t280 * t99;
          t289 = t1 + ph14;
          t290 = sin(t289);
          t291 = t290 * t279;
          t292 = cos(t119);
          t293 = t292 * omm;
          t294 = t291 * t293;
          t297 = 2 + t291 * t120;
          t298 = t297 * t90;
          t305 = t1 + ph15;
          t306 = sin(t305);
          t307 = t5 + ph8;
          t308 = cos(t307);
          t309 = t306 * t308;
          t310 = cos(t152);
          t311 = t310 * omm;
          t312 = t309 * t311;
          t315 = 2 + t309 * t153;
          t319 = sin(t1 + phm);
          t321 = cos(t5 + phm);
          t324 = sin(t8 + phm);
          t327 = amprho * (2 + t319 * t321 * t324);
          t329 = cv * cv;
          t330 = t33 * t329;
          t333 = -t262 * t88 - t238 * t131 - t280 * t281 * omm * t104 -
                 t286 * t90 * t235 - t294 * t124 + t298 * t93 - t267 * t137 -
                 t270 * t82 * t166 - t294 * t145 - t298 * t254 - t312 * t159 +
                 t315 * t90 * t141 - t327 * t31 * t330 * t24;
          threecomp[0] = t139 + t217 + t273 + t333;
          t337 = cos(t173);
          t340 = t337 * omm * t183 * t179;
          t344 = t199 * t16;
          t349 = cos(t189);
          t352 = t349 * omm * t183 * t193;
          t355 = cos(t207);
          t363 = cos(t225);
          t366 = t363 * omm * t219 * t10;
          t368 = t233 * t82;
          t372 = t186 * t90 * t93 - t340 * t88 + t52 * t31 * t35 - t344 * t77 -
                 t43 * t25 - t186 * t16 * t55 - t352 * t70 - t340 * t47 -
                 t355 * omm * t219 * t213 * t104 - 2 * t222 * t82 * t112 -
                 t366 * t124 - 2 * t368 * t131 - t352 * t137;
          t373 = t199 * t90;
          t377 = cos(t242);
          t382 = t252 * t82;
          t384 = sin(t61);
          t386 = omm * t99;
          t394 = sin(t278);
          t396 = t277 * t394 * t386;
          t402 = sin(t1 + ph16);
          t403 = sin(t307);
          t407 = t373 * t141 - t366 * t145 - t368 * t147 -
                 t377 * omm * t249 * t41 * t159 - t382 * t166 +
                 t106 * t384 * t386 * t25 - t109 * t16 * t55 + t215 * t47 +
                 t222 * t31 * t35 + t396 * t70 - t286 * t31 * t201 +
                 t215 * t88 + t402 * t403 * t180 * t104;
          t415 = sin(t1 + ph17);
          t416 = t5 + ph9;
          t417 = sin(t416);
          t419 = t415 * t417 * t194;
          t421 = cos(t416);
          t422 = t415 * t421;
          t425 = (2 + t422 * t193) * t90;
          t429 = t286 * t82;
          t434 = sin(t1 + ph18);
          t438 = t434 * t421;
          t441 = (2 + t438 * t213) * t90;
          t445 = t129 * t16;
          t447 = cos(t9);
          t449 = t231 * t447 * omm;
          t451 = t233 * t31;
          t453 = (10 + t402 * t308 * t179) * t90 * t93 + t419 * t124 -
                 2 * t425 * t235 + t396 * t137 - t429 * t147 + t419 * t145 +
                 t425 * t141 + t434 * t417 * t214 * t159 - t441 * t254 -
                 t127 * t293 * t25 - t445 * t77 - t449 * t47 - t451 * t201;
          t455 = t297 * t31;
          t461 = sin(t1 + ph19);
          t462 = t5 + ph10;
          t463 = cos(t462);
          t464 = t461 * t463;
          t465 = t8 + ph10;
          t466 = cos(t465);
          t468 = t464 * t466 * omm;
          t470 = sin(t465);
          t473 = (10 + t464 * t470) * t90;
          t476 = t297 * t82;
          t481 = sin(t1 + ph20);
          t482 = t481 * t463;
          t483 = t8 + ph11;
          t484 = cos(t483);
          t486 = t482 * t484 * omm;
          t488 = sin(t483);
          t491 = (2 + t482 * t488) * t90;
          t493 = t327 * t90;
          t497 = -t294 * t70 + t455 * t35 - t449 * t88 - t422 * t261 * t104 -
                 t468 * t124 + t473 * t93 - t294 * t137 - t476 * t166 -
                 t468 * t145 - t473 * t254 - t486 * t159 + t491 * t141 -
                 t493 * t92 * t329 * t23;
          threecomp[1] = t372 + t407 + t453 + t497;
          t499 = cos(t276);
          t512 = cos(t264);
          t515 = t512 * omm * t249 * t64;
          t519 = cos(t289);
          t522 = t519 * omm * t279 * t120;
          t526 = -t499 * omm * t279 * t99 * t104 - t429 * t112 +
                 t75 * t31 * t35 - t352 * t88 - t352 * t47 - t66 * t25 -
                 t270 * t16 * t77 - t515 * t70 - t344 * t55 + t373 * t93 -
                 t522 * t124 - t476 * t131 - t515 * t137;
          t532 = cos(t305);
          t540 = sin(t117);
          t542 = omm * t120;
          t549 = t290 * t394 * t542;
          t554 = t270 * t90 * t141 - t522 * t145 - 2 * t476 * t147 -
                 t532 * omm * t308 * t153 * t159 - 2 * t315 * t82 * t166 +
                 t126 * t540 * t542 * t25 - t445 * t55 + t229 * t47 +
                 t451 * t35 + t549 * t70 - t455 * t201 + t229 * t88 -
                 t368 * t112;
          t558 = sin(t462);
          t561 = t461 * t558 * omm * t470;
          t577 = cos(t40);
          t579 = t250 * t577 * omm;
          t583 = t419 * t104 + t425 * t93 + t561 * t124 - t473 * t235 +
                 t549 * t137 + t561 * t145 + t473 * t141 +
                 t481 * t558 * omm * t488 * t159 - 2 * t491 * t254 -
                 t162 * t311 * t25 - t164 * t16 * t77 - t579 * t47 -
                 t252 * t31 * t201;
          t589 = cos(t212);
          t599 = sin(t1 + ph21);
          t600 = t599 * t463;
          t601 = t8 + ph12;
          t602 = cos(t601);
          t606 = sin(t601);
          t613 = -t312 * t70 + t315 * t31 * t35 - t579 * t88 - t382 * t131 -
                 t438 * t589 * omm * t104 - t441 * t235 - t486 * t124 +
                 t491 * t93 - t312 * t137 - t486 * t145 -
                 t600 * t602 * omm * t159 + (10 + t600 * t606) * t90 * t141 -
                 t493 * t136 * t330;
          threecomp[2] = t526 + t554 + t583 + t613;
          size_t ind = base + i + ni * j + nij * k;
          fo[ind] = threecomp[0];
          fo[ind + nijk] = threecomp[1];
          fo[ind + 2 * nijk] = threecomp[2];
        }
  }
}

//-----------------------------------------------------------------------
void EW::tw_aniso_curvi_force_ci(int ifirst, int ilast, int jfirst, int jlast,
                                 int kfirst, int klast,
                                 float_sw4* __restrict__ fo, float_sw4 t,
                                 float_sw4 om, float_sw4 cv, float_sw4 ph,
                                 float_sw4 omm, float_sw4 phm, float_sw4 amprho,
                                 float_sw4 phc[21], float_sw4* __restrict__ xx,
                                 float_sw4* __restrict__ yy,
                                 float_sw4* __restrict__ zz)
// curvilinear case
{
  SW4_MARK_FUNCTION;
  const size_t ni = ilast - ifirst + 1;
  const size_t nij = ni * (jlast - jfirst + 1);
  const size_t nijk = nij * (klast - kfirst + 1);
  const size_t base = -(ifirst + ni * jfirst + nij * kfirst);
#pragma omp parallel
  {
    // local variables
    float_sw4 t1, t10, t102, t104, t106, t109, t111, t112, t114, t115, t117,
        t118, t119, t120, t122, t124, t126, t127, t129, t13, t130, t131, t134,
        t135, t136, t137, t139, t141, t145, t147, t149, t15, t150, t152, t153,
        t157, t159, t16, t161, t162, t164, t166, t168, t170, t173, t174, t175,
        t176, t179, t180, t181, t183, t186, t189, t19, t190, t192, t193, t194,
        t195, t197, t199, t2, t20, t201, t207, t208, t209, t210, t212, t213,
        t214, t215, t217, t219, t22, t222, t225, t226, t229, t23, t231, t233,
        t234, t235, t238, t24, t242, t243, t244, t245, t249, t25, t250, t252,
        t254, t256, t257, t260, t261, t262, t264, t265, t266, t267, t27, t270,
        t273, t276, t277, t278, t279, t280, t281, t286, t289, t290, t291, t292,
        t293, t294, t297, t298, t3, t305, t306, t307, t308, t309, t31, t310,
        t311, t312, t315, t319, t321, t324, t327, t329, t33, t330, t333, t337,
        t34, t340, t344, t349, t35, t352, t355, t363, t366, t368, t37, t372,
        t373, t377, t38, t382, t384, t386, t394, t396, t40, t402, t403, t407,
        t41, t415, t416, t417, t419, t421, t422, t425, t429, t43, t434, t438,
        t44, t441, t445, t447, t449, t451, t453, t455, t46, t461, t462, t463,
        t464, t465, t466, t468, t47, t470, t473, t476, t481, t482, t483, t484,
        t486, t488, t49, t491, t493, t497, t499, t5, t512, t515, t519, t52,
        t522, t526, t532, t54, t540, t542, t549, t55, t554, t558, t561, t577,
        t579, t58, t583, t589, t59, t599, t6, t600, t601, t602, t606, t61, t613,
        t62, t63, t64, t66, t68, t69, t7, t70, t72, t73, t75, t77, t8, t81, t82,
        t83, t85, t86, t88, t9, t90, t92, t93, t95, t96, t98, t99, threecomp[3],
        ph1, ph2, ph3, ph4, ph5, ph6, ph7, ph8, ph9, ph10, ph11, ph12, ph13,
        ph14, ph15, ph16, ph17, ph18, ph19, ph20, ph21;

    // extract all the phase angles for the stress matrix
    ph1 = phc[0];
    ph2 = phc[1];
    ph3 = phc[2];
    ph4 = phc[3];
    ph5 = phc[4];
    ph6 = phc[5];
    ph7 = phc[6];
    ph8 = phc[7];
    ph9 = phc[8];
    ph10 = phc[9];
    ph11 = phc[10];
    ph12 = phc[11];
    ph13 = phc[12];
    ph14 = phc[13];
    ph15 = phc[14];
    ph16 = phc[15];
    ph17 = phc[16];
    ph18 = phc[17];
    ph19 = phc[18];
    ph20 = phc[19];
    ph21 = phc[20];

#pragma omp for
    for (int k = kfirst; k <= klast; k++)
      for (int j = jfirst; j <= jlast; j++)
        for (int i = ifirst; i <= ilast; i++) {
          size_t ind = base + i + ni * j + nij * k;
          float_sw4 x = xx[ind];
          float_sw4 y = yy[ind];
          float_sw4 z = zz[ind];

          t1 = omm * x;
          t2 = t1 + ph1;
          t3 = cos(t2);
          t5 = omm * y;
          t6 = t5 + ph1;
          t7 = cos(t6);
          t8 = omm * z;
          t9 = t8 + ph1;
          t10 = sin(t9);
          t13 = cv * t;
          t15 = om * (x - t13);
          t16 = cos(t15);
          t19 = om * y + ph;
          t20 = sin(t19);
          t22 = om * z + ph;
          t23 = sin(t22);
          t24 = t20 * t23;
          t25 = t16 * om * t24;
          t27 = sin(t2);
          t31 = sin(t15);
          t33 = om * om;
          t34 = t33 * t20;
          t35 = t34 * t23;
          t37 = t1 + ph2;
          t38 = cos(t37);
          t40 = t8 + ph2;
          t41 = sin(t40);
          t43 = t38 * omm * t7 * t41;
          t44 = cos(t19);
          t46 = om * t23;
          t47 = t31 * t44 * t46;
          t49 = sin(t37);
          t52 = 2 + t49 * t7 * t41;
          t54 = t33 * t44;
          t55 = t54 * t23;
          t58 = t1 + ph3;
          t59 = cos(t58);
          t61 = t5 + ph2;
          t62 = cos(t61);
          t63 = t8 + ph3;
          t64 = sin(t63);
          t66 = t59 * omm * t62 * t64;
          t68 = cos(t22);
          t69 = t68 * om;
          t70 = t31 * t20 * t69;
          t72 = sin(t58);
          t73 = t72 * t62;
          t75 = 2 + t73 * t64;
          t77 = t34 * t68;
          t81 = om * x + ph;
          t82 = cos(t81);
          t83 = t82 * om;
          t85 = om * (y - t13);
          t86 = sin(t85);
          t88 = t83 * t86 * t23;
          t90 = sin(t81);
          t92 = t33 * t86;
          t93 = t92 * t23;
          t95 = t1 + ph4;
          t96 = cos(t95);
          t98 = t8 + ph4;
          t99 = sin(t98);
          t102 = cos(t85);
          t104 = t90 * t102 * t46;
          t106 = sin(t95);
          t109 = 2 + t106 * t62 * t99;
          t111 = t33 * t102;
          t112 = t111 * t23;
          t114 = t1 + ph5;
          t115 = cos(t114);
          t117 = t5 + ph3;
          t118 = cos(t117);
          t119 = t8 + ph5;
          t120 = sin(t119);
          t122 = t115 * omm * t118 * t120;
          t124 = t90 * t86 * t69;
          t126 = sin(t114);
          t127 = t126 * t118;
          t129 = 2 + t127 * t120;
          t130 = t129 * t82;
          t131 = t92 * t68;
          t134 = om * (z - t13);
          t135 = sin(t134);
          t136 = t20 * t135;
          t137 = t83 * t136;
          t139 = -t3 * omm * t7 * t10 * t25 +
                 (10 + t27 * t7 * t10) * t31 * t35 - t43 * t47 -
                 2 * t52 * t16 * t55 - t66 * t70 - 2 * t75 * t16 * t77 -
                 t43 * t88 + t52 * t90 * t93 - t96 * omm * t62 * t99 * t104 -
                 t109 * t82 * t112 - t122 * t124 - t130 * t131 - t66 * t137;
          t141 = t34 * t135;
          t145 = t90 * t44 * om * t135;
          t147 = t54 * t135;
          t149 = t1 + ph6;
          t150 = cos(t149);
          t152 = t8 + ph6;
          t153 = sin(t152);
          t157 = cos(t134);
          t159 = t90 * t20 * t157 * om;
          t161 = sin(t149);
          t162 = t161 * t118;
          t164 = 2 + t162 * t153;
          t166 = t34 * t157;
          t168 = sin(t6);
          t170 = omm * t41;
          t173 = t1 + ph7;
          t174 = sin(t173);
          t175 = t5 + ph4;
          t176 = sin(t175);
          t179 = sin(t8 + ph7);
          t180 = omm * t179;
          t181 = t174 * t176 * t180;
          t183 = cos(t175);
          t186 = 10 + t174 * t183 * t179;
          t189 = t1 + ph8;
          t190 = sin(t189);
          t192 = t8 + ph8;
          t193 = sin(t192);
          t194 = omm * t193;
          t195 = t190 * t176 * t194;
          t197 = t190 * t183;
          t199 = 2 + t197 * t193;
          t201 = t54 * t68;
          t207 = t1 + ph9;
          t208 = sin(t207);
          t209 = t5 + ph5;
          t210 = sin(t209);
          t212 = t8 + ph9;
          t213 = sin(t212);
          t214 = omm * t213;
          t215 = t208 * t210 * t214;
          t217 = t75 * t90 * t141 - t122 * t145 - t130 * t147 -
                 t150 * omm * t118 * t153 * t159 - t164 * t82 * t166 +
                 t49 * t168 * t170 * t25 + t181 * t47 + t186 * t31 * t35 +
                 t195 * t70 - 2 * t199 * t31 * t201 + t181 * t88 -
                 t186 * t82 * t112 + t215 * t104;
          t219 = cos(t209);
          t222 = 2 + t208 * t219 * t213;
          t225 = t1 + ph10;
          t226 = sin(t225);
          t229 = t226 * t210 * omm * t10;
          t231 = t226 * t219;
          t233 = 2 + t231 * t10;
          t234 = t233 * t90;
          t235 = t111 * t68;
          t238 = t199 * t82;
          t242 = t1 + ph11;
          t243 = sin(t242);
          t244 = t5 + ph6;
          t245 = sin(t244);
          t249 = cos(t244);
          t250 = t243 * t249;
          t252 = 2 + t250 * t41;
          t254 = t54 * t157;
          t256 = cos(t63);
          t257 = t256 * omm;
          t260 = cos(t192);
          t261 = t260 * omm;
          t262 = t197 * t261;
          t264 = t1 + ph12;
          t265 = sin(t264);
          t266 = t265 * t249;
          t267 = t266 * t257;
          t270 = 10 + t266 * t64;
          t273 = t222 * t90 * t93 + t229 * t124 - t234 * t235 + t195 * t137 -
                 t238 * t147 + t229 * t145 + t234 * t141 +
                 t243 * t245 * t170 * t159 - t252 * t90 * t254 -
                 t73 * t257 * t25 - t262 * t47 - t267 * t70 + t270 * t31 * t35;
          t276 = t1 + ph13;
          t277 = sin(t276);
          t278 = t5 + ph7;
          t279 = cos(t278);
          t280 = t277 * t279;
          t281 = cos(t98);
          t286 = 2 + t280 * t99;
          t289 = t1 + ph14;
          t290 = sin(t289);
          t291 = t290 * t279;
          t292 = cos(t119);
          t293 = t292 * omm;
          t294 = t291 * t293;
          t297 = 2 + t291 * t120;
          t298 = t297 * t90;
          t305 = t1 + ph15;
          t306 = sin(t305);
          t307 = t5 + ph8;
          t308 = cos(t307);
          t309 = t306 * t308;
          t310 = cos(t152);
          t311 = t310 * omm;
          t312 = t309 * t311;
          t315 = 2 + t309 * t153;
          t319 = sin(t1 + phm);
          t321 = cos(t5 + phm);
          t324 = sin(t8 + phm);
          t327 = amprho * (2 + t319 * t321 * t324);
          t329 = cv * cv;
          t330 = t33 * t329;
          t333 = -t262 * t88 - t238 * t131 - t280 * t281 * omm * t104 -
                 t286 * t90 * t235 - t294 * t124 + t298 * t93 - t267 * t137 -
                 t270 * t82 * t166 - t294 * t145 - t298 * t254 - t312 * t159 +
                 t315 * t90 * t141 - t327 * t31 * t330 * t24;
          threecomp[0] = t139 + t217 + t273 + t333;
          t337 = cos(t173);
          t340 = t337 * omm * t183 * t179;
          t344 = t199 * t16;
          t349 = cos(t189);
          t352 = t349 * omm * t183 * t193;
          t355 = cos(t207);
          t363 = cos(t225);
          t366 = t363 * omm * t219 * t10;
          t368 = t233 * t82;
          t372 = t186 * t90 * t93 - t340 * t88 + t52 * t31 * t35 - t344 * t77 -
                 t43 * t25 - t186 * t16 * t55 - t352 * t70 - t340 * t47 -
                 t355 * omm * t219 * t213 * t104 - 2 * t222 * t82 * t112 -
                 t366 * t124 - 2 * t368 * t131 - t352 * t137;
          t373 = t199 * t90;
          t377 = cos(t242);
          t382 = t252 * t82;
          t384 = sin(t61);
          t386 = omm * t99;
          t394 = sin(t278);
          t396 = t277 * t394 * t386;
          t402 = sin(t1 + ph16);
          t403 = sin(t307);
          t407 = t373 * t141 - t366 * t145 - t368 * t147 -
                 t377 * omm * t249 * t41 * t159 - t382 * t166 +
                 t106 * t384 * t386 * t25 - t109 * t16 * t55 + t215 * t47 +
                 t222 * t31 * t35 + t396 * t70 - t286 * t31 * t201 +
                 t215 * t88 + t402 * t403 * t180 * t104;
          t415 = sin(t1 + ph17);
          t416 = t5 + ph9;
          t417 = sin(t416);
          t419 = t415 * t417 * t194;
          t421 = cos(t416);
          t422 = t415 * t421;
          t425 = (2 + t422 * t193) * t90;
          t429 = t286 * t82;
          t434 = sin(t1 + ph18);
          t438 = t434 * t421;
          t441 = (2 + t438 * t213) * t90;
          t445 = t129 * t16;
          t447 = cos(t9);
          t449 = t231 * t447 * omm;
          t451 = t233 * t31;
          t453 = (10 + t402 * t308 * t179) * t90 * t93 + t419 * t124 -
                 2 * t425 * t235 + t396 * t137 - t429 * t147 + t419 * t145 +
                 t425 * t141 + t434 * t417 * t214 * t159 - t441 * t254 -
                 t127 * t293 * t25 - t445 * t77 - t449 * t47 - t451 * t201;
          t455 = t297 * t31;
          t461 = sin(t1 + ph19);
          t462 = t5 + ph10;
          t463 = cos(t462);
          t464 = t461 * t463;
          t465 = t8 + ph10;
          t466 = cos(t465);
          t468 = t464 * t466 * omm;
          t470 = sin(t465);
          t473 = (10 + t464 * t470) * t90;
          t476 = t297 * t82;
          t481 = sin(t1 + ph20);
          t482 = t481 * t463;
          t483 = t8 + ph11;
          t484 = cos(t483);
          t486 = t482 * t484 * omm;
          t488 = sin(t483);
          t491 = (2 + t482 * t488) * t90;
          t493 = t327 * t90;
          t497 = -t294 * t70 + t455 * t35 - t449 * t88 - t422 * t261 * t104 -
                 t468 * t124 + t473 * t93 - t294 * t137 - t476 * t166 -
                 t468 * t145 - t473 * t254 - t486 * t159 + t491 * t141 -
                 t493 * t92 * t329 * t23;
          threecomp[1] = t372 + t407 + t453 + t497;
          t499 = cos(t276);
          t512 = cos(t264);
          t515 = t512 * omm * t249 * t64;
          t519 = cos(t289);
          t522 = t519 * omm * t279 * t120;
          t526 = -t499 * omm * t279 * t99 * t104 - t429 * t112 +
                 t75 * t31 * t35 - t352 * t88 - t352 * t47 - t66 * t25 -
                 t270 * t16 * t77 - t515 * t70 - t344 * t55 + t373 * t93 -
                 t522 * t124 - t476 * t131 - t515 * t137;
          t532 = cos(t305);
          t540 = sin(t117);
          t542 = omm * t120;
          t549 = t290 * t394 * t542;
          t554 = t270 * t90 * t141 - t522 * t145 - 2 * t476 * t147 -
                 t532 * omm * t308 * t153 * t159 - 2 * t315 * t82 * t166 +
                 t126 * t540 * t542 * t25 - t445 * t55 + t229 * t47 +
                 t451 * t35 + t549 * t70 - t455 * t201 + t229 * t88 -
                 t368 * t112;
          t558 = sin(t462);
          t561 = t461 * t558 * omm * t470;
          t577 = cos(t40);
          t579 = t250 * t577 * omm;
          t583 = t419 * t104 + t425 * t93 + t561 * t124 - t473 * t235 +
                 t549 * t137 + t561 * t145 + t473 * t141 +
                 t481 * t558 * omm * t488 * t159 - 2 * t491 * t254 -
                 t162 * t311 * t25 - t164 * t16 * t77 - t579 * t47 -
                 t252 * t31 * t201;
          t589 = cos(t212);
          t599 = sin(t1 + ph21);
          t600 = t599 * t463;
          t601 = t8 + ph12;
          t602 = cos(t601);
          t606 = sin(t601);
          t613 = -t312 * t70 + t315 * t31 * t35 - t579 * t88 - t382 * t131 -
                 t438 * t589 * omm * t104 - t441 * t235 - t486 * t124 +
                 t491 * t93 - t312 * t137 - t486 * t145 -
                 t600 * t602 * omm * t159 + (10 + t600 * t606) * t90 * t141 -
                 t493 * t136 * t330;
          threecomp[2] = t526 + t554 + t583 + t613;
          fo[ind] = threecomp[0];
          fo[ind + nijk] = threecomp[1];
          fo[ind + 2 * nij] = threecomp[2];
        }
  }
}

//-----------------------------------------------------------------------
void EW::tw_aniso_free_surf_z_ci(int ifirst, int ilast, int jfirst, int jlast,
                                 int kfirst, int klast, int kz, float_sw4 t,
                                 float_sw4 om, float_sw4 cv, float_sw4 ph,
                                 float_sw4 omm, float_sw4 phc[21],
                                 float_sw4* __restrict__ bforce, float_sw4 h,
                                 float_sw4 zmin) {
  SW4_MARK_FUNCTION;
  //      real(dp), intent(out):: bforce(3,ifirst:ilast,jfirst:jlast)
  const size_t ni = ilast - ifirst + 1;
  // const size_t nij   = ni*(jlast-jfirst+1);
  const size_t base = -(ifirst + ni * jfirst);
#pragma omp parallel
  {
    float_sw4 ph1, ph2, ph3, ph4, ph5, ph6, ph7, ph8, ph9, ph10, ph11, ph12,
        ph13, ph14, ph15, ph17, ph18, ph19, ph20, ph21, forces[3], t1, t10,
        t101, t104, t106, t108, t110, t113, t115, t122, t124, t127, t129, t13,
        t137, t139, t146, t148, t15, t151, t154, t16, t160, t163, t166, t169,
        t176, t179, t181, t189, t19, t192, t20, t202, t205, t21, t23, t24, t25,
        t28, t3, t30, t33, t35, t36, t38, t39, t4, t40, t43, t45, t48, t50, t52,
        t55, t56, t59, t6, t60, t62, t65, t67, t70, t73, t75, t77, t8, t80, t83,
        t85, t86, t88, t92, t93, t94, t96, t99;

    // extract all the phase angles for the stress matrix
    ph1 = phc[0];
    ph2 = phc[1];
    ph3 = phc[2];
    ph4 = phc[3];
    ph5 = phc[4];
    ph6 = phc[5];
    ph7 = phc[6];
    ph8 = phc[7];
    ph9 = phc[8];
    ph10 = phc[9];
    ph11 = phc[10];
    ph12 = phc[11];
    ph13 = phc[12];
    ph14 = phc[13];
    ph15 = phc[14];
    //   ph16 = phc[15];
    ph17 = phc[16];
    ph18 = phc[17];
    ph19 = phc[18];
    ph20 = phc[19];
    ph21 = phc[20];

    float_sw4 z = (kz - 1) * h + zmin;
#pragma omp for
    for (int j = jfirst; j <= jlast; j++) {
      float_sw4 y = (j - 1) * h;
      for (int i = ifirst; i <= ilast; i++) {
        float_sw4 x = (i - 1) * h;
        t1 = omm * x;
        t3 = sin(t1 + ph3);
        t4 = omm * y;
        t6 = cos(t4 + ph2);
        t8 = omm * z;
        t10 = sin(t8 + ph3);
        t13 = cv * t;
        t15 = om * (x - t13);
        t16 = cos(t15);
        t19 = om * y + ph;
        t20 = sin(t19);
        t21 = om * t20;
        t23 = om * z + ph;
        t24 = sin(t23);
        t25 = t21 * t24;
        t28 = sin(t1 + ph8);
        t30 = cos(t4 + ph4);
        t33 = sin(t8 + ph8);
        t35 = 2 + t28 * t30 * t33;
        t36 = sin(t15);
        t38 = cos(t19);
        t39 = t38 * om;
        t40 = t39 * t24;
        t43 = sin(t1 + ph12);
        t45 = cos(t4 + ph6);
        t48 = 10 + t43 * t45 * t10;
        t50 = cos(t23);
        t52 = t20 * t50 * om;
        t55 = om * x + ph;
        t56 = cos(t55);
        t59 = om * (y - t13);
        t60 = sin(t59);
        t62 = om * t60 * t24;
        t65 = sin(t1 + ph13);
        t67 = cos(t4 + ph7);
        t70 = sin(t8 + ph4);
        t73 = sin(t55);
        t75 = cos(t59);
        t77 = t75 * om * t24;
        t80 = sin(t1 + ph14);
        t83 = sin(t8 + ph5);
        t85 = 2 + t80 * t67 * t83;
        t86 = t85 * t73;
        t88 = t60 * t50 * om;
        t92 = om * (z - t13);
        t93 = sin(t92);
        t94 = t21 * t93;
        t96 = t39 * t93;
        t99 = sin(t1 + ph15);
        t101 = cos(t4 + ph8);
        t104 = sin(t8 + ph6);
        t106 = 2 + t99 * t101 * t104;
        t108 = cos(t92);
        t110 = t20 * t108 * om;
        forces[0] = (2 + t3 * t6 * t10) * t16 * t25 + t35 * t36 * t40 +
                    t48 * t36 * t52 + t35 * t56 * t62 +
                    (2 + t65 * t67 * t70) * t73 * t77 + t86 * t88 +
                    t48 * t56 * t94 + t86 * t96 + t106 * t73 * t110;
        t113 = sin(t1 + ph5);
        t115 = cos(t4 + ph3);
        t122 = sin(t1 + ph10);
        t124 = cos(t4 + ph5);
        t127 = sin(t8 + ph1);
        t129 = 2 + t122 * t124 * t127;
        t137 = sin(t1 + ph17);
        t139 = cos(t4 + ph9);
        t146 = sin(t1 + ph19);
        t148 = cos(t4 + ph10);
        t151 = sin(t8 + ph10);
        t154 = (10 + t146 * t148 * t151) * t73;
        t160 = sin(t1 + ph20);
        t163 = sin(t8 + ph11);
        t166 = (2 + t160 * t148 * t163) * t73;
        forces[1] = (2 + t113 * t115 * t83) * t16 * t25 + t129 * t36 * t40 +
                    t85 * t36 * t52 + t129 * t56 * t62 +
                    (2 + t137 * t139 * t33) * t73 * t77 + t154 * t88 +
                    t85 * t56 * t94 + t154 * t96 + t166 * t110;
        t169 = sin(t1 + ph6);
        t176 = sin(t1 + ph11);
        t179 = sin(t8 + ph2);
        t181 = 2 + t176 * t45 * t179;
        t189 = sin(t1 + ph18);
        t192 = sin(t8 + ph9);
        t202 = sin(t1 + ph21);
        t205 = sin(t8 + ph12);
        forces[2] = (2 + t169 * t115 * t104) * t16 * t25 + t181 * t36 * t40 +
                    t106 * t36 * t52 + t181 * t56 * t62 +
                    (2 + t189 * t139 * t192) * t73 * t77 + t166 * t88 +
                    t106 * t56 * t94 + t166 * t96 +
                    (10 + t202 * t148 * t205) * t73 * t110;
        bforce[3 * (base + i + ni * j)] = forces[0];
        bforce[1 + 3 * (base + i + ni * j)] = forces[1];
        bforce[2 + 3 * (base + i + ni * j)] = forces[2];
      }
    }
  }
}
