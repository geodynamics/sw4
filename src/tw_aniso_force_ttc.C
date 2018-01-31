#include "sw4.h"

#include "EW.h"

void EW::tw_aniso_force_tt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       float_sw4* __restrict__ fo, float_sw4 t,float_sw4 om,float_sw4 cv,float_sw4 ph,
			       float_sw4 omm,float_sw4 phm,float_sw4 amprho,float_sw4 phc[21],float_sw4 h,
			       float_sw4 zmin)
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);

#pragma omp parallel
   {
float_sw4 forces[3],t1,t10,t100,t101,t102,t103,t104,t105,t106,t108,t109,t11,t110,t111,t112,t113,t114,t115,t117,t118,t121,t124,t125,t130,t132,t133,t134,t135,t136,t137,t139,t14,t142,t146,t147,t148,t149,t15,t151,t152,t154,t155,t156,t158,t159,t161,t163,t167,t168,t169,t17,t170,t173,t174,t175,t176,t177,t179,t18,t181,t184,t185,t187,t189,t193,t195,t196,t198,t199,t2,t20,t200,t201,t204,t206,t210,t214,t217,t218,t22,t221,t222,t224,t225,t227,t23,t233,t236,t24,t240,t243,t245,t25,t251,t252,t253,t254,t255,t257,t261,t267,t268,t269,t27,t270,t271,t273,t277,t278,t280,t282,t285,t288,t29,t296,t299,t3,t301,t303,t31,t310,t314,t315,t316,t317,t32,t322,t328,t33,t333,t335,t338,t34,t341,t343,t344,t352,t353,t354,t356,t36,t361,t366,t37,t373,t375,t378,t381,t384,t397,t398,t399,t4,t401,t404,t405,t406,t407,t408,t409,t41,t412,t416,t419,t42,t421,t425,t428,t430,t436,t438,t439,t44,t443,t444,t447,t450,t46,t464,t467,t469,t472,t474,t477,t48,t483,t488,t49,t490,t493,t496,t498,t5,t50,t511,t512,t513,t514,t516,t518,t52,t521,t524,t526,t53,t531,t533,t538,t540,t55,t552,t554,t556,t559,t56,t562,t563,t567,t57,t572,t573,t574,t575,t58,t581,t586,t589,t596,t6,t60,t602,t605,t616,t619,t62,t63,t630,t643,t65,t657,t659,t66,t662,t68,t684,t70,t71,t72,t73,t76,t78,t79,t8,t80,t81,t84,t85,t86,t87,t88,t89,t9,t90,t91,t93,t94,t95,t96,t97,t98,t99,ph1,ph2,ph3,ph4,ph5,ph6,ph7,ph8,ph9,ph10,ph11,ph12,ph13,ph14,ph15,ph16,ph17,ph18,ph19,ph20,ph21;
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
   for( int k=kfirst; k<= klast; k++ )
      for( int j=jfirst; j<= jlast; j++ )
	 for( int i=ifirst; i<= ilast; i++ )
	 {
            float_sw4 x=(i-1)*h;
            float_sw4 y=(j-1)*h;
            float_sw4 z=zmin+(k-1)*h;
        t1 = omm*x;
        t2 = t1+ph11;
        t3 = sin(t2);
        t4 = omm*y;
        t5 = t4+ph6;
        t6 = sin(t5);
        t8 = omm*z;
        t9 = t8+ph2;
        t10 = sin(t9);
        t11 = omm*t10;
        t14 = om*x+ph;
        t15 = sin(t14);
        t17 = om*y+ph;
        t18 = sin(t17);
        t20 = cv*t;
        t22 = om*(z-t20);
        t23 = cos(t22);
        t24 = om*om;
        t25 = t24*om;
        t27 = cv*cv;
        t29 = t15*t18*t23*t25*t27;
        t31 = t1+ph1;
        t32 = sin(t31);
        t33 = t4+ph1;
        t34 = cos(t33);
        t36 = t8+ph1;
        t37 = sin(t36);
        t41 = om*(x-t20);
        t42 = sin(t41);
        t44 = t24*t24;
        t46 = t27*t18;
        t48 = om*z+ph;
        t49 = sin(t48);
        t50 = t46*t49;
        t52 = t1+ph5;
        t53 = cos(t52);
        t55 = t4+ph3;
        t56 = cos(t55)   ;
        t57 = t8+ph5;
        t58 = sin(t57);
        t60 = t53*omm*t56*t58;
        t62 = om*(y-t20);
        t63 = sin(t62);
        t65 = t25*t27;
        t66 = cos(t48);
        t68 = t15*t63*t65*t66;
        t70 = t1+ph10;
        t71 = sin(t70);
        t72 = t4+ph5;
        t73 = sin(t72);
        t76 = t71*t73*omm*t37;
        t78 = t1+ph7;
        t79 = sin(t78);
        t80 = t4+ph4;
        t81 = sin(t80);
        t84 = sin(t8+ph7);
        t85 = omm*t84;
        t86 = t79*t81*t85;
        t87 = t42*t25;
        t88 = cos(t17);
        t89 = t27*t88;
        t90 = t89*t49;
        t91 = t87*t90;
        t93 = t1+ph12;
        t94 = sin(t93);
        t95 = cos(t5);
        t96 = t94*t95;
        t97 = t8+ph3;
        t98 = cos(t97);
        t99 = t98*omm;
        t100 = t96*t99;
        t101 = cos(t14);
        t102 = t101*t25;
        t103 = sin(t22);
        t104 = t18*t103;
        t105 = t104*t27;
        t106 = t102*t105;
        t108 = t1+ph14;
        t109 = sin(t108);
        t110 = t4+ph7;
        t111 = cos(t110);
        t112 = t109*t111;
        t113 = cos(t57);
        t114 = t113*omm;
        t115 = t112*t114;
        t117 = t1+ph2;
        t118 = sin(t117);
        t121 = 2+t118*t34*t10;
        t124 = t63*t27;
        t125 = t124*t49;
        t130 = t15*t88*t25*t103*t27;
        t132 = t1+ph8;
        t133 = sin(t132);
        t134 = cos(t80);
        t135 = t133*t134;
        t136 = t8+ph8;
        t137 = sin(t136);
        t139 = 2+t135*t137;
        t142 = t89*t66;
        t146 = omm*t137;
        t147 = t133*t81*t146;
        t148 = t46*t66;
        t149 = t87*t148;
        t151 = cos(t72);
        t152 = t71*t151;
        t154 = 2+t152*t37;
        t155 = t154*t15;
        t156 = cos(t62);
        t158 = t44*t27;
        t159 = t158*t66;
        t161 = sin(t97);
        t163 = 10+t96*t161;
        t167 = -t3*t6*t11*t29-(10+t32*t34*t37)*t42*t44*t50+t60*t68-t76*t68-t86*t91+
	   t100*t106+t115*t68-t121*t15*t44*t125+t115*t130+2*t139*t42*t44*t142-
	   t147*t149+t155*t156*t159-t163*t42*t44*t50;
        t168 = cos(t136);
        t169 = t168*omm;
        t170 = t135*t169;
        t173 = t1+ph3;
        t174 = sin(t173);
        t175 = t4+ph2;
        t176 = cos(t175);
        t177 = t174*t176;
        t179 = cos(t41);
        t181 = t179*t25*t50;
        t184 = t139*t101*t44;
        t185 = t124*t66;
        t187 = t3*t95;
        t189 = 2+t187*t10;
        t193 = t44*t23*t27;
        t195 = t1+ph9;
        t196 = sin(t195);
        t198 = t8+ph9;
        t199 = sin(t198);
        t200 = omm*t199;
        t201 = t196*t73*t200;
        t204 = t15*t156*t65*t49;
        t206 = t102*t125;
        t210 = 10+t79*t134*t84;
        t214 = t156*t27*t49;
        t217 = 2+t112*t58;
        t218 = t217*t15;
        t221 = t1+ph4;
        t222 = sin(t221);
        t224 = t8+ph4;
        t225 = sin(t224);
        t227 = 2+t222*t176*t225;
        t233 = t44*t103*t27;
        t236 = 2+t177*t161;
        t240 = cos(t173);
        t243 = t240*omm*t176*t161;
        t245 = t170*t91-t147*t106+t177*t99*t181+t184*t185+t189*t15*t88*t193-
	   t201*t204+t170*t206+t210*t101*t44*t214+t218*t88*t193+t227*t101
	   *t44*t214-t155*t18*t233-t236*t15*t44*t105+t243*t106;
        t251 = t1+ph6;
        t252 = sin(t251);
        t253 = t252*t56;
        t254 = t8+ph6;
        t255 = sin(t254);
        t257 = 2+t253*t255;
        t261 = t18*t23*t27;
        t267 = t1+ph15;
        t268 = sin(t267);
        t269 = t4+ph8;
        t270 = cos(t269);
        t271 = t268*t270;
        t273 = 2+t271*t255;
        t277 = sin(t52);
        t278 = t277*t56;
        t280 = 2+t278*t58;
        t282 = t280*t101*t44;
        t285 = t158*t49;
        t288 = t88*t103*t27;
        t296 = cos(t117);
        t299 = t296*omm*t34*t10;
        t301 = 2*t121*t179*t44*t90+t257*t101*t44*t261-t86*t206+t163*t101
     *t44*t261-t273*t15*t18*t233+t282*t185-t218*t63*t285+t184*t288-t76*
      t130-t210*t42*t44*t50+t282*t288+t60*t130+t299*t206;
        t303 = cos(t31);
        t310 = 2+t196*t151*t199;
        t314 = t1+ph13;
        t315 = sin(t314);
        t316 = t315*t111;
        t317 = cos(t224);
        t322 = cos(t221);
        t328 = 2+t316*t225;
        t333 = sin(t1+phm);
        t335 = cos(t4+phm);
        t338 = sin(t8+phm);
        t341 = amprho*(2+t333*t335*t338);
        t343 = t27*t27;
        t344 = t44*t343;
        t352 = cos(t254);
        t353 = t352*omm;
        t354 = t271*t353;
        t356 = sin(t33);
        t361 = cos(t251);
        t366 = t100*t149+t303*omm*t34*t37*t181-t310*t15*t63*t285+t316*t317*
	   omm*t204+t299*t91+t322*omm*t176*t225*t204+t328*t15*t156*t159+
	   t341*t42*t344*t18*t49+2*t236*t179*t44*t148+t354*t29-t118*t356*
	   t11*t181+t243*t149+t361*omm*t56*t255*t29;
        forces[0] = t167+t245+t301+t366;
        t373 = t280*t179*t44;
        t375 = cos(t70);
        t378 = t375*omm*t151*t37;
        t381 = t154*t101*t44;
        t384 = sin(t1+ph16);
        t397 = sin(t1+ph17);
        t398 = t4+ph9;
        t399 = sin(t398);
        t401 = t397*t399*t146;
        t404 = sin(t1+ph20);
        t405 = t4+ph10;
        t406 = cos(t405);
        t407 = t404*t406;
        t408 = t8+ph11;
        t409 = sin(t408);
        t412 = (2+t407*t409)*t15;
        t416 = t328*t101*t44;
        t419 = t139*t179*t44;
        t421 = -t201*t91-t210*t15*t44*t125+t373*t148+t378*t130+t381*t288
	   -(10+t384*t270*t84)*t15*t63*t285+t115*t149+t328*t42*t44*t142+t299*
	   t181-t401*t68-t412*t18*t233+t416*t288+t419*t148;
        t425 = cos(t132);
        t428 = t425*omm*t134*t137;
        t430 = cos(t2);
        t436 = t154*t42*t44;
        t438 = cos(t398);
        t439 = t397*t438;
        t443 = sin(t1+ph18);
        t444 = t443*t438;
        t447 = (2+t444*t199)*t15;
        t450 = cos(t195);
        t464 = (2+t439*t137)*t15;
        t467 = cos(t36);
        t469 = t152*t467*omm;
        t472 = t139*t15*t44;
        t474 = t210*t179*t44*t90+t428*t149+t430*omm*t95*t10*t29+t436*t142+
	   t439*t169*t204+t447*t88*t193+t450*omm*t151*t199*t204+t115*t106+
	   t227*t179*t44*t90-t310*t42*t44*t50-t464*t18*t233+t469*t206-t472*t105;
        t477 = t189*t101*t44;
        t483 = sin(t269);
        t488 = sin(t175);
        t490 = omm*t225;
        t493 = cos(t78);
        t496 = t493*omm*t134*t84;
        t498 = t341*t15;
        t511 = sin(t1+ph19);
        t512 = t511*t406;
        t513 = t8+ph10;
        t514 = cos(t513);
        t516 = t512*t514*omm;
        t518 = sin(t513);
        t521 = (10+t512*t518)*t15;
        t524 = t477*t261+2*t310*t101*t44*t214-t384*t483*t85*t204+t378*t68-
	   t222*t488*t490*t181+t496*t91+t498*t63*t44*t343*t49+t428*t106-
	   t121*t42*t44*t50+t278*t114*t181-t201*t206+t516*t130+t521*t88*t193;
        t526 = t217*t42*t44;
        t531 = cos(t408);
        t533 = t407*t531*omm;
        t538 = sin(t110);
        t540 = t315*t538*t490;
        t552 = t217*t101*t44;
        t554 = -t526*t50+2*t464*t156*t159+t533*t29-t521*t63*t285+t496*t206-
	   t540*t149-t443*t399*t200*t29-t401*t130-t540*t106+t469*t91+t516*
	   t68+2*t381*t185+t552*t261;
        forces[1] = t421+t474+t524+t554;
        t556 = cos(t108);
        t559 = t556*omm*t111*t58;
        t562 = omm*t58;
        t563 = t109*t538*t562;
        t567 = cos(t198);
        t572 = sin(t1+ph21);
        t573 = t572*t406;
        t574 = t8+ph12;
        t575 = sin(t574);
        t581 = sin(t55);
        t586 = sin(t405);
        t589 = t511*t586*omm*t518;
        t596 = cos(t267);
        t602 = t559*t68-t563*t149+t419*t90+t373*t90+t444*t567*omm*t204-(
     10+t573*t575)*t15*t18*t233-t277*t581*t562*t181+t243*t181-t589*t130
     +t354*t106+2*t273*t101*t44*t261+t596*omm*t270*t255*t29-t472*t125;
        t605 = cos(t314);
        t616 = cos(t93);
        t619 = t616*omm*t95*t161;
        t630 = t381*t214+t552*t185+t605*omm*t111*t225*t204-t76*t91+t526*
	   t142-t412*t63*t285-t521*t18*t233+t619*t149+2*t412*t88*t193+
	   t428*t206-t563*t106+t354*t149-t163*t15*t44*t105;
        t643 = cos(t574);
        t657 = cos(t9);
        t659 = t187*t657*omm;
        t662 = t163*t179*t44*t148+t257*t179*t44*t148+t416*t214-t404*t586
     *omm*t409*t29+t573*t643*omm*t29+t447*t156*t159+t533*t130-t76*t206-
	   t464*t63*t285+t533*t68-t236*t42*t44*t50+t659*t206+t619*t106;
        t684 = -t273*t42*t44*t50-t401*t204-t436*t50+t189*t42*t44*t142+
	   t521*t156*t159+t659*t91+t428*t91+t498*t104*t344-t589*t68+
	   t477*t185+t253*t353*t181+2*t552*t288+t559*t130;
        forces[2] = t602+t630+t662+t684;

	size_t ind= base+i+ni*j+nij*k;
	fo[ind]        = forces[0];
	fo[ind+nijk]   = forces[1];
	fo[ind+2*nijk] = forces[2];
	 }
   }
}

//-----------------------------------------------------------------------
void EW::tw_aniso_curvi_force_tt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
				     float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
				     float_sw4 cv, float_sw4 ph, float_sw4 omm, float_sw4 phm,
				     float_sw4 amprho, float_sw4 phc[21],
				     float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
				     float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);

#pragma omp parallel
   {
      float_sw4 forces[3],t1,t10,t100,t101,t102,t103,t104,t105,t106,t108,t109,t11,t110,t111,t112,t113,t114,t115,t117,t118,t121,t124,t125,t130,t132,t133,t134,t135,t136,t137,t139,t14,t142,t146,t147,t148,t149,t15,t151,t152,t154,t155,t156,t158,t159,t161,t163,t167,t168,t169,t17,t170,t173,t174,t175,t176,t177,t179,t18,t181,t184,t185,t187,t189,t193,t195,t196,t198,t199,t2,t20,t200,t201,t204,t206,t210,t214,t217,t218,t22,t221,t222,t224,t225,t227,t23,t233,t236,t24,t240,t243,t245,t25,t251,t252,t253,t254,t255,t257,t261,t267,t268,t269,t27,t270,t271,t273,t277,t278,t280,t282,t285,t288,t29,t296,t299,t3,t301,t303,t31,t310,t314,t315,t316,t317,t32,t322,t328,t33,t333,t335,t338,t34,t341,t343,t344,t352,t353,t354,t356,t36,t361,t366,t37,t373,t375,t378,t381,t384,t397,t398,t399,t4,t401,t404,t405,t406,t407,t408,t409,t41,t412,t416,t419,t42,t421,t425,t428,t430,t436,t438,t439,t44,t443,t444,t447,t450,t46,t464,t467,t469,t472,t474,t477,t48,t483,t488,t49,t490,t493,t496,t498,t5,t50,t511,t512,t513,t514,t516,t518,t52,t521,t524,t526,t53,t531,t533,t538,t540,t55,t552,t554,t556,t559,t56,t562,t563,t567,t57,t572,t573,t574,t575,t58,t581,t586,t589,t596,t6,t60,t602,t605,t616,t619,t62,t63,t630,t643,t65,t657,t659,t66,t662,t68,t684,t70,t71,t72,t73,t76,t78,t79,t8,t80,t81,t84,t85,t86,t87,t88,t89,t9,t90,t91,t93,t94,t95,t96,t97,t98,t99,ph1,ph2,ph3,ph4,ph5,ph6,ph7,ph8,ph9,ph10,ph11,ph12,ph13,ph14,ph15,ph16,ph17,ph18,ph19,ph20,ph21;
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
   for( int k=kfirst; k<= klast; k++ )
      for( int j=jfirst; j<= jlast; j++ )
	 for( int i=ifirst; i<= ilast; i++ )
	 {
	    size_t ind = base + i + ni*j + nij*k;
            float_sw4 x=xx[ind];
            float_sw4 y=yy[ind];
            float_sw4 z=zz[ind];
        t1 = omm*x;
        t2 = t1+ph11;
        t3 = sin(t2);
        t4 = omm*y;
        t5 = t4+ph6;
        t6 = sin(t5);
        t8 = omm*z;
        t9 = t8+ph2;
        t10 = sin(t9);
        t11 = omm*t10;
        t14 = om*x+ph;
        t15 = sin(t14);
        t17 = om*y+ph;
        t18 = sin(t17);
        t20 = cv*t;
        t22 = om*(z-t20);
        t23 = cos(t22);
        t24 = om*om;
        t25 = t24*om;
        t27 = cv*cv;
        t29 = t15*t18*t23*t25*t27;
        t31 = t1+ph1;
        t32 = sin(t31);
        t33 = t4+ph1;
        t34 = cos(t33);
        t36 = t8+ph1;
        t37 = sin(t36);
        t41 = om*(x-t20);
        t42 = sin(t41);
        t44 = t24*t24;
        t46 = t27*t18;
        t48 = om*z+ph;
        t49 = sin(t48);
        t50 = t46*t49;
        t52 = t1+ph5;
        t53 = cos(t52);
        t55 = t4+ph3;
        t56 = cos(t55)   ;
        t57 = t8+ph5;
        t58 = sin(t57);
        t60 = t53*omm*t56*t58;
        t62 = om*(y-t20);
        t63 = sin(t62);
        t65 = t25*t27;
        t66 = cos(t48);
        t68 = t15*t63*t65*t66;
        t70 = t1+ph10;
        t71 = sin(t70);
        t72 = t4+ph5;
        t73 = sin(t72);
        t76 = t71*t73*omm*t37;
        t78 = t1+ph7;
        t79 = sin(t78);
        t80 = t4+ph4;
        t81 = sin(t80);
        t84 = sin(t8+ph7);
        t85 = omm*t84;
        t86 = t79*t81*t85;
        t87 = t42*t25;
        t88 = cos(t17);
        t89 = t27*t88;
        t90 = t89*t49;
        t91 = t87*t90;
        t93 = t1+ph12;
        t94 = sin(t93);
        t95 = cos(t5);
        t96 = t94*t95;
        t97 = t8+ph3;
        t98 = cos(t97);
        t99 = t98*omm;
        t100 = t96*t99;
        t101 = cos(t14);
        t102 = t101*t25;
        t103 = sin(t22);
        t104 = t18*t103;
        t105 = t104*t27;
        t106 = t102*t105;
        t108 = t1+ph14;
        t109 = sin(t108);
        t110 = t4+ph7;
        t111 = cos(t110);
        t112 = t109*t111;
        t113 = cos(t57);
        t114 = t113*omm;
        t115 = t112*t114;
        t117 = t1+ph2;
        t118 = sin(t117);
        t121 = 2+t118*t34*t10;
        t124 = t63*t27;
        t125 = t124*t49;
        t130 = t15*t88*t25*t103*t27;
        t132 = t1+ph8;
        t133 = sin(t132);
        t134 = cos(t80);
        t135 = t133*t134;
        t136 = t8+ph8;
        t137 = sin(t136);
        t139 = 2+t135*t137;
        t142 = t89*t66;
        t146 = omm*t137;
        t147 = t133*t81*t146;
        t148 = t46*t66;
        t149 = t87*t148;
        t151 = cos(t72);
        t152 = t71*t151;
        t154 = 2+t152*t37;
        t155 = t154*t15;
        t156 = cos(t62);
        t158 = t44*t27;
        t159 = t158*t66;
        t161 = sin(t97);
        t163 = 10+t96*t161;
        t167 = -t3*t6*t11*t29-(10+t32*t34*t37)*t42*t44*t50+t60*t68-t76*t68-t86*t91+
	   t100*t106+t115*t68-t121*t15*t44*t125+t115*t130+2*t139*t42*t44*t142-
	   t147*t149+t155*t156*t159-t163*t42*t44*t50;
        t168 = cos(t136);
        t169 = t168*omm;
        t170 = t135*t169;
        t173 = t1+ph3;
        t174 = sin(t173);
        t175 = t4+ph2;
        t176 = cos(t175);
        t177 = t174*t176;
        t179 = cos(t41);
        t181 = t179*t25*t50;
        t184 = t139*t101*t44;
        t185 = t124*t66;
        t187 = t3*t95;
        t189 = 2+t187*t10;
        t193 = t44*t23*t27;
        t195 = t1+ph9;
        t196 = sin(t195);
        t198 = t8+ph9;
        t199 = sin(t198);
        t200 = omm*t199;
        t201 = t196*t73*t200;
        t204 = t15*t156*t65*t49;
        t206 = t102*t125;
        t210 = 10+t79*t134*t84;
        t214 = t156*t27*t49;
        t217 = 2+t112*t58;
        t218 = t217*t15;
        t221 = t1+ph4;
        t222 = sin(t221);
        t224 = t8+ph4;
        t225 = sin(t224);
        t227 = 2+t222*t176*t225;
        t233 = t44*t103*t27;
        t236 = 2+t177*t161;
        t240 = cos(t173);
        t243 = t240*omm*t176*t161;
        t245 = t170*t91-t147*t106+t177*t99*t181+t184*t185+t189*t15*t88*t193-
	   t201*t204+t170*t206+t210*t101*t44*t214+t218*t88*t193+t227*t101
	   *t44*t214-t155*t18*t233-t236*t15*t44*t105+t243*t106;
        t251 = t1+ph6;
        t252 = sin(t251);
        t253 = t252*t56;
        t254 = t8+ph6;
        t255 = sin(t254);
        t257 = 2+t253*t255;
        t261 = t18*t23*t27;
        t267 = t1+ph15;
        t268 = sin(t267);
        t269 = t4+ph8;
        t270 = cos(t269);
        t271 = t268*t270;
        t273 = 2+t271*t255;
        t277 = sin(t52);
        t278 = t277*t56;
        t280 = 2+t278*t58;
        t282 = t280*t101*t44;
        t285 = t158*t49;
        t288 = t88*t103*t27;
        t296 = cos(t117);
        t299 = t296*omm*t34*t10;
        t301 = 2*t121*t179*t44*t90+t257*t101*t44*t261-t86*t206+t163*t101
     *t44*t261-t273*t15*t18*t233+t282*t185-t218*t63*t285+t184*t288-t76*
      t130-t210*t42*t44*t50+t282*t288+t60*t130+t299*t206;
        t303 = cos(t31);
        t310 = 2+t196*t151*t199;
        t314 = t1+ph13;
        t315 = sin(t314);
        t316 = t315*t111;
        t317 = cos(t224);
        t322 = cos(t221);
        t328 = 2+t316*t225;
        t333 = sin(t1+phm);
        t335 = cos(t4+phm);
        t338 = sin(t8+phm);
        t341 = amprho*(2+t333*t335*t338);
        t343 = t27*t27;
        t344 = t44*t343;
        t352 = cos(t254);
        t353 = t352*omm;
        t354 = t271*t353;
        t356 = sin(t33);
        t361 = cos(t251);
        t366 = t100*t149+t303*omm*t34*t37*t181-t310*t15*t63*t285+t316*t317*
	   omm*t204+t299*t91+t322*omm*t176*t225*t204+t328*t15*t156*t159+
	   t341*t42*t344*t18*t49+2*t236*t179*t44*t148+t354*t29-t118*t356*
	   t11*t181+t243*t149+t361*omm*t56*t255*t29;
        forces[0] = t167+t245+t301+t366;
        t373 = t280*t179*t44;
        t375 = cos(t70);
        t378 = t375*omm*t151*t37;
        t381 = t154*t101*t44;
        t384 = sin(t1+ph16);
        t397 = sin(t1+ph17);
        t398 = t4+ph9;
        t399 = sin(t398);
        t401 = t397*t399*t146;
        t404 = sin(t1+ph20);
        t405 = t4+ph10;
        t406 = cos(t405);
        t407 = t404*t406;
        t408 = t8+ph11;
        t409 = sin(t408);
        t412 = (2+t407*t409)*t15;
        t416 = t328*t101*t44;
        t419 = t139*t179*t44;
        t421 = -t201*t91-t210*t15*t44*t125+t373*t148+t378*t130+t381*t288
	   -(10+t384*t270*t84)*t15*t63*t285+t115*t149+t328*t42*t44*t142+t299*
	   t181-t401*t68-t412*t18*t233+t416*t288+t419*t148;
        t425 = cos(t132);
        t428 = t425*omm*t134*t137;
        t430 = cos(t2);
        t436 = t154*t42*t44;
        t438 = cos(t398);
        t439 = t397*t438;
        t443 = sin(t1+ph18);
        t444 = t443*t438;
        t447 = (2+t444*t199)*t15;
        t450 = cos(t195);
        t464 = (2+t439*t137)*t15;
        t467 = cos(t36);
        t469 = t152*t467*omm;
        t472 = t139*t15*t44;
        t474 = t210*t179*t44*t90+t428*t149+t430*omm*t95*t10*t29+t436*t142+
	   t439*t169*t204+t447*t88*t193+t450*omm*t151*t199*t204+t115*t106+
	   t227*t179*t44*t90-t310*t42*t44*t50-t464*t18*t233+t469*t206-t472*t105;
        t477 = t189*t101*t44;
        t483 = sin(t269);
        t488 = sin(t175);
        t490 = omm*t225;
        t493 = cos(t78);
        t496 = t493*omm*t134*t84;
        t498 = t341*t15;
        t511 = sin(t1+ph19);
        t512 = t511*t406;
        t513 = t8+ph10;
        t514 = cos(t513);
        t516 = t512*t514*omm;
        t518 = sin(t513);
        t521 = (10+t512*t518)*t15;
        t524 = t477*t261+2*t310*t101*t44*t214-t384*t483*t85*t204+t378*t68-
	   t222*t488*t490*t181+t496*t91+t498*t63*t44*t343*t49+t428*t106-
	   t121*t42*t44*t50+t278*t114*t181-t201*t206+t516*t130+t521*t88*t193;
        t526 = t217*t42*t44;
        t531 = cos(t408);
        t533 = t407*t531*omm;
        t538 = sin(t110);
        t540 = t315*t538*t490;
        t552 = t217*t101*t44;
        t554 = -t526*t50+2*t464*t156*t159+t533*t29-t521*t63*t285+t496*t206-
	   t540*t149-t443*t399*t200*t29-t401*t130-t540*t106+t469*t91+t516*
	   t68+2*t381*t185+t552*t261;
        forces[1] = t421+t474+t524+t554;
        t556 = cos(t108);
        t559 = t556*omm*t111*t58;
        t562 = omm*t58;
        t563 = t109*t538*t562;
        t567 = cos(t198);
        t572 = sin(t1+ph21);
        t573 = t572*t406;
        t574 = t8+ph12;
        t575 = sin(t574);
        t581 = sin(t55);
        t586 = sin(t405);
        t589 = t511*t586*omm*t518;
        t596 = cos(t267);
        t602 = t559*t68-t563*t149+t419*t90+t373*t90+t444*t567*omm*t204-(
     10+t573*t575)*t15*t18*t233-t277*t581*t562*t181+t243*t181-t589*t130
     +t354*t106+2*t273*t101*t44*t261+t596*omm*t270*t255*t29-t472*t125;
        t605 = cos(t314);
        t616 = cos(t93);
        t619 = t616*omm*t95*t161;
        t630 = t381*t214+t552*t185+t605*omm*t111*t225*t204-t76*t91+t526*
	   t142-t412*t63*t285-t521*t18*t233+t619*t149+2*t412*t88*t193+
	   t428*t206-t563*t106+t354*t149-t163*t15*t44*t105;
        t643 = cos(t574);
        t657 = cos(t9);
        t659 = t187*t657*omm;
        t662 = t163*t179*t44*t148+t257*t179*t44*t148+t416*t214-t404*t586
     *omm*t409*t29+t573*t643*omm*t29+t447*t156*t159+t533*t130-t76*t206-
	   t464*t63*t285+t533*t68-t236*t42*t44*t50+t659*t206+t619*t106;
        t684 = -t273*t42*t44*t50-t401*t204-t436*t50+t189*t42*t44*t142+
	   t521*t156*t159+t659*t91+t428*t91+t498*t104*t344-t589*t68+
	   t477*t185+t253*t353*t181+2*t552*t288+t559*t130;
        forces[2] = t602+t630+t662+t684;
	 
	fo[ind]        = forces[0];
	fo[ind+nijk]   = forces[1];
	fo[ind+2*nijk] = forces[2];
	 }
   }
}

