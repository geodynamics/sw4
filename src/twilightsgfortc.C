#include "sw4.h"
#include "EW.h"
#include "caliper.h"
//#include <math.h>
//#include <sys/types.h>
//-----------------------------------------------------------------------
void EW::forcingfortsg_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4 h, float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
			   float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
   RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<DEFAULT_LOOP3>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t10,t100,t102,t103,t104,t105,t106,t113,t114,t115,t116,t117,t119,t121,t122,t123,t124,t126,t129,t13,t132,t133,t14,t142,t143,t144,t145,t147,t148,t150,t152,t157,t158,t16,t167,t168,t17,t172,t173,t177,t180,t181,t184,t185,t19,t191,t195,t2,t20,t21,t214,t215,t220,t226,t229,t23,t231,t24,t250,t253,t26,t265,t27,t28,t3,t31,t32,t34,t35,t36,t37,t38,t41,t43,t44,t45,t48,t5,t50,t51,t56,t6,t61,t62,t63,t72,t73,t75,t76,t78,t80,t81,t82,t84,t85,t87,t88,t9,t91,t92,t96,t97,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
//       {
// 	 float_sw4 z=(k-1)*h+zmin;
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	 {
// 	    float_sw4 y=(j-1)*h;
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       float_sw4 x=(i-1)*h;
	       float_sw4 y=(j-1)*h;
	       float_sw4 z=(k-1)*h+zmin;
        t2 = omm*x+phm;
        t3 = sin(t2);
        t5 = omm*y+phm;
        t6 = cos(t5);
        t9 = omm*z+phm;
        t10 = sin(t9);
        t13 = amprho*(2+t3*t6*t10);
        t14 = c*t;
        t16 = om*(x-t14);
        t17 = sin(t16);
        t19 = om*om;
        t20 = c*c;
        t21 = t19*t20;
        t23 = om*y+ph;
        t24 = sin(t23);
        t26 = om*z+ph;
        t27 = sin(t26);
        t28 = t24*t27;
        t31 = omstrx*x;
        t32 = sin(t31);
        t34 = 1+t32/2;
        t35 = ampmu*t3;
        t36 = sin(t5);
        t37 = omm*t36;
        t38 = t37*t10;
        t41 = cos(t2);
        t43 = cos(t9);
        t44 = t37*t43;
        t45 = amplambda*t41*t44;
        t48 = cos(t16);
        t50 = om*t24;
        t51 = t50*t27;
        t56 = ampmu*(3+t41*t36*t10);
        t61 = amplambda*(2+t3*t36*t43);
        t62 = 2*t56+t61;
        t63 = cos(t31);
        t72 = t19*t24;
        t73 = t72*t27;
        t75 = omstry*y;
        t76 = sin(t75);
        t78 = 1+t76/2;
        t80 = om*x+ph;
        t81 = sin(t80);
        t82 = t78*t81;
        t84 = om*(y-t14);
        t85 = cos(t84);
        t87 = t85*om*t27;
        t88 = t82*t87;
        t91 = cos(t80);
        t92 = t91*t19;
        t96 = omstrz*z;
        t97 = sin(t96);
        t99 = 1+t97/2;
        t100 = t99*t81;
        t102 = om*(z-t14);
        t103 = cos(t102);
        t104 = t24*t103;
        t105 = t104*om;
        t106 = t100*t105;
        t113 = ampmu*t41;
        t114 = t113*t6;
        t115 = omm*t10;
        t116 = t34*t91;
        t117 = sin(t84);
        t119 = om*t117*t27;
        t121 = t78*t17;
        t122 = cos(t23);
        t123 = t122*om;
        t124 = t123*t27;
        t126 = t116*t119+t121*t124;
        t129 = t19*t85;
        t132 = cos(t75);
        t133 = t132*omstry;
        t142 = t113*omm;
        t143 = t36*t43;
        t144 = sin(t102);
        t145 = t50*t144;
        t147 = t99*t17;
        t148 = cos(t26);
        t150 = t24*t148*om;
        t152 = t116*t145+t147*t150;
        t157 = cos(t96);
        t158 = t157*omstrz;
        forces[0] = -t13*t17*t21*t28-t34*((-2*t35*t38+t45)*t34*t48*t51+
	   t62*t63*omstrx*t48*om*t28/2-t62*t34*t17*t73+t45*t88+t61*t78*t92*t85
     *t27+t45*t106+t61*t99*t92*t104)-t78*(t114*t115*t126+t56*(t116*t129
     *t27+t133*t17*t124/2-t121*t73))-t99*(t142*t143*t152+t56*(t116*t72*
      t103+t158*t17*t150/2-t147*t73));
        t167 = t13*t81;
        t168 = t117*t19;
        t172 = t35*omm;
        t173 = t36*t10;
        t177 = t63*omstrx*t91;
        t180 = t34*t81;
        t181 = t168*t27;
        t184 = t19*t122;
        t185 = t184*t27;
        t191 = t6*omm;
        t195 = amplambda*t3;
        t214 = t34*t48;
        t215 = t214*t51;
        t220 = t184*t103;
        t226 = t123*t144;
        t229 = t117*t148*om;
        t231 = t82*t226+t100*t229;
        forces[1] = -t167*t168*t20*t27-t34*(-t172*t173*t126+t56*(t177*t119/2-
	   t180*t181+t78*t48*t185))-t78*((2*t113*t191*t10+t195*t191*t43)
     *t78*t81*t87+t62*t132*omstry*t81*t85*om*t27/2-t62*t78*t81*t181+
       t195*t6*omm*t43*(t215+t106)+t61*(t214*t185+t100*t220))-t99*(t142*t143
      *t231+t56*(t82*t220+t158*t81*t229/2-t100*t181));
        t250 = t72*t144;
        t253 = t72*t148;
        t265 = t129*t148;
        forces[2] = -t167*t24*t144*t21-t34*(-t172*t173*t152+t56*(t177*t145/2-
           t180*t250+t99*t48*t253))-t78*(t114*t115*t231+t56*(t133*t81*
          t226/2-t82*t250+t100*t265))-t99*((2*t113*t44-t195*t38)*t99*t81*t105+
     t62*t157*omstrz*t81*t24*t103*om/2-t62*t99*t81*t250-t195*omm*t173*(
     t215+t88)+t61*(t214*t253+t82*t265));

	size_t ind = base+i+ni*j+nij*k;
	fo[ind] = forces[0];
	fo[ind+nijk] = forces[1];
	fo[ind+2*nijk] = forces[2];
      // 	    }
      // 	 }
      // }
			       });
}

//-----------------------------------------------------------------------
void EW::forcingttfortsg_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4 h, float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
			   float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
 #else
 using LOCAL_POL = DEFAULT_LOOP3;
 #endif
     RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t10,t100,t102,t103,t106,t107,t110,t120,t121,t122,t123,t124,t125,t127,t129,t13,t131,t133,t135,t138,t14,t142,t143,t150,t151,t157,t158,t159,t16,t160,t161,t163,t165,t166,t168,t17,t171,t174,t175,t187,t188,t19,t192,t193,t197,t198,t2,t20,t203,t21,t212,t216,t22,t223,t23,t236,t238,t242,t246,t25,t252,t256,t259,t26,t279,t28,t29,t3,t30,t303,t315,t33,t34,t36,t37,t38,t39,t40,t43,t45,t46,t47,t5,t50,t52,t53,t54,t59,t6,t64,t65,t66,t74,t77,t78,t80,t82,t83,t84,t86,t87,t88,t9,t90,t92,t96,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
//       {
// 	 float_sw4 z=(k-1)*h+zmin;
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	 {
// 	    float_sw4 y=(j-1)*h;
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       float_sw4 x=(i-1)*h;
	       float_sw4 y=(j-1)*h;
	       float_sw4 z=(k-1)*h+zmin;
        t2 = omm*x+phm;
        t3 = sin(t2);
        t5 = omm*y+phm;
        t6 = cos(t5);
        t9 = omm*z+phm;
        t10 = sin(t9);
        t13 = amprho*(2+t3*t6*t10);
        t14 = c*t;
        t16 = om*(x-t14);
        t17 = sin(t16);
        t19 = om*om;
        t20 = t19*t19;
        t21 = c*c;
        t22 = t21*t21;
        t23 = t20*t22;
        t25 = om*y+ph;
        t26 = sin(t25);
        t28 = om*z+ph;
        t29 = sin(t28);
        t30 = t26*t29;
        t33 = omstrx*x;
        t34 = sin(t33);
        t36 = 1+t34/2;
        t37 = ampmu*t3;
        t38 = sin(t5);
        t39 = omm*t38;
        t40 = t39*t10;
        t43 = cos(t2);
        t45 = cos(t9);
        t46 = t39*t45;
        t47 = amplambda*t43*t46;
        t50 = cos(t16);
        t52 = t19*om;
        t53 = t52*t21;
        t54 = t53*t30;
        t59 = ampmu*(3+t43*t38*t10);
        t64 = amplambda*(2+t3*t38*t45);
        t65 = 2*t59+t64;
        t66 = cos(t33);
        t74 = t20*t21;
        t77 = omstry*y;
        t78 = sin(t77);
        t80 = 1+t78/2;
        t82 = om*x+ph;
        t83 = sin(t82);
        t84 = t80*t83;
        t86 = om*(y-t14);
        t87 = cos(t86);
        t88 = t84*t87;
        t90 = t88*t53*t29;
        t92 = cos(t82);
        t96 = t21*t29;
        t99 = omstrz*z;
        t100 = sin(t99);
        t102 = 1+t100/2;
        t103 = t102*t83;
        t106 = om*(z-t14);
        t107 = cos(t106);
        t110 = t103*t26*t107*t52*t21;
        t120 = ampmu*t43;
        t121 = t120*t6;
        t122 = omm*t10;
        t123 = t36*t92;
        t124 = t123*t52;
        t125 = sin(t86);
        t127 = t125*t21*t29;
        t129 = t80*t17;
        t131 = cos(t25);
        t133 = t21*t131*t29;
        t135 = -t124*t127-t129*t52*t133;
        t138 = t123*t20;
        t142 = cos(t77);
        t143 = t142*omstry;
        t150 = t21*t26;
        t151 = t150*t29;
        t157 = t120*omm;
        t158 = t38*t45;
        t159 = sin(t106);
        t160 = t26*t159;
        t161 = t160*t21;
        t163 = t102*t17;
        t165 = cos(t28);
        t166 = t150*t165;
        t168 = -t124*t161-t163*t52*t166;
        t171 = t26*t107;
        t174 = cos(t99);
        t175 = t174*omstrz;
        forces[0] = t13*t17*t23*t30-t36*(-(-2*t37*t40+t47)*t36*t50*t54-
       t65*t66*omstrx*t50*t54/2+t65*t36*t17*t74*t30-t47*t90-t64*t80*t92*
        t20*t87*t96-t47*t110-t64*t102*t92*t20*t26*t107*t21)-t80*(t121*t122*t135+
         t59*(-t138*t87*t21*t29-t143*t17*t53*t131*t29/2+t129*t20*t151))
     -t102*(t157*t158*t168+t59*(-t138*t171*t21-t175*t17*t53*t26*t165/2+
     t163*t20*t151));
        t187 = t13*t83;
        t188 = t125*t20;
        t192 = t37*omm;
        t193 = t38*t10;
        t197 = t66*omstrx*t92;
        t198 = t52*t125;
        t203 = t36*t83*t20;
        t212 = t6*omm;
        t216 = amplambda*t3;
        t223 = t87*t52*t96;
        t236 = t36*t50;
        t238 = t236*t52*t151;
        t242 = t236*t20;
        t246 = t20*t107*t21;
        t252 = t84*t131;
        t256 = t103*t125;
        t259 = -t252*t52*t159*t21-t256*t53*t165;
        forces[1] = t187*t188*t22*t29-t36*(-t192*t193*t135+t59*(
        -t197*t198*t96/2+t203*t127-t80*t50*t20*t133))-t80*(-(2*t120*t212*t10+t216*
     t212*t45)*t80*t83*t223-t65*t142*omstry*t83*t223/2+t65*t80*t83*t188
     *t96+t216*t6*omm*t45*(-t238-t110)+t64*(-t242*t133-t103*t131*t246))
     -t102*(t157*t158*t259+t59*(-t252*t246-t175*t83*t198*t21*t165/2+t256*t74*t29));
        t279 = t159*t21;
        t303 = t74*t165;
        t315 = t171*t53;
        forces[2] = t187*t160*t23-t36*(-t192*t193*t168+t59*(
     -t197*t52*t26*t279/2+t203*t161-t102*t50*t20*t166))-t80*(t121*t122*t259+t59*
       (-t143*t83*t131*t52*t279/2+t84*t26*t20*t159*t21-t103*t87*t303))-t102*
     (-(2*t120*t46-t216*t40)*t102*t83*t315-t65*t174*omstrz*t83*t315/2+
       t65*t102*t83*t160*t74-t216*omm*t193*(-t238-t90)+t64*(-t242*t166-t88
        *t303));

	size_t ind = base+i+ni*j+nij*k;
	fo[ind] = forces[0];
	fo[ind+nijk] = forces[1];
	fo[ind+2*nijk] = forces[2];
      // 	    }
      // 	 }
      // }
			       });
}

//-----------------------------------------------------------------------
void EW::forcingfortcsg_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			float_sw4* __restrict__ xx, float_sw4* __restrict yy, 
			float_sw4* __restrict__ zz, float_sw4 omstrx, float_sw4 omstry,
			float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
     using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
#else
     using LOCAL_POL = DEFAULT_LOOP3;
#endif
     RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t10,t100,t102,t103,t104,t105,t106,t113,t114,t115,t116,t117,t119,t121,t122,t123,t124,t126,t129,t13,t132,t133,t14,t142,t143,t144,t145,t147,t148,t150,t152,t157,t158,t16,t167,t168,t17,t172,t173,t177,t180,t181,t184,t185,t19,t191,t195,t2,t20,t21,t214,t215,t220,t226,t229,t23,t231,t24,t250,t253,t26,t265,t27,t28,t3,t31,t32,t34,t35,t36,t37,t38,t41,t43,t44,t45,t48,t5,t50,t51,t56,t6,t61,t62,t63,t72,t73,t75,t76,t78,t80,t81,t82,t84,t85,t87,t88,t9,t91,t92,t96,t97,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
// 	 for( int j=jfirst; j<=jlast; j++ )
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 x=xx[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 z=zz[ind];
        t2 = omm*x+phm;
        t3 = sin(t2);
        t5 = omm*y+phm;
        t6 = cos(t5);
        t9 = omm*z+phm;
        t10 = sin(t9);
        t13 = amprho*(2+t3*t6*t10);
        t14 = c*t;
        t16 = om*(x-t14);
        t17 = sin(t16);
        t19 = om*om;
        t20 = c*c;
        t21 = t19*t20;
        t23 = om*y+ph;
        t24 = sin(t23);
        t26 = om*z+ph;
        t27 = sin(t26);
        t28 = t24*t27;
        t31 = omstrx*x;
        t32 = sin(t31);
        t34 = 1+t32/2;
        t35 = ampmu*t3;
        t36 = sin(t5);
        t37 = omm*t36;
        t38 = t37*t10;
        t41 = cos(t2);
        t43 = cos(t9);
        t44 = t37*t43;
        t45 = amplambda*t41*t44;
        t48 = cos(t16);
        t50 = om*t24;
        t51 = t50*t27;
        t56 = ampmu*(3+t41*t36*t10);
        t61 = amplambda*(2+t3*t36*t43);
        t62 = 2*t56+t61;
        t63 = cos(t31);
        t72 = t19*t24;
        t73 = t72*t27;
        t75 = omstry*y;
        t76 = sin(t75);
        t78 = 1+t76/2;
        t80 = om*x+ph;
        t81 = sin(t80);
        t82 = t78*t81;
        t84 = om*(y-t14);
        t85 = cos(t84);
        t87 = t85*om*t27;
        t88 = t82*t87;
        t91 = cos(t80);
        t92 = t91*t19;
        t96 = omstrz*z;
        t97 = sin(t96);
        t99 = 1+t97/2;
        t100 = t99*t81;
        t102 = om*(z-t14);
        t103 = cos(t102);
        t104 = t24*t103;
        t105 = t104*om;
        t106 = t100*t105;
        t113 = ampmu*t41;
        t114 = t113*t6;
        t115 = omm*t10;
        t116 = t34*t91;
        t117 = sin(t84);
        t119 = om*t117*t27;
        t121 = t78*t17;
        t122 = cos(t23);
        t123 = t122*om;
        t124 = t123*t27;
        t126 = t116*t119+t121*t124;
        t129 = t19*t85;
        t132 = cos(t75);
        t133 = t132*omstry;
        t142 = t113*omm;
        t143 = t36*t43;
        t144 = sin(t102);
        t145 = t50*t144;
        t147 = t99*t17;
        t148 = cos(t26);
        t150 = t24*t148*om;
        t152 = t116*t145+t147*t150;
        t157 = cos(t96);
        t158 = t157*omstrz;
        forces[0] = -t13*t17*t21*t28-t34*((-2*t35*t38+t45)*t34*t48*t51+
	   t62*t63*omstrx*t48*om*t28/2-t62*t34*t17*t73+t45*t88+t61*t78*t92*t85
     *t27+t45*t106+t61*t99*t92*t104)-t78*(t114*t115*t126+t56*(t116*t129
     *t27+t133*t17*t124/2-t121*t73))-t99*(t142*t143*t152+t56*(t116*t72*
      t103+t158*t17*t150/2-t147*t73));
        t167 = t13*t81;
        t168 = t117*t19;
        t172 = t35*omm;
        t173 = t36*t10;
        t177 = t63*omstrx*t91;
        t180 = t34*t81;
        t181 = t168*t27;
        t184 = t19*t122;
        t185 = t184*t27;
        t191 = t6*omm;
        t195 = amplambda*t3;
        t214 = t34*t48;
        t215 = t214*t51;
        t220 = t184*t103;
        t226 = t123*t144;
        t229 = t117*t148*om;
        t231 = t82*t226+t100*t229;
        forces[1] = -t167*t168*t20*t27-t34*(-t172*t173*t126+t56*(t177*t119/2-
	   t180*t181+t78*t48*t185))-t78*((2*t113*t191*t10+t195*t191*t43)
     *t78*t81*t87+t62*t132*omstry*t81*t85*om*t27/2-t62*t78*t81*t181+
       t195*t6*omm*t43*(t215+t106)+t61*(t214*t185+t100*t220))-t99*(t142*t143
      *t231+t56*(t82*t220+t158*t81*t229/2-t100*t181));
        t250 = t72*t144;
        t253 = t72*t148;
        t265 = t129*t148;
        forces[2] = -t167*t24*t144*t21-t34*(-t172*t173*t152+t56*(t177*t145/2-
           t180*t250+t99*t48*t253))-t78*(t114*t115*t231+t56*(t133*t81*
          t226/2-t82*t250+t100*t265))-t99*((2*t113*t44-t195*t38)*t99*t81*t105+
     t62*t157*omstrz*t81*t24*t103*om/2-t62*t99*t81*t250-t195*omm*t173*(
     t215+t88)+t61*(t214*t253+t82*t265));

	fo[ind] = forces[0];
	fo[ind+nijk] = forces[1];
	fo[ind+2*nijk] = forces[2];
			       }); SYNC_STREAM;
   //}
}

//-----------------------------------------------------------------------
void EW::forcingttfortcsg_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			  float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			  float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			  float_sw4* __restrict__ xx, float_sw4* __restrict yy, 
			  float_sw4* __restrict__ zz, float_sw4 omstrx, float_sw4 omstry,
			  float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
 #else
 using LOCAL_POL = DEFAULT_LOOP3;
#endif
     RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t10,t100,t102,t103,t106,t107,t110,t120,t121,t122,t123,t124,t125,t127,t129,t13,t131,t133,t135,t138,t14,t142,t143,t150,t151,t157,t158,t159,t16,t160,t161,t163,t165,t166,t168,t17,t171,t174,t175,t187,t188,t19,t192,t193,t197,t198,t2,t20,t203,t21,t212,t216,t22,t223,t23,t236,t238,t242,t246,t25,t252,t256,t259,t26,t279,t28,t29,t3,t30,t303,t315,t33,t34,t36,t37,t38,t39,t40,t43,t45,t46,t47,t5,t50,t52,t53,t54,t59,t6,t64,t65,t66,t74,t77,t78,t80,t82,t83,t84,t86,t87,t88,t9,t90,t92,t96,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
// 	 for( int j=jfirst; j<=jlast; j++ )
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 x=xx[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 z=zz[ind];
        t2 = omm*x+phm;
        t3 = sin(t2);
        t5 = omm*y+phm;
        t6 = cos(t5);
        t9 = omm*z+phm;
        t10 = sin(t9);
        t13 = amprho*(2+t3*t6*t10);
        t14 = c*t;
        t16 = om*(x-t14);
        t17 = sin(t16);
        t19 = om*om;
        t20 = t19*t19;
        t21 = c*c;
        t22 = t21*t21;
        t23 = t20*t22;
        t25 = om*y+ph;
        t26 = sin(t25);
        t28 = om*z+ph;
        t29 = sin(t28);
        t30 = t26*t29;
        t33 = omstrx*x;
        t34 = sin(t33);
        t36 = 1+t34/2;
        t37 = ampmu*t3;
        t38 = sin(t5);
        t39 = omm*t38;
        t40 = t39*t10;
        t43 = cos(t2);
        t45 = cos(t9);
        t46 = t39*t45;
        t47 = amplambda*t43*t46;
        t50 = cos(t16);
        t52 = t19*om;
        t53 = t52*t21;
        t54 = t53*t30;
        t59 = ampmu*(3+t43*t38*t10);
        t64 = amplambda*(2+t3*t38*t45);
        t65 = 2*t59+t64;
        t66 = cos(t33);
        t74 = t20*t21;
        t77 = omstry*y;
        t78 = sin(t77);
        t80 = 1+t78/2;
        t82 = om*x+ph;
        t83 = sin(t82);
        t84 = t80*t83;
        t86 = om*(y-t14);
        t87 = cos(t86);
        t88 = t84*t87;
        t90 = t88*t53*t29;
        t92 = cos(t82);
        t96 = t21*t29;
        t99 = omstrz*z;
        t100 = sin(t99);
        t102 = 1+t100/2;
        t103 = t102*t83;
        t106 = om*(z-t14);
        t107 = cos(t106);
        t110 = t103*t26*t107*t52*t21;
        t120 = ampmu*t43;
        t121 = t120*t6;
        t122 = omm*t10;
        t123 = t36*t92;
        t124 = t123*t52;
        t125 = sin(t86);
        t127 = t125*t21*t29;
        t129 = t80*t17;
        t131 = cos(t25);
        t133 = t21*t131*t29;
        t135 = -t124*t127-t129*t52*t133;
        t138 = t123*t20;
        t142 = cos(t77);
        t143 = t142*omstry;
        t150 = t21*t26;
        t151 = t150*t29;
        t157 = t120*omm;
        t158 = t38*t45;
        t159 = sin(t106);
        t160 = t26*t159;
        t161 = t160*t21;
        t163 = t102*t17;
        t165 = cos(t28);
        t166 = t150*t165;
        t168 = -t124*t161-t163*t52*t166;
        t171 = t26*t107;
        t174 = cos(t99);
        t175 = t174*omstrz;
        forces[0] = t13*t17*t23*t30-t36*(-(-2*t37*t40+t47)*t36*t50*t54-
       t65*t66*omstrx*t50*t54/2+t65*t36*t17*t74*t30-t47*t90-t64*t80*t92*
        t20*t87*t96-t47*t110-t64*t102*t92*t20*t26*t107*t21)-t80*(t121*t122*t135+
         t59*(-t138*t87*t21*t29-t143*t17*t53*t131*t29/2+t129*t20*t151))
     -t102*(t157*t158*t168+t59*(-t138*t171*t21-t175*t17*t53*t26*t165/2+
     t163*t20*t151));
        t187 = t13*t83;
        t188 = t125*t20;
        t192 = t37*omm;
        t193 = t38*t10;
        t197 = t66*omstrx*t92;
        t198 = t52*t125;
        t203 = t36*t83*t20;
        t212 = t6*omm;
        t216 = amplambda*t3;
        t223 = t87*t52*t96;
        t236 = t36*t50;
        t238 = t236*t52*t151;
        t242 = t236*t20;
        t246 = t20*t107*t21;
        t252 = t84*t131;
        t256 = t103*t125;
        t259 = -t252*t52*t159*t21-t256*t53*t165;
        forces[1] = t187*t188*t22*t29-t36*(-t192*t193*t135+t59*(
        -t197*t198*t96/2+t203*t127-t80*t50*t20*t133))-t80*(-(2*t120*t212*t10+t216*
     t212*t45)*t80*t83*t223-t65*t142*omstry*t83*t223/2+t65*t80*t83*t188
     *t96+t216*t6*omm*t45*(-t238-t110)+t64*(-t242*t133-t103*t131*t246))
     -t102*(t157*t158*t259+t59*(-t252*t246-t175*t83*t198*t21*t165/2+t256*t74*t29));
        t279 = t159*t21;
        t303 = t74*t165;
        t315 = t171*t53;
        forces[2] = t187*t160*t23-t36*(-t192*t193*t168+t59*(
     -t197*t52*t26*t279/2+t203*t161-t102*t50*t20*t166))-t80*(t121*t122*t259+t59*
       (-t143*t83*t131*t52*t279/2+t84*t26*t20*t159*t21-t103*t87*t303))-t102*
     (-(2*t120*t46-t216*t40)*t102*t83*t315-t65*t174*omstrz*t83*t315/2+
       t65*t102*t83*t160*t74-t216*omm*t193*(-t238-t90)+t64*(-t242*t166-t88
        *t303));

	fo[ind] = forces[0];
	fo[ind+nijk] = forces[1];
	fo[ind+2*nijk] = forces[2];
			       });
   //}
}

//-----------------------------------------------------------------------
void EW::forcingfortsgatt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			  float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase, 
			  float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			  float_sw4 h, float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
			  float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
   //#pragma omp parallel
   RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<DEFAULT_LOOP3>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
				 float_sw4 forces[3],t1,t10,t100,t103,t104,t105,t106,t107,t11,t110,t111,t116,t117,t119,t120,t122,t125,t128,t132,t133,t134,t135,t137,t138,t14,t140,t141,t142,t144,t148,t15,t151,t152,t16,t164,t167,t17,t176,t18,t183,t194,t197,t198,t2,t20,t203,t210,t212,t214,t215,t217,t219,t24,t243,t26,t265,t27,t272,t28,t31,t32,t35,t36,t37,t39,t4,t40,t41,t42,t44,t50,t56,t57,t58,t6,t64,t67,t7,t73,t74,t76,t77,t78,t8,t81,t82,t85,t86,t87,t88,t92,t95,t97,t98;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
//       {
// 	 
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	 {
// 	    
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
				 float_sw4 y=(j-1)*h;
	       float_sw4 x=(i-1)*h;
	       float_sw4 z=(k-1)*h+zmin;
        t1 = omstrx*x;
        t2 = sin(t1);
        t4 = 1+t2/2;
        t6 = momega*x+mphase;
        t7 = sin(t6);
        t8 = ampmu*t7;
        t10 = momega*y+mphase;
        t11 = cos(t10);
        t14 = momega*z+mphase;
        t15 = sin(t14);
        t16 = momega*t11*t15;
        t17 = t8*t16;
        t18 = cos(t6);
        t20 = amplambda*t18*t16;
        t24 = c*t;
        t26 = omega*(x-t24);
        t27 = -t26-phase;
        t28 = sin(t27);
        t31 = omega*x+phase;
        t32 = sin(t31);
        t35 = -omega*(z-t24)-phase;
        t36 = cos(t35);
        t37 = t32*t36;
        t39 = cos(t27);
        t40 = cos(t31);
        t41 = t39*t40;
        t42 = omega*t36;
        t44 = t28*omega*t37+t41*t42;
        t50 = ampmu*(3./2.+t18*t11*t15/2);
        t56 = amplambda*(1./2.+t7*t11*t15/4);
        t57 = 2*t50+t56;
        t58 = cos(t1);
        t64 = omega*omega;
        t67 = t28*t64;
        t73 = omstry*y;
        t74 = sin(t73);
        t76 = 1+t74/2;
        t77 = sin(t26);
        t78 = t76*t77;
        t81 = -omega*(y-t24)-phase;
        t82 = sin(t81);
        t85 = omega*z+phase;
        t86 = cos(t85);
        t87 = t82*omega*t86;
        t88 = t78*t87;
        t92 = cos(t26);
        t95 = t92*t64*t82*t86;
        t97 = omstrz*z;
        t98 = sin(t97);
        t100 = 1+t98/2;
        t103 = omega*y+phase;
        t104 = cos(t103);
        t105 = t104*t36;
        t106 = t105*omega;
        t107 = t100*t40*t106;
        t110 = t56*t100;
        t111 = t32*t64;
        t116 = ampmu*t18;
        t117 = sin(t10);
        t119 = t117*momega*t15;
        t120 = t116*t119;
        t122 = cos(t81);
        t125 = t4*t92*omega*t122*t86;
        t128 = t50*t4;
        t132 = t116*t11;
        t133 = cos(t14);
        t134 = t133*momega;
        t135 = t4*t32;
        t137 = sin(t35);
        t138 = omega*t104*t137;
        t140 = t100*t39;
        t141 = t32*t137;
        t142 = t141*omega;
        t144 = t135*t138+t140*t142;
        t148 = t64*t104;
        t151 = cos(t97);
        t152 = t151*omstrz;
        forces[0] = t4*((-t17+t20/4)*t4*t44+t57*t58*omstrx*t44/2+t57*t4*
     (-2*t39*t64*t37+2*t67*t40*t36)+t20*t88/4+t56*t76*t95+t20*t107/4-
     t110*t111*t105)+t76*(-t120*t125/2+t128*t95)+t100*(t132*t134*t144/2+
       t50*(-t135*t148*t36+t152*t39*t142/2-t140*t37*t64));
        t164 = t58*omstrx;
        t167 = t122*t86;
        t176 = amplambda*t7;
        t183 = cos(t73);
        t194 = t122*t64*t86;
        t197 = momega*t15;
        t198 = t4*t44;
        t203 = sin(t103);
        t210 = t76*t40;
        t212 = t203*omega*t137;
        t214 = t100*t77;
        t215 = sin(t85);
        t217 = t122*t215*omega;
        t219 = t210*t212-t214*t217;
        forces[1] = t4*(-t17*t125/2+t50*t164*t92*omega*t167/2-t128*t77*
         t64*t167)+t76*((-t120-t176*t119/4)*t76*t77*t87+t57*t183*omstry*t77*
     t82*omega*t86/2-t57*t76*t77*t194-t176*t117*t197*(t198+t107)/4-t110
     *t40*t203*t64*t36)+t100*(t132*t134*t219/2+t50*(-t210*t203*t64*t36-
     t152*t77*t217/2-t214*t194));
        t243 = t148*t137;
        t265 = t82*t64*t215;
        t272 = t11*t133*momega;
        forces[2] = t4*(-t8*momega*t11*t15*t144/2+t50*(t164*t32*t138/2+
       t4*t40*t243+t100*t28*t111*t137+t140*t40*t64*t137))+t76*(-t116*t117*
     t197*t219/2+t50*(t183*omstry*t40*t212/2+t210*t243-t214*t265))+t100
     *((t116*t272+t176*t272/4)*t100*t40*t106+t57*t151*omstrz*t40*t104*
       t42/2+t57*t100*t40*t243+t176*t11*t134*(t198+t88)/4+t56*(t4*(t67*
       t141+t41*t64*t137)-t78*t265));

	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind]      += forces[0];
	       fo[ind+nijk] += forces[1];
	       fo[ind+2*nijk] += forces[2];
      // 	    }
      // 	 }
      // }
			       }); SYNC_STREAM;
}

//-----------------------------------------------------------------------
void EW::forcingttfortsgatt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			    int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			    float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase, 
			    float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			    float_sw4 h, float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
			    float_sw4 omstrz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
 #else
 using LOCAL_POL = DEFAULT_LOOP3;
 #endif
     RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t1,t10,t100,t101,t104,t105,t107,t108,t109,t11,t110,t114,t115,t116,t120,t121,t124,t125,t127,t128,t130,t131,t135,t14,t141,t146,t147,t149,t15,t150,t151,t153,t155,t157,t16,t161,t164,t167,t168,t169,t17,t170,t172,t173,t175,t179,t18,t183,t188,t191,t192,t194,t2,t20,t202,t204,t205,t207,t214,t215,t224,t228,t230,t234,t237,t238,t24,t245,t253,t254,t26,t260,t266,t267,t27,t271,t273,t274,t277,t279,t28,t282,t29,t297,t30,t300,t31,t314,t32,t34,t35,t352,t359,t36,t39,t4,t40,t41,t43,t44,t45,t46,t48,t49,t50,t52,t55,t6,t61,t67,t68,t69,t7,t75,t76,t78,t8,t82,t85,t86,t88,t89,t90,t94,t95,t96,t98,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
//       {
// 	 float_sw4 z=(k-1)*h+zmin;
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	 {
// 	    float_sw4 y=(j-1)*h;
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       float_sw4 x=(i-1)*h;
	       float_sw4 y=(j-1)*h;
	       float_sw4 z=(k-1)*h+zmin;
        t1 = omstrx*x;
        t2 = sin(t1);
        t4 = 1+t2/2;
        t6 = momega*x+mphase;
        t7 = sin(t6);
        t8 = ampmu*t7;
        t10 = momega*y+mphase;
        t11 = cos(t10);
        t14 = momega*z+mphase;
        t15 = sin(t14);
        t16 = momega*t11*t15;
        t17 = t8*t16;
        t18 = cos(t6);
        t20 = amplambda*t18*t16;
        t24 = c*t;
        t26 = omega*(x-t24);
        t27 = -t26-phase;
        t28 = sin(t27);
        t29 = omega*omega;
        t30 = t29*omega;
        t31 = t28*t30;
        t32 = c*c;
        t34 = omega*x+phase;
        t35 = sin(t34);
        t36 = t32*t35;
        t39 = -omega*(z-t24)-phase;
        t40 = cos(t39);
        t41 = t36*t40;
        t43 = cos(t27);
        t44 = t43*t30;
        t45 = sin(t39);
        t46 = t36*t45;
        t48 = cos(t34);
        t49 = t32*t48;
        t50 = t49*t40;
        t52 = t49*t45;
        t55 = -2*t31*t41-2*t44*t46-2*t44*t50+2*t31*t52;
        t61 = ampmu*(3./2.+t18*t11*t15/2);
        t67 = amplambda*(1./2.+t7*t11*t15/4);
        t68 = 2*t61+t67;
        t69 = cos(t1);
        t75 = t29*t29;
        t76 = t43*t75;
        t78 = t28*t75;
        t82 = t76*t41-t78*t46-t78*t50-t76*t52;
        t85 = omstry*y;
        t86 = sin(t85);
        t88 = 1+t86/2;
        t89 = sin(t26);
        t90 = t88*t89;
        t94 = -omega*(y-t24)-phase;
        t95 = sin(t94);
        t96 = t32*t95;
        t98 = omega*z+phase;
        t99 = cos(t98);
        t100 = t96*t99;
        t101 = t90*t30*t100;
        t104 = cos(t26);
        t105 = t88*t104;
        t107 = cos(t94);
        t108 = t32*t107;
        t109 = t108*t99;
        t110 = t105*t30*t109;
        t114 = t75*t32;
        t115 = t95*t99;
        t116 = t114*t115;
        t120 = t107*t99;
        t121 = t114*t120;
        t124 = omstrz*z;
        t125 = sin(t124);
        t127 = 1+t125/2;
        t128 = t127*t48;
        t130 = omega*y+phase;
        t131 = cos(t130);
        t135 = t128*t131*t40*t30*t32;
        t141 = t40*t32;
        t146 = ampmu*t18;
        t147 = sin(t10);
        t149 = t147*momega*t15;
        t150 = t146*t149;
        t151 = t4*t104;
        t153 = t151*t30*t109;
        t155 = t4*t89;
        t157 = t155*t30*t100;
        t161 = 2*t61*t151*t116;
        t164 = 2*t61*t155*t121;
        t167 = t146*t11;
        t168 = cos(t14);
        t169 = t168*momega;
        t170 = t4*t35;
        t172 = t131*t45;
        t173 = t172*t32;
        t175 = t127*t43;
        t179 = t127*t28;
        t183 = -t170*t30*t173-2*t175*t30*t46-2*t179*t30*t41;
        t188 = t131*t40;
        t191 = cos(t124);
        t192 = t191*omstrz;
        t194 = t30*t32;
        t202 = t175*t75;
        t204 = 2*t202*t41;
        t205 = t179*t75;
        t207 = 2*t205*t46;
        forces[0] = t4*((-t17+t20/4)*t4*t55+t68*t69*omstrx*t55/2+4*t68*t4*t82-
        t20*t101/2-t20*t110/2-2*t67*t105*t116+2*t67*t90*t121-t20*t135/4+
       t67*t127*t35*t75*t131*t141)+t88*(t150*t153+t150*t157-t161+t164
      )+t127*(t167*t169*t183/2+t61*(t170*t75*t188*t32-t192*t43*t194*t35*
     t45-t192*t28*t194*t35*t40+t204-t207));
        t214 = t69*omstrx;
        t215 = t61*t214;
        t224 = amplambda*t7;
        t228 = (-t150-t224*t149/4)*t88;
        t230 = t194*t115;
        t234 = t194*t120;
        t237 = cos(t85);
        t238 = t68*t237;
        t245 = t68*t88;
        t253 = momega*t15;
        t254 = t4*t55;
        t260 = sin(t130);
        t266 = t88*t48;
        t267 = t266*t260;
        t271 = t127*t89;
        t273 = sin(t98);
        t274 = t108*t273;
        t277 = t127*t104;
        t279 = t96*t273;
        t282 = -t267*t30*t45*t32+2*t271*t30*t274-2*t277*t30*t279;
        t297 = t271*t75;
        t300 = t277*t75;
        forces[1] = t4*(t17*t153+t17*t157-t215*t104*t30*t109-t215*t89*t30*t100-
        t161+t164)+t88*(-2*t228*t89*t230-2*t228*t104*t234-t238*omstry*t89*t230-
       t238*omstry*t104*t234+2*t245*t89*t121-2*t245*t104*t116
     -t224*t147*t253*(t254-t135)/4+t67*t128*t260*t75*t141)+t127*(t167*t169*t282/2+
        t61*(t267*t75*t40*t32+t192*t89*t194*t107*t273-t192*t104
        *t194*t95*t273+2*t297*t109-2*t300*t100));
        t314 = t45*t32;
        t352 = t11*t168*momega;
        t359 = t188*t194;
        forces[2] = t4*(-t8*momega*t11*t15*t183/2+t61*(-t214*t35*t30*t131*t314/2-
     t4*t48*t75*t173+t204-t207-2*t202*t52-2*t205*t50))+t88*(-t146*t147*t253*t282/2+
         t61*(-t237*omstry*t48*t260*t30*t314/2-t266*t131*t75*t45*t32+2*t297*t279+
       2*t300*t274))+t127*(-(t146*t352+t224*t352/4)*t127*t48*t359-t68*t191*omstrz*t48*t359/2-
		t68*t127*t48*t172*t114+t224*t11*t169*(t254-2*t101-2*t110)/4+
       t67*(2*t4*t82+2*t90*t75*t279+2*t105*t75*t274));

	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind]      += forces[0];
	       fo[ind+nijk] += forces[1];
	       fo[ind+2*nijk] += forces[2];

      // 	    }
      // 	 }
      // }
			       });
}

//-----------------------------------------------------------------------
void EW::forcingfortsgattc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
		 	   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			   float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			   float_sw4* __restrict__ zz, float_sw4 omstrx, float_sw4 omstry,
			   float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
     using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
#else
     using LOCAL_POL = DEFAULT_LOOP3;
#endif     
       RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t1,t10,t100,t103,t104,t105,t106,t107,t11,t110,t111,t116,t117,t119,t120,t122,t125,t128,t132,t133,t134,t135,t137,t138,t14,t140,t141,t142,t144,t148,t15,t151,t152,t16,t164,t167,t17,t176,t18,t183,t194,t197,t198,t2,t20,t203,t210,t212,t214,t215,t217,t219,t24,t243,t26,t265,t27,t272,t28,t31,t32,t35,t36,t37,t39,t4,t40,t41,t42,t44,t50,t56,t57,t58,t6,t64,t67,t7,t73,t74,t76,t77,t78,t8,t81,t82,t85,t86,t87,t88,t92,t95,t97,t98;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
// 	 for( int j=jfirst; j<=jlast; j++ )
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];

        t1 = omstrx*x;
        t2 = sin(t1);
        t4 = 1+t2/2;
        t6 = momega*x+mphase;
        t7 = sin(t6);
        t8 = ampmu*t7;
        t10 = momega*y+mphase;
        t11 = cos(t10);
        t14 = momega*z+mphase;
        t15 = sin(t14);
        t16 = momega*t11*t15;
        t17 = t8*t16;
        t18 = cos(t6);
        t20 = amplambda*t18*t16;
        t24 = c*t;
        t26 = omega*(x-t24);
        t27 = -t26-phase;
        t28 = sin(t27);
        t31 = omega*x+phase;
        t32 = sin(t31);
        t35 = -omega*(z-t24)-phase;
        t36 = cos(t35);
        t37 = t32*t36;
        t39 = cos(t27);
        t40 = cos(t31);
        t41 = t39*t40;
        t42 = omega*t36;
        t44 = t28*omega*t37+t41*t42;
        t50 = ampmu*(3./2.+t18*t11*t15/2);
        t56 = amplambda*(1./2.+t7*t11*t15/4);
        t57 = 2*t50+t56;
        t58 = cos(t1);
        t64 = omega*omega;
        t67 = t28*t64;
        t73 = omstry*y;
        t74 = sin(t73);
        t76 = 1+t74/2;
        t77 = sin(t26);
        t78 = t76*t77;
        t81 = -omega*(y-t24)-phase;
        t82 = sin(t81);
        t85 = omega*z+phase;
        t86 = cos(t85);
        t87 = t82*omega*t86;
        t88 = t78*t87;
        t92 = cos(t26);
        t95 = t92*t64*t82*t86;
        t97 = omstrz*z;
        t98 = sin(t97);
        t100 = 1+t98/2;
        t103 = omega*y+phase;
        t104 = cos(t103);
        t105 = t104*t36;
        t106 = t105*omega;
        t107 = t100*t40*t106;
        t110 = t56*t100;
        t111 = t32*t64;
        t116 = ampmu*t18;
        t117 = sin(t10);
        t119 = t117*momega*t15;
        t120 = t116*t119;
        t122 = cos(t81);
        t125 = t4*t92*omega*t122*t86;
        t128 = t50*t4;
        t132 = t116*t11;
        t133 = cos(t14);
        t134 = t133*momega;
        t135 = t4*t32;
        t137 = sin(t35);
        t138 = omega*t104*t137;
        t140 = t100*t39;
        t141 = t32*t137;
        t142 = t141*omega;
        t144 = t135*t138+t140*t142;
        t148 = t64*t104;
        t151 = cos(t97);
        t152 = t151*omstrz;
        forces[0] = t4*((-t17+t20/4)*t4*t44+t57*t58*omstrx*t44/2+t57*t4*
     (-2*t39*t64*t37+2*t67*t40*t36)+t20*t88/4+t56*t76*t95+t20*t107/4-
     t110*t111*t105)+t76*(-t120*t125/2+t128*t95)+t100*(t132*t134*t144/2+
       t50*(-t135*t148*t36+t152*t39*t142/2-t140*t37*t64));
        t164 = t58*omstrx;
        t167 = t122*t86;
        t176 = amplambda*t7;
        t183 = cos(t73);
        t194 = t122*t64*t86;
        t197 = momega*t15;
        t198 = t4*t44;
        t203 = sin(t103);
        t210 = t76*t40;
        t212 = t203*omega*t137;
        t214 = t100*t77;
        t215 = sin(t85);
        t217 = t122*t215*omega;
        t219 = t210*t212-t214*t217;
        forces[1] = t4*(-t17*t125/2+t50*t164*t92*omega*t167/2-t128*t77*
         t64*t167)+t76*((-t120-t176*t119/4)*t76*t77*t87+t57*t183*omstry*t77*
     t82*omega*t86/2-t57*t76*t77*t194-t176*t117*t197*(t198+t107)/4-t110
     *t40*t203*t64*t36)+t100*(t132*t134*t219/2+t50*(-t210*t203*t64*t36-
     t152*t77*t217/2-t214*t194));
        t243 = t148*t137;
        t265 = t82*t64*t215;
        t272 = t11*t133*momega;
        forces[2] = t4*(-t8*momega*t11*t15*t144/2+t50*(t164*t32*t138/2+
       t4*t40*t243+t100*t28*t111*t137+t140*t40*t64*t137))+t76*(-t116*t117*
     t197*t219/2+t50*(t183*omstry*t40*t212/2+t210*t243-t214*t265))+t100
     *((t116*t272+t176*t272/4)*t100*t40*t106+t57*t151*omstrz*t40*t104*
       t42/2+t57*t100*t40*t243+t176*t11*t134*(t198+t88)/4+t56*(t4*(t67*
       t141+t41*t64*t137)-t78*t265));

	       fo[ind]      += forces[0];
	       fo[ind+nijk] += forces[1];
	       fo[ind+2*nijk] += forces[2];
			       });
   //}
}

//-----------------------------------------------------------------------
void EW::forcingttfortsgattc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
		 	   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			   float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			   float_sw4* __restrict__ zz, float_sw4 omstrx, float_sw4 omstry,
			   float_sw4 omstrz )
{
  SW4_MARK_FUNCTION;
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
// #pragma omp parallel
//    {
#ifdef ENABLE_CUDA
using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
 #else
 using LOCAL_POL = DEFAULT_LOOP3;
 #endif
     RAJA::RangeSegment k_range(kfirst,klast+1);
   RAJA::RangeSegment j_range(jfirst,jlast+1);
   RAJA::RangeSegment i_range(ifirst,ilast+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
      float_sw4 forces[3],t1,t10,t100,t101,t104,t105,t107,t108,t109,t11,t110,t114,t115,t116,t120,t121,t124,t125,t127,t128,t130,t131,t135,t14,t141,t146,t147,t149,t15,t150,t151,t153,t155,t157,t16,t161,t164,t167,t168,t169,t17,t170,t172,t173,t175,t179,t18,t183,t188,t191,t192,t194,t2,t20,t202,t204,t205,t207,t214,t215,t224,t228,t230,t234,t237,t238,t24,t245,t253,t254,t26,t260,t266,t267,t27,t271,t273,t274,t277,t279,t28,t282,t29,t297,t30,t300,t31,t314,t32,t34,t35,t352,t359,t36,t39,t4,t40,t41,t43,t44,t45,t46,t48,t49,t50,t52,t55,t6,t61,t67,t68,t69,t7,t75,t76,t78,t8,t82,t85,t86,t88,t89,t90,t94,t95,t96,t98,t99;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
// 	 for( int j=jfirst; j<=jlast; j++ )
// #pragma ivdep
// #pragma simd
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];

        t1 = omstrx*x;
        t2 = sin(t1);
        t4 = 1+t2/2;
        t6 = momega*x+mphase;
        t7 = sin(t6);
        t8 = ampmu*t7;
        t10 = momega*y+mphase;
        t11 = cos(t10);
        t14 = momega*z+mphase;
        t15 = sin(t14);
        t16 = momega*t11*t15;
        t17 = t8*t16;
        t18 = cos(t6);
        t20 = amplambda*t18*t16;
        t24 = c*t;
        t26 = omega*(x-t24);
        t27 = -t26-phase;
        t28 = sin(t27);
        t29 = omega*omega;
        t30 = t29*omega;
        t31 = t28*t30;
        t32 = c*c;
        t34 = omega*x+phase;
        t35 = sin(t34);
        t36 = t32*t35;
        t39 = -omega*(z-t24)-phase;
        t40 = cos(t39);
        t41 = t36*t40;
        t43 = cos(t27);
        t44 = t43*t30;
        t45 = sin(t39);
        t46 = t36*t45;
        t48 = cos(t34);
        t49 = t32*t48;
        t50 = t49*t40;
        t52 = t49*t45;
        t55 = -2*t31*t41-2*t44*t46-2*t44*t50+2*t31*t52;
        t61 = ampmu*(3./2.+t18*t11*t15/2);
        t67 = amplambda*(1./2.+t7*t11*t15/4);
        t68 = 2*t61+t67;
        t69 = cos(t1);
        t75 = t29*t29;
        t76 = t43*t75;
        t78 = t28*t75;
        t82 = t76*t41-t78*t46-t78*t50-t76*t52;
        t85 = omstry*y;
        t86 = sin(t85);
        t88 = 1+t86/2;
        t89 = sin(t26);
        t90 = t88*t89;
        t94 = -omega*(y-t24)-phase;
        t95 = sin(t94);
        t96 = t32*t95;
        t98 = omega*z+phase;
        t99 = cos(t98);
        t100 = t96*t99;
        t101 = t90*t30*t100;
        t104 = cos(t26);
        t105 = t88*t104;
        t107 = cos(t94);
        t108 = t32*t107;
        t109 = t108*t99;
        t110 = t105*t30*t109;
        t114 = t75*t32;
        t115 = t95*t99;
        t116 = t114*t115;
        t120 = t107*t99;
        t121 = t114*t120;
        t124 = omstrz*z;
        t125 = sin(t124);
        t127 = 1+t125/2;
        t128 = t127*t48;
        t130 = omega*y+phase;
        t131 = cos(t130);
        t135 = t128*t131*t40*t30*t32;
        t141 = t40*t32;
        t146 = ampmu*t18;
        t147 = sin(t10);
        t149 = t147*momega*t15;
        t150 = t146*t149;
        t151 = t4*t104;
        t153 = t151*t30*t109;
        t155 = t4*t89;
        t157 = t155*t30*t100;
        t161 = 2*t61*t151*t116;
        t164 = 2*t61*t155*t121;
        t167 = t146*t11;
        t168 = cos(t14);
        t169 = t168*momega;
        t170 = t4*t35;
        t172 = t131*t45;
        t173 = t172*t32;
        t175 = t127*t43;
        t179 = t127*t28;
        t183 = -t170*t30*t173-2*t175*t30*t46-2*t179*t30*t41;
        t188 = t131*t40;
        t191 = cos(t124);
        t192 = t191*omstrz;
        t194 = t30*t32;
        t202 = t175*t75;
        t204 = 2*t202*t41;
        t205 = t179*t75;
        t207 = 2*t205*t46;
        forces[0] = t4*((-t17+t20/4)*t4*t55+t68*t69*omstrx*t55/2+4*t68*t4*t82-
        t20*t101/2-t20*t110/2-2*t67*t105*t116+2*t67*t90*t121-t20*t135/4+
       t67*t127*t35*t75*t131*t141)+t88*(t150*t153+t150*t157-t161+t164
      )+t127*(t167*t169*t183/2+t61*(t170*t75*t188*t32-t192*t43*t194*t35*
     t45-t192*t28*t194*t35*t40+t204-t207));
        t214 = t69*omstrx;
        t215 = t61*t214;
        t224 = amplambda*t7;
        t228 = (-t150-t224*t149/4)*t88;
        t230 = t194*t115;
        t234 = t194*t120;
        t237 = cos(t85);
        t238 = t68*t237;
        t245 = t68*t88;
        t253 = momega*t15;
        t254 = t4*t55;
        t260 = sin(t130);
        t266 = t88*t48;
        t267 = t266*t260;
        t271 = t127*t89;
        t273 = sin(t98);
        t274 = t108*t273;
        t277 = t127*t104;
        t279 = t96*t273;
        t282 = -t267*t30*t45*t32+2*t271*t30*t274-2*t277*t30*t279;
        t297 = t271*t75;
        t300 = t277*t75;
        forces[1] = t4*(t17*t153+t17*t157-t215*t104*t30*t109-t215*t89*t30*t100-
        t161+t164)+t88*(-2*t228*t89*t230-2*t228*t104*t234-t238*omstry*t89*t230-
       t238*omstry*t104*t234+2*t245*t89*t121-2*t245*t104*t116
     -t224*t147*t253*(t254-t135)/4+t67*t128*t260*t75*t141)+t127*(t167*t169*t282/2+
        t61*(t267*t75*t40*t32+t192*t89*t194*t107*t273-t192*t104
        *t194*t95*t273+2*t297*t109-2*t300*t100));
        t314 = t45*t32;
        t352 = t11*t168*momega;
        t359 = t188*t194;
        forces[2] = t4*(-t8*momega*t11*t15*t183/2+t61*(-t214*t35*t30*t131*t314/2-
     t4*t48*t75*t173+t204-t207-2*t202*t52-2*t205*t50))+t88*(-t146*t147*t253*t282/2+
         t61*(-t237*omstry*t48*t260*t30*t314/2-t266*t131*t75*t45*t32+2*t297*t279+
       2*t300*t274))+t127*(-(t146*t352+t224*t352/4)*t127*t48*t359-t68*t191*omstrz*t48*t359/2-
		t68*t127*t48*t172*t114+t224*t11*t169*(t254-2*t101-2*t110)/4+
       t67*(2*t4*t82+2*t90*t75*t279+2*t105*t75*t274));

	       fo[ind]      += forces[0];
	       fo[ind+nijk] += forces[1];
	       fo[ind+2*nijk] += forces[2];
			       });
   //}
}
