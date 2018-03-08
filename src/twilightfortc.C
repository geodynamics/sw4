#include "sw4.h"
#include "EW.h"
#include "Mspace.h"
#include "policies.h"
void EW::twilightfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ u, float_sw4 t, float_sw4 om,
			  float_sw4 cv, float_sw4 ph, float_sw4 h, float_sw4 zmin )
// new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
// u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
// v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
// w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const ssize_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
   {
      float_sw4 z = (k-1)*h + zmin;
      for( int j=jfirst; j<=jlast; j++ )
      {
         float_sw4 y = (j-1)*h;
#pragma ivdep
#pragma simd
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    float_sw4 x = (i-1)*h;
	    size_t ind = base + i + ni*j + nij*k;
            u[ind] =       sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
            u[ind+nijk]   = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
            u[ind+2*nijk] = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::twilightfortwind_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			     int klast, float_sw4* __restrict__ u, float_sw4 t, float_sw4 om, 
			     float_sw4 cv, float_sw4 ph, float_sw4 h, float_sw4 zmin,
			     int i1, int i2, int j1, int j2, int k1, int k2 )
{
// new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
// u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
// v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
// w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));

   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
   for( int k=k1; k<=k2; k++ )
   {
      float_sw4 z = (k-1)*h + zmin;
      for( int j=j1; j<=j2; j++ )
      {
         float_sw4 y = (j-1)*h;
	 for( int i=i1; i<=i2; i++ )
	 {
	    float_sw4 x = (i-1)*h;
	    size_t ind = base + i + ni*j + nij*k;
            u[ind]        = sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
            u[ind+nijk]   = sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
            u[ind+2*nijk] = sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::twilightfortc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
			  int klast, float_sw4* __restrict__ u, float_sw4 t, float_sw4 om,
			  float_sw4 cv, float_sw4 ph, float_sw4* __restrict__ x,
			  float_sw4* __restrict__ y, float_sw4* __restrict__ z )
{
// new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
// u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
// v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
// w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
      for( int j=jfirst; j<=jlast; j++ )
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    size_t ind = base + i + ni*j + nij*k;
            u[ind]        = sin(om*(x[ind]-cv*t))*sin(om*y[ind]+ph)*
	       sin(om*z[ind]+ph);
            u[ind+nijk]   = sin(om*x[ind]+ph)*sin(om*(y[ind]-cv*t))*
	       sin(om*z[ind]+ph);
            u[ind+2*nijk] = sin(om*x[ind]+ph)*sin(om*y[ind]+ph)*
	       sin(om*(z[ind]-cv*t));
	 }
}


//-----------------------------------------------------------------------
void EW::twilightfortatt_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
			    int klast, float_sw4* __restrict__ alpha, float_sw4 t, float_sw4 om,
			    float_sw4 cv,float_sw4 ph,float_sw4 h,float_sw4 zmin )
{
// new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
// Attenuation memory variables
// att1 := cos(omega*(x-c*t)+phase)*sin(omega*x+phase)*cos(omega*(z-c*t)+phase);
// att2 := sin(omega*(x-c*t))*cos(omega*(y-c*t)+phase)*cos(omega*z+phase);
// att3 := cos(omega*x+phase)*cos(omega*y+phase)*sin(omega*(z-c*t)+phase);

   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
   {
      float_sw4 z = (k-1)*h + zmin;
      for( int j=jfirst; j<=jlast; j++ )
      {
         float_sw4 y = (j-1)*h;
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    float_sw4 x = (i-1)*h;
	    size_t ind = base + i + ni*j + nij*k;
            alpha[ind] = cos(om*(x-cv*t)+ph)*sin(om*x+ph)*
	       cos(om*(z-cv*t)+ph);
            alpha[ind+nijk] = sin(om*(x-cv*t))*cos(om*(y-cv*t)+ph)*
	       cos(om*z+ph);
            alpha[ind+2*nijk] = cos(om*x+ph)*cos(om*y+ph)*
	       sin(om*(z-cv*t)+ph);
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::twilightfortattc_ci(int ifirst, int ilast, int jfirst, int jlast, int kfirst,
			     int klast, float_sw4* __restrict__ alpha, float_sw4 t,
			     float_sw4 om, float_sw4 cv, float_sw4 ph, 
			     float_sw4* __restrict__ x, float_sw4* __restrict__ y,
			     float_sw4* __restrict__ z )
{
// new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
// u      := sin(om*(x-cv*t))*sin(om*y+ph)*sin(om*z+ph);
// v      := sin(om*x+ph)*sin(om*(y-cv*t))*sin(om*z+ph);
// w      := sin(om*x+ph)*sin(om*y+ph)*sin(om*(z-cv*t));
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
      for( int j=jfirst; j<=jlast; j++ )
#pragma ivdep
#pragma simd
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    size_t ind = base + i + ni*j + nij*k;
            alpha[ind]        = cos(om*(x[ind]-cv*t)+ph)*
	       sin(om*x[ind]+ph)*cos(om*(z[ind]-cv*t)+ph);
            alpha[ind+nijk]   = sin(om*(x[ind]-cv*t))*
	       cos(om*(y[ind]-cv*t)+ph)*cos(om*z[ind]+ph);
            alpha[ind+2*nijk] = cos(om*x[ind]+ph)*cos(om*y[ind]+ph)*
	       sin(om*(z[ind]-cv*t)+ph);
	 }
}

//-----------------------------------------------------------------------
void EW::exactrhsfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			  float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			  float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			  float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t100,t101,t102,t103,t105,t106,t108,t11,t114,t115,t118,t119,t12,t120,t122,t125,t129,t135,t143,t15,t150,t152,t159,t16,t166,t168,t17,t173,t175,t18,t2,t21,t23,t24,t27,t28,t3,t31,t32,t38,t4,t43,t44,t45,t47,t48,t51,t53,t54,t55,t57,t58,t6,t60,t63,t64,t69,t7,t70,t71,t72,t77,t78,t79,t8,t80,t81,t82,t83,t85,t87,t89,t92,t95,t97;
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
   {
      float_sw4 z=(k-1)*h+zmin;
      for( int j=jfirst; j<=jlast; j++ )
      {
	 float_sw4 y=(j-1)*h;
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    float_sw4 x=(i-1)*h;
	    t2 = omm*x+phm;
	    t3 = sin(t2);
	    t4 = ampmu*t3;
	    t6 = omm*y+phm;
	    t7 = sin(t6);
	    t8 = omm*t7;
	    t10 = omm*z+phm;
	    t11 = sin(t10);
	    t12 = t8*t11;
	    t15 = cos(t2);
	    t16 = amplambda*t15;
	    t17 = cos(t10);
	    t18 = t8*t17;
	    t21 = c*t;
	    t23 = om*(x-t21);
	    t24 = cos(t23);
	    t27 = om*y+ph;
	    t28 = sin(t27);
	    t31 = om*z+ph;
	    t32 = sin(t31);
	    t38 = ampmu*(3+t15*t7*t11);
	    t43 = amplambda*(2+t3*t7*t17);
	    t44 = 2*t38+t43;
	    t45 = sin(t23);
	    t47 = om*om;
	    t48 = t47*t28;
	    t51 = t16*t8;
	    t53 = om*x+ph;
	    t54 = sin(t53);
	    t55 = t17*t54;
	    t57 = om*(y-t21);
	    t58 = cos(t57);
	    t60 = t58*om*t32;
	    t63 = cos(t53);
	    t64 = t43*t63;
	    t69 = om*(z-t21);
	    t70 = cos(t69);
	    t71 = t28*t70;
	    t72 = t71*om;
	    t77 = ampmu*t15;
	    t78 = cos(t6);
	    t79 = t77*t78;
	    t80 = omm*t11;
	    t81 = t63*om;
	    t82 = sin(t57);
	    t83 = t82*t32;
	    t85 = cos(t27);
	    t87 = om*t32;
	    t89 = t81*t83+t45*t85*t87;
	    t92 = t63*t47;
	    t95 = t45*t28;
	    t97 = t95*t47*t32;
	    t100 = t77*omm;
	    t101 = t7*t17;
	    t102 = sin(t69);
	    t103 = t28*t102;
	    t105 = cos(t31);
	    t106 = t105*om;
	    t108 = t81*t103+t95*t106;
	    forces[0] = (-2*t4*t12+t16*t18)*t24*om*t28*t32-t44*t45*t48*t32+
	       t51*t55*t60+t64*t47*t58*t32+t51*t55*t72+t64*t48*t70+t79*t80*t89+
	       t38*(t92*t58*t32-t97)+t100*t101*t108+t38*(t92*t71-t97);
	    t114 = t4*omm;
	    t115 = t7*t11;
	    t118 = t54*t47;
	    t119 = t118*t83;
	    t120 = t24*t47;
	    t122 = t120*t85*t32;
	    t125 = t78*omm;
	    t129 = amplambda*t3;
	    t135 = t44*t54;
	    t143 = t24*om*t28*t32;
	    t150 = t54*t85;
	    t152 = t150*t47*t70;
	    t159 = t150*om*t102+t54*t82*t106;
	    forces[1] = -t114*t115*t89+t38*(-t119+t122)+(2*t77*t125*t11+t129*t125*t17)*
	       t54*t60-t135*t82*t47*t32+t129*t78*omm*t17*(t143+t54*t28*t70*om)+
	       t43*(t122+t152)+t100*t101*t159+t38*(t152-t119);
	    t166 = t118*t103;
	    t168 = t120*t28*t105;
	    t173 = t54*t58;
	    t175 = t173*t47*t105;
	    forces[2] = -t114*t115*t108+t38*(-t166+t168)+t79*t80*t159+t38*(-t166+t175)+
	       (2*t77*t18-t129*t12)*t54*t72-t135*t103*t47
	       -t129*omm*t115*(t143+t173*t87)+t43*(t168+t175);
	    size_t ind = base+i+ni*j+nij*k;
	    fo[ind]       = forces[0];
	    fo[ind+nijk]  = forces[1];
	    fo[ind+2*nijk]= forces[2];
	 }
      }
   }
   }
}

//-----------------------------------------------------------------------
void EW::exactrhsfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			   float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t100,t101,t102,t103,t105,t106,t108,t11,t114,t115,t118,t119,t12,t120,t122,t125,t129,t135,t143,t15,t150,t152,t159,t16,t166,t168,t17,t173,t175,t18,t2,t21,t23,t24,t27,t28,t3,t31,t32,t38,t4,t43,t44,t45,t47,t48,t51,t53,t54,t55,t57,t58,t6,t60,t63,t64,t69,t7,t70,t71,t72,t77,t78,t79,t8,t80,t81,t82,t83,t85,t87,t89,t92,t95,t97;
#pragma omp for
   for( int k=kfirst; k<=klast; k++ )
      for( int j=jfirst; j<=jlast; j++ )
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    size_t ind = base+i+ni*j+nij*k;
	    float_sw4 z=zz[ind];
	    float_sw4 y=yy[ind];
	    float_sw4 x=xx[ind];
	    t2 = omm*x+phm;
	    t3 = sin(t2);
	    t4 = ampmu*t3;
	    t6 = omm*y+phm;
	    t7 = sin(t6);
	    t8 = omm*t7;
	    t10 = omm*z+phm;
	    t11 = sin(t10);
	    t12 = t8*t11;
	    t15 = cos(t2);
	    t16 = amplambda*t15;
	    t17 = cos(t10);
	    t18 = t8*t17;
	    t21 = c*t;
	    t23 = om*(x-t21);
	    t24 = cos(t23);
	    t27 = om*y+ph;
	    t28 = sin(t27);
	    t31 = om*z+ph;
	    t32 = sin(t31);
	    t38 = ampmu*(3+t15*t7*t11);
	    t43 = amplambda*(2+t3*t7*t17);
	    t44 = 2*t38+t43;
	    t45 = sin(t23);
	    t47 = om*om;
	    t48 = t47*t28;
	    t51 = t16*t8;
	    t53 = om*x+ph;
	    t54 = sin(t53);
	    t55 = t17*t54;
	    t57 = om*(y-t21);
	    t58 = cos(t57);
	    t60 = t58*om*t32;
	    t63 = cos(t53);
	    t64 = t43*t63;
	    t69 = om*(z-t21);
	    t70 = cos(t69);
	    t71 = t28*t70;
	    t72 = t71*om;
	    t77 = ampmu*t15;
	    t78 = cos(t6);
	    t79 = t77*t78;
	    t80 = omm*t11;
	    t81 = t63*om;
	    t82 = sin(t57);
	    t83 = t82*t32;
	    t85 = cos(t27);
	    t87 = om*t32;
	    t89 = t81*t83+t45*t85*t87;
	    t92 = t63*t47;
	    t95 = t45*t28;
	    t97 = t95*t47*t32;
	    t100 = t77*omm;
	    t101 = t7*t17;
	    t102 = sin(t69);
	    t103 = t28*t102;
	    t105 = cos(t31);
	    t106 = t105*om;
	    t108 = t81*t103+t95*t106;
	    forces[0] = (-2*t4*t12+t16*t18)*t24*om*t28*t32-t44*t45*t48*t32+
	       t51*t55*t60+t64*t47*t58*t32+t51*t55*t72+t64*t48*t70+t79*t80*t89+
	       t38*(t92*t58*t32-t97)+t100*t101*t108+t38*(t92*t71-t97);
	    t114 = t4*omm;
	    t115 = t7*t11;
	    t118 = t54*t47;
	    t119 = t118*t83;
	    t120 = t24*t47;
	    t122 = t120*t85*t32;
	    t125 = t78*omm;
	    t129 = amplambda*t3;
	    t135 = t44*t54;
	    t143 = t24*om*t28*t32;
	    t150 = t54*t85;
	    t152 = t150*t47*t70;
	    t159 = t150*om*t102+t54*t82*t106;
	    forces[1] = -t114*t115*t89+t38*(-t119+t122)+(2*t77*t125*t11+t129*t125*t17)*
	       t54*t60-t135*t82*t47*t32+t129*t78*omm*t17*(t143+t54*t28*t70*om)+
	       t43*(t122+t152)+t100*t101*t159+t38*(t152-t119);
	    t166 = t118*t103;
	    t168 = t120*t28*t105;
	    t173 = t54*t58;
	    t175 = t173*t47*t105;
	    forces[2] = -t114*t115*t108+t38*(-t166+t168)+t79*t80*t159+t38*(-t166+t175)+
	       (2*t77*t18-t129*t12)*t54*t72-t135*t103*t47
	       -t129*omm*t115*(t143+t173*t87)+t43*(t168+t175);
	    fo[ind]       = forces[0];
	    fo[ind+nijk]  = forces[1];
	    fo[ind+2*nijk]= forces[2];
	 }
   }
}

//-----------------------------------------------------------------------
void EW::exactaccfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ utt, float_sw4 t, float_sw4 om,
			  float_sw4 c, float_sw4 ph, float_sw4 h, float_sw4 zmin )
{ 
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 acc[3], t1,t4,t5,t7,t10,t14,t19,t22,t30;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       t1 = c*t;
	       t4 = sin(om*(x-t1));
	       t5 = om*om;
	       t7 = c*c;
	       t10 = sin(om*y+ph);
	       t14 = sin(om*z+ph);
	       acc[0] = -t4*t5*t7*t10*t14;
	       t19 = sin(om*x+ph);
	       t22 = sin(om*(y-t1));
	       acc[1] = -t19*t22*t5*t7*t14;
	       t30 = sin(om*(z-t1));
	       acc[2] = -t19*t10*t30*t5*t7;
	       size_t ind=base+i+ni*j+nij*k;
	       utt[ind] = acc[0];
	       utt[ind+nijk] = acc[1];
	       utt[ind+2*nijk] = acc[2];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::exactaccfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ utt, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4* __restrict__ x,
			   float_sw4* __restrict__ y, float_sw4* __restrict__ z )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 acc[3], t1,t4,t5,t7,t10,t14,t19,t22,t30;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       t1 = c*t;
	       t4 = sin(om*(x[ind]-t1));
	       t5 = om*om;
	       t7 = c*c;
	       t10 = sin(om*y[ind]+ph);
	       t14 = sin(om*z[ind]+ph);
	       acc[0] = -t4*t5*t7*t10*t14;
	       t19 = sin(om*x[ind]+ph);
	       t22 = sin(om*(y[ind]-t1));
	       acc[1] = -t19*t22*t5*t7*t14;
	       t30 = sin(om*(z[ind]-t1));
	       acc[2] = -t19*t10*t30*t5*t7;
	       utt[ind] = acc[0];
	       utt[ind+nijk] = acc[1];
	       utt[ind+2*nijk] = acc[2];
	    }
   }
}

//-----------------------------------------------------------------------
void EW::forcingfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			 int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			 float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			 float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			 float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t102,t105,t107,t110,t111,t112,t113,t115,t116,t118,t124,t125,t129,t13,t130,t133,t134,t135,t137,t14,t140,t144,t150,t156,t16,t163,t165,t17,t172,t181,t183,t188,t19,t190,t2,t20,t21,t23,t24,t26,t27,t28,t3,t31,t32,t33,t34,t37,t38,t39,t40,t43,t5,t51,t56,t57,t59,t6,t62,t64,t65,t66,t68,t69,t71,t74,t75,t80,t81,t82,t83,t88,t89,t9,t90,t91,t92,t93,t95,t97,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
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
        t31 = ampmu*t3;
        t32 = sin(t5);
        t33 = omm*t32;
        t34 = t33*t10;
        t37 = cos(t2);
        t38 = amplambda*t37;
        t39 = cos(t9);
        t40 = t33*t39;
        t43 = cos(t16);
        t51 = ampmu*(3+t37*t32*t10);
        t56 = amplambda*(2+t3*t32*t39);
        t57 = 2*t51+t56;
        t59 = t19*t24;
        t62 = t38*t33;
        t64 = om*x+ph;
        t65 = sin(t64);
        t66 = t39*t65;
        t68 = om*(y-t14);
        t69 = cos(t68);
        t71 = t69*om*t27;
        t74 = cos(t64);
        t75 = t56*t74;
        t80 = om*(z-t14);
        t81 = cos(t80);
        t82 = t24*t81;
        t83 = t82*om;
        t88 = ampmu*t37;
        t89 = t88*t6;
        t90 = omm*t10;
        t91 = t74*om;
        t92 = sin(t68);
        t93 = t92*t27;
        t95 = cos(t23);
        t97 = om*t27;
        t99 = t91*t93+t17*t95*t97;
        t102 = t74*t19;
        t105 = t17*t24;
        t107 = t105*t19*t27;
        t110 = t88*omm;
        t111 = t32*t39;
        t112 = sin(t80);
        t113 = t24*t112;
        t115 = cos(t26);
        t116 = t115*om;
        t118 = t91*t113+t105*t116;
        forces[0] = -t13*t17*t21*t28-(-2*t31*t34+t38*t40)*t43*om*t24*t27
	   +t57*t17*t59*t27-t62*t66*t71-t75*t19*t69*t27-t62*t66*t83-t75*t59*t81-
	   t89*t90*t99-t51*(t102*t69*t27-t107)-t110*t111*t118-t51*(t102*t82-t107);
        t124 = t13*t65;
        t125 = t92*t19;
        t129 = t31*omm;
        t130 = t32*t10;
        t133 = t65*t19;
        t134 = t133*t93;
        t135 = t43*t19;
        t137 = t135*t95*t27;
        t140 = t6*omm;
        t144 = amplambda*t3;
        t150 = t57*t65;
        t156 = t43*om*t28;
        t163 = t65*t95;
        t165 = t163*t19*t81;
        t172 = t163*om*t112+t65*t92*t116;
        forces[1] = -t124*t125*t20*t27+t129*t130*t99-t51*(-t134+t137)-(2
	   *t88*t140*t10+t144*t140*t39)*t65*t71+t150*t125*t27-t144*t6*omm*t39
	   *(t156+t65*t24*t81*om)-t56*(t137+t165)-t110*t111*t172-t51*(t165-t134);
        t181 = t133*t113;
        t183 = t135*t24*t115;
        t188 = t65*t69;
        t190 = t188*t19*t115;
        forces[2] = -t124*t113*t21+t129*t130*t118-t51*(-t181+t183)-
	   t89*t90*t172-t51*(-t181+t190)-(2*t88*t40-t144*t34)*t65*t83+t150*t113*t19+
	   t144*omm*t130*(t156+t188*t97)-t56*(t183+t190);
	size_t ind = base+i+ni*j+nij*k;
	fo[ind] = forces[0];
	fo[ind+nijk] = forces[1];
	fo[ind+2*nijk] = forces[2];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::forcingttfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
   ASSERT_MANAGED(fo);
   //#pragma omp parallel
   {
      // float_sw4 forces[3],t10,t100,t102,t103,t105,t107,t110,t115,t118,t119,t120,t121,t122,t124,t125,t127,t13,t135,t14,t140,t141,t144,t145,t146,t147,t150,t154,t16,t161,t163,t169,t17,t173,t176,t185,t19,t194,t195,t2,t20,t201,t21,t22,t23,t25,t26,t28,t29,t3,t33,t34,t35,t36,t39,t41,t42,t43,t45,t47,t49,t5,t50,t55,t6,t60,t61,t66,t67,t69,t70,t71,t72,t73,t74,t76,t77,t84,t85,t87,t88,t9,t94,t95,t96,t97,t98;
      RAJA::RangeSegment k_range(kfirst,klast+1);
      RAJA::RangeSegment j_range(jfirst,jlast+1);
      RAJA::RangeSegment i_range(ifirst,ilast+1);
      RAJA::nested::forall(RHS4_EXEC_POL{},
			  RAJA::make_tuple(k_range, j_range,i_range),
			  [=]RAJA_DEVICE (int k,int j,int i) {
			    float_sw4 forces[3],t10,t100,t102,t103,t105,t107,t110,t115,t118,t119,t120,t121,t122,t124,t125,t127,t13,t135,t14,t140,t141,t144,t145,t146,t147,t150,t154,t16,t161,t163,t169,t17,t173,t176,t185,t19,t194,t195,t2,t20,t201,t21,t22,t23,t25,t26,t28,t29,t3,t33,t34,t35,t36,t39,t41,t42,t43,t45,t47,t49,t5,t50,t55,t6,t60,t61,t66,t67,t69,t70,t71,t72,t73,t74,t76,t77,t84,t85,t87,t88,t9,t94,t95,t96,t97,t98;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
//       {
// 	 float_sw4 z=(k-1)*h+zmin;
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	 {
// 	    float_sw4 y=(j-1)*h;
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
			    float_sw4 z=(k-1)*h+zmin;
			    float_sw4 y=(j-1)*h;
	       float_sw4 x=(i-1)*h;
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
	       t33 = ampmu*t3;
	       t34 = sin(t5);
	       t35 = omm*t34;
	       t36 = t35*t10;
	       t39 = cos(t2);
	       t41 = cos(t9);
	       t42 = t35*t41;
	       t43 = amplambda*t39*t42;
	       t45 = cos(t16);
	       t47 = t19*om;
	       t49 = t21*t26;
	       t50 = t49*t29;
	       t55 = ampmu*(3+t39*t34*t10);
	       t60 = amplambda*(2+t3*t34*t41);
	       t61 = 2*t55+t60;
	       t66 = om*x+ph;
	       t67 = sin(t66);
	       t69 = om*(y-t14);
	       t70 = cos(t69);
	       t71 = t67*t70;
	       t72 = t47*t21;
	       t73 = t72*t29;
	       t74 = t71*t73;
	       t76 = cos(t66);
	       t77 = t60*t76;
	       t84 = om*(z-t14);
	       t85 = cos(t84);
	       t87 = t85*t47*t21;
	       t88 = t67*t26*t87;
	       t94 = ampmu*t39;
	       t95 = t94*t6;
	       t96 = omm*t10;
	       t97 = t76*t47;
	       t98 = sin(t69);
	       t100 = t98*t21*t29;
	       t102 = t17*t47;
	       t103 = cos(t25);
	       t105 = t21*t103*t29;
	       t107 = -t97*t100-t102*t105;
	       t110 = t76*t20;
	       t115 = t17*t20*t50;
	       t118 = t94*omm;
	       t119 = t34*t41;
	       t120 = sin(t84);
	       t121 = t26*t120;
	       t122 = t121*t21;
	       t124 = cos(t28);
	       t125 = t49*t124;
	       t127 = -t97*t122-t102*t125;
	       forces[0] = t13*t17*t23*t26*t29+(-2*t33*t36+t43)*t45*t47*t50-t61
		  *t17*t20*t50+t43*t74+t77*t20*t70*t21*t29+t43*t88+t77*t20*t26*t85*t21-
		  t95*t96*t107-t55*(-t110*t70*t21*t29+t115)-t118*t119*t127-t55*(-
										t110*t26*t85*t21+t115);
	       t135 = t13*t67;
	       t140 = t33*omm;
	       t141 = t34*t10;
	       t144 = t67*t20;
	       t145 = t144*t100;
	       t146 = t45*t20;
	       t147 = t146*t105;
	       t150 = t6*omm;
	       t154 = amplambda*t3;
	       t161 = t61*t67;
	       t163 = t20*t21;
	       t169 = t45*t47*t50;
	       t173 = t67*t103;
	       t176 = t173*t20*t85*t21;
	       t185 = -t173*t47*t120*t21-t67*t98*t72*t124;
	       forces[1] = t135*t98*t20*t22*t29+t140*t141*t107-t55*(t145-t147)+
		  (2*t94*t150*t10+t154*t150*t41)*t67*t70*t73-t161*t98*t163*t29-t154*
		  t6*omm*t41*(-t169-t88)-t60*(-t147-t176)-t118*t119*t185-t55*(-t176+
									      t145);
	       t194 = t144*t122;
	       t195 = t146*t125;
	       t201 = t71*t163*t124;
	       forces[2] = t135*t121*t23+t140*t141*t127-t55*(t194-t195)-t95*t96
		  *t185-t55*(t194-t201)+(2*t94*t42-t154*t36)*t67*t26*t87-t161*t26*t120*t20*t21+
		  t154*omm*t141*(-t169-t74)-t60*(-t195-t201);

	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind] = forces[0];
	       fo[ind+nijk] = forces[1];
	       fo[ind+2*nijk] = forces[2];
	       //}
	       // }
			  }); // End of RAJA loop
   }
}

//-----------------------------------------------------------------------
void EW::forcingfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			  float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			  float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			  float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			  float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t102,t105,t107,t110,t111,t112,t113,t115,t116,t118,t124,t125,t129,t13,t130,t133,t134,t135,t137,t14,t140,t144,t150,t156,t16,t163,t165,t17,t172,t181,t183,t188,t19,t190,t2,t20,t21,t23,t24,t26,t27,t28,t3,t31,t32,t33,t34,t37,t38,t39,t40,t43,t5,t51,t56,t57,t59,t6,t62,t64,t65,t66,t68,t69,t71,t74,t75,t80,t81,t82,t83,t88,t89,t9,t90,t91,t92,t93,t95,t97,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];
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
	       t31 = ampmu*t3;
	       t32 = sin(t5);
	       t33 = omm*t32;
	       t34 = t33*t10;
	       t37 = cos(t2);
	       t38 = amplambda*t37;
	       t39 = cos(t9);
	       t40 = t33*t39;
	       t43 = cos(t16);
	       t51 = ampmu*(3+t37*t32*t10);
	       t56 = amplambda*(2+t3*t32*t39);
	       t57 = 2*t51+t56;
	       t59 = t19*t24;
	       t62 = t38*t33;
	       t64 = om*x+ph;
	       t65 = sin(t64);
	       t66 = t39*t65;
	       t68 = om*(y-t14);
	       t69 = cos(t68);
	       t71 = t69*om*t27;
	       t74 = cos(t64);
	       t75 = t56*t74;
	       t80 = om*(z-t14);
	       t81 = cos(t80);
	       t82 = t24*t81;
	       t83 = t82*om;
	       t88 = ampmu*t37;
	       t89 = t88*t6;
	       t90 = omm*t10;
	       t91 = t74*om;
	       t92 = sin(t68);
	       t93 = t92*t27;
	       t95 = cos(t23);
	       t97 = om*t27;
	       t99 = t91*t93+t17*t95*t97;
	       t102 = t74*t19;
	       t105 = t17*t24;
	       t107 = t105*t19*t27;
	       t110 = t88*omm;
	       t111 = t32*t39;
	       t112 = sin(t80);
	       t113 = t24*t112;
	       t115 = cos(t26);
	       t116 = t115*om;
	       t118 = t91*t113+t105*t116;
	       forces[0] = -t13*t17*t21*t28-(-2*t31*t34+t38*t40)*t43*om*t24*t27
		  +t57*t17*t59*t27-t62*t66*t71-t75*t19*t69*t27-t62*t66*t83-t75*t59*t81-
		  t89*t90*t99-t51*(t102*t69*t27-t107)-t110*t111*t118-t51*(t102*t82-t107);
	       t124 = t13*t65;
	       t125 = t92*t19;
	       t129 = t31*omm;
	       t130 = t32*t10;
	       t133 = t65*t19;
	       t134 = t133*t93;
	       t135 = t43*t19;
	       t137 = t135*t95*t27;
	       t140 = t6*omm;
	       t144 = amplambda*t3;
	       t150 = t57*t65;
	       t156 = t43*om*t28;
	       t163 = t65*t95;
	       t165 = t163*t19*t81;
	       t172 = t163*om*t112+t65*t92*t116;
	       forces[1] = -t124*t125*t20*t27+t129*t130*t99-t51*(-t134+t137)-(2
 	         *t88*t140*t10+t144*t140*t39)*t65*t71+t150*t125*t27-t144*t6*omm*t39
		  *(t156+t65*t24*t81*om)-t56*(t137+t165)-t110*t111*t172-t51*(t165-t134);
	       t181 = t133*t113;
	       t183 = t135*t24*t115;
	       t188 = t65*t69;
	       t190 = t188*t19*t115;
	       forces[2] = -t124*t113*t21+t129*t130*t118-t51*(-t181+t183)-
		  t89*t90*t172-t51*(-t181+t190)-(2*t88*t40-t144*t34)*t65*t83+t150*t113*t19+
		  t144*omm*t130*(t156+t188*t97)-t56*(t183+t190);

	       fo[ind] = forces[0];
	       fo[ind+nijk] = forces[1];
	       fo[ind+2*nijk] = forces[2];
	    }
   }
}

//-----------------------------------------------------------------------
void EW::forcingttfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			   float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			    float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			    float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
   //#pragma omp parallel
   {
     ASSERT_MANAGED(fo);
     RAJA::RangeSegment k_range(kfirst,klast+1);
     RAJA::RangeSegment j_range(jfirst,jlast+1);
     RAJA::RangeSegment i_range(ifirst,ilast+1);
     RAJA::nested::forall(RHS4_EXEC_POL{},
			  RAJA::make_tuple(k_range, j_range,i_range),
			  [=]RAJA_DEVICE (int k,int j,int i) {
      float_sw4 forces[3],t10,t100,t102,t103,t105,t107,t110,t115,t118,t119,t120,t121,t122,t124,t125,t127,t13,t135,t14,t140,t141,t144,t145,t146,t147,t150,t154,t16,t161,t163,t169,t17,t173,t176,t185,t19,t194,t195,t2,t20,t201,t21,t22,t23,t25,t26,t28,t29,t3,t33,t34,t35,t36,t39,t41,t42,t43,t45,t47,t49,t5,t50,t55,t6,t60,t61,t66,t67,t69,t70,t71,t72,t73,t74,t76,t77,t84,t85,t87,t88,t9,t94,t95,t96,t97,t98;
// #pragma omp for
//       for( int k=kfirst; k<=klast; k++ )
// 	 for( int j=jfirst; j<=jlast; j++ )
// 	    for( int i=ifirst; i<=ilast; i++ )
// 	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];
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
	       t33 = ampmu*t3;
	       t34 = sin(t5);
	       t35 = omm*t34;
	       t36 = t35*t10;
	       t39 = cos(t2);
	       t41 = cos(t9);
	       t42 = t35*t41;
	       t43 = amplambda*t39*t42;
	       t45 = cos(t16);
	       t47 = t19*om;
	       t49 = t21*t26;
	       t50 = t49*t29;
	       t55 = ampmu*(3+t39*t34*t10);
	       t60 = amplambda*(2+t3*t34*t41);
	       t61 = 2*t55+t60;
	       t66 = om*x+ph;
	       t67 = sin(t66);
	       t69 = om*(y-t14);
	       t70 = cos(t69);
	       t71 = t67*t70;
	       t72 = t47*t21;
	       t73 = t72*t29;
	       t74 = t71*t73;
	       t76 = cos(t66);
	       t77 = t60*t76;
	       t84 = om*(z-t14);
	       t85 = cos(t84);
	       t87 = t85*t47*t21;
	       t88 = t67*t26*t87;
	       t94 = ampmu*t39;
	       t95 = t94*t6;
	       t96 = omm*t10;
	       t97 = t76*t47;
	       t98 = sin(t69);
	       t100 = t98*t21*t29;
	       t102 = t17*t47;
	       t103 = cos(t25);
	       t105 = t21*t103*t29;
	       t107 = -t97*t100-t102*t105;
	       t110 = t76*t20;
	       t115 = t17*t20*t50;
	       t118 = t94*omm;
	       t119 = t34*t41;
	       t120 = sin(t84);
	       t121 = t26*t120;
	       t122 = t121*t21;
	       t124 = cos(t28);
	       t125 = t49*t124;
	       t127 = -t97*t122-t102*t125;
	       forces[0] = t13*t17*t23*t26*t29+(-2*t33*t36+t43)*t45*t47*t50-t61
		  *t17*t20*t50+t43*t74+t77*t20*t70*t21*t29+t43*t88+t77*t20*t26*t85*t21-
		  t95*t96*t107-t55*(-t110*t70*t21*t29+t115)-t118*t119*t127-t55*(-
										t110*t26*t85*t21+t115);
	       t135 = t13*t67;
	       t140 = t33*omm;
	       t141 = t34*t10;
	       t144 = t67*t20;
	       t145 = t144*t100;
	       t146 = t45*t20;
	       t147 = t146*t105;
	       t150 = t6*omm;
	       t154 = amplambda*t3;
	       t161 = t61*t67;
	       t163 = t20*t21;
	       t169 = t45*t47*t50;
	       t173 = t67*t103;
	       t176 = t173*t20*t85*t21;
	       t185 = -t173*t47*t120*t21-t67*t98*t72*t124;
	       forces[1] = t135*t98*t20*t22*t29+t140*t141*t107-t55*(t145-t147)+
		  (2*t94*t150*t10+t154*t150*t41)*t67*t70*t73-t161*t98*t163*t29-t154*
		  t6*omm*t41*(-t169-t88)-t60*(-t147-t176)-t118*t119*t185-t55*(-t176+
									      t145);
	       t194 = t144*t122;
	       t195 = t146*t125;
	       t201 = t71*t163*t124;
	       forces[2] = t135*t121*t23+t140*t141*t127-t55*(t194-t195)-t95*t96
		  *t185-t55*(t194-t201)+(2*t94*t42-t154*t36)*t67*t26*t87-t161*t26*t120*t20*t21+
		  t154*omm*t141*(-t169-t74)-t60*(-t195-t201);

	       fo[ind] = forces[0];
	       fo[ind+nijk] = forces[1];
	       fo[ind+2*nijk] = forces[2];
			  }); // End of RAJA loop
   }
}

//-----------------------------------------------------------------------
void EW::exactmatfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			  int klast, float_sw4* __restrict__ rho, float_sw4* __restrict__ mu,
			  float_sw4* __restrict__ la, float_sw4 omm, float_sw4 phm,
			  float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			  float_sw4 h, float_sw4 zmin )
{
   // new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
   // rho    := amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
   // mu     := ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*sin(omm*z+phm) );
   // lambda := amplambda*(2 + sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) );
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
  //   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       size_t ind = base+i+ni*j+nij*k;
	       rho[ind] = amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*
                sin(omm*z+phm) );
	       mu[ind]  = ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*
				sin(omm*z+phm) );
	       la[ind]  = amplambda*(2 + 
		  sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) );
	    }
	 }
      }
}

//-----------------------------------------------------------------------
void EW::exactmatfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			   int klast, float_sw4* __restrict__ rho, float_sw4* __restrict__ mu,
			   float_sw4* __restrict__ la, float_sw4 omm, float_sw4 phm,
			   float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			   float_sw4* __restrict__ x, float_sw4* __restrict__ y,
			   float_sw4* __restrict__ z )
{
   // new 3d twilight functions (corresponding to subroutines fg, fgt and twrfsurz, see below
   // rho    := amprho*(2 + sin(omm*x+phm)*cos(omm*y+phm)*sin(omm*z+phm) );
   // mu     := ampmu*(3 + cos(omm*x+phm)*sin(omm*y+phm)*sin(omm*z+phm) );
   // lambda := amplambda*(2 + sin(omm*x+phm)*sin(omm*y+phm)*cos(omm*z+phm) );
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
  //   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
            rho[ind] = amprho*(2 + sin(omm*x[ind]+phm)*
				 cos(omm*y[ind]+phm)*sin(omm*z[ind]+phm) );
            mu[ind]  = ampmu*(3 + cos(omm*x[ind]+phm)*
				sin(omm*y[ind]+phm)*sin(omm*z[ind]+phm) );
            la[ind]  = amplambda*(2 + sin(omm*x[ind]+phm)*
				    sin(omm*y[ind]+phm)*cos(omm*z[ind]+phm) );
	       }
}

//-----------------------------------------------------------------------
void EW::exactrhsfortsg_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			    int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			    float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			    float_sw4 amprho, float_sw4 ampmu, float_sw4 ampla,
			    float_sw4 h, float_sw4 zmin, float_sw4 omstrx, float_sw4 omstry,
			    float_sw4 omstrz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t1,t10,t103,t104,t105,t106,t107,t108,t11,t110,t112,t113,t114,t115,t117,t12,t120,t123,t124,t133,t134,t135,t136,t138,t139,t14,t141,t143,t148,t149,t15,t158,t159,t16,t163,t166,t168,t171,t172,t178,t182,t19,t2,t201,t202,t207,t21,t213,t216,t218,t22,t23,t234,t237,t249,t26,t28,t29,t32,t33,t34,t36,t37,t38,t4,t43,t48,t49,t50,t59,t6,t61,t62,t63,t65,t66,t68,t7,t70,t71,t72,t74,t75,t77,t78,t8,t81,t82,t86,t87,t89,t90,t92,t93,t94,t95,t96;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       t1 = omstrx*x;
	       t2 = sin(t1);
	       t4 = 1+t2/2;
	       t6 = omm*x+phm;
	       t7 = sin(t6);
	       t8 = ampmu*t7;
	       t10 = omm*y+phm;
	       t11 = sin(t10);
	       t12 = omm*t11;
	       t14 = omm*z+phm;
	       t15 = sin(t14);
	       t16 = t12*t15;
	       t19 = cos(t6);
	       t21 = cos(t14);
	       t22 = t12*t21;
	       t23 = ampla*t19*t22;
	       t26 = c*t;
	       t28 = om*(x-t26);
	       t29 = cos(t28);
	       t32 = om*y+ph;
	       t33 = sin(t32);
	       t34 = om*t33;
	       t36 = om*z+ph;
	       t37 = sin(t36);
	       t38 = t34*t37;
	       t43 = ampmu*(3+t19*t11*t15);
	       t48 = ampla*(2+t7*t11*t21);
	       t49 = 2*t43+t48;
	       t50 = cos(t1);
	       t59 = sin(t28);
	       t61 = om*om;
	       t62 = t61*t33;
	       t63 = t62*t37;
	       t65 = omstry*y;
	       t66 = sin(t65);
	       t68 = 1+t66/2;
	       t70 = om*x+ph;
	       t71 = sin(t70);
	       t72 = t68*t71;
	       t74 = om*(y-t26);
	       t75 = cos(t74);
	       t77 = t75*om*t37;
	       t78 = t72*t77;
	       t81 = cos(t70);
	       t82 = t81*t61;
	       t86 = omstrz*z;
	       t87 = sin(t86);
	       t89 = 1+t87/2;
	       t90 = t89*t71;
	       t92 = om*(z-t26);
	       t93 = cos(t92);
	       t94 = t33*t93;
	       t95 = t94*om;
	       t96 = t90*t95;
	       t103 = ampmu*t19;
	       t104 = cos(t10);
	       t105 = t103*t104;
	       t106 = omm*t15;
	       t107 = t4*t81;
	       t108 = sin(t74);
	       t110 = om*t108*t37;
	       t112 = t68*t59;
	       t113 = cos(t32);
	       t114 = t113*om;
	       t115 = t114*t37;
	       t117 = t107*t110+t112*t115;
	       t120 = t61*t75;
	       t123 = cos(t65);
	       t124 = t123*omstry;
	       t133 = t103*omm;
	       t134 = t11*t21;
	       t135 = sin(t92);
	       t136 = t34*t135;
	       t138 = t89*t59;
	       t139 = cos(t36);
	       t141 = t33*t139*om;
	       t143 = t107*t136+t138*t141;
	       t148 = cos(t86);
	       t149 = t148*omstrz;
	       forces[0] = t4*((-2*t8*t16+t23)*t4*t29*t38+t49*t50*omstrx*t29*om
			       *t33*t37/2-t49*t4*t59*t63+t23*t78+t48*t68*t82*t75*t37+t23*t96+t48*
			       t89*t82*t94)+t68*(t105*t106*t117+t43*(t107*t120*t37+t124*t59*t115/
			      2-t112*t63))+t89*(t133*t134*t143+t43*(t107*t62*t93+t149*t59*t141/2-t138*t63));
	       t158 = t8*omm;
	       t159 = t11*t15;
	       t163 = t50*omstrx*t81;
	       t166 = t4*t71;
	       t168 = t61*t108*t37;
	       t171 = t61*t113;
	       t172 = t171*t37;
	       t178 = t104*omm;
	       t182 = ampla*t7;
	       t201 = t4*t29;
	       t202 = t201*t38;
	       t207 = t171*t93;
	       t213 = t114*t135;
	       t216 = t108*t139*om;
	       t218 = t72*t213+t90*t216;
	       forces[1] = t4*(-t158*t159*t117+t43*(t163*t110/2-t166*t168+t68*t29*t172))+
		  t68*((2*t103*t178*t15+t182*t178*t21)*t68*t71*t77+t49*t123*omstry*t71*t75*om*t37/2
		       -t49*t68*t71*t168+t182*t104*omm*t21*(t202+t96)+t48*(t201*t172+t90*t207))+
		  t89*(t133*t134*t218+t43*(t72*t207+t149*t71*t216/2-t90*t168));
	       t234 = t62*t135;
	       t237 = t62*t139;
	       t249 = t120*t139;
	       forces[2] = t4*(-t158*t159*t143+t43*(t163*t136/2-t166*t234+t89*t29*t237))+
		  t68*(t105*t106*t218+t43*(t124*t71*t213/2-t72*t234+t90*t249))+
		  t89*((2*t103*t22-t182*t16)*t89*t71*t95+t49*t148*omstrz*t71*t33*t93*om/2-
		       t49*t89*t71*t234-t182*omm*t159*(t202+t78)+t48*(t201*t237+t72*t249));
	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind] = forces[0];
	       fo[ind+nijk] = forces[1];
	       fo[ind+2*nijk] = forces[2];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::exactrhsfortsgc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			     int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 om,
			     float_sw4 c, float_sw4 ph, float_sw4 omm, float_sw4 phm, 
			     float_sw4 amprho, float_sw4 ampmu, float_sw4 ampla,
			     float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			     float_sw4* __restrict__ zz,
			     float_sw4 omstrx, float_sw4 omstry, float_sw4 omstrz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t1,t10,t103,t104,t105,t106,t107,t108,t11,t110,t112,t113,t114,t115,t117,t12,t120,t123,t124,t133,t134,t135,t136,t138,t139,t14,t141,t143,t148,t149,t15,t158,t159,t16,t163,t166,t168,t171,t172,t178,t182,t19,t2,t201,t202,t207,t21,t213,t216,t218,t22,t23,t234,t237,t249,t26,t28,t29,t32,t33,t34,t36,t37,t38,t4,t43,t48,t49,t50,t59,t6,t61,t62,t63,t65,t66,t68,t7,t70,t71,t72,t74,t75,t77,t78,t8,t81,t82,t86,t87,t89,t90,t92,t93,t94,t95,t96;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];
	       t1 = omstrx*x;
	       t2 = sin(t1);
	       t4 = 1+t2/2;
	       t6 = omm*x+phm;
	       t7 = sin(t6);
	       t8 = ampmu*t7;
	       t10 = omm*y+phm;
	       t11 = sin(t10);
	       t12 = omm*t11;
	       t14 = omm*z+phm;
	       t15 = sin(t14);
	       t16 = t12*t15;
	       t19 = cos(t6);
	       t21 = cos(t14);
	       t22 = t12*t21;
	       t23 = ampla*t19*t22;
	       t26 = c*t;
	       t28 = om*(x-t26);
	       t29 = cos(t28);
	       t32 = om*y+ph;
	       t33 = sin(t32);
	       t34 = om*t33;
	       t36 = om*z+ph;
	       t37 = sin(t36);
	       t38 = t34*t37;
	       t43 = ampmu*(3+t19*t11*t15);
	       t48 = ampla*(2+t7*t11*t21);
	       t49 = 2*t43+t48;
	       t50 = cos(t1);
	       t59 = sin(t28);
	       t61 = om*om;
	       t62 = t61*t33;
	       t63 = t62*t37;
	       t65 = omstry*y;
	       t66 = sin(t65);
	       t68 = 1+t66/2;
	       t70 = om*x+ph;
	       t71 = sin(t70);
	       t72 = t68*t71;
	       t74 = om*(y-t26);
	       t75 = cos(t74);
	       t77 = t75*om*t37;
	       t78 = t72*t77;
	       t81 = cos(t70);
	       t82 = t81*t61;
	       t86 = omstrz*z;
	       t87 = sin(t86);
	       t89 = 1+t87/2;
	       t90 = t89*t71;
	       t92 = om*(z-t26);
	       t93 = cos(t92);
	       t94 = t33*t93;
	       t95 = t94*om;
	       t96 = t90*t95;
	       t103 = ampmu*t19;
	       t104 = cos(t10);
	       t105 = t103*t104;
	       t106 = omm*t15;
	       t107 = t4*t81;
	       t108 = sin(t74);
	       t110 = om*t108*t37;
	       t112 = t68*t59;
	       t113 = cos(t32);
	       t114 = t113*om;
	       t115 = t114*t37;
	       t117 = t107*t110+t112*t115;
	       t120 = t61*t75;
	       t123 = cos(t65);
	       t124 = t123*omstry;
	       t133 = t103*omm;
	       t134 = t11*t21;
	       t135 = sin(t92);
	       t136 = t34*t135;
	       t138 = t89*t59;
	       t139 = cos(t36);
	       t141 = t33*t139*om;
	       t143 = t107*t136+t138*t141;
	       t148 = cos(t86);
	       t149 = t148*omstrz;
	       forces[0] = t4*((-2*t8*t16+t23)*t4*t29*t38+t49*t50*omstrx*t29*om
			       *t33*t37/2-t49*t4*t59*t63+t23*t78+t48*t68*t82*t75*t37+t23*t96+t48*
			       t89*t82*t94)+t68*(t105*t106*t117+t43*(t107*t120*t37+t124*t59*t115/
			      2-t112*t63))+t89*(t133*t134*t143+t43*(t107*t62*t93+t149*t59*t141/2-t138*t63));
	       t158 = t8*omm;
	       t159 = t11*t15;
	       t163 = t50*omstrx*t81;
	       t166 = t4*t71;
	       t168 = t61*t108*t37;
	       t171 = t61*t113;
	       t172 = t171*t37;
	       t178 = t104*omm;
	       t182 = ampla*t7;
	       t201 = t4*t29;
	       t202 = t201*t38;
	       t207 = t171*t93;
	       t213 = t114*t135;
	       t216 = t108*t139*om;
	       t218 = t72*t213+t90*t216;
	       forces[1] = t4*(-t158*t159*t117+t43*(t163*t110/2-t166*t168+t68*t29*t172))+
		  t68*((2*t103*t178*t15+t182*t178*t21)*t68*t71*t77+t49*t123*omstry*t71*t75*om*t37/2
		       -t49*t68*t71*t168+t182*t104*omm*t21*(t202+t96)+t48*(t201*t172+t90*t207))+
		  t89*(t133*t134*t218+t43*(t72*t207+t149*t71*t216/2-t90*t168));
	       t234 = t62*t135;
	       t237 = t62*t139;
	       t249 = t120*t139;
	       forces[2] = t4*(-t158*t159*t143+t43*(t163*t136/2-t166*t234+t89*t29*t237))+
		  t68*(t105*t106*t218+t43*(t124*t71*t213/2-t72*t234+t90*t249))+
		  t89*((2*t103*t22-t182*t16)*t89*t71*t95+t49*t148*omstrz*t71*t33*t93*om/2-
		       t49*t89*t71*t234-t182*omm*t159*(t202+t78)+t48*(t201*t237+t72*t249));
	       fo[ind] = forces[0];
	       fo[ind+nijk] = forces[1];
	       fo[ind+2*nijk] = forces[2];
	    }
   }
}

//-----------------------------------------------------------------------
void EW::exactmatfortatt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			     int klast, float_sw4* __restrict__ mu,
			     float_sw4* __restrict__ la, float_sw4 momega, float_sw4 mphase,
			     float_sw4 ampmu, float_sw4 amplambda,
			     float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
  //   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       size_t ind = base+i+ni*j+nij*k;
	       mu[ind] = ampmu*(1.5 + 0.5*cos(momega*x+mphase)*
				cos(momega*y+mphase)*sin(momega*z+mphase) );
	       la[ind] = amplambda*(0.5 + 0.25*sin(momega*x+mphase)
				    *cos(momega*y+mphase)*sin(momega*z+mphase) );
	    }
	 }
      }
}

//-----------------------------------------------------------------------
void EW::exactmatfortattc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			      int klast, float_sw4* __restrict__ mu,
			      float_sw4* __restrict__ la, float_sw4 momega, float_sw4 mphase,
			      float_sw4 ampmu, float_sw4 amplambda,
			      float_sw4* __restrict__ x, float_sw4* __restrict__ y,
			      float_sw4* __restrict__ z )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
  //   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       mu[ind] = ampmu*(1.5 + 0.5*cos(momega*x[ind]+mphase)*
			       cos(momega*y[ind]+mphase)*sin(momega*z[ind]+mphase) );
	       la[ind] = amplambda*(0.5 + 0.25*sin(momega*x[ind]+mphase)
				 *cos(momega*y[ind]+mphase)*sin(momega*z[ind]+mphase) );
	    }
}

//-----------------------------------------------------------------------
void EW::forcingfortatt_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			    int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			    float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase,
			    float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			    float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t100,t103,t108,t109,t11,t110,t112,t113,t116,t118,t12,t131,t133,t135,t14,t144,t15,t152,t156,t158,t159,t162,t178,t180,t182,t189,t19,t191,t195,t2,t21,t22,t23,t26,t27,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t4,t45,t51,t52,t53,t55,t56,t6,t62,t63,t67,t68,t7,t71,t72,t73,t77,t8,t80,t84,t85,t86,t87,t95,t96,t97;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       t2 = momega*x+mphase;
	       t3 = sin(t2);
	       t4 = ampmu*t3;
	       t6 = momega*y+mphase;
	       t7 = cos(t6);
	       t8 = momega*t7;
	       t10 = momega*z+mphase;
	       t11 = sin(t10);
	       t12 = t8*t11;
	       t14 = cos(t2);
	       t15 = amplambda*t14;
	       t19 = c*t;
	       t21 = omega*(x-t19);
	       t22 = -t21-phase;
	       t23 = sin(t22);
	       t26 = omega*x+phase;
	       t27 = sin(t26);
	       t30 = -omega*(z-t19)-phase;
	       t31 = cos(t30);
	       t32 = t27*t31;
	       t33 = t23*omega*t32;
	       t34 = cos(t22);
	       t35 = cos(t26);
	       t36 = t34*t35;
	       t37 = omega*t31;
	       t38 = t36*t37;
	       t45 = ampmu*(3.0/2.0+t14*t7*t11/2);
	       t51 = amplambda*(1.0/2.0+t3*t7*t11/4);
	       t52 = 2*t45+t51;
	       t53 = omega*omega;
	       t55 = t34*t53*t32;
	       t56 = t23*t53;
	       t62 = t15*t8;
	       t63 = sin(t21);
	       t67 = -omega*(y-t19)-phase;
	       t68 = sin(t67);
	       t71 = omega*z+phase;
	       t72 = cos(t71);
	       t73 = t68*omega*t72;
	       t77 = cos(t21);
	       t80 = t53*t68*t72;
	       t84 = omega*y+phase;
	       t85 = cos(t84);
	       t86 = t85*t31;
	       t87 = t86*omega;
	       t95 = ampmu*t14;
	       t96 = sin(t6);
	       t97 = t96*momega;
	       t100 = cos(t67);
	       t103 = t11*t77*omega*t100*t72;
	       t108 = t95*t7;
	       t109 = cos(t10);
	       t110 = t109*momega;
	       t112 = sin(t30);
	       t113 = t85*t112;
	       t116 = t112*omega;
	       t118 = t27*omega*t113+t34*t27*t116;
	       forces[0] = (-t4*t12+t15*t12/4)*(t33+t38)+t52*(-2*t55+2*t56*t35*
		      t31)+t62*t11*t63*t73/4+t51*t77*t80+t62*t11*t35*t87/4-t51*t27*t53*t85*t31
		  -t95*t97*t103/2+t45*t77*t80+t108*t110*t118/2+t45*(-t27*t53*t86-t55);
	       t131 = t53*t100*t72;
	       t133 = t97*t11;
	       t135 = amplambda*t3;
	       t144 = momega*t11;
	       t152 = sin(t84);
	       t156 = t35*t152;
	       t158 = t63*t100;
	       t159 = sin(t71);
	       t162 = t156*t116-t158*t159*omega;
	       forces[1] = -t4*t8*t103/2-t45*t63*t131+(-t95*t133-t135*t133/4)*t63*t73-
		  t52*t63*t131-t135*t96*t144*(t33+t38+t35*t85*t37)/4-t51*t35*
		  t152*t53*t31+t108*t110*t162/2+t45*(-t156*t53*t31-t158*t72*t53);
	       t178 = t35*t53*t113;
	       t180 = t56*t27*t112;
	       t182 = t36*t53*t112;
	       t189 = t63*t68;
	       t191 = t189*t53*t159;
	       t195 = t7*t109*momega;
	       forces[2] = -t4*momega*t7*t11*t118/2+t45*(t178+t180+t182)-
		  t95*t96*t144*t162/2+t45*(t178-t191)+(t95*t195+t135*t195/4)*t35*t87+
		  t52*t35*t113*t53+t135*t7*t110*(t33+t38+t189*omega*t72)/4+t51*(t180+t182-t191);
	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind]        += forces[0];
	       fo[ind+nijk]   += forces[1];
	       fo[ind+2*nijk] += forces[2];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::forcingttattfort_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			      int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			      float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase,
			      float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			      float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t104,t108,t109,t11,t112,t113,t118,t12,t121,t122,t124,t125,t13,t130,t133,t134,t135,t136,t139,t14,t141,t142,t143,t151,t152,t157,t16,t160,t178,t186,t19,t190,t194,t195,t198,t2,t201,t208,t21,t211,t22,t222,t23,t231,t233,t238,t24,t25,t26,t27,t29,t3,t30,t31,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t56,t6,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t74,t75,t78,t79,t80,t82,t83,t84,t85,t88,t89,t90,t91,t92,t93,t97,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       t2 = momega*x+mphase;
	       t3 = sin(t2);
	       t4 = ampmu*t3;
	       t6 = momega*y+mphase;
	       t7 = cos(t6);
	       t10 = momega*z+mphase;
	       t11 = sin(t10);
	       t12 = momega*t7*t11;
	       t13 = t4*t12;
	       t14 = cos(t2);
	       t16 = amplambda*t14*t12;
	       t19 = c*t;
	       t21 = omega*(x-t19);
	       t22 = -t21-phase;
	       t23 = sin(t22);
	       t24 = omega*omega;
	       t25 = t24*omega;
	       t26 = t23*t25;
	       t27 = c*c;
	       t29 = omega*x+phase;
	       t30 = sin(t29);
	       t31 = t27*t30;
	       t34 = -omega*(z-t19)-phase;
	       t35 = cos(t34);
	       t36 = t31*t35;
	       t37 = t26*t36;
	       t38 = cos(t22);
	       t39 = t38*t25;
	       t40 = sin(t34);
	       t41 = t31*t40;
	       t42 = t39*t41;
	       t43 = cos(t29);
	       t44 = t27*t43;
	       t45 = t44*t35;
	       t46 = t39*t45;
	       t47 = t44*t40;
	       t48 = t26*t47;
	       t56 = ampmu*(3.0/2.0+t14*t7*t11/2);
	       t62 = amplambda*(1.0/2.0+t3*t7*t11/4);
	       t63 = 2*t56+t62;
	       t64 = t24*t24;
	       t65 = t38*t64;
	       t66 = t65*t36;
	       t67 = t23*t64;
	       t68 = t67*t41;
	       t69 = t67*t45;
	       t70 = t65*t47;
	       t74 = sin(t21);
	       t75 = t74*t25;
	       t78 = -omega*(y-t19)-phase;
	       t79 = sin(t78);
	       t80 = t27*t79;
	       t82 = omega*z+phase;
	       t83 = cos(t82);
	       t84 = t80*t83;
	       t85 = t75*t84;
	       t88 = cos(t21);
	       t89 = t88*t25;
	       t90 = cos(t78);
	       t91 = t27*t90;
	       t92 = t91*t83;
	       t93 = t89*t92;
	       t97 = t64*t27;
	       t99 = t97*t79*t83;
	       t104 = t97*t90*t83;
	       t108 = omega*y+phase;
	       t109 = cos(t108);
	       t112 = t35*t25*t27;
	       t113 = t43*t109*t112;
	       t118 = t35*t27;
	       t121 = ampmu*t14;
	       t122 = sin(t6);
	       t124 = t122*momega*t11;
	       t125 = t121*t124;
	       t130 = 2*t56*t88*t99;
	       t133 = 2*t56*t74*t104;
	       t134 = t121*t7;
	       t135 = cos(t10);
	       t136 = t135*momega;
	       t139 = t109*t40*t27;
	       t141 = 2*t37;
	       t142 = 2*t42;
	       t143 = -t30*t25*t139-t141-t142;
	       t151 = 2*t66;
	       t152 = 2*t68;
	       forces[0] = (-t13+t16/4)*(-2*t37-2*t42-2*t46+2*t48)+t63*
		  (4*t66-4*t68-4*t69-4*t70)-t16*t85/2-t16*t93/2-2*t62*t88*t99+2*t62*t74*t104
		  -t16*t113/4+t62*t30*t64*t109*t118+t125*t93+t125*t85-t130+t133+t134
		  *t136*t143/2+t56*(t30*t64*t109*t35*t27+t151-t152);
	       t157 = amplambda*t3;
	       t160 = -t125-t157*t124/4;
	       t178 = momega*t11;
	       t186 = sin(t108);
	       t190 = t43*t186;
	       t194 = sin(t82);
	       t195 = t91*t194;
	       t198 = t80*t194;
	       t201 = -t190*t25*t40*t27+2*t75*t195-2*t89*t198;
	       t208 = t74*t64;
	       t211 = t88*t64;
	       forces[1] = t13*t93+t13*t85-t130+t133-2*t160*t74*t25*t84-2*t160*
		  t88*t25*t92+2*t63*t74*t64*t92-2*t63*t88*t64*t84-t157*t122*t178*
		  (-t141-t142-2*t46+2*t48-t113)/4+t62*t43*t186*t64*t118+t134*t136*t201/
		  2+t56*(t190*t64*t35*t27+2*t208*t92-2*t211*t84);
	       t222 = t43*t64*t139;
	       t231 = t208*t198;
	       t233 = t211*t195;
	       t238 = t7*t135*momega;
	       forces[2] = -t4*momega*t7*t11*t143/2+t56*(-t222+t151-t152-2*t69-2*t70)-
		  t121*t122*t178*t201/2+t56*(-t222+2*t231+2*t233)-(t121*t238+t157*t238/4)*
		  t43*t109*t112-t63*t43*t109*t40*t64*t27+t157*t7*t136*
		  (-2*t37-2*t42-2*t46+2*t48-2*t85-2*t93)/4+
		  t62*(2*t66-2*t68-2*t69-2*t70+2*t231+2*t233);
	       size_t ind = base+i+ni*j+nij*k;
	       fo[ind]        += forces[0];
	       fo[ind+nijk]   += forces[1];
	       fo[ind+2*nijk] += forces[2];
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::addmemvarforcing_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			      int klast, float_sw4* __restrict__ alpha, float_sw4 t, float_sw4 omega,
			      float_sw4 c, float_sw4 phase, float_sw4 omegaVE, float_sw4 dt,
			      float_sw4 h, float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const ssize_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 dto = dt*omegaVE;
      float_sw4 icp = 12*dto/( 6 + 6*dto + 3*dto*dto + dto*dto*dto );
      float_sw4 forces[3],t1,t10,t100,t105,t109,t113,t13,t139,t14,t141,t142,t144,t145,t146,t147,t148,t149,t15,t150,t154,t157,t17,t18,t19,t2,t20,t23,t25,t26,t27,t29,t30,t33,t34,t35,t36,t37,t4,t40,t42,t43,t45,t48,t5,t52,t53,t6,t60,t62,t63,t74,t82,t83,t84,t86,t88,t89,t9,t91,t93,t95,t97,t98,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
      {
	 float_sw4 z=(k-1)*h+zmin;
	 for( int j=jfirst; j<=jlast; j++ )
	 {
	    float_sw4 y=(j-1)*h;
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       float_sw4 x=(i-1)*h;
	       t1 = 1/omegaVE;
	       t2 = c*t;
	       t4 = omega*(x-t2);
	       t5 = -t4-phase;
	       t6 = sin(t5);
	       t9 = omega*x+phase;
	       t10 = sin(t9);
	       t13 = omega*(z-t2);
	       t14 = -t13-phase;
	       t15 = cos(t14);
	       t17 = t6*omega*c*t10*t15;
	       t18 = cos(t5);
	       t19 = t18*t10;
	       t20 = sin(t14);
	       t23 = t19*t20*omega*c;
	       t25 = t1*(-t17-t23);
	       t26 = t19*t15;
	       t27 = sin(t4);
	       t29 = omega*y+phase;
	       t30 = sin(t29);
	       t33 = omega*z+phase;
	       t34 = sin(t33);
	       t35 = t27*t30*t34;
	       t36 = dt*dt;
	       t37 = omegaVE*omegaVE;
	       t40 = omega*omega;
	       t42 = c*c;
	       t43 = t42*t10;
	       t45 = t18*t40*t43*t15;
	       t48 = t6*t40*t43*t20;
	       t52 = cos(t4);
	       t53 = t52*omega;
	       t60 = t40*omega;
	       t62 = t42*c;
	       t63 = t62*t10;
	       t74 = t27*t40;
	       forces[0] = t25+t26-t35+t36*(t37*(t25+t26-t35)+2*omegaVE*
			    (t1*(-2*t45+2*t48)-t17-t23+t53*c*t30*t34)+
			    t1*(4*t6*t60*t63*t15+4*t18*t60*t63*t20)-2*t45+2*t48+t74*t42*t30*t34)/6;
	       t82 = omega*(y-t2);
	       t83 = -t82-phase;
	       t84 = cos(t83);
	       t86 = cos(t33);
	       t88 = t53*c*t84*t86;
	       t89 = sin(t83);
	       t91 = omega*c;
	       t93 = t27*t89*t91*t86;
	       t95 = t1*(-t88-t93);
	       t97 = t27*t84*t86;
	       t98 = sin(t82);
	       t99 = t10*t98;
	       t100 = t99*t34;
	       t105 = t74*t42*t84*t86;
	       t109 = t52*t40*t42*t89*t86;
	       t113 = cos(t82);
	       forces[1] = t95+t97-t100+t36*(t37*(t95+t97-t100)+2*omegaVE*(t1*(
		   -2*t105+2*t109)-t88-t93+t10*t113*t91*t34)+t1*(4*t52*t60*t62*t84*t86+	
                	 4*t27*t60*t62*t89*t86)-2*t105+2*t109+t99*t40*t42*t34)/6;
	       t139 = cos(t9);
	       t141 = cos(t29);
	       t142 = t1*t139*t141;
	       t144 = t15*omega*c;
	       t145 = t142*t144;
	       t146 = t139*t141;
	       t147 = t146*t20;
	       t148 = t10*t30;
	       t149 = sin(t13);
	       t150 = t148*t149;
	       t154 = t20*t40*t42;
	       t157 = cos(t13);
	       forces[2] = -t145-t147-t150+t36*(t37*(-t145-t147-t150)+2*omegaVE
						*(t142*t154-t146*t144+t148*t157*omega*c)+
						t142*t15*t60*t62+t146*t154+t148*t149*t40*t42)/6;
	       size_t ind = base+i+ni*j+nij*k;
	       alpha[ind]        += forces[0]*icp;
	       alpha[ind+nijk]   += forces[1]*icp;
	       alpha[ind+2*nijk] += forces[2]*icp;
	    }
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::memvarforcesurf_ci( int ifirst, int ilast, int jfirst, int jlast,
			     int k, float_sw4* __restrict__ fo, float_sw4 t,
			     float_sw4 omega, float_sw4 c, float_sw4 phase,
			     float_sw4 omegaVE, float_sw4 dt, float_sw4 h,
			     float_sw4 zmin )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
//   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t1,t10,t100,t105,t109,t113,t13,t139,t14,t141,t142,t144,t145,t146,t147,t148,t149,t15,t150,t154,t157,t17,t18,t19,t2,t20,t23,t25,t26,t27,t29,t30,t33,t34,t35,t36,t37,t4,t40,t42,t43,t45,t48,t5,t52,t53,t6,t60,t62,t63,t74,t82,t83,t84,t86,t88,t89,t9,t91,t93,t95,t97,t98,t99;
      float_sw4 z = (k-1)*h+zmin;
#pragma omp for
      for( int j=jfirst; j<=jlast; j++ )
      {
	 float_sw4 y=(j-1)*h;
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    float_sw4 x=(i-1)*h;
        t1 = 1/omegaVE;
        t2 = c*t;
        t4 = omega*(x-t2);
        t5 = -t4-phase;
        t6 = sin(t5);
        t9 = omega*x+phase;
        t10 = sin(t9);
        t13 = omega*(z-t2);
        t14 = -t13-phase;
        t15 = cos(t14);
        t17 = t6*omega*c*t10*t15;
        t18 = cos(t5);
        t19 = t18*t10;
        t20 = sin(t14);
        t23 = t19*t20*omega*c;
        t25 = t1*(-t17-t23);
        t26 = t19*t15;
        t27 = sin(t4);
        t29 = omega*y+phase;
        t30 = sin(t29);
        t33 = omega*z+phase;
        t34 = sin(t33);
        t35 = t27*t30*t34;
        t36 = dt*dt;
        t37 = omegaVE*omegaVE;
        t40 = omega*omega;
        t42 = c*c;
        t43 = t42*t10;
        t45 = t18*t40*t43*t15;
        t48 = t6*t40*t43*t20;
        t52 = cos(t4);
        t53 = t52*omega;
        t60 = t40*omega;
        t62 = t42*c;
        t63 = t62*t10;
        t74 = t27*t40;
        forces[0] = t25+t26-t35+t36*(t37*(t25+t26-t35)+2*omegaVE*(t1*(-2
         *t45+2*t48)-t17-t23+t53*c*t30*t34)+t1*(4*t6*t60*t63*t15+4*t18*t60*
         t63*t20)-2*t45+2*t48+t74*t42*t30*t34)/6;
        t82 = omega*(y-t2);
        t83 = -t82-phase;
        t84 = cos(t83);
        t86 = cos(t33);
        t88 = t53*c*t84*t86;
        t89 = sin(t83);
        t91 = omega*c;
        t93 = t27*t89*t91*t86;
        t95 = t1*(-t88-t93);
        t97 = t27*t84*t86;
        t98 = sin(t82);
        t99 = t10*t98;
        t100 = t99*t34;
        t105 = t74*t42*t84*t86;
        t109 = t52*t40*t42*t89*t86;
        t113 = cos(t82);
        forces[1] = t95+t97-t100+t36*(t37*(t95+t97-t100)+2*omegaVE*(t1*(
               -2*t105+2*t109)-t88-t93+t10*t113*t91*t34)+
           t1*(4*t52*t60*t62*t84*t86+4*t27*t60*t62*t89*t86)-
             2*t105+2*t109+t99*t40*t42*t34)/6;
        t139 = cos(t9);
        t141 = cos(t29);
        t142 = t1*t139*t141;
        t144 = t15*omega*c;
        t145 = t142*t144;
        t146 = t139*t141;
        t147 = t146*t20;
        t148 = t10*t30;
        t149 = sin(t13);
        t150 = t148*t149;
        t154 = t20*t40*t42;
        t157 = cos(t13);
        forces[2] = -t145-t147-t150+t36*(t37*(-t145-t147-t150)+2*omegaVE
             *(t142*t154-t146*t144+t148*t157*omega*c)+t142*t15*t60*t62+
               t146*t154+t148*t149*t40*t42)/6;
	    size_t ind = base+i+ni*j;
	    fo[ind] = forces[0];
	    fo[ind+nij] = forces[1];
	    fo[ind+2*nij] = forces[2];
	 }
      }
   }
}

//-----------------------------------------------------------------------
void EW::forcingfortattc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			     int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			     float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase,
			     float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			     float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			     float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t100,t103,t108,t109,t11,t110,t112,t113,t116,t118,t12,t131,t133,t135,t14,t144,t15,t152,t156,t158,t159,t162,t178,t180,t182,t189,t19,t191,t195,t2,t21,t22,t23,t26,t27,t3,t30,t31,t32,t33,t34,t35,t36,t37,t38,t4,t45,t51,t52,t53,t55,t56,t6,t62,t63,t67,t68,t7,t71,t72,t73,t77,t8,t80,t84,t85,t86,t87,t95,t96,t97;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];
	       t2 = momega*x+mphase;
	       t3 = sin(t2);
	       t4 = ampmu*t3;
	       t6 = momega*y+mphase;
	       t7 = cos(t6);
	       t8 = momega*t7;
	       t10 = momega*z+mphase;
	       t11 = sin(t10);
	       t12 = t8*t11;
	       t14 = cos(t2);
	       t15 = amplambda*t14;
	       t19 = c*t;
	       t21 = omega*(x-t19);
	       t22 = -t21-phase;
	       t23 = sin(t22);
	       t26 = omega*x+phase;
	       t27 = sin(t26);
	       t30 = -omega*(z-t19)-phase;
	       t31 = cos(t30);
	       t32 = t27*t31;
	       t33 = t23*omega*t32;
	       t34 = cos(t22);
	       t35 = cos(t26);
	       t36 = t34*t35;
	       t37 = omega*t31;
	       t38 = t36*t37;
	       t45 = ampmu*(3.0/2.0+t14*t7*t11/2);
	       t51 = amplambda*(1.0/2.0+t3*t7*t11/4);
	       t52 = 2*t45+t51;
	       t53 = omega*omega;
	       t55 = t34*t53*t32;
	       t56 = t23*t53;
	       t62 = t15*t8;
	       t63 = sin(t21);
	       t67 = -omega*(y-t19)-phase;
	       t68 = sin(t67);
	       t71 = omega*z+phase;
	       t72 = cos(t71);
	       t73 = t68*omega*t72;
	       t77 = cos(t21);
	       t80 = t53*t68*t72;
	       t84 = omega*y+phase;
	       t85 = cos(t84);
	       t86 = t85*t31;
	       t87 = t86*omega;
	       t95 = ampmu*t14;
	       t96 = sin(t6);
	       t97 = t96*momega;
	       t100 = cos(t67);
	       t103 = t11*t77*omega*t100*t72;
	       t108 = t95*t7;
	       t109 = cos(t10);
	       t110 = t109*momega;
	       t112 = sin(t30);
	       t113 = t85*t112;
	       t116 = t112*omega;
	       t118 = t27*omega*t113+t34*t27*t116;
	       forces[0] = (-t4*t12+t15*t12/4)*(t33+t38)+t52*(-2*t55+2*t56*t35*
		      t31)+t62*t11*t63*t73/4+t51*t77*t80+t62*t11*t35*t87/4-t51*t27*t53*t85*t31
		  -t95*t97*t103/2+t45*t77*t80+t108*t110*t118/2+t45*(-t27*t53*t86-t55);
	       t131 = t53*t100*t72;
	       t133 = t97*t11;
	       t135 = amplambda*t3;
	       t144 = momega*t11;
	       t152 = sin(t84);
	       t156 = t35*t152;
	       t158 = t63*t100;
	       t159 = sin(t71);
	       t162 = t156*t116-t158*t159*omega;
	       forces[1] = -t4*t8*t103/2-t45*t63*t131+(-t95*t133-t135*t133/4)*t63*t73-
		  t52*t63*t131-t135*t96*t144*(t33+t38+t35*t85*t37)/4-t51*t35*
		  t152*t53*t31+t108*t110*t162/2+t45*(-t156*t53*t31-t158*t72*t53);
	       t178 = t35*t53*t113;
	       t180 = t56*t27*t112;
	       t182 = t36*t53*t112;
	       t189 = t63*t68;
	       t191 = t189*t53*t159;
	       t195 = t7*t109*momega;
	       forces[2] = -t4*momega*t7*t11*t118/2+t45*(t178+t180+t182)-
		  t95*t96*t144*t162/2+t45*(t178-t191)+(t95*t195+t135*t195/4)*t35*t87+
		  t52*t35*t113*t53+t135*t7*t110*(t33+t38+t189*omega*t72)/4+t51*(t180+t182-t191);
	       fo[ind]        += forces[0];
	       fo[ind+nijk]   += forces[1];
	       fo[ind+2*nijk] += forces[2];
	    }
   }
}

//-----------------------------------------------------------------------
void EW::forcingttattfortc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			       int klast, float_sw4* __restrict__ fo, float_sw4 t, float_sw4 omega,
			       float_sw4 c, float_sw4 phase, float_sw4 momega, float_sw4 mphase,
			       float_sw4 amprho, float_sw4 ampmu, float_sw4 amplambda,
			       float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			       float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t10,t104,t108,t109,t11,t112,t113,t118,t12,t121,t122,t124,t125,t13,t130,t133,t134,t135,t136,t139,t14,t141,t142,t143,t151,t152,t157,t16,t160,t178,t186,t19,t190,t194,t195,t198,t2,t201,t208,t21,t211,t22,t222,t23,t231,t233,t238,t24,t25,t26,t27,t29,t3,t30,t31,t34,t35,t36,t37,t38,t39,t4,t40,t41,t42,t43,t44,t45,t46,t47,t48,t56,t6,t62,t63,t64,t65,t66,t67,t68,t69,t7,t70,t74,t75,t78,t79,t80,t82,t83,t84,t85,t88,t89,t90,t91,t92,t93,t97,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];

	       t2 = momega*x+mphase;
	       t3 = sin(t2);
	       t4 = ampmu*t3;
	       t6 = momega*y+mphase;
	       t7 = cos(t6);
	       t10 = momega*z+mphase;
	       t11 = sin(t10);
	       t12 = momega*t7*t11;
	       t13 = t4*t12;
	       t14 = cos(t2);
	       t16 = amplambda*t14*t12;
	       t19 = c*t;
	       t21 = omega*(x-t19);
	       t22 = -t21-phase;
	       t23 = sin(t22);
	       t24 = omega*omega;
	       t25 = t24*omega;
	       t26 = t23*t25;
	       t27 = c*c;
	       t29 = omega*x+phase;
	       t30 = sin(t29);
	       t31 = t27*t30;
	       t34 = -omega*(z-t19)-phase;
	       t35 = cos(t34);
	       t36 = t31*t35;
	       t37 = t26*t36;
	       t38 = cos(t22);
	       t39 = t38*t25;
	       t40 = sin(t34);
	       t41 = t31*t40;
	       t42 = t39*t41;
	       t43 = cos(t29);
	       t44 = t27*t43;
	       t45 = t44*t35;
	       t46 = t39*t45;
	       t47 = t44*t40;
	       t48 = t26*t47;
	       t56 = ampmu*(3.0/2.0+t14*t7*t11/2);
	       t62 = amplambda*(1.0/2.0+t3*t7*t11/4);
	       t63 = 2*t56+t62;
	       t64 = t24*t24;
	       t65 = t38*t64;
	       t66 = t65*t36;
	       t67 = t23*t64;
	       t68 = t67*t41;
	       t69 = t67*t45;
	       t70 = t65*t47;
	       t74 = sin(t21);
	       t75 = t74*t25;
	       t78 = -omega*(y-t19)-phase;
	       t79 = sin(t78);
	       t80 = t27*t79;
	       t82 = omega*z+phase;
	       t83 = cos(t82);
	       t84 = t80*t83;
	       t85 = t75*t84;
	       t88 = cos(t21);
	       t89 = t88*t25;
	       t90 = cos(t78);
	       t91 = t27*t90;
	       t92 = t91*t83;
	       t93 = t89*t92;
	       t97 = t64*t27;
	       t99 = t97*t79*t83;
	       t104 = t97*t90*t83;
	       t108 = omega*y+phase;
	       t109 = cos(t108);
	       t112 = t35*t25*t27;
	       t113 = t43*t109*t112;
	       t118 = t35*t27;
	       t121 = ampmu*t14;
	       t122 = sin(t6);
	       t124 = t122*momega*t11;
	       t125 = t121*t124;
	       t130 = 2*t56*t88*t99;
	       t133 = 2*t56*t74*t104;
	       t134 = t121*t7;
	       t135 = cos(t10);
	       t136 = t135*momega;
	       t139 = t109*t40*t27;
	       t141 = 2*t37;
	       t142 = 2*t42;
	       t143 = -t30*t25*t139-t141-t142;
	       t151 = 2*t66;
	       t152 = 2*t68;
	       forces[0] = (-t13+t16/4)*(-2*t37-2*t42-2*t46+2*t48)+t63*
		  (4*t66-4*t68-4*t69-4*t70)-t16*t85/2-t16*t93/2-2*t62*t88*t99+2*t62*t74*t104
		  -t16*t113/4+t62*t30*t64*t109*t118+t125*t93+t125*t85-t130+t133+t134
		  *t136*t143/2+t56*(t30*t64*t109*t35*t27+t151-t152);
	       t157 = amplambda*t3;
	       t160 = -t125-t157*t124/4;
	       t178 = momega*t11;
	       t186 = sin(t108);
	       t190 = t43*t186;
	       t194 = sin(t82);
	       t195 = t91*t194;
	       t198 = t80*t194;
	       t201 = -t190*t25*t40*t27+2*t75*t195-2*t89*t198;
	       t208 = t74*t64;
	       t211 = t88*t64;
	       forces[1] = t13*t93+t13*t85-t130+t133-2*t160*t74*t25*t84-2*t160*
		  t88*t25*t92+2*t63*t74*t64*t92-2*t63*t88*t64*t84-t157*t122*t178*
		  (-t141-t142-2*t46+2*t48-t113)/4+t62*t43*t186*t64*t118+t134*t136*t201/
		  2+t56*(t190*t64*t35*t27+2*t208*t92-2*t211*t84);
	       t222 = t43*t64*t139;
	       t231 = t208*t198;
	       t233 = t211*t195;
	       t238 = t7*t135*momega;
	       forces[2] = -t4*momega*t7*t11*t143/2+t56*(-t222+t151-t152-2*t69-2*t70)-
		  t121*t122*t178*t201/2+t56*(-t222+2*t231+2*t233)-(t121*t238+t157*t238/4)*
		  t43*t109*t112-t63*t43*t109*t40*t64*t27+t157*t7*t136*
		  (-2*t37-2*t42-2*t46+2*t48-2*t85-2*t93)/4+
		  t62*(2*t66-2*t68-2*t69-2*t70+2*t231+2*t233);
	       fo[ind]        += forces[0];
	       fo[ind+nijk]   += forces[1];
	       fo[ind+2*nijk] += forces[2];
	    }
   }
}

//-----------------------------------------------------------------------
void EW::addmemvarforcingc_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, 
			       int klast, float_sw4* __restrict__ alpha, float_sw4 t, float_sw4 omega,
			       float_sw4 c, float_sw4 phase, float_sw4 omegaVE, float_sw4 dt,
			       float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			       float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base  = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 dto = dt*omegaVE;
      float_sw4 icp = 12*dto/( 6 + 6*dto + 3*dto*dto + dto*dto*dto );
      float_sw4 forces[3],t1,t10,t100,t105,t109,t113,t13,t139,t14,t141,t142,t144,t145,t146,t147,t148,t149,t15,t150,t154,t157,t17,t18,t19,t2,t20,t23,t25,t26,t27,t29,t30,t33,t34,t35,t36,t37,t4,t40,t42,t43,t45,t48,t5,t52,t53,t6,t60,t62,t63,t74,t82,t83,t84,t86,t88,t89,t9,t91,t93,t95,t97,t98,t99;
#pragma omp for
      for( int k=kfirst; k<=klast; k++ )
	 for( int j=jfirst; j<=jlast; j++ )
	    for( int i=ifirst; i<=ilast; i++ )
	    {
	       size_t ind = base+i+ni*j+nij*k;
	       float_sw4 z=zz[ind];
	       float_sw4 y=yy[ind];
	       float_sw4 x=xx[ind];
	       t1 = 1/omegaVE;
	       t2 = c*t;
	       t4 = omega*(x-t2);
	       t5 = -t4-phase;
	       t6 = sin(t5);
	       t9 = omega*x+phase;
	       t10 = sin(t9);
	       t13 = omega*(z-t2);
	       t14 = -t13-phase;
	       t15 = cos(t14);
	       t17 = t6*omega*c*t10*t15;
	       t18 = cos(t5);
	       t19 = t18*t10;
	       t20 = sin(t14);
	       t23 = t19*t20*omega*c;
	       t25 = t1*(-t17-t23);
	       t26 = t19*t15;
	       t27 = sin(t4);
	       t29 = omega*y+phase;
	       t30 = sin(t29);
	       t33 = omega*z+phase;
	       t34 = sin(t33);
	       t35 = t27*t30*t34;
	       t36 = dt*dt;
	       t37 = omegaVE*omegaVE;
	       t40 = omega*omega;
	       t42 = c*c;
	       t43 = t42*t10;
	       t45 = t18*t40*t43*t15;
	       t48 = t6*t40*t43*t20;
	       t52 = cos(t4);
	       t53 = t52*omega;
	       t60 = t40*omega;
	       t62 = t42*c;
	       t63 = t62*t10;
	       t74 = t27*t40;
	       forces[0] = t25+t26-t35+t36*(t37*(t25+t26-t35)+2*omegaVE*
			    (t1*(-2*t45+2*t48)-t17-t23+t53*c*t30*t34)+
			    t1*(4*t6*t60*t63*t15+4*t18*t60*t63*t20)-2*t45+2*t48+t74*t42*t30*t34)/6;
	       t82 = omega*(y-t2);
	       t83 = -t82-phase;
	       t84 = cos(t83);
	       t86 = cos(t33);
	       t88 = t53*c*t84*t86;
	       t89 = sin(t83);
	       t91 = omega*c;
	       t93 = t27*t89*t91*t86;
	       t95 = t1*(-t88-t93);
	       t97 = t27*t84*t86;
	       t98 = sin(t82);
	       t99 = t10*t98;
	       t100 = t99*t34;
	       t105 = t74*t42*t84*t86;
	       t109 = t52*t40*t42*t89*t86;
	       t113 = cos(t82);
	       forces[1] = t95+t97-t100+t36*(t37*(t95+t97-t100)+2*omegaVE*(t1*(
		   -2*t105+2*t109)-t88-t93+t10*t113*t91*t34)+t1*(4*t52*t60*t62*t84*t86+	
                	 4*t27*t60*t62*t89*t86)-2*t105+2*t109+t99*t40*t42*t34)/6;
	       t139 = cos(t9);
	       t141 = cos(t29);
	       t142 = t1*t139*t141;
	       t144 = t15*omega*c;
	       t145 = t142*t144;
	       t146 = t139*t141;
	       t147 = t146*t20;
	       t148 = t10*t30;
	       t149 = sin(t13);
	       t150 = t148*t149;
	       t154 = t20*t40*t42;
	       t157 = cos(t13);
	       forces[2] = -t145-t147-t150+t36*(t37*(-t145-t147-t150)+2*omegaVE
						*(t142*t154-t146*t144+t148*t157*omega*c)+
						t142*t15*t60*t62+t146*t154+t148*t149*t40*t42)/6;

	       alpha[ind]        += forces[0]*icp;
	       alpha[ind+nijk]   += forces[1]*icp;
	       alpha[ind+2*nijk] += forces[2]*icp;
	    }
   }
}

//-----------------------------------------------------------------------
void EW::memvarforcesurfc_ci( int ifirst, int ilast, int jfirst, int jlast,
			      int kfirst, int klast, int k, float_sw4* __restrict__ fo, 
			      float_sw4 t, float_sw4 omega, float_sw4 c, float_sw4 phase,
			      float_sw4 omegaVE, float_sw4 dt, 
			      float_sw4* __restrict__ xx, float_sw4* __restrict__ yy,
			      float_sw4* __restrict__ zz )
{
   const size_t ni    = ilast-ifirst+1;
   const size_t nij   = ni*(jlast-jfirst+1);
   //const size_t nijk  = nij*(klast-kfirst+1);
   const size_t base2 = -(ifirst+ni*jfirst);
   const size_t base3 = -(ifirst+ni*jfirst+nij*kfirst);
#pragma omp parallel
   {
      float_sw4 forces[3],t1,t10,t100,t105,t109,t113,t13,t139,t14,t141,t142,t144,t145,t146,t147,t148,t149,t15,t150,t154,t157,t17,t18,t19,t2,t20,t23,t25,t26,t27,t29,t30,t33,t34,t35,t36,t37,t4,t40,t42,t43,t45,t48,t5,t52,t53,t6,t60,t62,t63,t74,t82,t83,t84,t86,t88,t89,t9,t91,t93,t95,t97,t98,t99;
#pragma omp for
      for( int j=jfirst; j<=jlast; j++ )
	 for( int i=ifirst; i<=ilast; i++ )
	 {
	    size_t ind = base3+i+ni*j+nij*k;
	    float_sw4 z=zz[ind];
	    float_sw4 y=yy[ind];
	    float_sw4 x=xx[ind];
        t1 = 1/omegaVE;
        t2 = c*t;
        t4 = omega*(x-t2);
        t5 = -t4-phase;
        t6 = sin(t5);
        t9 = omega*x+phase;
        t10 = sin(t9);
        t13 = omega*(z-t2);
        t14 = -t13-phase;
        t15 = cos(t14);
        t17 = t6*omega*c*t10*t15;
        t18 = cos(t5);
        t19 = t18*t10;
        t20 = sin(t14);
        t23 = t19*t20*omega*c;
        t25 = t1*(-t17-t23);
        t26 = t19*t15;
        t27 = sin(t4);
        t29 = omega*y+phase;
        t30 = sin(t29);
        t33 = omega*z+phase;
        t34 = sin(t33);
        t35 = t27*t30*t34;
        t36 = dt*dt;
        t37 = omegaVE*omegaVE;
        t40 = omega*omega;
        t42 = c*c;
        t43 = t42*t10;
        t45 = t18*t40*t43*t15;
        t48 = t6*t40*t43*t20;
        t52 = cos(t4);
        t53 = t52*omega;
        t60 = t40*omega;
        t62 = t42*c;
        t63 = t62*t10;
        t74 = t27*t40;
        forces[0] = t25+t26-t35+t36*(t37*(t25+t26-t35)+2*omegaVE*(t1*(-2
         *t45+2*t48)-t17-t23+t53*c*t30*t34)+t1*(4*t6*t60*t63*t15+4*t18*t60*
         t63*t20)-2*t45+2*t48+t74*t42*t30*t34)/6;
        t82 = omega*(y-t2);
        t83 = -t82-phase;
        t84 = cos(t83);
        t86 = cos(t33);
        t88 = t53*c*t84*t86;
        t89 = sin(t83);
        t91 = omega*c;
        t93 = t27*t89*t91*t86;
        t95 = t1*(-t88-t93);
        t97 = t27*t84*t86;
        t98 = sin(t82);
        t99 = t10*t98;
        t100 = t99*t34;
        t105 = t74*t42*t84*t86;
        t109 = t52*t40*t42*t89*t86;
        t113 = cos(t82);
        forces[1] = t95+t97-t100+t36*(t37*(t95+t97-t100)+2*omegaVE*(t1*(
               -2*t105+2*t109)-t88-t93+t10*t113*t91*t34)+
           t1*(4*t52*t60*t62*t84*t86+4*t27*t60*t62*t89*t86)-
             2*t105+2*t109+t99*t40*t42*t34)/6;
        t139 = cos(t9);
        t141 = cos(t29);
        t142 = t1*t139*t141;
        t144 = t15*omega*c;
        t145 = t142*t144;
        t146 = t139*t141;
        t147 = t146*t20;
        t148 = t10*t30;
        t149 = sin(t13);
        t150 = t148*t149;
        t154 = t20*t40*t42;
        t157 = cos(t13);
        forces[2] = -t145-t147-t150+t36*(t37*(-t145-t147-t150)+2*omegaVE
             *(t142*t154-t146*t144+t148*t157*omega*c)+t142*t15*t60*t62+
               t146*t154+t148*t149*t40*t42)/6;
	    size_t ind2 = base2+i+ni*j;
	    fo[ind2]   = forces[0];
	    fo[ind2+nij] = forces[1];
	    fo[ind2+2*nij] = forces[2];
	 }
   }
}
