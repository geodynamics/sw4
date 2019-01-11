#include "Sarray.h"
#include "Mspace.h"
#include "caliper.h"
#include "policies.h"
//-----------------From WPP------------------------------------------------------
void  addMemVarPredCart( float_sw4 zMin, float_sw4 h, float_sw4 t, Sarray &alpha,
                         float_sw4 omegaVE, float_sw4 dt ,float_sw4 omega, float_sw4 phase, float_sw4 c)
{
  SW4_MARK_FUNCTION;
   //  float_sw4 x, y, z;
  // double omega =m_omega;
  // double phase =m_phase;
  // double c=m_c;
  float_sw4 cof = 2*dt/(1+omegaVE*dt);
// #pragma omp parallel for
//   for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
//     for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
// #pragma ivdep
// #pragma simd
//       for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
//       {
  SView &alphaV = alpha.getview();
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION==6
  using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
#elif SW4_RAJA_VERSION==7

  using LOCAL_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			    RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
						  RAJA::statement::Tile<2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
									RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
											     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
														  RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
																       RAJA::statement::Lambda<0> >>>>>>>>;

#endif

#else
  using LOCAL_POL = DEFAULT_LOOP3;
#endif
	   RAJA::RangeSegment k_range(alpha.m_kb,alpha.m_ke+1);
	   RAJA::RangeSegment j_range(alpha.m_jb,alpha.m_je+1);
	   RAJA::RangeSegment i_range(alpha.m_ib,alpha.m_ie+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
	float_sw4 x = (i-1)*h;
	float_sw4 y = (j-1)*h;
	float_sw4 z = zMin + (k-1)*h;
	{
	   float_sw4 t1;
	   float_sw4 t12;
	   float_sw4 t13;
	   float_sw4 t14;
	   float_sw4 t17;
	   float_sw4 t18;
	   float_sw4 t19;
	   float_sw4 t24;
	   float_sw4 t26;
	   float_sw4 t27;
	   float_sw4 t3;
	   float_sw4 t30;
	   float_sw4 t31;
	   float_sw4 t35;
	   float_sw4 t38;
	   float_sw4 t39;
	   float_sw4 t4;
	   float_sw4 t40;
	   float_sw4 t42;
	   float_sw4 t45;
	   float_sw4 t5;
	   float_sw4 t52;
	   float_sw4 t57;
	   float_sw4 t58;
	   float_sw4 t59;
	   float_sw4 t65;
	   float_sw4 t8;
	   float_sw4 t9;
           float_sw4 forces[3];
	   {
	      t1 = c*t;
	      t3 = omega*(x-t1);
	      t4 = -t3-phase;
	      t5 = sin(t4);
	      t8 = omega*x+phase;
	      t9 = sin(t8);
	      t12 = omega*(z-t1);
	      t13 = -t12-phase;
	      t14 = cos(t13);
	      t17 = cos(t4);
	      t18 = t17*t9;
	      t19 = sin(t13);
	      t24 = sin(t3);
	      t26 = omega*y+phase;
	      t27 = sin(t26);
	      t30 = omega*z+phase;
	      t31 = sin(t30);
	      forces[0] = -t5*omega*c*t9*t14-t18*t19*omega*c+omegaVE*(t18*t14-t24*t27*t31
								      );
              alphaV(1,i,j,k) += cof*forces[0];
	      t35 = cos(t3);
	      t38 = omega*(y-t1);
	      t39 = -t38-phase;
	      t40 = cos(t39);
	      t42 = cos(t30);
	      t45 = sin(t39);
	      t52 = sin(t38);
	      forces[1] = -t35*omega*c*t40*t42-t24*t45*omega*c*t42+omegaVE*(t24*t40*t42-
									    t9*t52*t31);
              alphaV(2,i,j,k) += cof*forces[1];
	      t57 = cos(t8);
	      t58 = cos(t26);
	      t59 = t57*t58;
	      t65 = sin(t12);
	      forces[2] = -t59*t14*omega*c+omegaVE*(-t59*t19-t9*t27*t65);
              alphaV(3,i,j,k) += cof*forces[2];
	   }
	}
			       }); SYNC_STREAM;
}

//-----------------------FROM WPP------------------------------------------------
void addMemVarPredCurvilinear( Sarray& a_X, Sarray& a_Y, Sarray& a_Z, float_sw4 t,
			       Sarray& alpha, float_sw4 omegaVE, float_sw4 dt, float_sw4 omega,
			       float_sw4 phase, float_sw4 c )
{
  SW4_MARK_FUNCTION;
  // double omega =m_omega;
  // double phase =m_phase;
  // double c=m_c;
  float_sw4 cof = 2*dt/(1+omegaVE*dt);
// #pragma omp parallel for
//   for( int k=a_X.m_kb ; k<=a_X.m_ke; k++ )
//     for( int j=a_X.m_jb ; j<=a_X.m_je; j++ )
// #pragma ivdep
// #pragma simd
//       for( int i=a_X.m_ib ; i<=a_X.m_ie; i++ )
//       {
  SView &alphaV = alpha.getview();
  SView &a_XV = a_X.getview();
   SView &a_YV = a_Y.getview();
   SView &a_ZV = a_Z.getview();
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION==6
   using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernel<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;

#elif SW4_RAJA_VERSION==7

  using LOCAL_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			    RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
						  RAJA::statement::Tile<2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
									RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
											     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
														  RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
																       RAJA::statement::Lambda<0> >>>>>>>>;

#endif


#else
   using LOCAL_POL = DEFAULT_LOOP3;
#endif
	   RAJA::RangeSegment k_range(a_X.m_kb,a_X.m_ke+1);
	   RAJA::RangeSegment j_range(a_X.m_jb,a_X.m_je+1);
	   RAJA::RangeSegment i_range(a_X.m_ib,a_X.m_ie+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {

	float_sw4 x = a_XV(i,j,k);
	float_sw4 y = a_YV(i,j,k);
	float_sw4 z = a_ZV(i,j,k);
	{
	   float_sw4 t1;
	   float_sw4 t12;
	   float_sw4 t13;
	   float_sw4 t14;
	   float_sw4 t17;
	   float_sw4 t18;
	   float_sw4 t19;
	   float_sw4 t24;
	   float_sw4 t26;
	   float_sw4 t27;
	   float_sw4 t3;
	   float_sw4 t30;
	   float_sw4 t31;
	   float_sw4 t35;
	   float_sw4 t38;
	   float_sw4 t39;
	   float_sw4 t4;
	   float_sw4 t40;
	   float_sw4 t42;
	   float_sw4 t45;
	   float_sw4 t5;
	   float_sw4 t52;
	   float_sw4 t57;
	   float_sw4 t58;
	   float_sw4 t59;
	   float_sw4 t65;
	   float_sw4 t8;
	   float_sw4 t9;
           float_sw4 forces[3];
	   {
	      t1 = c*t;
	      t3 = omega*(x-t1);
	      t4 = -t3-phase;
	      t5 = sin(t4);
	      t8 = omega*x+phase;
	      t9 = sin(t8);
	      t12 = omega*(z-t1);
	      t13 = -t12-phase;
	      t14 = cos(t13);
	      t17 = cos(t4);
	      t18 = t17*t9;
	      t19 = sin(t13);
	      t24 = sin(t3);
	      t26 = omega*y+phase;
	      t27 = sin(t26);
	      t30 = omega*z+phase;
	      t31 = sin(t30);
	      forces[0] = -t5*omega*c*t9*t14-t18*t19*omega*c+omegaVE*(t18*t14-t24*t27*t31
								      );
              alphaV(1,i,j,k) += cof*forces[0];
	      t35 = cos(t3);
	      t38 = omega*(y-t1);
	      t39 = -t38-phase;
	      t40 = cos(t39);
	      t42 = cos(t30);
	      t45 = sin(t39);
	      t52 = sin(t38);
	      forces[1] = -t35*omega*c*t40*t42-t24*t45*omega*c*t42+omegaVE*(t24*t40*t42-
									    t9*t52*t31);
              alphaV(2,i,j,k) += cof*forces[1];
	      t57 = cos(t8);
	      t58 = cos(t26);
	      t59 = t57*t58;
	      t65 = sin(t12);
	      forces[2] = -t59*t14*omega*c+omegaVE*(-t59*t19-t9*t27*t65);
              alphaV(3,i,j,k) += cof*forces[2];
	   }
	}
			       }); SYNC_STREAM;
}

// //---------------------------- NEW MAPLE GENERATED CODE ----------------------------------
// void addMemVarCorrCart(double zMin, double h, double t, Sarray &alpha,
//                        double omegaVE, double dt, double omega, double phase, double c )
// {
//    double x, y, z, dto=omegaVE*dt;
//    double cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

//    double forces[3];
//   double t1;
//   double t11;
//   double t12;
//   double t121;
//   double t123;
//   double t124;
//   double t125;
//   double t127;
//   double t129;
//   double t13;
//   double t130;
//   double t134;
//   double t137;
//   double t15;
//   double t16;
//   double t17;
//   double t18;
//   double t19;
//   double t2;
//   double t20;
//   double t21;
//   double t22;
//   double t24;
//   double t27;
//   double t28;
//   double t3;
//   double t30;
//   double t31;
//   double t34;
//   double t35;
//   double t36;
//   double t39;
//   double t4;
//   double t40;
//   double t41;
//   double t42;
//   double t44;
//   double t48;
//   double t5;
//   double t57;
//   double t7;
//   double t70;
//   double t71;
//   double t72;
//   double t74;
//   double t75;
//   double t76;
//   double t78;
//   double t8;
//   double t81;
//   double t85;
//   double t86;
//   double t9;
//   double t91;
//   double t94;
//   double t99;

//   for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
//     for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
//       for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
//       {
// 	x = (i-1)*h;
// 	y = (j-1)*h;
// 	z = zMin + (k-1)*h;

//         t1 = dt*dt;
//         t2 = omegaVE*omegaVE;
//         t3 = 1/omegaVE;
//         t4 = omega*c;
//         t5 = c*t;
//         t7 = omega*(-t5+x);
//         t8 = -t7-phase;
//         t9 = sin(t8);
//         t11 = omega*x+phase;
//         t12 = sin(t11);
//         t13 = t9*t12;
//         t15 = omega*(-t5+z);
//         t16 = -t15-phase;
//         t17 = cos(t16);
//         t18 = t13*t17;
//         t19 = t4*t18;
//         t20 = cos(t8);
//         t21 = t20*t12;
//         t22 = sin(t16);
//         t24 = t21*t4*t22;
//         t27 = t21*t17;
//         t28 = sin(t7);
//         t30 = omega*y+phase;
//         t31 = sin(t30);
//         t34 = omega*z+phase;
//         t35 = sin(t34);
//         t36 = t28*t31*t35;
//         t39 = omega*omega;
//         t40 = c*c;
//         t41 = t39*t40;
//         t42 = t41*t27;
//         t44 = t41*t13*t22;
//         t48 = cos(t7);
//         t57 = t39*omega*t40*c;
//         forces[0] = t1*(t2*(t3*(-t19-t24)+t27-t36)+2.0*omegaVE*(t3*(-2.0*t42+2.0*
//                                                                     t44)-t19-t24+t4*t48*t31*t35)+t3*(4.0*t21*t22*t57+4.0*t18*t57)-2.0*t42+2.0*t44+
//                         t41*t36)/6.0;
//         t70 = omega*(-t5+y);
//         t71 = -t70-phase;
//         t72 = cos(t71);
//         t74 = cos(t34);
//         t75 = t48*t72*t74;
//         t76 = t4*t75;
//         t78 = sin(t71);
//         t81 = t28*omega*c*t78*t74;
//         t85 = t28*t72*t74;
//         t86 = sin(t70);
//         t91 = t41*t85;
//         t94 = t41*t48*t78*t74;
//         t99 = cos(t70);
//         forces[1] = t1*(t2*(t3*(-t76-t81)+t85-t12*t86*t35)+2.0*omegaVE*(t3*(-2.0*
//                                                                             t91+2.0*t94)-t76-t81+t12*omega*c*t99*t35)+t3*(4.0*t28*t57*t74*t78+4.0*t57*t75)
//                         -2.0*t91+2.0*t94+t12*t39*t40*t86*t35)/6.0;
//         t121 = cos(t11);
//         t123 = cos(t30);
//         t124 = t3*t121*t123;
//         t125 = t4*t17;
//         t127 = t121*t123;
//         t129 = t12*t31;
//         t130 = sin(t15);
//         t134 = t41*t22;
//         t137 = cos(t15);
//         forces[2] = t1*(t2*(-t124*t125-t127*t22-t129*t130)+2.0*omegaVE*(t129*t137*
//                                                                         t4+t124*t134-t125*t127)+t124*t57*t17+t127*t134+t129*t41*t130)/6.0;
//         alpha(1,i,j,k) += cof*forces[0];
//         alpha(2,i,j,k) += cof*forces[1];
//         alpha(3,i,j,k) += cof*forces[2];
//       } // end for*3
//   return;
// } // end function

// //----------------------------------------------------------------------
// void addMemVarCorrCurvilinear( Sarray& a_X, Sarray& a_Y, Sarray& a_Z, double t,
//                                                       Sarray& alpha, double omegaVE, double dt, double omega, double phase, double c )
// {
//    double x, y, z, dto=omegaVE*dt;
//    double cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

//    double forces[3];
//   double t1;
//   double t11;
//   double t12;
//   double t121;
//   double t123;
//   double t124;
//   double t125;
//   double t127;
//   double t129;
//   double t13;
//   double t130;
//   double t134;
//   double t137;
//   double t15;
//   double t16;
//   double t17;
//   double t18;
//   double t19;
//   double t2;
//   double t20;
//   double t21;
//   double t22;
//   double t24;
//   double t27;
//   double t28;
//   double t3;
//   double t30;
//   double t31;
//   double t34;
//   double t35;
//   double t36;
//   double t39;
//   double t4;
//   double t40;
//   double t41;
//   double t42;
//   double t44;
//   double t48;
//   double t5;
//   double t57;
//   double t7;
//   double t70;
//   double t71;
//   double t72;
//   double t74;
//   double t75;
//   double t76;
//   double t78;
//   double t8;
//   double t81;
//   double t85;
//   double t86;
//   double t9;
//   double t91;
//   double t94;
//   double t99;

//   for( int k=a_X.m_kb ; k<=a_X.m_ke; k++ )
//     for( int j=a_X.m_jb ; j<=a_X.m_je; j++ )
//       for( int i=a_X.m_ib ; i<=a_X.m_ie; i++ )
//       {
// 	x = a_X(i,j,k);
// 	y = a_Y(i,j,k);
// 	z = a_Z(i,j,k);

//         t1 = dt*dt;
//         t2 = omegaVE*omegaVE;
//         t3 = 1/omegaVE;
//         t4 = omega*c;
//         t5 = c*t;
//         t7 = omega*(-t5+x);
//         t8 = -t7-phase;
//         t9 = sin(t8);
//         t11 = omega*x+phase;
//         t12 = sin(t11);
//         t13 = t9*t12;
//         t15 = omega*(-t5+z);
//         t16 = -t15-phase;
//         t17 = cos(t16);
//         t18 = t13*t17;
//         t19 = t4*t18;
//         t20 = cos(t8);
//         t21 = t20*t12;
//         t22 = sin(t16);
//         t24 = t21*t4*t22;
//         t27 = t21*t17;
//         t28 = sin(t7);
//         t30 = omega*y+phase;
//         t31 = sin(t30);
//         t34 = omega*z+phase;
//         t35 = sin(t34);
//         t36 = t28*t31*t35;
//         t39 = omega*omega;
//         t40 = c*c;
//         t41 = t39*t40;
//         t42 = t41*t27;
//         t44 = t41*t13*t22;
//         t48 = cos(t7);
//         t57 = t39*omega*t40*c;
//         forces[0] = t1*(t2*(t3*(-t19-t24)+t27-t36)+2.0*omegaVE*(t3*(-2.0*t42+2.0*
//                                                                     t44)-t19-t24+t4*t48*t31*t35)+t3*(4.0*t21*t22*t57+4.0*t18*t57)-2.0*t42+2.0*t44+
//                         t41*t36)/6.0;
//         t70 = omega*(-t5+y);
//         t71 = -t70-phase;
//         t72 = cos(t71);
//         t74 = cos(t34);
//         t75 = t48*t72*t74;
//         t76 = t4*t75;
//         t78 = sin(t71);
//         t81 = t28*omega*c*t78*t74;
//         t85 = t28*t72*t74;
//         t86 = sin(t70);
//         t91 = t41*t85;
//         t94 = t41*t48*t78*t74;
//         t99 = cos(t70);
//         forces[1] = t1*(t2*(t3*(-t76-t81)+t85-t12*t86*t35)+2.0*omegaVE*(t3*(-2.0*
//                                                                             t91+2.0*t94)-t76-t81+t12*omega*c*t99*t35)+t3*(4.0*t28*t57*t74*t78+4.0*t57*t75)
//                         -2.0*t91+2.0*t94+t12*t39*t40*t86*t35)/6.0;
//         t121 = cos(t11);
//         t123 = cos(t30);
//         t124 = t3*t121*t123;
//         t125 = t4*t17;
//         t127 = t121*t123;
//         t129 = t12*t31;
//         t130 = sin(t15);
//         t134 = t41*t22;
//         t137 = cos(t15);
//         forces[2] = t1*(t2*(-t124*t125-t127*t22-t129*t130)+2.0*omegaVE*(t129*t137*
//                                                                         t4+t124*t134-t125*t127)+t124*t57*t17+t127*t134+t129*t41*t130)/6.0;
//         alpha(1,i,j,k) += cof*forces[0];
//         alpha(2,i,j,k) += cof*forces[1];
//         alpha(3,i,j,k) += cof*forces[2];
//       } // end for*3
//   return;
// } // end function

// NEW June 14, 2017: Full forcing for the memory variables correction stage
//---------------------------- NEW MAPLE GENERATED CODE ----------------------------------
void addMemVarCorr2Cart(float_sw4 zMin, float_sw4 h, float_sw4 t, Sarray &alpha,
                        float_sw4 omegaVE, float_sw4 dt, float_sw4 omega, float_sw4 phase, float_sw4 c )
{
  SW4_MARK_FUNCTION;
   float_sw4 dto=omegaVE*dt;
   float_sw4 cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

// #pragma omp parallel for
//    for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
//       for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
// #pragma ivdep
// #pragma simd
//          for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
//          {
   SView &alphaV = alpha.getview();
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION==6
   using LOCAL_POL = 
  RAJA::KernelPolicy< 
  RAJA::statement::CudaKernelAsync<
    RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			 RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
					      RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								   RAJA::statement::Lambda<0> >>>>>;
#elif SW4_RAJA_VERSION==7

  using LOCAL_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			    RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
						  RAJA::statement::Tile<2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
									RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
											     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
														  RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
																       RAJA::statement::Lambda<0> >>>>>>>>;

#endif


#else
   using LOCAL_POL = DEFAULT_LOOP3;
#endif
	   RAJA::RangeSegment k_range(alpha.m_kb,alpha.m_ke+1);
	   RAJA::RangeSegment j_range(alpha.m_jb,alpha.m_je+1);
	   RAJA::RangeSegment i_range(alpha.m_ib,alpha.m_ie+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
	    float_sw4 x, y, z;
	    float_sw4 forces[3];
	    float_sw4 t1;
	    float_sw4 t10;
	    float_sw4 t100;
	    float_sw4 t11;
	    float_sw4 t123;
	    float_sw4 t125;
	    float_sw4 t126;
	    float_sw4 t127;
	    float_sw4 t128;
	    float_sw4 t129;
	    float_sw4 t13;
	    float_sw4 t130;
	    float_sw4 t131;
	    float_sw4 t132;
	    float_sw4 t133;
	    float_sw4 t136;
	    float_sw4 t139;
	    float_sw4 t14;
	    float_sw4 t15;
	    float_sw4 t16;
	    float_sw4 t17;
	    float_sw4 t18;
	    float_sw4 t19;
	    float_sw4 t2;
	    float_sw4 t20;
	    float_sw4 t22;
	    float_sw4 t24;
	    float_sw4 t25;
	    float_sw4 t26;
	    float_sw4 t28;
	    float_sw4 t29;
	    float_sw4 t3;
	    float_sw4 t32;
	    float_sw4 t33;
	    float_sw4 t34;
	    float_sw4 t35;
	    float_sw4 t36;
	    float_sw4 t39;
	    float_sw4 t40;
	    float_sw4 t41;
	    float_sw4 t42;
	    float_sw4 t44;
	    float_sw4 t48;
	    float_sw4 t5;
	    float_sw4 t57;
	    float_sw4 t6;
	    float_sw4 t7;
	    float_sw4 t71;
	    float_sw4 t72;
	    float_sw4 t73;
	    float_sw4 t75;
	    float_sw4 t76;
	    float_sw4 t77;
	    float_sw4 t79;
	    float_sw4 t82;
	    float_sw4 t84;
	    float_sw4 t86;
	    float_sw4 t87;
	    float_sw4 t89;
	    float_sw4 t9;
	    float_sw4 t92;
	    float_sw4 t95;
            x = (i-1)*h;
            y = (j-1)*h;
            z = zMin + (k-1)*h;

            t1 = 1/omegaVE;
            t2 = omega*c;
            t3 = c*t;
            t5 = omega*(-t3+x);
            t6 = -t5-phase;
            t7 = sin(t6);
            t9 = omega*x+phase;
            t10 = sin(t9);
            t11 = t7*t10;
            t13 = omega*(-t3+z);
            t14 = -t13-phase;
            t15 = cos(t14);
            t16 = t11*t15;
            t17 = t2*t16;
            t18 = cos(t6);
            t19 = t18*t10;
            t20 = sin(t14);
            t22 = t19*t2*t20;
            t24 = t1*(-t17-t22);
            t25 = t19*t15;
            t26 = sin(t5);
            t28 = omega*y+phase;
            t29 = sin(t28);
            t32 = omega*z+phase;
            t33 = sin(t32);
            t34 = t26*t29*t33;
            t35 = dt*dt;
            t36 = omegaVE*omegaVE;
            t39 = omega*omega;
            t40 = c*c;
            t41 = t39*t40;
            t42 = t41*t25;
            t44 = t41*t11*t20;
            t48 = cos(t5);
            t57 = t39*omega*t40*c;
            forces[0] = t24+t25-t34+t35*(t36*(t24+t25-t34)+2.0*omegaVE*(t1*(-2.0*t42+
                                                                            2.0*t44)-t17-t22+t2*t48*t29*t33)+t1*(4.0*t57*t19*t20+4.0*t57*t16)-2.0*t42+2.0*
                                         t44+t41*t34)/6.0;
            t71 = omega*(-t3+y);
            t72 = -t71-phase;
            t73 = cos(t72);
            t75 = cos(t32);
            t76 = t48*t73*t75;
            t77 = t2*t76;
            t79 = sin(t72);
            t82 = t26*omega*c*t79*t75;
            t84 = t1*(-t77-t82);
            t86 = t26*t73*t75;
            t87 = sin(t71);
            t89 = t10*t87*t33;
            t92 = t41*t86;
            t95 = t41*t48*t79*t75;
            t100 = cos(t71);
            forces[1] = t84+t86-t89+t35*(t36*(t84+t86-t89)+2.0*omegaVE*(t1*(-2.0*t92+
                2.0*t95)-t77-t82+t10*omega*c*t100*t33)+t1*(4.0*t57*t26*t79*t75+4.0*t57*t76)-2.0
                *t92+2.0*t95+t10*t39*t40*t87*t33)/6.0;
            t123 = cos(t9);
            t125 = cos(t28);
            t126 = t1*t123*t125;
            t127 = t2*t15;
            t128 = t126*t127;
            t129 = t123*t125;
            t130 = t129*t20;
            t131 = t10*t29;
            t132 = sin(t13);
            t133 = t131*t132;
            t136 = t41*t20;
            t139 = cos(t13);
            forces[2] = -t128-t130-t133+t35*(t36*(-t128-t130-t133)+2.0*omegaVE*(t131*t2
                 *t139+t126*t136-t129*t127)+t126*t57*t15+t129*t136+t131*t41*t132)/6.0;

            alphaV(1,i,j,k) += cof*forces[0];
            alphaV(2,i,j,k) += cof*forces[1];
            alphaV(3,i,j,k) += cof*forces[2];
			       } );// end for*3
   return;
} // end function

//----------------------------------------------------------------------
void addMemVarCorr2Curvilinear( Sarray& a_X, Sarray& a_Y, Sarray& a_Z, float_sw4 t,
				Sarray& alpha, float_sw4 omegaVE, float_sw4 dt, float_sw4 omega,
				float_sw4 phase, float_sw4 c )
{
  SW4_MARK_FUNCTION;
   float_sw4 dto=omegaVE*dt;
   float_sw4 cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

// #pragma omp parallel for
//    for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
//       for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
// #pragma ivdep
// #pragma simd
//          for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
// 	 {
   SView &alphaV = alpha.getview();
   SView &a_XV = a_X.getview();
   SView &a_YV = a_Y.getview();
   SView &a_ZV = a_Z.getview();
#ifdef ENABLE_CUDA

#if SW4_RAJA_VERSION==6

   using LOCAL_POL = 
     RAJA::KernelPolicy< 
     RAJA::statement::CudaKernelAsync<
       RAJA::statement::For<0, RAJA::cuda_threadblock_exec<4>, 
			    RAJA::statement::For<1, RAJA::cuda_threadblock_exec<4>, 
						 RAJA::statement::For<2, RAJA::cuda_threadblock_exec<64>,
								      RAJA::statement::Lambda<0> >>>>>;

#elif SW4_RAJA_VERSION==7

   using LOCAL_POL = RAJA::KernelPolicy<
    RAJA::statement::CudaKernel<
      RAJA::statement::Tile<0, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_y_loop,
			    RAJA::statement::Tile<1, RAJA::statement::tile_fixed<4>, RAJA::cuda_block_x_loop,
						  RAJA::statement::Tile<2, RAJA::statement::tile_fixed<64>, RAJA::cuda_block_z_loop,
									RAJA::statement::For<0, RAJA::cuda_thread_y_direct,
											     RAJA::statement::For<1, RAJA::cuda_thread_x_direct,
														  RAJA::statement::For<2, RAJA::cuda_thread_z_direct,
																       RAJA::statement::Lambda<0> >>>>>>>>;

#endif

#else
   using LOCAL_POL = DEFAULT_LOOP3;
#endif
	   RAJA::RangeSegment k_range(alpha.m_kb,alpha.m_ke+1);
	   RAJA::RangeSegment j_range(alpha.m_jb,alpha.m_je+1);
	   RAJA::RangeSegment i_range(alpha.m_ib,alpha.m_ie+1);
   RAJA::kernel<LOCAL_POL>(
			       RAJA::make_tuple(k_range,j_range,i_range),
			       [=]RAJA_DEVICE (int k, int j,int i) {
	    float_sw4 x, y, z;
	    float_sw4 forces[3];
	    float_sw4 t1;
	    float_sw4 t10;
	    float_sw4 t100;
	    float_sw4 t11;
	    float_sw4 t123;
	    float_sw4 t125;
	    float_sw4 t126;
	    float_sw4 t127;
	    float_sw4 t128;
	    float_sw4 t129;
	    float_sw4 t13;
	    float_sw4 t130;
	    float_sw4 t131;
	    float_sw4 t132;
	    float_sw4 t133;
	    float_sw4 t136;
	    float_sw4 t139;
	    float_sw4 t14;
	    float_sw4 t15;
	    float_sw4 t16;
	    float_sw4 t17;
	    float_sw4 t18;
	    float_sw4 t19;
	    float_sw4 t2;
	    float_sw4 t20;
	    float_sw4 t22;
	    float_sw4 t24;
	    float_sw4 t25;
	    float_sw4 t26;
	    float_sw4 t28;
	    float_sw4 t29;
	    float_sw4 t3;
	    float_sw4 t32;
	    float_sw4 t33;
	    float_sw4 t34;
	    float_sw4 t35;
	    float_sw4 t36;
	    float_sw4 t39;
	    float_sw4 t40;
	    float_sw4 t41;
	    float_sw4 t42;
	    float_sw4 t44;
	    float_sw4 t48;
	    float_sw4 t5;
	    float_sw4 t57;
	    float_sw4 t6;
	    float_sw4 t7;
	    float_sw4 t71;
	    float_sw4 t72;
	    float_sw4 t73;
	    float_sw4 t75;
	    float_sw4 t76;
	    float_sw4 t77;
	    float_sw4 t79;
	    float_sw4 t82;
	    float_sw4 t84;
	    float_sw4 t86;
	    float_sw4 t87;
	    float_sw4 t89;
	    float_sw4 t9;
	    float_sw4 t92;
	    float_sw4 t95;
	    x = a_XV(i,j,k);
	    y = a_YV(i,j,k);
	    z = a_ZV(i,j,k);
            t1 = 1/omegaVE;
            t2 = omega*c;
            t3 = c*t;
            t5 = omega*(-t3+x);
            t6 = -t5-phase;
            t7 = sin(t6);
            t9 = omega*x+phase;
            t10 = sin(t9);
            t11 = t7*t10;
            t13 = omega*(-t3+z);
            t14 = -t13-phase;
            t15 = cos(t14);
            t16 = t11*t15;
            t17 = t2*t16;
            t18 = cos(t6);
            t19 = t18*t10;
            t20 = sin(t14);
            t22 = t19*t2*t20;
            t24 = t1*(-t17-t22);
            t25 = t19*t15;
            t26 = sin(t5);
            t28 = omega*y+phase;
            t29 = sin(t28);
            t32 = omega*z+phase;
            t33 = sin(t32);
            t34 = t26*t29*t33;
            t35 = dt*dt;
            t36 = omegaVE*omegaVE;
            t39 = omega*omega;
            t40 = c*c;
            t41 = t39*t40;
            t42 = t41*t25;
            t44 = t41*t11*t20;
            t48 = cos(t5);
            t57 = t39*omega*t40*c;
            forces[0] = t24+t25-t34+t35*(t36*(t24+t25-t34)+2.0*omegaVE*(t1*(-2.0*t42+
                                                                            2.0*t44)-t17-t22+t2*t48*t29*t33)+t1*(4.0*t57*t19*t20+4.0*t57*t16)-2.0*t42+2.0*
                                         t44+t41*t34)/6.0;
            t71 = omega*(-t3+y);
            t72 = -t71-phase;
            t73 = cos(t72);
            t75 = cos(t32);
            t76 = t48*t73*t75;
            t77 = t2*t76;
            t79 = sin(t72);
            t82 = t26*omega*c*t79*t75;
            t84 = t1*(-t77-t82);
            t86 = t26*t73*t75;
            t87 = sin(t71);
            t89 = t10*t87*t33;
            t92 = t41*t86;
            t95 = t41*t48*t79*t75;
            t100 = cos(t71);
            forces[1] = t84+t86-t89+t35*(t36*(t84+t86-t89)+2.0*omegaVE*(t1*(-2.0*t92+
                                                                            2.0*t95)-t77-t82+t10*omega*c*t100*t33)+t1*(4.0*t57*t26*t79*t75+4.0*t57*t76)-2.0
                                         *t92+2.0*t95+t10*t39*t40*t87*t33)/6.0;
            t123 = cos(t9);
            t125 = cos(t28);
            t126 = t1*t123*t125;
            t127 = t2*t15;
            t128 = t126*t127;
            t129 = t123*t125;
            t130 = t129*t20;
            t131 = t10*t29;
            t132 = sin(t13);
            t133 = t131*t132;
            t136 = t41*t20;
            t139 = cos(t13);
            forces[2] = -t128-t130-t133+t35*(t36*(-t128-t130-t133)+2.0*omegaVE*(t131*t2
                                                                                *t139+t126*t136-t129*t127)+t126*t57*t15+t129*t136+t131*t41*t132)/6.0;

            alphaV(1,i,j,k) += cof*forces[0];
            alphaV(2,i,j,k) += cof*forces[1];
            alphaV(3,i,j,k) += cof*forces[2];
			       } );// end for*3
  return;
} // end function
