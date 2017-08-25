#include "Sarray.h"
// from WPP
//-----------------------------------------------------------------------
void  addmemvarforcing2( float_sw4 zMin, float_sw4 h, float_sw4 t, Sarray &alpha,
                         float_sw4 omegaVE, float_sw4 dt, float_sw4 omega, float_sw4 phase, float_sw4 c )
{
   //  float_sw4 x, y, z;
  // double omega =m_omega;
  // double phase =m_phase;
  // double c=m_c;
  float_sw4 cof = 2*dt/(1+omegaVE*dt);
#pragma omp parallel for
  for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
    for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
#pragma ivdep
#pragma simd
      for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
      {
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
              alpha(1,i,j,k) += cof*forces[0];
	      t35 = cos(t3);
	      t38 = omega*(y-t1);
	      t39 = -t38-phase;
	      t40 = cos(t39);
	      t42 = cos(t30);
	      t45 = sin(t39);
	      t52 = sin(t38);
	      forces[1] = -t35*omega*c*t40*t42-t24*t45*omega*c*t42+omegaVE*(t24*t40*t42-
									    t9*t52*t31);
              alpha(2,i,j,k) += cof*forces[1];
	      t57 = cos(t8);
	      t58 = cos(t26);
	      t59 = t57*t58;
	      t65 = sin(t12);
	      forces[2] = -t59*t14*omega*c+omegaVE*(-t59*t19-t9*t27*t65);
              alpha(3,i,j,k) += cof*forces[2];
	   }
	}
      }
}
