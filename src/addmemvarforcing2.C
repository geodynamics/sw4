#include "Sarray.h"
//-----------------From WPP------------------------------------------------------
void  addmemvarforcing2( double zMin, double h, double t, Sarray &alpha,
                                                  double omegaVE, double dt ,double omega, double phase, double c)
{
  double x, y, z;
  // double omega =m_omega;
  // double phase =m_phase;
  // double c=m_c;
  double cof = 2*dt/(1+omegaVE*dt);
  for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
    for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
      for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
      {
	x = (i-1)*h;
	y = (j-1)*h;
	z = zMin + (k-1)*h;
	{
	   double t1;
	   double t12;
	   double t13;
	   double t14;
	   double t17;
	   double t18;
	   double t19;
	   double t24;
	   double t26;
	   double t27;
	   double t3;
	   double t30;
	   double t31;
	   double t35;
	   double t38;
	   double t39;
	   double t4;
	   double t40;
	   double t42;
	   double t45;
	   double t5;
	   double t52;
	   double t57;
	   double t58;
	   double t59;
	   double t65;
	   double t8;
	   double t9;
           double forces[3];
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
