#include "Sarray.h"
//-----------------From WPP------------------------------------------------------
void  addMemVarPredCart( double zMin, double h, double t, Sarray &alpha,
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

//-----------------------FROM WPP------------------------------------------------
void addMemVarPredCurvilinear( Sarray& a_X, Sarray& a_Y, Sarray& a_Z, double t,
                                                      Sarray& alpha, double omegaVE, double dt, double omega, double phase, double c )
{
  double x, y, z;
  // double omega =m_omega;
  // double phase =m_phase;
  // double c=m_c;
  double cof = 2*dt/(1+omegaVE*dt);
  for( int k=a_X.m_kb ; k<=a_X.m_ke; k++ )
    for( int j=a_X.m_jb ; j<=a_X.m_je; j++ )
      for( int i=a_X.m_ib ; i<=a_X.m_ie; i++ )
      {
	x = a_X(i,j,k);
	y = a_Y(i,j,k);
	z = a_Z(i,j,k);
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

//---------------------------- NEW MAPLE GENERATED CODE ----------------------------------
void addMemVarCorrCart(double zMin, double h, double t, Sarray &alpha,
                       double omegaVE, double dt, double omega, double phase, double c )
{
   double x, y, z, dto=omegaVE*dt;
   double cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

   double forces[3];
  double t1;
  double t11;
  double t12;
  double t121;
  double t123;
  double t124;
  double t125;
  double t127;
  double t129;
  double t13;
  double t130;
  double t134;
  double t137;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t24;
  double t27;
  double t28;
  double t3;
  double t30;
  double t31;
  double t34;
  double t35;
  double t36;
  double t39;
  double t4;
  double t40;
  double t41;
  double t42;
  double t44;
  double t48;
  double t5;
  double t57;
  double t7;
  double t70;
  double t71;
  double t72;
  double t74;
  double t75;
  double t76;
  double t78;
  double t8;
  double t81;
  double t85;
  double t86;
  double t9;
  double t91;
  double t94;
  double t99;

  for( int k=alpha.m_kb ; k<= alpha.m_ke; k++ )
    for( int j=alpha.m_jb ; j<= alpha.m_je; j++ )
      for( int i=alpha.m_ib ; i<= alpha.m_ie; i++ )
      {
	x = (i-1)*h;
	y = (j-1)*h;
	z = zMin + (k-1)*h;

        t1 = dt*dt;
        t2 = omegaVE*omegaVE;
        t3 = 1/omegaVE;
        t4 = omega*c;
        t5 = c*t;
        t7 = omega*(-t5+x);
        t8 = -t7-phase;
        t9 = sin(t8);
        t11 = omega*x+phase;
        t12 = sin(t11);
        t13 = t9*t12;
        t15 = omega*(-t5+z);
        t16 = -t15-phase;
        t17 = cos(t16);
        t18 = t13*t17;
        t19 = t4*t18;
        t20 = cos(t8);
        t21 = t20*t12;
        t22 = sin(t16);
        t24 = t21*t4*t22;
        t27 = t21*t17;
        t28 = sin(t7);
        t30 = omega*y+phase;
        t31 = sin(t30);
        t34 = omega*z+phase;
        t35 = sin(t34);
        t36 = t28*t31*t35;
        t39 = omega*omega;
        t40 = c*c;
        t41 = t39*t40;
        t42 = t41*t27;
        t44 = t41*t13*t22;
        t48 = cos(t7);
        t57 = t39*omega*t40*c;
        forces[0] = t1*(t2*(t3*(-t19-t24)+t27-t36)+2.0*omegaVE*(t3*(-2.0*t42+2.0*
                                                                    t44)-t19-t24+t4*t48*t31*t35)+t3*(4.0*t21*t22*t57+4.0*t18*t57)-2.0*t42+2.0*t44+
                        t41*t36)/6.0;
        t70 = omega*(-t5+y);
        t71 = -t70-phase;
        t72 = cos(t71);
        t74 = cos(t34);
        t75 = t48*t72*t74;
        t76 = t4*t75;
        t78 = sin(t71);
        t81 = t28*omega*c*t78*t74;
        t85 = t28*t72*t74;
        t86 = sin(t70);
        t91 = t41*t85;
        t94 = t41*t48*t78*t74;
        t99 = cos(t70);
        forces[1] = t1*(t2*(t3*(-t76-t81)+t85-t12*t86*t35)+2.0*omegaVE*(t3*(-2.0*
                                                                            t91+2.0*t94)-t76-t81+t12*omega*c*t99*t35)+t3*(4.0*t28*t57*t74*t78+4.0*t57*t75)
                        -2.0*t91+2.0*t94+t12*t39*t40*t86*t35)/6.0;
        t121 = cos(t11);
        t123 = cos(t30);
        t124 = t3*t121*t123;
        t125 = t4*t17;
        t127 = t121*t123;
        t129 = t12*t31;
        t130 = sin(t15);
        t134 = t41*t22;
        t137 = cos(t15);
        forces[2] = t1*(t2*(-t124*t125-t127*t22-t129*t130)+2.0*omegaVE*(t129*t137*
                                                                        t4+t124*t134-t125*t127)+t124*t57*t17+t127*t134+t129*t41*t130)/6.0;
        alpha(1,i,j,k) += cof*forces[0];
        alpha(2,i,j,k) += cof*forces[1];
        alpha(3,i,j,k) += cof*forces[2];
      } // end for*3
  return;
} // end function

//----------------------------------------------------------------------
void addMemVarCorrCurvilinear( Sarray& a_X, Sarray& a_Y, Sarray& a_Z, double t,
                                                      Sarray& alpha, double omegaVE, double dt, double omega, double phase, double c )
{
   double x, y, z, dto=omegaVE*dt;
   double cof = 1.0/( 1.0/2 + 1.0/(2*dto) + dto/4 + dto*dto/12 );

   double forces[3];
  double t1;
  double t11;
  double t12;
  double t121;
  double t123;
  double t124;
  double t125;
  double t127;
  double t129;
  double t13;
  double t130;
  double t134;
  double t137;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t24;
  double t27;
  double t28;
  double t3;
  double t30;
  double t31;
  double t34;
  double t35;
  double t36;
  double t39;
  double t4;
  double t40;
  double t41;
  double t42;
  double t44;
  double t48;
  double t5;
  double t57;
  double t7;
  double t70;
  double t71;
  double t72;
  double t74;
  double t75;
  double t76;
  double t78;
  double t8;
  double t81;
  double t85;
  double t86;
  double t9;
  double t91;
  double t94;
  double t99;

  for( int k=a_X.m_kb ; k<=a_X.m_ke; k++ )
    for( int j=a_X.m_jb ; j<=a_X.m_je; j++ )
      for( int i=a_X.m_ib ; i<=a_X.m_ie; i++ )
      {
	x = a_X(i,j,k);
	y = a_Y(i,j,k);
	z = a_Z(i,j,k);

        t1 = dt*dt;
        t2 = omegaVE*omegaVE;
        t3 = 1/omegaVE;
        t4 = omega*c;
        t5 = c*t;
        t7 = omega*(-t5+x);
        t8 = -t7-phase;
        t9 = sin(t8);
        t11 = omega*x+phase;
        t12 = sin(t11);
        t13 = t9*t12;
        t15 = omega*(-t5+z);
        t16 = -t15-phase;
        t17 = cos(t16);
        t18 = t13*t17;
        t19 = t4*t18;
        t20 = cos(t8);
        t21 = t20*t12;
        t22 = sin(t16);
        t24 = t21*t4*t22;
        t27 = t21*t17;
        t28 = sin(t7);
        t30 = omega*y+phase;
        t31 = sin(t30);
        t34 = omega*z+phase;
        t35 = sin(t34);
        t36 = t28*t31*t35;
        t39 = omega*omega;
        t40 = c*c;
        t41 = t39*t40;
        t42 = t41*t27;
        t44 = t41*t13*t22;
        t48 = cos(t7);
        t57 = t39*omega*t40*c;
        forces[0] = t1*(t2*(t3*(-t19-t24)+t27-t36)+2.0*omegaVE*(t3*(-2.0*t42+2.0*
                                                                    t44)-t19-t24+t4*t48*t31*t35)+t3*(4.0*t21*t22*t57+4.0*t18*t57)-2.0*t42+2.0*t44+
                        t41*t36)/6.0;
        t70 = omega*(-t5+y);
        t71 = -t70-phase;
        t72 = cos(t71);
        t74 = cos(t34);
        t75 = t48*t72*t74;
        t76 = t4*t75;
        t78 = sin(t71);
        t81 = t28*omega*c*t78*t74;
        t85 = t28*t72*t74;
        t86 = sin(t70);
        t91 = t41*t85;
        t94 = t41*t48*t78*t74;
        t99 = cos(t70);
        forces[1] = t1*(t2*(t3*(-t76-t81)+t85-t12*t86*t35)+2.0*omegaVE*(t3*(-2.0*
                                                                            t91+2.0*t94)-t76-t81+t12*omega*c*t99*t35)+t3*(4.0*t28*t57*t74*t78+4.0*t57*t75)
                        -2.0*t91+2.0*t94+t12*t39*t40*t86*t35)/6.0;
        t121 = cos(t11);
        t123 = cos(t30);
        t124 = t3*t121*t123;
        t125 = t4*t17;
        t127 = t121*t123;
        t129 = t12*t31;
        t130 = sin(t15);
        t134 = t41*t22;
        t137 = cos(t15);
        forces[2] = t1*(t2*(-t124*t125-t127*t22-t129*t130)+2.0*omegaVE*(t129*t137*
                                                                        t4+t124*t134-t125*t127)+t124*t57*t17+t127*t134+t129*t41*t130)/6.0;
        alpha(1,i,j,k) += cof*forces[0];
        alpha(2,i,j,k) += cof*forces[1];
        alpha(3,i,j,k) += cof*forces[2];
      } // end for*3
  return;
} // end function
