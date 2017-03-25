#include "EW.h"
#include "F77_FUNC.h"

extern "C" {
   void F77_FUNC(dgesv,DGESV)(int*,int*,double*,int*,int*,double*,int*,int*);
}

void EW::bcfreesurfcurvani_ci( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			       int nz, float_sw4* u, float_sw4* c, int side, float_sw4 sbop[5], 
			       float_sw4* bforce5, float_sw4* bforce6, float_sw4* strx, float_sw4* stry )
{
   const float_sw4 d4a = 2.0/3.0;
   const float_sw4 d4b = -1.0/12.0;
   const size_t ni  = ilast-ifirst+1;
   const size_t nij = ni*(jlast-jfirst+1);
   const size_t npts = nij*(klast-kfirst+1);
   double s0i = 1/sbop[0];
   // side is fortran enumeration 1,2,3,4,5,6
   if( side == 5 )
   {
      int k=1;
#pragma omp parallel for 
      for( int j=jfirst+2 ; j <= jlast-2 ; j++ )
	 for( int i=ifirst+2 ; i <= ilast-2 ; i++ )
	 {
	    size_t qq = i-ifirst+ni*(j-jfirst);
	    size_t ind = i-ifirst+ni*(j-jfirst)+nij*(k-kfirst);
	    float_sw4 du = strx[i-ifirst]*(d4a*(u[ind+1]-u[ind-1])+d4b*(u[ind+2]-u[ind-2]));
	    float_sw4 dv = strx[i-ifirst]*(d4a*(u[npts+ind+1]-u[npts+ind-1])+d4b*(u[npts+ind+2]-u[npts+ind-2]));
	    float_sw4 dw = strx[i-ifirst]*(d4a*(u[2*npts+ind+1]-u[2*npts+ind-1])+d4b*(u[2*npts+ind+2]-u[2*npts+ind-2]));

	    float_sw4 rhs1 = c[ind+27*npts]*du+c[ind+30*npts]*dv+c[ind+33*npts]*dw;
	    float_sw4 rhs2 = c[ind+28*npts]*du+c[ind+31*npts]*dv+c[ind+34*npts]*dw;
	    float_sw4 rhs3 = c[ind+29*npts]*du+c[ind+32*npts]*dv+c[ind+35*npts]*dw;

	    du = stry[j-jfirst]*(d4a*(u[ind+ni]  -u[ind-ni])+
					     d4b*(u[ind+2*ni]-u[ind-2*ni]));
	    dv = stry[j-jfirst]*(d4a*(u[npts+ind+ni]  -u[npts+ind-ni])+
					     d4b*(u[npts+ind+2*ni]-u[npts+ind-2*ni]));
	    dw = stry[j-jfirst]*(d4a*(u[2*npts+ind+ni  ]-u[2*npts+ind-ni  ])+
					     d4b*(u[2*npts+ind+2*ni]-u[2*npts+ind-2*ni]));
	    rhs1 += c[ind+36*npts]*du+c[ind+39*npts]*dv+c[ind+42*npts]*dw;
	    rhs2 += c[ind+37*npts]*du+c[ind+40*npts]*dv+c[ind+43*npts]*dw;
	    rhs3 += c[ind+38*npts]*du+c[ind+41*npts]*dv+c[ind+44*npts]*dw;

	    du=dv=dw=0;
	    for( int w=1 ; w <= 4 ; w++ )
	    {
	       du += sbop[w]*u[       ind+nij*(w-1)];
	       dv += sbop[w]*u[npts  +ind+nij*(w-1)];
	       dw += sbop[w]*u[2*npts+ind+nij*(w-1)];
	    }
	    rhs1 += c[ind+12*npts]*du+c[ind+13*npts]*dv+c[ind+14*npts]*dw-bforce5[  3*qq];
	    rhs2 += c[ind+13*npts]*du+c[ind+15*npts]*dv+c[ind+16*npts]*dw-bforce5[1+3*qq];
	    rhs3 += c[ind+14*npts]*du+c[ind+16*npts]*dv+c[ind+17*npts]*dw-bforce5[2+3*qq];
 // Solve system for ghost point values
	    float_sw4 x[3] = {rhs1,rhs2,rhs3};
	    float_sw4 a[9];
	    a[0] = c[ind+12*npts];
	    a[1] = c[ind+13*npts];
	    a[2] = c[ind+14*npts];
	    a[3] = c[ind+13*npts];
	    a[4] = c[ind+15*npts];
	    a[5] = c[ind+16*npts];
	    a[6] = c[ind+14*npts];
	    a[7] = c[ind+16*npts];
	    a[8] = c[ind+17*npts];
	    int dim=3, one=1, info=0, ipiv[3];
	    F77_FUNC(dgesv,DGESV)(&dim, &one, a, &dim, ipiv, x, &dim, &info );
	    if( info != 0 )
	       cout << "ERROR in bcfreesurfcurvanic_ci, call to DGESV returned info " << info << endl;
	    u[       ind-nij] = -s0i*x[0];
	    u[npts  +ind-nij] = -s0i*x[1];
	    u[2*npts+ind-nij] = -s0i*x[2];
	 }
   }
   else if( side == 6 )
   {
      int k=nz;
#pragma omp parallel for 
      for( int j=jfirst+2 ; j <= jlast-2 ; j++ )
	 for( int i=ifirst+2 ; i <= ilast-2 ; i++ )
	 {
	    size_t qq = i-ifirst+ni*(j-jfirst);
	    size_t ind = i-ifirst+ni*(j-jfirst)+nij*(k-kfirst);
	    float_sw4 du = strx[i-ifirst]*(d4a*(u[ind+1]-u[ind-1])+d4b*(u[ind+2]-u[ind-2]));
	    float_sw4 dv = strx[i-ifirst]*(d4a*(u[npts+ind+1]-u[npts+ind-1])+d4b*(u[npts+ind+2]-u[npts+ind-2]));
	    float_sw4 dw = strx[i-ifirst]*(d4a*(u[2*npts+ind+1]-u[2*npts+ind-1])+d4b*(u[2*npts+ind+2]-u[2*npts+ind-2]));

	    float_sw4 rhs1 = c[ind+27*npts]*du+c[ind+30*npts]*dv+c[ind+33*npts]*dw;
	    float_sw4 rhs2 = c[ind+28*npts]*du+c[ind+31*npts]*dv+c[ind+34*npts]*dw;
	    float_sw4 rhs3 = c[ind+29*npts]*du+c[ind+32*npts]*dv+c[ind+35*npts]*dw;

	    du = stry[j-jfirst]*(d4a*(u[ind+ni]  -u[ind-ni])+
			     d4b*(u[ind+2*ni]-u[ind-2*ni]));
	    dv = stry[j-jfirst]*(d4a*(u[npts+ind+ni]  -u[npts+ind-ni])+
			     d4b*(u[npts+ind+2*ni]-u[npts+ind-2*ni]));
	    dw = stry[j-jfirst]*(d4a*(u[2*npts+ind+ni  ]-u[2*npts+ind-ni  ])+
			     d4b*(u[2*npts+ind+2*ni]-u[2*npts+ind-2*ni]));
	    rhs1 += c[ind+36*npts]*du+c[ind+39*npts]*dv+c[ind+42*npts]*dw;
	    rhs2 += c[ind+37*npts]*du+c[ind+40*npts]*dv+c[ind+43*npts]*dw;
	    rhs3 += c[ind+38*npts]*du+c[ind+41*npts]*dv+c[ind+44*npts]*dw;

	    du=dv=dw=0;
	    for( int w=1 ; w <= 4 ; w++ )
	    {
	       du -= sbop[w]*u[       ind-nij*(w-1)];
	       dv -= sbop[w]*u[npts  +ind-nij*(w-1)];
	       dw -= sbop[w]*u[2*npts+ind-nij*(w-1)];
	    }
	    rhs1 += c[ind+12*npts]*du+c[ind+13*npts]*dv+c[ind+14*npts]*dw-bforce6[  3*qq];
	    rhs2 += c[ind+13*npts]*du+c[ind+15*npts]*dv+c[ind+16*npts]*dw-bforce6[1+3*qq];
	    rhs3 += c[ind+14*npts]*du+c[ind+16*npts]*dv+c[ind+17*npts]*dw-bforce6[2+3*qq];

 // Solve system for ghost point values
	    float_sw4 x[3] = {rhs1,rhs2,rhs3};
	    float_sw4 a[9];
	    a[0] = c[ind+12*npts];
	    a[1] = c[ind+13*npts];
	    a[2] = c[ind+14*npts];
	    a[3] = c[ind+13*npts];
	    a[4] = c[ind+15*npts];
	    a[5] = c[ind+16*npts];
	    a[6] = c[ind+14*npts];
	    a[7] = c[ind+16*npts];
	    a[8] = c[ind+17*npts];
	    int dim=3, one=1, info=0, ipiv[3];
	    F77_FUNC(dgesv,DGESV)(&dim, &one, a, &dim, ipiv, x, &dim, &info );
	    if( info != 0 )
	       cout << "ERROR in bcfreesurfcurvani_ci, call to DGESV returned info " << info << endl;
	    u[       ind+nij] = s0i*x[0];
	    u[npts  +ind+nij] = s0i*x[1];
	    u[2*npts+ind+nij] = s0i*x[2];
	 }
   }
   else
      cout << "ERROR in bcfreesurfcurvani_ci, b.c. not implemented for side= " << side << endl;
}
