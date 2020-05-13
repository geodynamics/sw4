#include <iostream>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <random>
#include <mpi.h>
#ifdef ENABLE_FFTW
#include <fftw3-mpi.h>
#endif

#include "EW.h"
#include "RandomizedMaterial.h"
#include "AllDims.h"

MPI_Datatype get_mpi_datatype( double* var ) {return MPI_DOUBLE;}
MPI_Datatype get_mpi_datatype( float* var ) {return MPI_FLOAT;}
MPI_Datatype get_mpi_datatype( int* var ) {return MPI_INT;}
MPI_Datatype get_mpi_datatype( std::complex<double>* var ) {return MPI_CXX_DOUBLE_COMPLEX;}
MPI_Datatype get_mpi_datatype( std::complex<float>* var ) {return MPI_CXX_FLOAT_COMPLEX;}


//-----------------------------------------------------------------------
RandomizedMaterial::RandomizedMaterial( EW * a_ew, float_sw4 zmin, float_sw4 zmax,
					float_sw4 corrlen, float_sw4 corrlenz, 
					float_sw4 hurst, float_sw4 sigma,
					unsigned int seed )
{
   mEW = a_ew;
   float_sw4 bbox[6];
   a_ew->getGlobalBoundingBox( bbox );

   float_sw4 global_xmax = bbox[1];
   float_sw4 global_ymax = bbox[3];
   float_sw4 global_zmin = bbox[4];
   float_sw4 global_zmax = bbox[5];
   if( zmin < global_zmin )
      zmin = global_zmin;
   if( zmax > global_zmax )
      zmax = global_zmax;

   m_corrlen  = corrlen;
   m_corrlenz = corrlenz;
   m_zmin  = zmin;
   m_zmax  = zmax;
   m_hurst = hurst;
   m_sigma = sigma;
   m_seed  = seed;
   m_vsmax = 1e38;
   
// Determine discretization based on correlation length.
   float_sw4 ppcl = 20; // grid points per correlation length
   m_nig = ppcl*(global_xmax)/corrlen;
   m_njg = ppcl*(global_ymax)/corrlen;
   m_nkg = ppcl*(zmax-zmin)/corrlenz;

   m_hh= global_xmax/(m_nig-1);
   m_hv= (zmax-zmin)/(m_nkg-1);
   
// Find finest grid intersecting this block, and limit the grid spacing to
// not be finer than spacing in the computational grid.
   int g=a_ew->mNumberOfGrids-1;

// Force spacings to be multiple of h to make the domain size exactly the same
// as for the computational grid.
   float_sw4 htmp = m_hh/a_ew->mGridSize[g];
   int ratio = static_cast<int>(round(htmp));
   m_hh = ratio*a_ew->mGridSize[g];
   m_nig = static_cast<int>(round(global_xmax/m_hh)+1);
   m_njg = static_cast<int>(round(global_ymax/m_hh)+1);

   bool tmptest=false;
   if( tmptest )
   {
   // temp testing, 1 grid case with equal spacing
      m_hh  = a_ew->mGridSize[g];
      m_nig = a_ew->m_global_nx[g];
      m_njg = a_ew->m_global_ny[g];
      m_hv = a_ew->mGridSize[g];
      m_nkg = a_ew->m_global_nz[g];
   }
   else
   {
      bool found = false;
      while( g >= 0 && !found )
      {
	 if( zmax > a_ew->m_zmin[g] )
	    found = true;
	 else
	    g--;
      }
      CHECK_INPUT( g >=0,"Error in RandomizeMaterial::RandomizeMaterial, g = "<< g  << "\n");

      if( m_hh < a_ew->mGridSize[g] )
      {
	 m_hh  = a_ew->mGridSize[g];
	 m_nig = a_ew->m_global_nx[g];
	 m_njg = a_ew->m_global_ny[g];
      }
      if( m_hv < a_ew->mGridSize[g] )
      {
	 m_hv = a_ew->mGridSize[g];
	 m_nkg = static_cast<int>(round((zmax-zmin)/m_hv+1));
	 m_hv = (zmax-zmin)/(m_nkg-1);
      }
      if( a_ew->getRank() == 0 )
      {
	 cout << "RANDMTRL spacing hh " << m_hh << " hv " << m_hv << endl;
	 cout << "RANDMTRL npts  ni " << m_nig << " nj " << m_njg << " nk " << m_nkg << endl;
      }

   }   
   int per2d[2], coord2d[2];
   MPI_Cart_get( a_ew->m_cartesian_communicator, 2, m_nproc2d, per2d, coord2d );
 //   cout << "RANDMTRL myrank " << a_ew->getRank() << " " << m_nproc2d[0] << " " << m_nproc2d[1] << endl;
   // Scale correlation lengths with domain lengths
   float_sw4 scaledcorrlenx = m_corrlen/global_xmax;
   float_sw4 scaledcorrleny = m_corrlen/global_ymax;
   float_sw4 scaledcorrlenz = m_corrlenz/(zmax-zmin);

   gen_random_mtrl_fft3d_fftw( m_nig, m_njg, m_nkg, scaledcorrlenx, scaledcorrleny, scaledcorrlenz, m_hurst );
   rescale_perturbation();

   // with linear interpolation between fd grid and the random material grid, one ghost point
   // should be enough to allow interpolation without communication.
   repad_sarray( mRndMaterial, 0, 2 );
   //   cout << "RANDMTRL dims padded array " << mRndMaterial.m_jb << " " << mRndMaterial.m_je << endl;
   //   MPI_Barrier(MPI_COMM_WORLD);
}

//-----------------------------------------------------------------------
void RandomizedMaterial::perturb_velocities( int g, Sarray& cs, Sarray& cp, 
					     float_sw4 h, float_sw4 zmin, float_sw4 zmax ) 
{
   //-----------------------------------------------------------------------
   // Input: g - Grid number.
   //        h - Grid spacing on grid g
   //        cs, cp - Unperturbated shear and pressure wave speeds arrays on grid g
   //        zmin, zmax - z limits for grid g
   // Output: cs, cp - Shear and pressure wave speeds on grid g with random perturbation.
   //-----------------------------------------------------------------------         

   if( m_zmax > zmin && zmax > m_zmin )
   {
      //      cout << "intersection, z lims grid block " << zmin << " " << zmax << endl;
      //      cout << "intersection, z lims rand block " << m_zmin << " " << m_zmax << endl;
  // Grid block intersects random material block
      //      bool curvilinear = g == mEW->mNumberOfGrids-1 && mEW->topographyExists(); // NOT verified for several curvilinear grids
      bool curvilinear = g >= mEW->mNumberOfCartesianGrids;
      // Interpolate to sw4 grid
      for( int k=mEW->m_kStartInt[g] ; k <= mEW->m_kEndInt[g] ; k++ )
	 for( int j=mEW->m_jStartInt[g] ; j <= mEW->m_jEndInt[g] ; j++ )
	    for( int i=mEW->m_iStartInt[g] ; i <= mEW->m_iEndInt[g] ; i++ )
	    {
	       float_sw4 x = (i-1)*h, y=(j-1)*h, z= zmin + (k-1)*h;
	       if( curvilinear )
	       {
		  x = mEW->mX[g](i,j,k);
		  y = mEW->mY[g](i,j,k);
		  z = mEW->mZ[g](i,j,k);
	       }
	       if( m_zmin <= z && z <= m_zmax )
	       {
	       int ip = x/m_hh, jp=y/m_hh, kp=(z-m_zmin)/m_hv;
	       if( ip >= mRndMaterial.m_ib && ip <= mRndMaterial.m_ie-1 &&
		   jp >= mRndMaterial.m_jb && jp <= mRndMaterial.m_je-1 &&
		   kp >= mRndMaterial.m_kb && kp <= mRndMaterial.m_ke-1 )
	       {
		  float_sw4 wghi= (x-ip*m_hh)/m_hh, wghj=(y-jp*m_hh)/m_hh, wghk=(z-(m_zmin+kp*m_hv))/m_hv;
		  float_sw4 rndpert =(1-wghk)*((1-wghj)*((1-wghi)*mRndMaterial(ip,jp,  kp)  + wghi*mRndMaterial(ip+1,jp,  kp))  +
					    (wghj) *((1-wghi)*mRndMaterial(ip,jp+1,kp)  + wghi*mRndMaterial(ip+1,jp+1,kp))) +
		     (wghk)*((1-wghj)*((1-wghi)*mRndMaterial(ip,jp,  kp+1)+ wghi*mRndMaterial(ip+1,jp,  kp+1))+
			     (wghj) *((1-wghi)*mRndMaterial(ip,jp+1,kp+1)+ wghi*mRndMaterial(ip+1,jp+1,kp+1)));
		  if( cs(i,j,k) <= m_vsmax )
		  {
		     cs(i,j,k) *= rndpert;
		     cp(i,j,k) *= rndpert;
		  }
	       }
	       else if( ip >= mRndMaterial.m_ib && ip <= mRndMaterial.m_ie &&
			jp >= mRndMaterial.m_jb && jp <= mRndMaterial.m_je &&
			kp >= mRndMaterial.m_kb && kp <= mRndMaterial.m_ke )
	       {
		  float_sw4 rndpert = mRndMaterial(ip,jp,kp);
		  if( cs(i,j,k) <= m_vsmax )
		  {
		     cs(i,j,k) *= rndpert;
		     cp(i,j,k) *= rndpert;
		  }
	       }
	       else
		  CHECK_INPUT(false,"ERROR: index " << ip << " " << jp << " " << kp << " not in material array bounds " <<
			      mRndMaterial.m_ib << " <= ip <= " << mRndMaterial.m_ie << "  " <<
			      mRndMaterial.m_jb << " <= jp <= " << mRndMaterial.m_je << "  " <<
			      mRndMaterial.m_kb << " <= kp <= " << mRndMaterial.m_ke << " y= " << y << " j= " << j<<endl );
                  }
	    }
   }
}

//-----------------------------------------------------------------------
void RandomizedMaterial::perturb_velocities( std::vector<Sarray> & cs, 
					     std::vector<Sarray> & cp ) 
{
   for( int g=0 ; g < cs.size() ; g++ )
   {
      //tmp testing the 1 grid case, with equal grid spacing

    // Note: mRndMaterial is a zero-based array -->  add 1 here 
      if( cs[0].m_ib<=mRndMaterial.m_ib+1 && cs[0].m_ie>=mRndMaterial.m_ie+1 &&
	  cs[0].m_jb<=mRndMaterial.m_jb+1 && cs[0].m_je>=mRndMaterial.m_je+1 &&
	  cs[0].m_kb<=mRndMaterial.m_kb+1 && cs[0].m_ke>=mRndMaterial.m_ke+1 )
      {
	 cout << "DOING MTRL " << endl;
	 for( int k=mRndMaterial.m_kb ; k <= mRndMaterial.m_ke ; k++ )
	    for( int j=mRndMaterial.m_jb ; j <= mRndMaterial.m_je ; j++ )
	       for( int i=mRndMaterial.m_ib ; i <= mRndMaterial.m_ie ; i++ )
	       {
		  if( cs[0](i+1,j+1,k+1) <= m_vsmax )
		  {
		     cs[0](i+1,j+1,k+1) *= mRndMaterial(i,j,k);
		     cp[0](i+1,j+1,k+1) *= mRndMaterial(i,j,k);
		  }
	       }
      }
   }
}

//-----------------------------------------------------------------------
void RandomizedMaterial::gen_random_mtrl_fft3d_fftw( int n1g, int n2g, int n3g, 
						     float_sw4 Lx, float_sw4 Ly, float_sw4 Lz, 
						     float_sw4 hurst )
{
#ifdef ENABLE_FFTW
   int nprocs;
   MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
   AllDims* dimobj = new AllDims( nprocs, 0, n1g-1, 0, n2g-1, 0, n3g-1, 0 );
   int dims[6];
   dimobj->getdims_nopad(dims);
   ptrdiff_t ib1 = dims[0], n1=dims[1]-dims[0]+1;
   complex<float_sw4>* uc = new complex<float_sw4>[dimobj->m_fftw_alloc_local];

   if( m_seed == 0 )
   {
      int fd=open("/dev/urandom",O_RDONLY);
      read(fd,&m_seed,sizeof(unsigned int));
      close(fd);
   }

// 1. Generate Fourier modes and setup FFTW plan 

   get_fourier_modes( uc, n1, ib1, n1g, n2g, n3g, Lx, Ly, Lz, hurst, m_seed );
   fftw_plan plan = fftw_mpi_plan_dft_3d( n1g, n2g, n3g, (fftw_complex*)uc, (fftw_complex*)uc,
					  MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE );
// 2. Enforce symmetries
   int r1=(n1g-1)/2, r2=(n2g-1)/2, r3=(n3g-1)/2;
   int n1h=n1g/2, n2h=n2g/2, n3h=n3g/2;

#define u(k1,k2,k3) uc[k3+n3g*(k2)+n2g*n3g*(k1-ib1)]
   if( ib1 == 0 )
      u(0,0,0) = 0;
   if( (n1g%2 == 0) && (ib1 <= n1h && n1h <= n1-1+ib1) )
      //      uim(n1h,0,0) = 0;
      u(n1h,0,0) = real(u(n1h,0,0));

   // Local symmetries (in processor)
   int ulim=0;
   if( n1g%2 == 0 )
      ulim = n1h;

   for( int k1=0 ; k1 <= ulim ; k1+= n1h )
   {
      if( ib1 <= k1 && k1 <= n1-1+ib1 )
      {
	 if( n2g%2 == 0 )
	    u(k1,n2h,0) = real(u(k1,n2h,0));
	 if( n3g%2 == 0 )
	    u(k1,0,n3h) = real(u(k1,0,n3h));
	 if( n2g%2 == 0 && n3g%2 == 0 )
	    u(k1,n2h,n3h) = real(u(k1,n2h,n3h));
	 for( int k2 = 1 ; k2 <= r2 ; k2++ )
	 {
	    u(k1,n2g-k2,0)   = conj(u(k1,k2,0));
	    if( n3g%2 == 0 )
	       u(k1,n2g-k2,n3h) = conj(u(k1,k2,n3h));
	 }
	 for( int k3 = 1 ; k3 <= r3 ; k3++ )
	 {
	    u(k1,0,  n3g-k3) = conj(u(k1,0,  k3));
	    if( n2g%2 == 0 )
	       u(k1,n2h,n3g-k3) = conj(u(k1,n2h,k3));
	 }
	 for( int k2 = 1 ; k2 <= r2 ; k2++ )
	    for( int k3 = 1 ; k3 <= r3 ; k3++ )
	    {
	       u(k1,n2g-k2,n3g-k3) = conj(u(k1,k2,    k3));
	       u(k1,k2,    n3g-k3) = conj(u(k1,n2g-k2,k3));
	    }
      }
   }

  // Processor interaction
   // get ucc from other proc

   MPI_Request* req = new MPI_Request[n1];
   int tag=349;
   for( int k1=1 ; k1 <= r1 ; k1++ )
   {
      if( ib1 <= k1 && k1 <= ib1+n1-1 )
      {
      // Isend plane k1 to owner of n1g-k1
	 int proc = dimobj->owner_i(n1g-k1);
	 if( proc == -1 )
	 {
	    std::cout << "k1 " << k1 << " will send to " << proc << " who owns "  << n1g-k1 << endl;
	    MPI_Abort(MPI_COMM_WORLD,-1);
	 }
	 if( proc != -1 )
	    MPI_Isend( &uc[n2g*n3g*(k1-ib1)], n2g*n3g, MPI_CXX_DOUBLE_COMPLEX, proc, tag,
		    MPI_COMM_WORLD, &req[k1-ib1] );
	 else
	    cout << "Error finding owner of " << n1g-k1 << " for send" << endl;
      }
   }

   complex<double>* ucc_ = new complex<double>[n2g*n3g];
#define ucc(k2,k3) ucc_[k3+n3g*(k2)]

   for( int k1=n1g-1 ; k1 >= n1g-r1 ; k1--)
   {
      if( ib1 <= k1 && k1 <= ib1+n1-1 )
      {
      //  receive plane to ucc from owner of n1g-k1
	 int proc = dimobj->owner_i(n1g-k1);
	 if( proc == -1 )
	 {
	    std::cout<< " k1 " << k1 << " will receive from " << proc << " who owns "  << n1g-k1 << endl;
	    MPI_Abort(MPI_COMM_WORLD,-1);
	 }
 //	 std::cout << "k1 " << k1 << " will receive from " << proc << " who owns "  << n1g-k1 << endl;
	 MPI_Status status;
	 if( proc != -1 )
	    MPI_Recv( ucc_, n2g*n3g, MPI_CXX_DOUBLE_COMPLEX, proc, tag, MPI_COMM_WORLD, &status );
	 else
	    cout << "Error finding owner of " << n1g-k1 << " for receive" << endl;	    

      // enforce symmetry
	 u(k1,0,0  ) = conj(ucc(0,0));
	 if( n2g%2 == 0 )
	    u(k1,n2h,0) = conj(ucc(n2h,0));
	 if( n3g%2 == 0 )
	    u(k1,0,n3h) = conj(ucc(0,n3h));
	 if( n2g%2 == 0 && n3g%2 == 0 )
	    u(k1,n2h,n3h) = conj(ucc(n2h,n3h));
	 for( int k2=1 ; k2 <= r2 ; k2++ )
	 {
	    u(k1,k2,0)     = conj(ucc(n2g-k2,0));
	    u(k1,n2g-k2,0) = conj(ucc(k2,0));
	    if( n3g%2 == 0 )
	    {
	       u(k1,k2,n3h)     = conj(ucc(n2g-k2,n3h));
	       u(k1,n2g-k2,n3h) = conj(ucc(k2,n3h));
	    }
	 }
	 for( int k3=1 ; k3 <= r3 ; k3++ )
	 {
	    u(k1,0,k3)     = conj(ucc(0,n3g-k3));
	    u(k1,0,n3g-k3) = conj(ucc(0,k3));
	    if( n2g%2 == 0 )
	    {
	       u(k1,n2h,k3)     = conj(ucc(n2h,n3g-k3));
	       u(k1,n2h,n3g-k3) = conj(ucc(n2h,k3));
	    }
	 }
	 for( int k2 = 1 ; k2 <= r2 ; k2++ )
	    for( int k3 = 1 ; k3 <= r3 ; k3++ )
	    {
	       u(k1,k2,    k3    ) = conj(ucc(n2g-k2,n3g-k3));
	       u(k1,n2g-k2,n3g-k3) = conj(ucc(k2,    k3    ));
	       u(k1,k2,    n3g-k3) = conj(ucc(n2g-k2,k3    ));
	       u(k1,n2g-k2,k3    ) = conj(ucc(k2,    n3g-k3));
	    }
      }
   }
   for( int k1=1 ; k1 <= r1 ; k1++ )
   {
      // Wait for Isend to complete 
      if( ib1 <= k1 && k1 <= ib1+n1-1 )
      {
	 MPI_Status status;	 
	 MPI_Wait( &req[k1-ib1], &status );
      }
   }
#undef u
#undef ucc
   delete[] ucc_;
   delete[] req;

// 3. Transform back
   fftw_execute(plan);
   fftw_destroy_plan(plan);

   size_t size=static_cast<size_t>(n1)*n2g*n3g;
   float_sw4* u = new float_sw4[size];
   double imnrm = 0;
   for( size_t i=0 ; i < size ; i++ )
   {
      u[i] = real(uc[i]);
      imnrm = imnrm>fabs(imag(uc[i]))?imnrm:fabs(imag(uc[i]));
   }
   if( imnrm > 1e-12 )
      std::cout << "imnrm = "  << imnrm << std::endl;
   delete[] uc;


   AllDims* sw4dims= mEW->get_fine_alldimobject( );
   AllDims sarobj( sw4dims, 0, n1g-1, 0, n2g-1, 0, n3g-1, 0, 0 );
   //   AllDims sarobj(m_nproc2d[0], m_nproc2d[1], 1, 0, n1g-1, 0, n2g-1, 0, n3g-1, 0, 0 );

   sarobj.getdims_nopad(dims);
   mRndMaterial.define(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
   redistribute_array<float_sw4>( *dimobj, sarobj, u, mRndMaterial.c_ptr() );

   delete[] u;
#else
   cout << "ERROR: Can not generate random material without FFTW" << endl;
#endif
}

//-----------------------------------------------------------------------
void RandomizedMaterial::get_fourier_modes( complex<float_sw4>* uhat, int n1, int ib1, int n1g,
					    int n2, int n3, float_sw4 l1, float_sw4 l2, float_sw4 l3, 
					    float_sw4 hurst, unsigned int seed )
{
   const complex<float_sw4> I(0.0,1.0);
   float_sw4 A0=1; // Amplitude
   int D=3; // For three space dimensions
   float_sw4 A0isq2 = A0/sqrt(2.0);
   float_sw4 hhalf = 0.5*(hurst+D*0.5);
   float_sw4 ll1=l1*l1, ll2=l2*l2, ll3=l3*l3;
   
   default_random_engine generator(seed);
   normal_distribution<float_sw4> ndist(0.0,1.0);

   int r1=(n1g-1)/2, r2=(n2-1)/2, r3=(n3-1)/2;

   for( int k1=ib1 ; k1 <= n1-1+ib1 ; k1++ )
   {
      float_sw4 k1eff=k1;
      if( k1 > r1 )
	 k1eff = k1-n1g;
      for( int k2=0 ; k2 <= n2-1 ; k2++ )
      {
	 float_sw4 k2eff=k2;
	 if( k2 > r2 )
	    k2eff = k2-n2;
	 for( int k3=0 ; k3 <= n3-1 ; k3++ )
	 {
	    float_sw4 k3eff=k3;
	    if( k3 > r3 )
	       k3eff = k3-n3;
	    uhat[k3+n3*k2+n2*n3*(k1-ib1)] = (A0isq2/pow(1+k1eff*k1eff*ll1+k2eff*k2eff*ll2+k3eff*k3eff*ll3,hhalf))
	       *(ndist(generator)+I*ndist(generator));
	    //	       *(1+I);
	 }
      }
   }
}
//-----------------------------------------------------------------------
void RandomizedMaterial::rescale_perturbation()
{
   // Average should be zero by construction. Verify it, and compute s.d.
   float_sw4 locavg=0, locsum2=0;
   size_t n=static_cast<size_t>(mRndMaterial.m_nc)*mRndMaterial.m_ni*
                                mRndMaterial.m_nj*mRndMaterial.m_nk;

   float_sw4* mtrl = mRndMaterial.c_ptr();
   for( size_t ind = 0 ; ind < n ;ind++ )
   {
      locavg  += mtrl[ind];
      locsum2 += mtrl[ind]*mtrl[ind];
   }
   float_sw4 avg;
   MPI_Datatype mpifloat = get_mpi_datatype(mtrl);
   MPI_Allreduce(&locavg,&avg,1,mpifloat,MPI_SUM,MPI_COMM_WORLD);
   float_sw4 sum2;
   MPI_Allreduce(&locsum2,&sum2,1,mpifloat,MPI_SUM,MPI_COMM_WORLD);
   size_t nelem = static_cast<size_t>(m_nig)*m_njg*m_nkg;
   avg = avg/(nelem);
   float_sw4 stdev = sqrt(sum2/nelem);
   if( abs(avg) > 1e-10*stdev )
      cout << "Warning, average random perturbation is " << avg << endl;
   //   cout << "stdev = " << stdev << " avg = " << avg << " sigma= " << m_sigma << endl;
   // Rescale to desired sigma, add 1 for later multiplication with given material
   float_sw4 istdev = 1/stdev;
   for( size_t ind = 0 ; ind < n ;ind++ )
      mtrl[ind] = 1 + m_sigma*mtrl[ind]*istdev;
}

//#include <complex>
//#include <vector>
//#include <iostream>
//
//#include <mpi.h>

//#include "AllDims.h"

//-----------------------------------------------------------------------
#include "Patch.h"

//-----------------------------------------------------------------------
template<class T>
void RandomizedMaterial::redistribute_array( AllDims& src, AllDims& dest, T* src_array, T* dest_array )
{
   // Copies array src_array to array dest_array, where the dimensions of the
   // arrays on the processors are described by the AllDims objects src and dest respectively.

   std::vector<Patch*> sendlist, recvlist;
   Patch* selfsend;
   int dims[6], selfsnr=0, selfrnr=0;
   MPI_Datatype mpi_datatype = get_mpi_datatype( src_array );

   for( int p3 = 0 ; p3 < dest.m_nprock ; p3++)
      for( int p2 = 0 ; p2 < dest.m_nprocj ; p2++)
	 for( int p1 = 0 ; p1 < dest.m_nproci ; p1++)
	    if( dest.intersect(p1,p2,p3,src,dims) )
	    {
	       Patch* intersection = new Patch( dims, dest.proc1d(p1,p2,p3) );
	       if( !dest.owner(p1,p2,p3) )
		  sendlist.push_back(intersection);
	       else
	       {
		  selfsend = intersection;
		  selfsnr++;
	       }
	    }
   
   if( selfsnr > 1 )
      std::cout << "ERROR: found more than one self intersection" << std::endl;
   
   for( int p3 = 0 ; p3 < src.m_nprock ; p3++)
      for( int p2 = 0 ; p2 < src.m_nprocj ; p2++)
	 for( int p1 = 0 ; p1 < src.m_nproci ; p1++)
	    if( src.intersect(p1,p2,p3,dest,dims) )
	    {
	       Patch* intersection = new Patch( dims, src.proc1d(p1, p2, p3) );
	       if( !src.owner(p1,p2,p3) )
		  recvlist.push_back(intersection);
	       else
		  selfrnr++;
	    }
   if( selfrnr != selfsnr )
      std::cout << "ERROR: different number of self send and self recieve patches"<< std::endl;

   T** recv_array_patch;
   MPI_Request* req;
   int tag=667;

   //      if( src.m_myid1d == 0 )
   //      {
   //         cout << "sendlist      " << sendlist.size() << endl;
   //         cout << "recvlist      " << recvlist.size() << endl;
   //         cout << "selfintersect " << selfsnr << endl;
   //      }
   if( selfsnr == 1 )
   {
     // Copy patches between my own parts of the arrays
      selfsend->selfcopy( src, src_array, dest, dest_array );
   }

   // Communication between tasks
   if( recvlist.size() > 0 )
   {
      recv_array_patch = new T*[recvlist.size()];
      req = new MPI_Request[recvlist.size()];
   }

   // Post receives
   for( int r=0 ; r < recvlist.size() ; r++ )
   {
      size_t size=recvlist[r]->size();
      recv_array_patch[r] = new T[size];
      MPI_Irecv( recv_array_patch[r], size, mpi_datatype, recvlist[r]->m_procid, tag,
		 MPI_COMM_WORLD, &req[r] );
   }
   // Pack data and send
   for( int s=0 ; s < sendlist.size() ; s++ )
   {
      size_t size=sendlist[s]->size();
      T* array_patch = new T[size];
      sendlist[s]->pack( src_array, src, array_patch );
      MPI_Send( array_patch, size, mpi_datatype, sendlist[s]->m_procid, tag, MPI_COMM_WORLD );
      delete[] array_patch;
   }

   // Do the receieve and unpack
   for( int r=0 ; r < recvlist.size() ; r++ )
   {
      MPI_Status status;
      MPI_Wait( &req[r], &status );
      recvlist[r]->unpack( dest_array, dest, recv_array_patch[r] );
      delete[] recv_array_patch[r];
   }

   if( recvlist.size() > 0 )
   {
      delete[] recv_array_patch;
      delete[] req;
   }
}


// -----------------------------------------------------------------------
void RandomizedMaterial::repad_sarray( Sarray& sar, int old_padding, int new_padding )
{
   // Assuming 2D processor decomp.
   if( old_padding == new_padding )
      return;

   // Shorter names
   int neigh[4]={mEW->m_neighbor[0],mEW->m_neighbor[1],mEW->m_neighbor[2],mEW->m_neighbor[3]};

   // Leave physical boundaries unchanged.
   int ib = sar.m_ib, ie=sar.m_ie;
   if( neigh[0] != MPI_PROC_NULL )
      ib = ib + old_padding - new_padding;
   if( neigh[1] != MPI_PROC_NULL )
      ie = ie - old_padding + new_padding;

   int jb = sar.m_jb, je=sar.m_je;
   if( neigh[2] != MPI_PROC_NULL )
      jb = jb + old_padding - new_padding;
   if( neigh[3] != MPI_PROC_NULL )
      je = je - old_padding + new_padding;

   int kb=sar.m_kb, ke=sar.m_ke;
      
   if( old_padding < new_padding )
   {
      Sarray tmp(sar.m_nc,ib,ie,jb,je,kb,ke);
   //   cout << "REPAD before " << sar.m_ib << " " << sar.m_ie << " " << sar.m_jb << " " << sar.m_je
   //	<< " " << sar.m_kb << " " << sar.m_ke << endl;
   //   cout << "REPAD after " << ib << " " << ie << " " << jb << " " << je
   //	<< " " << kb << " " << ke << endl;
      tmp.insert_subarray(sar.m_ib,sar.m_ie,sar.m_jb,sar.m_je,sar.m_kb,sar.m_ke,sar.c_ptr());
      comm_sarray( tmp, neigh, new_padding );
      sar.copy( tmp );
   }
   if( old_padding > new_padding )
   {
      cout << "WARNING: RandomizedMaterial::repad_sarray: reduction of pad points NYI" << endl;
   }
}

// -----------------------------------------------------------------------
void RandomizedMaterial::comm_sarray( Sarray& sar, int neigh[4], int padding )
{
   const int ib=sar.m_ib, ie=sar.m_ie, jb=sar.m_jb, je=sar.m_je, kb=sar.m_kb, ke=sar.m_ke;

   size_t npts1 = static_cast<size_t>(padding)*(je-jb+1)*(ke-kb+1)*sar.m_nc;
   size_t npts2 = static_cast<size_t>(padding)*(ie-ib+1)*(ke-kb+1)*sar.m_nc;
   size_t ptsmax = npts1>npts2?npts1:npts2;

   float_sw4* sbuf = new float_sw4[ptsmax];
   float_sw4* rbuf = new float_sw4[ptsmax];
   MPI_Datatype mpi_float = get_mpi_datatype(sbuf);
   MPI_Comm comm = MPI_COMM_WORLD;

   MPI_Status status;
   int xtag1 = 3343, xtag2 = 3344;
   int ytag1 = 3345, ytag2 = 3346;

   //   memset(sbuf,0,ptsmax*sizeof(float_sw4)); // To shut up memory checker
// I-direction communication
   if( neigh[0] != MPI_PROC_NULL )
      sar.extract_subarray( ib+padding, ib+2*padding-1, jb, je, kb, ke, sbuf );

   MPI_Sendrecv( sbuf, npts1, mpi_float, neigh[0], xtag1, 
                 rbuf, npts1, mpi_float, neigh[1], xtag1, comm, &status );
   if( neigh[1] != MPI_PROC_NULL )
   {
      sar.insert_subarray( ie-padding+1, ie, jb, je, kb, ke, rbuf );
      sar.extract_subarray( ie-2*padding+1, ie-padding, jb, je, kb, ke, sbuf );
   }
   MPI_Sendrecv( sbuf, npts1, mpi_float, neigh[1], xtag2, 
                 rbuf, npts1, mpi_float, neigh[0], xtag2, comm, &status );
   if( neigh[0] != MPI_PROC_NULL )
      sar.insert_subarray( ib, ib+padding-1, jb, je, kb, ke, rbuf );

// J-direction communication
   if( neigh[2] != MPI_PROC_NULL )
      sar.extract_subarray( ib, ie, jb+padding, jb+2*padding-1, kb, ke, sbuf );
   MPI_Sendrecv( sbuf, npts2, mpi_float, neigh[2], ytag1, 
                 rbuf, npts2, mpi_float, neigh[3], ytag1, comm, &status);
   if( neigh[3] != MPI_PROC_NULL )
   {
      sar.insert_subarray( ib, ie, je-padding+1, je, kb, ke, rbuf );
      sar.extract_subarray( ib, ie, je-2*padding+1, je-padding, kb, ke, sbuf );
   }
   MPI_Sendrecv( sbuf, npts2, mpi_float, neigh[3], ytag2, 
                 rbuf, npts2, mpi_float, neigh[2], ytag2, comm, &status);
   if( neigh[2] != MPI_PROC_NULL )
      sar.insert_subarray( ib, ie, jb, jb+padding-1, kb, ke, rbuf );

   delete[] sbuf;
   delete[] rbuf;
}

//-----------------------------------------------------------------------
void RandomizedMaterial::set_vsmax( float_sw4 vsmax )
{
   m_vsmax = vsmax;
}

//-----------------------------------------------------------------------
// Instantiations
template void RandomizedMaterial::redistribute_array<float_sw4>( AllDims& src, AllDims& dest, float_sw4* src_array, float_sw4* dest_array );
