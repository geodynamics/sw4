#include <iostream>
#include "Policies.h"
#include "RAJA/RAJA.hpp"
int main() {
#ifdef POLICY_TEST
  std::cout << "Hello policies\n";
  const int M = 500;
  const int N = M*M*M;
  //double *d = new double[N];
  double *d,*a,*b;

  cudaMallocManaged(&d,N*sizeof(double));
  cudaMallocManaged(&a,N*sizeof(double));
  cudaMallocManaged(&b,N*sizeof(double));
  //using Pol = Policies::Policy<Policies::Cuda,Policies::Default>::type;
  //  using Pol2 = Policies::Policy<Policies::Cuda, Policies::Loop1<32>>::type;
  // RAJA::forall<Pol>(
//       // RAJA::forall<NEWP>(
//       RAJA::RangeSegment(0, N), [=] RAJA_DEVICE(size_t i) { d[i] = i; });
//   for (int i = 0; i < N; i += N / 10) std::cout << d[i] << ",";
//   std::cout << "\n";
//   delete[] d;
  forallasync(0,N,[=]__host__ __device__ (int i){
      d[i]=i;});
  
  forallasync(0,N,[=]__host__ __device__ (int i){
      double x = static_cast<double>(i)/N;
      d[i]=i*i+sin(x)+cos(x)+tan(x)+i%10+i%999;
    });


  // COmpare GS vs normal for compute intensive kernels
  // Here GS takes 18% less time most likely due
  // ILP within the thread
  for( int i=0;i<20;i++){
    forallasync(0,N,[=]__host__ __device__ (int i){
	double x = static_cast<double>(i)/N;
	d[i]=i*i+sin(x)+cos(x)+tan(x)+i%10+i%999;
      });
    forallX<1024>(0,N,[=]__host__ __device__ (int i){
	double x = static_cast<double>(i)/N;
	d[i]=i*i+sin(x)+cos(x)+tan(x)+i%10+i%999;
      });
  }
 
  //  Compare GS vs plain for bandwidth limited kernels
  // Here GS takes 0.5% more time
  forallasync(0,N,[=]__host__ __device__ (int i){
      a[i]=3*b[i]+d[i];});

  for( int i=0;i<25;i++){

    forallasync(0,N,[=]__host__ __device__ (int i){
	a[i]=3*b[i]+d[i];});
    
    
    forallX<1024>(0,N,[=]__host__ __device__ (int i){
	a[i]=3*b[i]+d[i];});
}


  // Do the same using 3D kernels. GS is 28% less
  // AVerage time is less for 1D kernel
#define D(i,j,k) d[i+M*j+M*M*k]
  for( int i=0;i<30;i++){
    forall3(0,M,0,M,0,M,[=]__host__ __device__ (int i,int j, int k){
	double x = static_cast<double>(i)/N;
	D(i,j,k)=i*i+sin(x)+cos(x)+tan(x)+i%10+i%999;
      });

    forall3X<1024>(0,M,0,M,0,M,[=]__host__ __device__ (int i,int j, int k){
	double x = static_cast<double>(i)/N;
	D(i,j,k)=i*i+sin(x)+cos(x)+tan(x)+i%10+i%999;
      });
  }
 

  // for Bandwidth limited kernels,GS takes 2.1X more time
#define A(i,j,k) a[i+M*j+M*M*k]
#define B(i,j,k) b[i+M*j+M*M*k]
  for( int i=0;i<35;i++){

    forall3(0,M,0,M,0,M,[=]__host__ __device__ (int i,int j, int k){
	A(i,j,k) = 3*B(i,j,k)+D(i,j,k);});
    
    
    forall3X<1024>(0,M,0,M,0,M,[=]__host__ __device__ (int i,int j, int k){
	A(i,j,k) = 3*B(i,j,k)+D(i,j,k);});
}


  cudaStreamSynchronize(0);
#endif
  // INtrinsics tests using Parray;

  const int N=10;

  Parray Aa(N,N,N,N);
  Parray Ba(N,N,N,N);
  Parray Ca(N,N,N,N);
  forall3(0,N,0,N,0,N,[=]__device__(int i, int j, int k){
      for (int c=0;c<N;c++) Aa(i,j,k,c)=i+j*N+k*N*N+c*N*N*N;
      for (int c=0;c<N;c++) Ba(i,j,k,c)=-(i+j*N+k*N*N+c*N*N*N);
      for (int c=0;c<N;c++) Ca(i,j,k,c)=-1.0;
  });

forall3(1,N-1,1,N-1,1,N-1,[=]__device__(int i, int j, int k){
    for (int c=0;c<N;c++) Ca(i,j,k,c)+=Aa(i,j,k,c)+Ba(i,j,k,c)+Ba(i-1,j,k,c)+Ba(i,j-1,k,c)+Ba(i,j,k-1,c);
  });
    


 cudaStreamSynchronize(0);
 const int V=Ca.offset(N/2,N/2,N/2,N/2);
 std::cout<<" Exact value is "<<Ca.data[V]<<" A = "<<Aa.data[V]<<" B = "<<Ba.data[V]<<"\n";
 unsigned int I;
 for (int i=0;i<N;i++) for(int j=0;j<N;j++){
 I = (i<<8)^j;
 I=i;
 //std::cout<<sizeof(int)<<" "<<i<<" "<<j<<" "<<" In hex "<<std::hex<<I<<"\n";
 I=i<<8;
 //std::cout<<std::dec<<sizeof(int)<<" "<<i<<" "<<j<<" "<<" In hex "<<std::hex<<I<<"\n";
 I=I^j;
 //std::cout<<std::dec<<sizeof(int)<<" "<<i<<" "<<j<<" "<<" In hex "<<std::hex<<I<<"\n";
 int jj = I & 0x00ff;
 int ii = (I>>8) & 0xff;
 if ((i!=ii) || (j!=jj)){
     std::cout<<std::dec<<sizeof(int)<<" OUT  "<<ii<<" "<<jj<<" "<<" In hex "<<std::hex<<I<<"\n";
   }
   }
 for (int l=0;l<5;l++){
   Ca.set(-1.0);
forall3(1,N-1,1,N-1,1,N-1,[=]__device__(int i, int j, int k){
    unsigned int I;
    unsigned int J,JJ;
    const unsigned int p10 = 0x0100, p01=0x0001;
    I = (i<<8)^j;
    JJ = (k<<8);
    unsigned int IM = __vsub2(I,p10);
    unsigned int JM = __vsub2(I,p01);
    unsigned int KM = __vsub2(JJ,p10);
      for (int c=0;c<N;c++) {
	J=JJ^c;
	
	Ca(I,J)+=Aa(I,J)+Ba(I,J)+Ba(IM,J)+Ba(JM,J)+Ba(I,KM^c);
      }
  });
 cudaStreamSynchronize(0);
 std::cout<<" Mid value is "<<Ca.data[V]<<" A = "<<Aa.data[V]<<" B = "<<Ba.data[V]<<"\n";
 }


 for (int l=0;l<7;l++){
   //Ca.set(-1.0);
forall3(1,N-1,1,N-1,1,N-1,[=]__device__(int i, int j, int k){
    unsigned int I;
    unsigned int J,JJ;
    const unsigned int p10 = 0x0100, p01=0x0001;
    I = (i<<8)^j;
    JJ = (k<<8);
      for (int c=0;c<N;c++) {
	J=JJ^c;
	Ca(I,J)+=Aa(i,j,k,c)+Ba(i,j,k,c)+Ba(i-1,j,k,c)+Ba(i,j-1,k,c)+Ba(i,j,k-1,c);
      }
  });
//cudaStreamSynchronize(0);
//std::cout<<" Baseline value is "<<Ca.data[V]<<" A = "<<Aa.data[V]<<" B = "<<Ba.data[V]<<"\n";

 }
  cudaStreamSynchronize(0);


  // int lc=0;
  // for (int i=1;i<N-1;i++) for(int j=1;j<N-1;j++) for( int k=1;k<N-1;k++) for(int c=1;c<N-1;c++){

  // 	  unsigned int I = (i<<8)^j;
  // 	  unsigned int J = (k<<8);
  // 	    J=J^c;
  // 	  int o1 = Ca.offset(i,j,k,c);
  // 	  int o2 = Ca.offset2(I,J);
  // 	  if (o1!=o2){
  // 	    std::cout<<"MISMATCH "<<o1<<" "<<o2<<"\n";
  // 	    break;
  // 	  } else lc++;
  // 	}
  // std::cout<<" lc is "<<lc<<"\n";
}
