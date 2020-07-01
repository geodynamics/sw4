#include <iostream>
#include "Policies.h"
#include "RAJA/RAJA.hpp"
int main() {
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
}
