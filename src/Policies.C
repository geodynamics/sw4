#include <iostream>

#include "Policies.h"
#include "RAJA/RAJA.hpp"
sw4_type main() {
#ifdef POLICY_TEST
  std::cout << "Hello policies\n";
  const sw4_type M = 500;
  const sw4_type N = M * M * M;
  // double *d = new double[N];
  double *d, *a, *b;

  cudaMallocManaged(&d, N * sizeof(double));
  cudaMallocManaged(&a, N * sizeof(double));
  cudaMallocManaged(&b, N * sizeof(double));
  // using Pol = Policies::Policy<Policies::Cuda,Policies::Default>::type;
  //  using Pol2 = Policies::Policy<Policies::Cuda, Policies::Loop1<32>>::type;
  // RAJA::forall<Pol>(
  //       // RAJA::forall<NEWP>(
  //       RAJA::RangeSegment(0, N), [=] RAJA_DEVICE(size_t i) { d[i] = i; });
  //   for (sw4_type i = 0; i < N; i += N / 10) std::cout << d[i] << ",";
  //   std::cout << "\n";
  //   delete[] d;
  forallasync(0, N, [=] __host__ __device__(sw4_type i) { d[i] = i; });

  forallasync(0, N, [=] __host__ __device__(sw4_type i) {
    double x = static_cast<double>(i) / N;
    d[i] = i * i + sin(x) + cos(x) + tan(x) + i % 10 + i % 999;
  });

  // COmpare GS vs normal for compute sw4_typeensive kernels
  // Here GS takes 18% less time most likely due
  // ILP within the thread
  for (sw4_type i = 0; i < 20; i++) {
    forallasync(0, N, [=] __host__ __device__(sw4_type i) {
      double x = static_cast<double>(i) / N;
      d[i] = i * i + sin(x) + cos(x) + tan(x) + i % 10 + i % 999;
    });
    forallX<1024>(0, N, [=] __host__ __device__(sw4_type i) {
      double x = static_cast<double>(i) / N;
      d[i] = i * i + sin(x) + cos(x) + tan(x) + i % 10 + i % 999;
    });
  }

  //  Compare GS vs plain for bandwidth limited kernels
  // Here GS takes 0.5% more time
  forallasync(0, N, [=] __host__ __device__(sw4_type i) { a[i] = 3 * b[i] + d[i]; });

  for (sw4_type i = 0; i < 25; i++) {
    forallasync(0, N,
                [=] __host__ __device__(sw4_type i) { a[i] = 3 * b[i] + d[i]; });

    forallX<1024>(0, N,
                  [=] __host__ __device__(sw4_type i) { a[i] = 3 * b[i] + d[i]; });
  }

  // Do the same using 3D kernels. GS is 28% less
  // AVerage time is less for 1D kernel
#define D(i, j, k) d[i + M * j + M * M * k]
  for (sw4_type i = 0; i < 30; i++) {
    forall3(0, M, 0, M, 0, M, [=] __host__ __device__(sw4_type i, sw4_type j, sw4_type k) {
      double x = static_cast<double>(i) / N;
      D(i, j, k) = i * i + sin(x) + cos(x) + tan(x) + i % 10 + i % 999;
    });

    forall3X<1024>(
        0, M, 0, M, 0, M, [=] __host__ __device__(sw4_type i, sw4_type j, sw4_type k) {
          double x = static_cast<double>(i) / N;
          D(i, j, k) = i * i + sin(x) + cos(x) + tan(x) + i % 10 + i % 999;
        });
  }

  // for Bandwidth limited kernels,GS takes 2.1X more time
#define A(i, j, k) a[i + M * j + M * M * k]
#define B(i, j, k) b[i + M * j + M * M * k]
  for (sw4_type i = 0; i < 35; i++) {
    forall3(0, M, 0, M, 0, M, [=] __host__ __device__(sw4_type i, sw4_type j, sw4_type k) {
      A(i, j, k) = 3 * B(i, j, k) + D(i, j, k);
    });

    forall3X<1024>(0, M, 0, M, 0, M,
                   [=] __host__ __device__(sw4_type i, sw4_type j, sw4_type k) {
                     A(i, j, k) = 3 * B(i, j, k) + D(i, j, k);
                   });
  }

  cudaStreamSynchronize(0);
#endif
  // Sw4_Typerinsics tests using Parray;

  const sw4_type N = 10;

  Parray Aa(N, N, N, N);
  Parray Ba(N, N, N, N);
  Parray Ca(N, N, N, N);
  forall3(0, N, 0, N, 0, N, [=] __device__(sw4_type i, sw4_type j, sw4_type k) {
    for (sw4_type c = 0; c < N; c++)
      Aa(i, j, k, c) = i + j * N + k * N * N + c * N * N * N;
    for (sw4_type c = 0; c < N; c++)
      Ba(i, j, k, c) = -(i + j * N + k * N * N + c * N * N * N);
    for (sw4_type c = 0; c < N; c++) Ca(i, j, k, c) = -1.0;
  });

  forall3(1, N - 1, 1, N - 1, 1, N - 1, [=] __device__(sw4_type i, sw4_type j, sw4_type k) {
    for (sw4_type c = 0; c < N; c++)
      Ca(i, j, k, c) += Aa(i, j, k, c) + Ba(i, j, k, c) + Ba(i - 1, j, k, c) +
                        Ba(i, j - 1, k, c) + Ba(i, j, k - 1, c);
  });

  cudaStreamSynchronize(0);
  const sw4_type V = Ca.offset(N / 2, N / 2, N / 2, N / 2);
  std::cout << " Exact value is " << Ca.data[V] << " A = " << Aa.data[V]
            << " B = " << Ba.data[V] << "\n";
  unsigned sw4_type I;
  for (sw4_type i = 0; i < N; i++)
    for (sw4_type j = 0; j < N; j++) {
      I = (i << 8) ^ j;
      I = i;
      // std::cout<<sizeof(sw4_type)<<" "<<i<<" "<<j<<" "<<" In hex
      // "<<std::hex<<I<<"\n";
      I = i << 8;
      // std::cout<<std::dec<<sizeof(sw4_type)<<" "<<i<<" "<<j<<" "<<" In hex
      // "<<std::hex<<I<<"\n";
      I = I ^ j;
      // std::cout<<std::dec<<sizeof(sw4_type)<<" "<<i<<" "<<j<<" "<<" In hex
      // "<<std::hex<<I<<"\n";
      sw4_type jj = I & 0x00ff;
      sw4_type ii = (I >> 8) & 0xff;
      if ((i != ii) || (j != jj)) {
        std::cout << std::dec << sizeof(sw4_type) << " OUT  " << ii << " " << jj
                  << " "
                  << " In hex " << std::hex << I << "\n";
      }
    }
  for (sw4_type l = 0; l < 5; l++) {
    Ca.set(-1.0);
    forall3(1, N - 1, 1, N - 1, 1, N - 1, [=] __device__(sw4_type i, sw4_type j, sw4_type k) {
      unsigned sw4_type I;
      unsigned sw4_type J, JJ;
      const unsigned sw4_type p10 = 0x0100, p01 = 0x0001;
      I = (i << 8) ^ j;
      JJ = (k << 8);
      unsigned sw4_type IM = __vsub2(I, p10);
      unsigned sw4_type JM = __vsub2(I, p01);
      unsigned sw4_type KM = __vsub2(JJ, p10);
      for (sw4_type c = 0; c < N; c++) {
        J = JJ ^ c;

        Ca(I, J) += Aa(I, J) + Ba(I, J) + Ba(IM, J) + Ba(JM, J) + Ba(I, KM ^ c);
      }
    });
    cudaStreamSynchronize(0);
    std::cout << " Mid value is " << Ca.data[V] << " A = " << Aa.data[V]
              << " B = " << Ba.data[V] << "\n";
  }

  for (sw4_type l = 0; l < 7; l++) {
    // Ca.set(-1.0);
    forall3(1, N - 1, 1, N - 1, 1, N - 1, [=] __device__(sw4_type i, sw4_type j, sw4_type k) {
      unsigned sw4_type I;
      unsigned sw4_type J, JJ;
      const unsigned sw4_type p10 = 0x0100, p01 = 0x0001;
      I = (i << 8) ^ j;
      JJ = (k << 8);
      for (sw4_type c = 0; c < N; c++) {
        J = JJ ^ c;
        Ca(I, J) += Aa(i, j, k, c) + Ba(i, j, k, c) + Ba(i - 1, j, k, c) +
                    Ba(i, j - 1, k, c) + Ba(i, j, k - 1, c);
      }
    });
    // cudaStreamSynchronize(0);
    // std::cout<<" Baseline value is "<<Ca.data[V]<<" A = "<<Aa.data[V]<<" B =
    // "<<Ba.data[V]<<"\n";
  }
  cudaStreamSynchronize(0);

  // sw4_type lc=0;
  // for (sw4_type i=1;i<N-1;i++) for(sw4_type j=1;j<N-1;j++) for( sw4_type k=1;k<N-1;k++)
  // for(sw4_type c=1;c<N-1;c++){

  // 	  unsigned sw4_type I = (i<<8)^j;
  // 	  unsigned sw4_type J = (k<<8);
  // 	    J=J^c;
  // 	  sw4_type o1 = Ca.offset(i,j,k,c);
  // 	  sw4_type o2 = Ca.offset2(I,J);
  // 	  if (o1!=o2){
  // 	    std::cout<<"MISMATCH "<<o1<<" "<<o2<<"\n";
  // 	    break;
  // 	  } else lc++;
  // 	}
  // std::cout<<" lc is "<<lc<<"\n";
}
