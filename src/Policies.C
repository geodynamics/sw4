#include <iostream>

#include "Policies.h"
#include "RAJA/RAJA.hpp"
int main() {
  std::cout << "Hello policies\n";
  const int N = 1000;
  double *d = new double[N];
  // using Pol = Policies::Policy<Policies::Cuda,Policies::Default>::type;
  using Pol2 = Policies::Policy<Policies::Cuda, Policies::Loop1<32>>::type;
  RAJA::forall<Pol2>(
      // RAJA::forall<NEWP>(
      RAJA::RangeSegment(0, N), [=] RAJA_DEVICE(size_t i) { d[i] = i; });
  for (int i = 0; i < N; i += N / 10) std::cout << d[i] << ",";
  std::cout << "\n";
  delete[] d;
}
