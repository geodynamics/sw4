#ifndef __POLICIES_FOR_SW4__
#define __POLICIES_FOR_SW4__
namespace Policies {
template <int T>
class Loop1 {
 public:
  const static int value = T;
};
struct Cuda {};
struct OpenMP {};
struct Default {};
template <typename T1, typename T2>
class Policy {
 public:
  typedef int type;  // to make compilation fail
};
template <typename T2>
class Policy<Cuda, T2> {
 public:
  typedef RAJA::cuda_exec<T2::value> type;
};

template <typename T1, typename... args>
class Policy2 {
 public:
  typedef int type;  // to make compilation fail
};
template <typename... args>
class Policy2<Cuda, args...> {
 public:
  typedef RAJA::cuda_exec<1024> type;
};
}  // namespace Policies
#endif
