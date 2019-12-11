#ifndef __FORALLSDECL_H__
#define __FORALLSDECL_H__ 1
template <int N, int M, int L>
std::map<std::tuple<int, int, int>, std::tuple<dim3, dim3>>
    RangeAT<N, M, L>::map;
template <int N, int M, int L>
std::map<std::tuple<int, int, int>, float> RangeAT<N, M, L>::mintime;
template <int N, int M, int L>
std::map<std::tuple<int, int, int>, int> RangeAT<N, M, L>::vloca;
template <int N, int M, int L>
dim3 RangeAT<N, M, L>::mint;
template <int N, int M, int L>
dim3 RangeAT<N, M, L>::minb;
template <int N, int M, int L>
std::map<std::tuple<int, int, int>, std::vector<std::tuple<dim3, dim3>>>
    RangeAT<N, M, L>::confs;
template <int N, int M, int L>
std::map<std::tuple<int, int, int>, std::ofstream> RangeAT<N, M, L>::files;
#endif
