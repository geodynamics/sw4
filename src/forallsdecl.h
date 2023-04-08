#ifndef __FORALLSDECL_H__
#define __FORALLSDECL_H__ 1

#ifdef ENABLE_GPU
template <sw4_type N, sw4_type M, sw4_type L>
std::map<std::tuple<sw4_type, sw4_type, sw4_type>, std::tuple<dim3, dim3>>
    RangeAT<N, M, L>::map;

template <sw4_type N, sw4_type M, sw4_type L>
std::map<std::tuple<sw4_type, sw4_type, sw4_type>, float> RangeAT<N, M, L>::msw4_typeime;

template <sw4_type N, sw4_type M, sw4_type L>
std::map<std::tuple<sw4_type, sw4_type, sw4_type>, sw4_type> RangeAT<N, M, L>::vloca;

template <sw4_type N, sw4_type M, sw4_type L>
dim3 RangeAT<N, M, L>::msw4_type;

template <sw4_type N, sw4_type M, sw4_type L>
dim3 RangeAT<N, M, L>::minb;

template <sw4_type N, sw4_type M, sw4_type L>
std::map<std::tuple<sw4_type, sw4_type, sw4_type>, std::vector<std::tuple<dim3, dim3>>>
    RangeAT<N, M, L>::confs;

template <sw4_type N, sw4_type M, sw4_type L>
std::map<std::tuple<sw4_type, sw4_type, sw4_type>, std::ofstream> RangeAT<N, M, L>::files;

#endif

#endif
