#include <sys/types.h>

#include "AllDims.h"
#include "Patch.h"
#include "sw4.h"

//-----------------------------------------------------------------------
Patch::Patch(sw4_type dims[6], sw4_type procid) {
  m_ib = dims[0];
  m_ie = dims[1];
  m_jb = dims[2];
  m_je = dims[3];
  m_kb = dims[4];
  m_ke = dims[5];
  m_procid = procid;
}

//-----------------------------------------------------------------------
template <class T>
void Patch::pack(T* array, AllDims& dims, T* array_patch) {
  sw4_type myib = dims.m_ib[dims.m_myid3di];
  sw4_type myie = dims.m_ie[dims.m_myid3di];
  sw4_type myjb = dims.m_jb[dims.m_myid3dj];
  sw4_type myje = dims.m_je[dims.m_myid3dj];
  sw4_type mykb = dims.m_kb[dims.m_myid3dk];
  sw4_type myke = dims.m_ke[dims.m_myid3dk];
  sw4_type ni = myie - myib + 1;
  sw4_type nj = myje - myjb + 1;
  sw4_type nk = myke - mykb + 1;
  size_t ind = 0;
  size_t ai = 1, aj = ni, ak = ni * nj;
  if (dims.m_indrev) {
    ai = nk * nj;
    aj = nk;
    ak = 1;
  }
  for (sw4_type k = m_kb; k <= m_ke; k++)
    for (sw4_type j = m_jb; j <= m_je; j++)
      for (sw4_type i = m_ib; i <= m_ie; i++)
        array_patch[ind++] =
            array[ai * (i - myib) + aj * (j - myjb) + ak * (k - mykb)];
}

//-----------------------------------------------------------------------
template <class T>
void Patch::unpack(T* array, AllDims& dims, T* array_patch) {
  sw4_type myib = dims.m_ib[dims.m_myid3di];
  sw4_type myie = dims.m_ie[dims.m_myid3di];
  sw4_type myjb = dims.m_jb[dims.m_myid3dj];
  sw4_type myje = dims.m_je[dims.m_myid3dj];
  sw4_type mykb = dims.m_kb[dims.m_myid3dk];
  sw4_type myke = dims.m_ke[dims.m_myid3dk];
  sw4_type ni = myie - myib + 1;
  sw4_type nj = myje - myjb + 1;
  sw4_type nk = myke - mykb + 1;
  size_t ind = 0;
  size_t ai = 1, aj = ni, ak = ni * nj;
  if (dims.m_indrev) {
    ai = nk * nj;
    aj = nk;
    ak = 1;
  }
  for (sw4_type k = m_kb; k <= m_ke; k++)
    for (sw4_type j = m_jb; j <= m_je; j++)
      for (sw4_type i = m_ib; i <= m_ie; i++)
        array[ai * (i - myib) + aj * (j - myjb) + ak * (k - mykb)] =
            array_patch[ind++];
}

//-----------------------------------------------------------------------
template <class T>
void Patch::selfcopy(AllDims& src, T* src_array, AllDims& dest, T* dest_array) {
  sw4_type myib = src.m_ib[src.m_myid3di];
  sw4_type myie = src.m_ie[src.m_myid3di];
  sw4_type myjb = src.m_jb[src.m_myid3dj];
  sw4_type myje = src.m_je[src.m_myid3dj];
  sw4_type mykb = src.m_kb[src.m_myid3dk];
  sw4_type myke = src.m_ke[src.m_myid3dk];
  sw4_type ni = myie - myib + 1;
  sw4_type nj = myje - myjb + 1;
  sw4_type nk = myke - mykb + 1;

  sw4_type myibd = dest.m_ib[dest.m_myid3di];
  sw4_type myied = dest.m_ie[dest.m_myid3di];
  sw4_type myjbd = dest.m_jb[dest.m_myid3dj];
  sw4_type myjed = dest.m_je[dest.m_myid3dj];
  sw4_type mykbd = dest.m_kb[dest.m_myid3dk];
  sw4_type myked = dest.m_ke[dest.m_myid3dk];
  sw4_type nid = myied - myibd + 1;
  sw4_type njd = myjed - myjbd + 1;
  sw4_type nkd = myked - mykbd + 1;

  size_t ai = 1, aj = ni, ak = ni * nj;
  if (src.m_indrev) {
    ai = nk * nj;
    aj = nk;
    ak = 1;
  }
  size_t aid = 1, ajd = nid, akd = nid * njd;
  if (dest.m_indrev) {
    aid = nkd * njd;
    ajd = nkd;
    akd = 1;
  }
  for (sw4_type k = m_kb; k <= m_ke; k++)
    for (sw4_type j = m_jb; j <= m_je; j++)
      for (sw4_type i = m_ib; i <= m_ie; i++)
        dest_array[aid * (i - myibd) + ajd * (j - myjbd) + akd * (k - mykbd)] =
            src_array[ai * (i - myib) + aj * (j - myjb) + ak * (k - mykb)];
}

//-----------------------------------------------------------------------
size_t Patch::size() {
  return static_cast<size_t>((m_ie - m_ib + 1)) * (m_je - m_jb + 1) *
         (m_ke - m_kb + 1);
}

//-----------------------------------------------------------------------
// Instantiations
template void Patch::pack<float_sw4>(float_sw4* array, AllDims& dims,
                                     float_sw4* array_patch);
template void Patch::unpack<float_sw4>(float_sw4* array, AllDims& dims,
                                       float_sw4* array_patch);
template void Patch::selfcopy<float_sw4>(AllDims& src, float_sw4* src_array,
                                         AllDims& dest, float_sw4* dest_array);
