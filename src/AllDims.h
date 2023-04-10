#ifndef SW4_ALLDIMS_H
#define SW4_ALLDIMS_H
#include "sw4.h"
#include <mpi.h>
#include <sys/types.h>

#include <vector>

class AllDims {
 public:
  sw4_type m_myid1d;
  sw4_type m_myid3di, m_myid3dj, m_myid3dk;

  sw4_type m_nproci, m_nprocj, m_nprock;
  sw4_type m_ibg, m_ieg, m_jbg, m_jeg, m_kbg,
      m_keg;  // Global dims, includes ghost points

  std::vector<sw4_type> m_ib, m_ie, m_jb, m_je, m_kb,
      m_ke;  // Local for each processor, includes padding and ghost points
  sw4_type m_npad, m_nghost;
  bool m_indrev;
  size_t m_fftw_alloc_local;
  MPI_Comm m_communicator;

  AllDims(sw4_type nproci, sw4_type nprocj, sw4_type nprock, sw4_type ibg, sw4_type ieg, sw4_type jbg,
          sw4_type jeg, sw4_type kbg, sw4_type keg, sw4_type nghost, sw4_type npad, MPI_Comm ewcomm);

  AllDims(sw4_type nprocs, sw4_type ibg, sw4_type ieg, sw4_type jbg, sw4_type jeg, sw4_type kbg, sw4_type keg,
          sw4_type nghost, MPI_Comm ewcomm);

  AllDims(sw4_type nproci, sw4_type nprocj, sw4_type nprock, sw4_type ibg, sw4_type ieg, sw4_type jbg,
          sw4_type jeg, sw4_type kbg, sw4_type keg, sw4_type nghost, sw4_type npad);
  AllDims(sw4_type nprocs, sw4_type ibg, sw4_type ieg, sw4_type jbg, sw4_type jeg, sw4_type kbg, sw4_type keg,
          sw4_type nghost);
  AllDims(AllDims* fine, sw4_type ibg, sw4_type ieg, sw4_type jbg, sw4_type jeg, sw4_type kbg, sw4_type keg,
          sw4_type nghost, sw4_type npad);

  void getdims_nopad(sw4_type dims[6], sw4_type p1 = -1, sw4_type p2 = -1, sw4_type p3 = -1);
  bool intersect(sw4_type p1, sw4_type p2, sw4_type p3, AllDims& other, sw4_type dims[6]);
  void compute_myid3d();
  sw4_type proc1d(sw4_type p1, sw4_type p2, sw4_type p3);
  bool owner(sw4_type p1, sw4_type p2, sw4_type p3);
  sw4_type owner_i(sw4_type i);
  void decomp1d(sw4_type nglobal, sw4_type myid, sw4_type nproc, sw4_type& s, sw4_type& e, sw4_type nghost,
                sw4_type npad);
  void decomp1d_2(sw4_type N, sw4_type myid, sw4_type nproc, sw4_type& s, sw4_type& e, sw4_type nghost,
                  sw4_type npad);
  void decomp1d_frf(sw4_type N, sw4_type myid, sw4_type nproc, sw4_type& s, sw4_type& e, sw4_type nghost,
                    sw4_type npad, sw4_type Nf, sw4_type sf, sw4_type ef, sw4_type nghostf, sw4_type npadf);
};

#endif
