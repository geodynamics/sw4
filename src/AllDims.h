#ifndef SW4_ALLDIMS_H
#define SW4_ALLDIMS_H

#include <mpi.h>
#include <sys/types.h>

#include <vector>

class AllDims {
 public:
  int m_myid1d;
  int m_myid3di, m_myid3dj, m_myid3dk;

  int m_nproci, m_nprocj, m_nprock;
  int m_ibg, m_ieg, m_jbg, m_jeg, m_kbg,
      m_keg;  // Global dims, includes ghost points

  std::vector<int> m_ib, m_ie, m_jb, m_je, m_kb,
      m_ke;  // Local for each processor, includes padding and ghost points
  int m_npad, m_nghost;
  bool m_indrev;
  size_t m_fftw_alloc_local;
  MPI_Comm m_communicator;

  AllDims(int nproci, int nprocj, int nprock, int ibg, int ieg, int jbg,
          int jeg, int kbg, int keg, int nghost, int npad, MPI_Comm ewcomm);

  AllDims(int nprocs, int ibg, int ieg, int jbg, int jeg, int kbg, int keg,
          int nghost, MPI_Comm ewcomm);

  AllDims(int nproci, int nprocj, int nprock, int ibg, int ieg, int jbg,
          int jeg, int kbg, int keg, int nghost, int npad);
  AllDims(int nprocs, int ibg, int ieg, int jbg, int jeg, int kbg, int keg,
          int nghost);
  AllDims(AllDims* fine, int ibg, int ieg, int jbg, int jeg, int kbg, int keg,
          int nghost, int npad);

  void getdims_nopad(int dims[6], int p1 = -1, int p2 = -1, int p3 = -1);
  bool intersect(int p1, int p2, int p3, AllDims& other, int dims[6]);
  void compute_myid3d();
  int proc1d(int p1, int p2, int p3);
  bool owner(int p1, int p2, int p3);
  int owner_i(int i);
  void decomp1d(int nglobal, int myid, int nproc, int& s, int& e, int nghost,
                int npad);
  void decomp1d_2(int N, int myid, int nproc, int& s, int& e, int nghost,
                  int npad);
  void decomp1d_frf(int N, int myid, int nproc, int& s, int& e, int nghost,
                    int npad, int Nf, int sf, int ef, int nghostf, int npadf);
};

#endif
