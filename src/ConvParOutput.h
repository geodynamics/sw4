#ifndef CONVPAROUTPUT_H
#define CONVPAROUTPUT_H
#include <vector>
#include "EW.h"

using namespace std;

class ConvParOutput {
 private:
  FILE* m_fd;
  FILE* m_fdx;
  bool m_stepfileio;
  int m_myrank;
  int m_verbose;
  vector<int> m_iconvdata;
  vector<double> m_convdata;
  vector<double> m_paradata;
  bool m_cg;
  int m_ind, m_n, m_nvar, m_nsize;
  string m_convfile, m_parafile;

 public:
  ConvParOutput(EW& ew, int n, int nvar, int myrank, int verbose, bool cg);
  ~ConvParOutput();
  void print_dfmsg(double f, double* df, double* sf);
  void print_xmsg(double* x, double rnorm, double dxnorm, int j, int k = 0);
  void print_vector(double* v, const char* msg, int m);
  void print_scalar(double s, const char* msg, int m);
  void save_step(double f, double* df, double* sf, double* x, double rnorm,
                 double dxnorm, int nreductions, int j, int k = 0);
  void finish();
};
#endif
