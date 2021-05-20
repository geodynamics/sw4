#ifndef __FARRAY_H__
#define __FARRAY_H__
#include <initializer_list>
#include <memory>
#include "Sarray.h"

class Farray {
 public:
  Farray() {}

  Farray(int ibeg, int iend) {
    len[0] = iend - ibeg + 1;
    size = len[0];
    data = new float_sw4[size];
    start[0] = ibeg;
    end[0] = iend;
    base = -start[0];
    count++;
    id = count;
    dims = 1;
    owns_data = true;
  }

  Farray(int ibeg, int iend, int jbeg, int jend) {
    len[0] = iend - ibeg + 1;
    len[1] = jend - jbeg + 1;
    size = len[0] * len[1];
    data = new float_sw4[size];
    start[0] = ibeg;
    end[0] = iend;
    start[1] = jbeg;
    end[1] = jend;
    base = -start[0] - start[1] * len[0];
    count++;
    id = count;
    dims = 2;
    owns_data = true;
  }

  Farray(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend) {
    start[0] = ibeg;
    end[0] = iend;

    start[1] = jbeg;
    end[1] = jend;

    start[2] = kbeg;
    end[2] = kend;

    for (int i = 0; i < 3; i++) len[i] = (end[i] - start[i] + 1);
    size = 1;
    for (int i = 0; i < 3; i++) size *= len[i];
    data = new float_sw4[size];

    base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1];
    count++;
    id = count;
    dims = 3;
    owns_data = true;
  }

  Farray(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg,
         int lend) {
    start[0] = ibeg;
    end[0] = iend;

    start[1] = jbeg;
    end[1] = jend;

    start[2] = kbeg;
    end[2] = kend;

    start[3] = lbeg;
    end[3] = lend;

    for (int i = 0; i < 4; i++) len[i] = (end[i] - start[i] + 1);
    size = 1;
    for (int i = 0; i < 4; i++) size *= len[i];
    data = new float_sw4[size];

    base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1] -
           start[3] * len[0] * len[1] * len[2];
    count++;
    id = count;
    dims = 4;
    owns_data = true;
  }

  Farray(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int lbeg,
         int lend, int mbeg, int mend) {
    start[0] = ibeg;
    end[0] = iend;

    start[1] = jbeg;
    end[1] = jend;

    start[2] = kbeg;
    end[2] = kend;

    start[3] = lbeg;
    end[3] = lend;

    start[4] = mbeg;
    end[4] = mend;

    for (int i = 0; i < 5; i++) len[i] = (end[i] - start[i] + 1);
    size = 1;
    for (int i = 0; i < 5; i++) size *= len[i];
    data = new float_sw4[size];

    base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1] -
           start[3] * len[0] * len[1] * len[2] -
           start[4] * len[0] * len[1] * len[2] * len[3];
    count++;
    id = count;
    dims = 5;
    owns_data = true;
  }

  void define(int ibeg, int iend) {
    len[0] = iend - ibeg + 1;
    size = len[0];
    data = new float_sw4[size];
    start[0] = ibeg;
    end[0] = iend;
    base = -start[0];
    count++;
    id = count;
    dims = 1;
    owns_data = true;
  }

  void define(int ibeg, int iend, int jbeg, int jend) {
    len[0] = iend - ibeg + 1;
    len[1] = jend - jbeg + 1;
    size = len[0] * len[1];
    data = new float_sw4[size];
    start[0] = ibeg;
    end[0] = iend;
    start[1] = jbeg;
    end[1] = jend;
    base = -start[0] - start[1] * len[0];
    count++;
    id = count;
    dims = 2;
    owns_data = true;
  }

  void define(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend) {
    start[0] = ibeg;
    end[0] = iend;

    start[1] = jbeg;
    end[1] = jend;

    start[2] = kbeg;
    end[2] = kend;

    for (int i = 0; i < 3; i++) len[i] = (end[i] - start[i] + 1);
    size = 1;
    for (int i = 0; i < 3; i++) size *= len[i];
    data = new float_sw4[size];

    base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1];
    count++;
    id = count;
    dims = 3;
    owns_data = true;
  }
  void define(int ibeg, int iend, int jbeg, int jend, int kbeg, int kend,
              int lbeg, int lend) {
    start[0] = ibeg;
    end[0] = iend;

    start[1] = jbeg;
    end[1] = jend;

    start[2] = kbeg;
    end[2] = kend;

    start[3] = lbeg;
    end[3] = lend;

    for (int i = 0; i < 4; i++) len[i] = (end[i] - start[i] + 1);
    size = 1;
    for (int i = 0; i < 4; i++) size *= len[i];
    data = new float_sw4[size];

    base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1] -
           start[3] * len[0] * len[1] * len[2];
    count++;
    id = count;
    dims = 4;
    owns_data = true;
  }

  ~Farray() {
    // std::cout<<"Deleting "<<id<<"\n"<<std::flush;
    if (owns_data)
      delete[] data;
    else
      std::cout << "Null op for subset data" << id << "\n";
  }

  inline float_sw4& operator()(int i) { return data[base + i]; }
  inline float_sw4& operator[](int i) { return data[i]; }

  inline void operator=(Farray& in) {
    if (size != in.size) {
      std::cerr << "Assignement of incompatible sizes " << size
                << "!=" << in.size << "\n";
      abort();
    }
    for (size_t i = 0; i < size; i++) data[i] = in.data[i];
  }

  inline Farray& operator*=(Farray& in) {
#ifdef DEBUG_MODE
    if (size != in.size) {
      std::cerr << "operator *= of incompatible sizes " << size
                << "!=" << in.size << "\n";
      abort();
    }
#endif
    for (size_t i = 0; i < size; i++) data[i] *= in.data[i];
    return *this;
  }

  inline void operator=(double in) {
    for (size_t i = 0; i < size; i++) data[i] = in;
  }

  inline float_sw4& operator()(int i, int j) {
#ifdef DEBUG_MODE
    if (dims != 2) {
      std::cerr << "ERROR OPERATOR2 " << dims << "\n";
      abort();
    }
#endif
    return data[base + i + j * len[0]];
  }

  inline float_sw4& operator()(int i, int j, int k) {
#ifdef DEBUG_MODE
    if (dims != 3) {
      std::cerr << "ERROR OPERATOR3 " << dims << "\n";
      abort();
    }
#endif
    return data[base + i + j * len[0] + k * len[0] * len[1]];
  }

  inline float_sw4& operator()(int i, int j, int k, int l) {
#ifdef DEBUG_MODE
    if (dims != 4) {
      std::cerr << "ERROR OPERATOR4 " << dims << "\n";
      abort();
    }
#endif
    return data[base + i + j * len[0] + k * len[0] * len[1] +
                l * len[0] * len[1] * len[2]];
  }

  inline float_sw4& operator()(int i, int j, int k, int l, int m) {
#ifdef DEBUG_MODE
    if (dims != 5) {
      std::cerr << "ERROR OPERATOR5 " << dims << "\n";
      abort();
    }
#endif
    return data[base + i + j * len[0] + k * len[0] * len[1] +
                l * len[0] * len[1] * len[2] +
                m * len[0] * len[1] * len[2] * len[3]];
  }

  inline float_sw4 max() {
    float_sw4 maxval = data[0];
    for (size_t i = 1; i < size; i++)
      maxval = (maxval < data[i]) ? data[i] : maxval;
    return maxval;
  }

  inline float_sw4* get() { return data; }
  inline float_sw4* c_ptr() { return data; }
  inline int debut(int n = 0) { return start[n]; }
  inline int fin(int n = 0) { return end[n]; }

  void print() {
    std::cout << "============================================================="
                 "==========\n";
    std::cout << "Farray ID " << id << " :: \n";
    std::cout << "============================================================="
                 "==========\n";
    if (dims == 1) {
      for (int i = start[0]; i <= end[0]; i++)
        std::cout << this->operator()(i) << ",";
      std::cout << "\n" << std::flush;
    } else {
      for (int i = start[0]; i <= end[0]; i++) {
        for (int j = start[1]; j <= end[1]; j++)
          std::cout << this->operator()(i, j) << ",";
        std::cout << "\n" << std::flush;
      }
    }
    std::cout << "============================================================="
                 "==========\n";
  }

  double maxabs() {
    double maxval = -1.0;
    for (size_t i = 0; i < size; i++)
      maxval = (fabs(data[i]) > maxval) ? fabs(data[i]) : maxval;
    return maxval;
  }

  // Hardwired for going from 5 dims to 3 dims
  std::shared_ptr<Farray> subset(int l, int m) {
    std::shared_ptr<Farray> ret = std::make_shared<Farray>();
    ret->dims = dims - 2;
    ret->owns_data = false;
    ret->data = &(this->operator()(start[0], start[1], start[2], l, m));
    ret->data = data + base + start[0] + start[1] * len[0] +
                start[2] * len[0] * len[1] + l * len[0] * len[1] * len[2] +
                m * len[0] * len[1] * len[2] * len[3];
    ret->size = 1;
    for (int i = 0; i < ret->dims; i++) {
      ret->start[i] = start[i];
      ret->end[i] = end[i];
      ret->len[i] = len[i];
      ret->size *= len[i];
    }
    ret->base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1];
    count++;
    ret->id = count;
    return ret;
  }

  std::shared_ptr<Farray> subset(int m) {
    std::shared_ptr<Farray> ret = std::make_shared<Farray>();
    ret->dims = dims - 1;
    ret->owns_data = false;
    ret->data = &(this->operator()(start[0], start[1], start[2], start[3], m));
    ret->data = data + base + start[0] + start[1] * len[0] +
                start[2] * len[0] * len[1] +
                start[3] * len[0] * len[1] * len[2] +
                m * len[0] * len[1] * len[2] * len[3];
    ret->size = 1;
    for (int i = 0; i < ret->dims; i++) {
      ret->start[i] = start[i];
      ret->end[i] = end[i];
      ret->len[i] = len[i];
      ret->size *= len[i];
    }
    ret->base = -start[0] - start[1] * len[0] - start[2] * len[0] * len[1] -
                start[3] * len[0] * len[1] * len[2];
    count++;
    ret->id = count;
    return ret;
  }

  void compare(Farray& b) {
    for (int i = start[0]; i <= end[0]; i++) {
      for (int j = start[1]; j <= end[1]; j++) {
        for (int k = start[2]; k <= end[2]; k++) {
          for (int l = start[3]; l <= end[3]; l++) {
            if (this->operator()(i, j, k, l, 3) != b(i, j, k, l, 3)) {
              std::cout << "NOPE" << i << " " << j << " " << k << " " << l
                        << " " << this->operator()(i, j, k, l, 3) << " "
                        << b(i, j, k, l, 3) << "\n";
            }
          }
        }
      }
    }
  }

  bool owns_data;

  int start[5];
  int end[5];
  int len[5];
  int base;
  float_sw4* data;
  static int count;
  int id;
  int dims;
  size_t size;
};

template <int N>
class Garray {
 public:
  Garray(std::initializer_list<int> l) {
    for (auto t = l.begin(); t != l.end(); t += 2) {
      std::cout << *t << " " << *(t + 1) << "\n";
    }
  }
  int start[N];
  int end[N];
  int len[N];
  double* data;
};

class PackArgs {
 public:
  int dim, nrg, n1_c, n2_c, n3_c, n1_f, n2_f, n3_f;
  float_sw4 tn, pi, l1, l2, l3, int_pos, h1phy_c, h1phy_f, h2phy_c, h2phy_f,
      h1_c, h2_c, h3_c, h1_f, h2_f, h3_f, amp, peak;
};

#endif
