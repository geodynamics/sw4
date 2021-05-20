#include "mpi.h"

#include <unistd.h>
#include <algorithm>
#include <iostream>
#include <sstream>

// mking directories
#include <errno.h>
#include <sys/stat.h>
#include <list>
#include <sstream>

#include "EW.h"
#include "F77_FUNC.h"

using namespace std;

#define SQR(x) ((x) * (x))

template <int iu, int il, int ju, int jl, int ku, int kl>
void evalLuCurv(int ib, int ie, int jb, int je, int kb, int ke, Sarray& u,
                Sarray& lu, float_sw4* a_mu, float_sw4* a_la, Sarray& met,
                Sarray& jac, int ilb, int ile, int jlb, int jle, int klb,
                int kle);

void evalLu_Dip(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_Dim(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_Djp(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_Djm(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_Dkp(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_Dkm(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_DkpDip(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_DkpDim(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_DkpDjp(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);
void evalLu_DkpDjm(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& a_u, Sarray& a_lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle);

//-----------------------------------------------------------------------
void EW::set_geodyn_data(string file, int nx, int nz, float_sw4 h,
                         float_sw4 origin[3], float_sw4 dt, int nsteps,
                         int faces) {
  m_do_geodynbc = true;
  m_geodyn_past_end = false;
  m_geodyn_filename = file;
  m_geodyn_ni = m_geodyn_nj = nx;
  m_geodyn_nk = nz;
  m_geodyn_h = h;
  m_geodyn_dt = dt;
  m_geodyn_origin[0] = origin[0];
  m_geodyn_origin[1] = origin[1];
  m_geodyn_origin[2] = origin[2];
  m_geodyn_faces = faces;
  m_geodyn_maxsteps = nsteps;

  int ny = nx;
  if (m_geodyn_faces == 6)
    m_geodyn_blocksize = 2 * (nx * ny + nx * nz + ny * nz);
  else if (m_geodyn_faces == 5)
    m_geodyn_blocksize = (nx * ny + 2 * nx * nz + 2 * ny * nz);

  m_geodyn_data1.resize(6);
  m_geodyn_data2.resize(6);

  m_geodyn_data1[0].define(3, ny, nz, 1);
  m_geodyn_data1[1].define(3, ny, nz, 1);
  m_geodyn_data2[0].define(3, ny, nz, 1);
  m_geodyn_data2[1].define(3, ny, nz, 1);

  m_geodyn_data1[2].define(3, nx, nz, 1);
  m_geodyn_data1[3].define(3, nx, nz, 1);
  m_geodyn_data2[2].define(3, nx, nz, 1);
  m_geodyn_data2[3].define(3, nx, nz, 1);

  m_geodyn_data1[4].define(3, nx, ny, 1);
  m_geodyn_data1[5].define(3, nx, ny, 1);
  m_geodyn_data2[4].define(3, nx, ny, 1);
  m_geodyn_data2[5].define(3, nx, ny, 1);

  m_geo_usgh.resize(4);

  int i0, i1, j0, j1, k0, k1;
  float_sw4 cubelen = (nx - 1) * m_geodyn_h;
  float_sw4 zcubelen = (nz - 1) * m_geodyn_h;

  //   m_geodyn_dims.resize(mNumberOfCartesianGrids);
  m_geodyn_dims.resize(mNumberOfGrids);
  m_geodyn_iwillread = false;

  for (int g = 0; g < mNumberOfGrids; g++) {
    i0 = static_cast<int>(round(m_geodyn_origin[0] / mGridSize[g] + 1));
    i1 = static_cast<int>(
        round((m_geodyn_origin[0] + cubelen) / mGridSize[g] + 1));
    j0 = static_cast<int>(round(m_geodyn_origin[1] / mGridSize[g] + 1));
    j1 = static_cast<int>(
        round((m_geodyn_origin[1] + cubelen) / mGridSize[g] + 1));
    k0 = static_cast<int>(
        round((m_geodyn_origin[2] - m_zmin[g]) / mGridSize[g] + 1));
    k1 = static_cast<int>(
        round((m_geodyn_origin[2] - m_zmin[g] + zcubelen) / mGridSize[g] + 1));

    if (g == mNumberOfGrids - 1 &&
        topographyExists())  // NOT verified for several curvilinear grids
    {
      // Curvilinear grid, inaccurate quick fix.
      // Assume upper side of cube is at the free surface
      k0 = 1;

      // Find k=const grid surface with smallest distance to the plane
      // z=cubelen.
      int icmin = m_iStartInt[g], icmax = m_iEndInt[g], jcmin = m_jStartInt[g],
          jcmax = m_jEndInt[g];
      if (i0 > icmin) icmin = i0;
      if (i1 < icmax) icmax = i1;
      if (j0 > jcmin) jcmin = j0;
      if (j1 < jcmax) jcmax = j1;

      float_sw4 kavgm = 0, kavgp = 0;
      int nptsij = 0;
      for (int j = jcmin; j <= jcmax; j++)
        for (int i = icmin; i <= icmax; i++) {
          int k = 1;
          while (k <= m_kEnd[g] && mZ[g](i, j, k) - mZ[g](i, j, 1) < zcubelen)
            k++;
          kavgm += k - 1;
          kavgp += k;
          nptsij++;
        }
      float_sw4 ktmp = kavgm;
      MPI_Allreduce(&ktmp, &kavgm, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
      ktmp = kavgp;
      MPI_Allreduce(&ktmp, &kavgp, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
      ktmp = nptsij;
      float_sw4 nptstot;
      MPI_Allreduce(&ktmp, &nptstot, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
      int km = static_cast<int>(round(kavgm / nptstot));
      int kp = static_cast<int>(round(kavgp / nptstot));

      float_sw4 deperrp = 0, deperrm = 0;
      if (nptsij > 0) {
        if (kp > m_kEnd[g]) kp = m_kEnd[g];
        if (km > m_kEnd[g]) km = m_kEnd[g];
        for (int j = jcmin; j <= jcmax; j++)
          for (int i = icmin; i <= icmax; i++) {
            deperrp += (mZ[g](i, j, kp) - mZ[g](i, j, 1) - zcubelen) *
                       (mZ[g](i, j, kp) - mZ[g](i, j, 1) - zcubelen);
            deperrm += (mZ[g](i, j, km) - mZ[g](i, j, 1) - zcubelen) *
                       (mZ[g](i, j, km) - mZ[g](i, j, 1) - zcubelen);
          }
      }
      float_sw4 errp, errm;
      MPI_Allreduce(&deperrp, &errp, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&deperrm, &errm, 1, m_mpifloat, MPI_SUM, MPI_COMM_WORLD);
      if (nptsij > 0) {
        if (errp < errm)
          k1 = kp;
        else
          k1 = km;
      } else
        k1 = 0;  // k1 =0 if cube not in my processor
      int k1tmp = k1;
      MPI_Allreduce(&k1tmp, &k1, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    bool cubeok = true;
    if (g == mNumberOfGrids - 1 && m_geodyn_faces == 5) {
      if (k0 != 1) {
        cubeok = false;
        if (m_myRank == 0)
          printf("Error in geodyn setup, k0 = %d, but number of faces is 5\n",
                 k0);
      }
    }
    if (g == mNumberOfGrids - 1) {
      if (1 < k0 && k0 < 5) {
        cubeok = false;
        if (m_myRank == 0)
          cout << "Error: Geodyn cube upper boundary at k = " << k0
               << ", is too close to the surface \n";
      }
      if (k0 == 5 || k0 == 6)
        if (m_myRank == 0)
          cout << "Warning: Geodyn cube upper boundary is close to surface, "
                  "simulation will run but results might be unreliable"
               << endl;
      if (k1 < 6) {
        cubeok = false;
        if (m_myRank == 0)
          cout << "Error: Geodyn cube bottom boundary at k = " << k1
               << ", is too close to the surface \n";
      }
    }
    if (m_myRank == 0)
      cout << "Geodyn cube dims [" << i0 << ".." << i1 << "]x[" << j0 << ".."
           << j1 << "]x[" << k0 << ".." << k1 << "]" << endl;
    VERIFY2(cubeok,
            "SW4 exiting because Geodyn cube placement is too close to the "
            "surface");

    m_geodyn_dims[g] = new int[6];
    int imin = m_iStart[g], jmin = m_jStart[g], kmin = m_kStart[g];
    int imax = m_iEnd[g], jmax = m_jEnd[g], kmax = m_kEnd[g];
    if ((i0 <= imax && i1 >= imin) && (j0 <= jmax && j1 >= jmin) &&
        (k0 <= kmax && k1 >= kmin)) {
      if (i0 < imin) i0 = imin;
      if (i1 > imax) i1 = imax;
      if (j0 < jmin) j0 = jmin;
      if (j1 > jmax) j1 = jmax;
      if (k0 < kmin) k0 = kmin;
      if (k1 > kmax) k1 = kmax;

      // Store index bounds for cube in proc
      m_geodyn_dims[g][0] = i0;
      m_geodyn_dims[g][1] = i1;
      m_geodyn_dims[g][2] = j0;
      m_geodyn_dims[g][3] = j1;
      m_geodyn_dims[g][4] = k0;
      m_geodyn_dims[g][5] = k1;
      m_geodyn_iwillread = true;
      if (g == mNumberOfGrids - 1 && k0 == 1) {
        m_geo_usgh[0] = new float_sw4[3 * (j1 - j0 + 1)];
        m_geo_usgh[1] = new float_sw4[3 * (j1 - j0 + 1)];
        m_geo_usgh[2] = new float_sw4[3 * (i1 - i0 + 1)];
        m_geo_usgh[3] = new float_sw4[3 * (i1 - i0 + 1)];
      }
    } else {
      // Empty cube
      m_geodyn_dims[g][0] = 0;
      m_geodyn_dims[g][1] = -1;
      m_geodyn_dims[g][2] = 0;
      m_geodyn_dims[g][3] = -1;
      m_geodyn_dims[g][4] = 0;
      m_geodyn_dims[g][5] = -1;
      //	 m_geodyn_iwillread = false;
    }
  }
  if (m_geodyn_iwillread) {
    m_geodynfile.open(m_geodyn_filename.c_str());
    VERIFY2(m_geodynfile.is_open(), "Error opening Geodyn input file " << file);
    bool done = false;
    char buffer[256];
    while (!m_geodynfile.eof() && !done) {
      m_geodynfile.getline(buffer, 256);
      if (startswith("begindata", buffer)) done = true;
    }
    m_geodyn_step = -2;
  }
}

//-----------------------------------------------------------------------
void EW::impose_geodyn_ibcdata(vector<Sarray>& u, vector<Sarray>& um,
                               float_sw4 t, vector<float_sw4**>& bforcing) {
  //   int n1 = static_cast<int>(floor(t/m_geodyn_dt));
  int i0, i1, j0, j1, k0, k1;

  // Zero out values inside the cube. Use -100000 for debug.
  for (int g = 0; g < mNumberOfGrids; g++) {
    i0 = m_geodyn_dims[g][0];
    i1 = m_geodyn_dims[g][1];
    j0 = m_geodyn_dims[g][2];
    j1 = m_geodyn_dims[g][3];
    k0 = m_geodyn_dims[g][4];
    k1 = m_geodyn_dims[g][5];
    for (int k = k0 + 1; k <= k1 - 1; k++)
      for (int j = j0 + 1; j <= j1 - 1; j++)
        for (int i = i0 + 1; i <= i1 - 1; i++)
          for (int c = 1; c <= 3; c++) {
            //		      u[g](c,i,j,k) = -100000;
            u[g](c, i, j, k) = 0;
          }
  }
  if (m_geodyn_past_end)
  //   if( n1 > m_geodyn_maxsteps-1 )
  {
    // Past geodyn end time, switch off geodyn boundary, impose zero values.
    // A volume reset is necessary in the debug -100000 case above.
    // Note: need to reset both n+1 and n levels.
    m_do_geodynbc = false;
    if (m_myRank == 0) cout << "Switching off Geodyn boundary " << endl;
    for (int g = 0; g < mNumberOfGrids; g++) {
      i0 = m_geodyn_dims[g][0];
      i1 = m_geodyn_dims[g][1];
      j0 = m_geodyn_dims[g][2];
      j1 = m_geodyn_dims[g][3];
      k0 = m_geodyn_dims[g][4];
      k1 = m_geodyn_dims[g][5];
#pragma omp parallel for
      for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++)
          for (int i = i0; i <= i1; i++)
            for (int c = 1; c <= 3; c++)
              u[g](c, i, j, k) = um[g](c, i, j, k) = 0;
    }
  } else {
    // // Find the right time step
    // if( m_geodyn_step+1 == n1 )
    // {
    // 	 // Advance one step
    //    copy_geodyn_timelevel( m_geodyn_data1, m_geodyn_data2 );
    //    get_geodyn_timelevel( m_geodyn_data2 );
    // }
    // else if( m_geodyn_step+1 < n1 )
    // {
    // 	 // Advance many steps
    //    if( m_geodyn_iwillread )
    // 	 {
    // 	    char buf[256];
    // 	    for( int i=0 ; i < m_geodyn_blocksize*(n1-m_geodyn_step-2) ; i++ )
    // 	       m_geodynfile.getline(buf,256);
    // 	 }
    //    get_geodyn_timelevel( m_geodyn_data1 );
    //    get_geodyn_timelevel( m_geodyn_data2 );
    // }
    // m_geodyn_step = n1;

    //      double twgh = ((n1+1)*m_geodyn_dt-t)/m_geodyn_dt;
    //      double cext1 = 2, cext2=-1, cext3=0;

    float_sw4 twgh = ((m_geodyn_step + 1) * m_geodyn_dt - t) / m_geodyn_dt;

    for (int g = 0; g < mNumberOfCartesianGrids; g++) {
      i0 = m_geodyn_dims[g][0];
      i1 = m_geodyn_dims[g][1];
      j0 = m_geodyn_dims[g][2];
      j1 = m_geodyn_dims[g][3];
      k0 = m_geodyn_dims[g][4];
      k1 = m_geodyn_dims[g][5];
      bool at_surface = g == mNumberOfGrids - 1 && k0 <= 1;

      float_sw4 h = mGridSize[g];
#pragma omp parallel for
      for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++) {
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(
              floor(((k - 1) * h - m_geodyn_origin[2]) / m_geodyn_h + 1));
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghk =
              ((k - 1) * h - (m_geodyn_origin[2] + (kg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;

          for (int c = 1; c <= 3; c++) {
            u[g](c, i0, j, k) =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[0](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[0](c, jg0 + 1, kg0 + 1, 1));
            u[g](c, i1, j, k) =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[1](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[1](c, jg0 + 1, kg0 + 1, 1));
          }
        }
#pragma omp parallel for
      for (int k = k0; k <= k1; k++)
        for (int i = i0; i <= i1; i++) {
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(
              floor(((k - 1) * h - m_geodyn_origin[2]) / m_geodyn_h + 1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghk =
              ((k - 1) * h - (m_geodyn_origin[2] + (kg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;

          for (int c = 1; c <= 3; c++) {
            u[g](c, i, j0, k) =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[2](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[2](c, ig0 + 1, kg0 + 1, 1));
            u[g](c, i, j1, k) =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[3](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[3](c, ig0 + 1, kg0 + 1, 1));
          }
        }
      //	 if( k0==1 && m_geodyn_faces == 5)
      //	    for( int i=i0 ;  i<=i1 ; i++ )
      //	       for( int c=1 ; c <= 3 ; c++ )
      //	       {
      //		  u[g](c,i,j0,k0-1) =
      //cext1*u[g](c,i,j0,k0)+cext2*u[g](c,i,j0,k0+1)+cext3*u[g](c,i,j0,k0+2);
      //		  u[g](c,i,j1,k0-1) =
      //cext1*u[g](c,i,j1,k0)+cext2*u[g](c,i,j1,k0+1)+cext3*u[g](c,i,j1,k0+2);
      //	       }

      Sarray& gd14 = m_geodyn_data1[4];
      Sarray& gd24 = m_geodyn_data2[4];
#pragma omp parallel for
      for (int j = j0; j <= j1; j++)
        for (int i = i0; i <= i1; i++) {
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          for (int c = 1; c <= 3; c++) {
            if (!at_surface && m_geodyn_faces == 6) {
              //		     u[g](c,i,j,k0) = twgh*(
              //(1-wghi)*(1-wghj)*m_geodyn_data1[4](c,ig0,jg0,1)+
              //					  wghi*(1-wghj)*m_geodyn_data1[4](c,ig0+1,jg0,1)
              //					  (1-wghi)*wghj*m_geodyn_data1[4](c,ig0,jg0+1,1)+
              //					     wghi*wghj*m_geodyn_data1[4](c,ig0+1,jg0+1,1)
              //)
              //		     +
              //(1-twgh)*((1-wghi)*(1-wghj)*m_geodyn_data2[4](c,ig0,jg0,1)+
              //				 wghi*(1-wghj)*m_geodyn_data2[4](c,ig0+1,jg0,1)
              //				 (1-wghi)*wghj*m_geodyn_data2[4](c,ig0,jg0+1,1)+
              //				 wghi*wghj*m_geodyn_data2[4](c,ig0+1,jg0+1,1)
              //));
              u[g](c, i, j, k0) =
                  twgh * ((1 - wghi) * (1 - wghj) * gd14(c, ig0, jg0, 1) +
                          wghi * (1 - wghj) * gd14(c, ig0 + 1, jg0, 1) +
                          (1 - wghi) * wghj * gd14(c, ig0, jg0 + 1, 1) +
                          wghi * wghj * gd14(c, ig0 + 1, jg0 + 1, 1)) +
                  (1 - twgh) * ((1 - wghi) * (1 - wghj) * gd24(c, ig0, jg0, 1) +
                                wghi * (1 - wghj) * gd24(c, ig0 + 1, jg0, 1) +
                                (1 - wghi) * wghj * gd24(c, ig0, jg0 + 1, 1) +
                                wghi * wghj * gd24(c, ig0 + 1, jg0 + 1, 1));
            }
            u[g](c, i, j, k1) =
                twgh *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data1[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data1[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data1[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data1[5](c, ig0 + 1, jg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data2[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data2[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data2[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data2[5](c, ig0 + 1, jg0 + 1, 1));
          }
        }
      // Impose ghost point at corner to free surface
      if (at_surface) {
        int ib = m_iStart[g];
        int jb = m_jStart[g];
        int ni = m_iEnd[g] - m_iStart[g] + 1;
        // Free surface condition
        for (int j = j0; j <= j1; j++) {
          int i = i0;
          size_t qq = i - ib + ni * (j - jb);
          // One sided x-derivatives
          float_sw4 wx = u[g](3, i, j, k0) - u[g](3, i - 1, j, k0);
          float_sw4 ux = u[g](1, i, j, k0) - u[g](1, i - 1, j, k0);
          float_sw4 wy = 0.5 * (u[g](3, i, j + 1, k0) - u[g](3, i, j - 1, k0));
          float_sw4 vy = 0.5 * (u[g](2, i, j + 1, k0) - u[g](2, i, j - 1, k0));
          float_sw4 mup = 0.5 * (mMu[g](i, j, k0) + mMu[g](i, j, k0 + 1));
          float_sw4 mum = 0.5 * (mMu[g](i, j, k0 - 1) + mMu[g](i, j, k0));

          u[g](1, i, j, k0 - 1) =
              u[g](1, i, j, k0) +
              (mup * (u[g](1, i, j, k0 + 1) - u[g](1, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wx - 2 * h * bforcing[g][4][3 * qq]) /
                  mum;
          u[g](2, i, j, k0 - 1) =
              u[g](2, i, j, k0) +
              (mup * (u[g](2, i, j, k0 + 1) - u[g](2, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wy - 2 * h * bforcing[g][4][1 + 3 * qq]) /
                  mum;
          float_sw4 lap =
              0.5 * (mLambda[g](i, j, k0) + mLambda[g](i, j, k0 + 1));
          float_sw4 lam =
              0.5 * (mLambda[g](i, j, k0 - 1) + mLambda[g](i, j, k0));
          u[g](3, i, j, k0 - 1) =
              u[g](3, i, j, k0) +
              ((2 * mup + lap) * (u[g](3, i, j, k0 + 1) - u[g](3, i, j, k0)) +
               2 * mLambda[g](i, j, k0) * (ux + vy) -
               2 * h * bforcing[g][4][2 + 3 * qq]) /
                  (2 * mum + lam);
          i = i1;
          qq = i - ib + ni * (j - jb);
          // One sided x-derivatives
          wx = u[g](3, i + 1, j, k0) - u[g](3, i, j, k0);
          ux = u[g](1, i + 1, j, k0) - u[g](1, i, j, k0);
          wy = 0.5 * (u[g](3, i, j + 1, k0) - u[g](3, i, j - 1, k0));
          vy = 0.5 * (u[g](2, i, j + 1, k0) - u[g](2, i, j - 1, k0));
          mup = 0.5 * (mMu[g](i, j, k0) + mMu[g](i, j, k0 + 1));
          mum = 0.5 * (mMu[g](i, j, k0 - 1) + mMu[g](i, j, k0));

          u[g](1, i, j, k0 - 1) =
              u[g](1, i, j, k0) +
              (mup * (u[g](1, i, j, k0 + 1) - u[g](1, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wx - 2 * h * bforcing[g][4][3 * qq]) /
                  mum;
          u[g](2, i, j, k0 - 1) =
              u[g](2, i, j, k0) +
              (mup * (u[g](2, i, j, k0 + 1) - u[g](2, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wy - 2 * h * bforcing[g][4][1 + 3 * qq]) /
                  mum;
          lap = 0.5 * (mLambda[g](i, j, k0) + mLambda[g](i, j, k0 + 1));
          lam = 0.5 * (mLambda[g](i, j, k0 - 1) + mLambda[g](i, j, k0));
          u[g](3, i, j, k0 - 1) =
              u[g](3, i, j, k0) +
              ((2 * mup + lap) * (u[g](3, i, j, k0 + 1) - u[g](3, i, j, k0)) +
               2 * mLambda[g](i, j, k0) * (ux + vy) -
               2 * h * bforcing[g][4][2 + 3 * qq]) /
                  (2 * mum + lam);
        }
        for (int i = i0; i <= i1; i++) {
          int j = j0;
          size_t qq = i - ib + ni * (j - jb);
          // One sided y-derivatives
          float_sw4 wx = 0.5 * (u[g](3, i + 1, j, k0) - u[g](3, i - 1, j, k0));
          float_sw4 ux = 0.5 * (u[g](1, i + 1, j, k0) - u[g](1, i - 1, j, k0));
          float_sw4 wy = (u[g](3, i, j, k0) - u[g](3, i, j - 1, k0));
          float_sw4 vy = (u[g](2, i, j, k0) - u[g](2, i, j - 1, k0));
          float_sw4 mup = 0.5 * (mMu[g](i, j, k0) + mMu[g](i, j, k0 + 1));
          float_sw4 mum = 0.5 * (mMu[g](i, j, k0 - 1) + mMu[g](i, j, k0));

          u[g](1, i, j, k0 - 1) =
              u[g](1, i, j, k0) +
              (mup * (u[g](1, i, j, k0 + 1) - u[g](1, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wx - 2 * h * bforcing[g][4][3 * qq]) /
                  mum;
          u[g](2, i, j, k0 - 1) =
              u[g](2, i, j, k0) +
              (mup * (u[g](2, i, j, k0 + 1) - u[g](2, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wy - 2 * h * bforcing[g][4][1 + 3 * qq]) /
                  mum;
          float_sw4 lap =
              0.5 * (mLambda[g](i, j, k0) + mLambda[g](i, j, k0 + 1));
          float_sw4 lam =
              0.5 * (mLambda[g](i, j, k0 - 1) + mLambda[g](i, j, k0));
          u[g](3, i, j, k0 - 1) =
              u[g](3, i, j, k0) +
              ((2 * mup + lap) * (u[g](3, i, j, k0 + 1) - u[g](3, i, j, k0)) +
               2 * mLambda[g](i, j, k0) * (ux + vy) -
               2 * h * bforcing[g][4][2 + 3 * qq]) /
                  (2 * mum + lam);
          j = j1;
          qq = i - ib + ni * (j - jb);
          // One sided y-derivatives
          wx = 0.5 * (u[g](3, i + 1, j, k0) - u[g](3, i - 1, j, k0));
          ux = 0.5 * (u[g](1, i + 1, j, k0) - u[g](1, i - 1, j, k0));
          wy = (u[g](3, i, j + 1, k0) - u[g](3, i, j, k0));
          vy = (u[g](2, i, j + 1, k0) - u[g](2, i, j, k0));
          mup = 0.5 * (mMu[g](i, j, k0) + mMu[g](i, j, k0 + 1));
          mum = 0.5 * (mMu[g](i, j, k0 - 1) + mMu[g](i, j, k0));

          u[g](1, i, j, k0 - 1) =
              u[g](1, i, j, k0) +
              (mup * (u[g](1, i, j, k0 + 1) - u[g](1, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wx - 2 * h * bforcing[g][4][3 * qq]) /
                  mum;
          u[g](2, i, j, k0 - 1) =
              u[g](2, i, j, k0) +
              (mup * (u[g](2, i, j, k0 + 1) - u[g](2, i, j, k0)) +
               2 * mMu[g](i, j, k0) * wy - 2 * h * bforcing[g][4][1 + 3 * qq]) /
                  mum;
          lap = 0.5 * (mLambda[g](i, j, k0) + mLambda[g](i, j, k0 + 1));
          lam = 0.5 * (mLambda[g](i, j, k0 - 1) + mLambda[g](i, j, k0));
          u[g](3, i, j, k0 - 1) =
              u[g](3, i, j, k0) +
              ((2 * mup + lap) * (u[g](3, i, j, k0 + 1) - u[g](3, i, j, k0)) +
               2 * mLambda[g](i, j, k0) * (ux + vy) -
               2 * h * bforcing[g][4][2 + 3 * qq]) /
                  (2 * mum + lam);
        }
      }
    }  // for g=0, mNumberOfCartesianGrids-1
    if (topographyExists()) {
      // Curvilinear, inaccurate quick fix
      int g = mNumberOfGrids - 1;
      i0 = m_geodyn_dims[g][0];
      i1 = m_geodyn_dims[g][1];
      j0 = m_geodyn_dims[g][2];
      j1 = m_geodyn_dims[g][3];
      k0 = m_geodyn_dims[g][4];
      k1 = m_geodyn_dims[g][5];
      float_sw4 h = mGridSize[g];
      float_sw4 zcubelen = (m_geodyn_nk - 1) * m_geodyn_h;
      bool at_surface = k0 == 1;
#pragma omp parallel for
      for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++) {
          float_sw4 strfact = (mZ[g](i0, j, k1) - mZ[g](i0, j, 1)) / zcubelen;
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(floor((mZ[g](i0, j, k) - mZ[g](i0, j, 1)) /
                                               (strfact * m_geodyn_h) +
                                           1));
          //               int kg0 = static_cast<int>(floor(((k-1)*h -
          //               m_geodyn_origin[2])/m_geodyn_h+1));
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghk = ((mZ[g](i0, j, k) - mZ[g](i0, j, 1)) / strfact -
                            ((kg0 - 1) * m_geodyn_h)) /
                           m_geodyn_h;

          for (int c = 1; c <= 3; c++) {
            u[g](c, i0, j, k) =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[0](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[0](c, jg0 + 1, kg0 + 1, 1));
          }
          strfact = (mZ[g](i1, j, k1) - mZ[g](i1, j, 1)) / zcubelen;
          kg0 = static_cast<int>(floor((mZ[g](i1, j, k) - mZ[g](i1, j, 1)) /
                                           (strfact * m_geodyn_h) +
                                       1));
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          wghk = ((mZ[g](i1, j, k) - mZ[g](i1, j, 1)) / strfact -
                  ((kg0 - 1) * m_geodyn_h)) /
                 m_geodyn_h;

          for (int c = 1; c <= 3; c++) {
            u[g](c, i1, j, k) =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[1](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[1](c, jg0 + 1, kg0 + 1, 1));
          }
        }
#pragma omp parallel for
      for (int k = k0; k <= k1; k++)
        for (int i = i0; i <= i1; i++) {
          float_sw4 strfact = (mZ[g](i, j0, k1) - mZ[g](i, j0, 1)) / zcubelen;
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(floor((mZ[g](i, j0, k) - mZ[g](i, j0, 1)) /
                                               (strfact * m_geodyn_h) +
                                           1));
          //	       int kg0 = static_cast<int>(floor(((k-1)*h -
          //m_geodyn_origin[2])/m_geodyn_h+1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          //	       float_sw4 wghk = ((k-1)*h -
          //(m_geodyn_origin[2]+(kg0-1)*m_geodyn_h))/m_geodyn_h;
          float_sw4 wghk = ((mZ[g](i, j0, k) - mZ[g](i, j0, 1)) / strfact -
                            ((kg0 - 1) * m_geodyn_h)) /
                           m_geodyn_h;

          for (int c = 1; c <= 3; c++) {
            u[g](c, i, j0, k) =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[2](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[2](c, ig0 + 1, kg0 + 1, 1));
          }
          strfact = (mZ[g](i, j1, k1) - mZ[g](i, j1, 1)) / zcubelen;
          kg0 = static_cast<int>(floor((mZ[g](i, j1, k) - mZ[g](i, j1, 1)) /
                                           (strfact * m_geodyn_h) +
                                       1));
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          wghk = ((mZ[g](i, j1, k) - mZ[g](i, j1, 1)) / strfact -
                  ((kg0 - 1) * m_geodyn_h)) /
                 m_geodyn_h;
          for (int c = 1; c <= 3; c++) {
            u[g](c, i, j1, k) =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[3](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[3](c, ig0 + 1, kg0 + 1, 1));
          }
        }
      Sarray& gd14 = m_geodyn_data1[4];
      Sarray& gd24 = m_geodyn_data2[4];
#pragma omp parallel for
      for (int j = j0; j <= j1; j++)
        for (int i = i0; i <= i1; i++) {
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          for (int c = 1; c <= 3; c++) {
            if (!at_surface && m_geodyn_faces == 6) {
              //		     u[g](c,i,j,k0) = twgh*(
              //(1-wghi)*(1-wghj)*m_geodyn_data1[4](c,ig0,jg0,1)+
              //					  wghi*(1-wghj)*m_geodyn_data1[4](c,ig0+1,jg0,1)
              //					  (1-wghi)*wghj*m_geodyn_data1[4](c,ig0,jg0+1,1)+
              //					     wghi*wghj*m_geodyn_data1[4](c,ig0+1,jg0+1,1)
              //)
              //		     +
              //(1-twgh)*((1-wghi)*(1-wghj)*m_geodyn_data2[4](c,ig0,jg0,1)+
              //				 wghi*(1-wghj)*m_geodyn_data2[4](c,ig0+1,jg0,1)
              //				 (1-wghi)*wghj*m_geodyn_data2[4](c,ig0,jg0+1,1)+
              //				 wghi*wghj*m_geodyn_data2[4](c,ig0+1,jg0+1,1)
              //));
              u[g](c, i, j, k0) =
                  twgh * ((1 - wghi) * (1 - wghj) * gd14(c, ig0, jg0, 1) +
                          wghi * (1 - wghj) * gd14(c, ig0 + 1, jg0, 1) +
                          (1 - wghi) * wghj * gd14(c, ig0, jg0 + 1, 1) +
                          wghi * wghj * gd14(c, ig0 + 1, jg0 + 1, 1)) +
                  (1 - twgh) * ((1 - wghi) * (1 - wghj) * gd24(c, ig0, jg0, 1) +
                                wghi * (1 - wghj) * gd24(c, ig0 + 1, jg0, 1) +
                                (1 - wghi) * wghj * gd24(c, ig0, jg0 + 1, 1) +
                                wghi * wghj * gd24(c, ig0 + 1, jg0 + 1, 1));
            }
            u[g](c, i, j, k1) =
                twgh *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data1[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data1[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data1[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data1[5](c, ig0 + 1, jg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data2[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data2[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data2[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data2[5](c, ig0 + 1, jg0 + 1, 1));
          }
        }
      if (at_surface) {
        // Ghost point at corner, curvilinear case
        bcsurf_curvilinear_2nd_order(0, i0, i1, j0, j1, k0, g, u[g],
                                     bforcing[g][4]);
        bcsurf_curvilinear_2nd_order(1, i0, i1, j0, j1, k0, g, u[g],
                                     bforcing[g][4]);
        bcsurf_curvilinear_2nd_order(2, i0, i1, j0, j1, k0, g, u[g],
                                     bforcing[g][4]);
        bcsurf_curvilinear_2nd_order(3, i0, i1, j0, j1, k0, g, u[g],
                                     bforcing[g][4]);
      }
    }  // topographyExists
  }
}

//-----------------------------------------------------------------------
void EW::advance_geodyn_time(float_sw4 t) {
  // Before calling this routine geodyn_step should be such that
  // geodyn_data1 contains data at time step geodyn_step, and
  // geodyn_data2 contains data at time step geodyn_step+1.
  //
  // This routine advances geodyn_step such that
  //     geodyn_step*geodyn_dt <= t < (geodyn_step+1)*geodyn_dt
  // and updates geodyn data such that geodyn_data1 and geodyn_data2 are
  // as described above. It is assumed that t >= geodyn_step*geodyn_dt,
  // if t < geodyn_step*geodyn_dt, this routine does not do anything.
  //
  int n1 = static_cast<int>(floor(t / m_geodyn_dt));
  if (n1 > m_geodyn_maxsteps - 2)
    m_geodyn_past_end = true;
  else {
    // Find the right time step
    if (m_geodyn_step + 1 == n1) {
      // Advance one step
      copy_geodyn_timelevel(m_geodyn_data1, m_geodyn_data2);
      get_geodyn_timelevel(m_geodyn_data2);
    } else if (m_geodyn_step + 1 < n1) {
      // Advance many steps
      if (m_geodyn_iwillread) {
        char buf[256];
        for (int i = 0; i < m_geodyn_blocksize * (n1 - m_geodyn_step - 2); i++)
          m_geodynfile.getline(buf, 256);
      }
      get_geodyn_timelevel(m_geodyn_data1);
      get_geodyn_timelevel(m_geodyn_data2);
    }
    m_geodyn_step = n1;
  }
}

//-----------------------------------------------------------------------
void EW::get_geodyn_timelevel(vector<Sarray>& geodyndata) {
  if (m_geodyn_iwillread) {
    VERIFY2(!m_geodynfile.eof(),
            "Error: trying to read past end of Geodyn file");
    for (int side = 0; side <= 1; side++) {
      for (int k = 1; k <= m_geodyn_nk; k++)
        for (int j = 1; j <= m_geodyn_nj; j++)
          m_geodynfile >> geodyndata[side](1, j, k, 1) >>
              geodyndata[side](2, j, k, 1) >> geodyndata[side](3, j, k, 1);
    }
    for (int side = 2; side <= 3; side++) {
      for (int k = 1; k <= m_geodyn_nk; k++)
        for (int i = 1; i <= m_geodyn_ni; i++)
          m_geodynfile >> geodyndata[side](1, i, k, 1) >>
              geodyndata[side](2, i, k, 1) >> geodyndata[side](3, i, k, 1);
    }
    int llim = 4;
    if (m_geodyn_faces == 5) llim = 5;
    for (int side = llim; side <= 5; side++) {
      for (int j = 1; j <= m_geodyn_nj; j++)
        for (int i = 1; i <= m_geodyn_ni; i++)
          m_geodynfile >> geodyndata[side](1, i, j, 1) >>
              geodyndata[side](2, i, j, 1) >> geodyndata[side](3, i, j, 1);
    }
    // Read past remaining eol
    char buf[10];
    m_geodynfile.getline(buf, 10);
  }
}

//-----------------------------------------------------------------------
void EW::copy_geodyn_timelevel(vector<Sarray>& geodyndata1,
                               vector<Sarray>& geodyndata2) {
  for (int side = 0; side <= 1; side++) {
#pragma omp parallel for
    for (int k = 1; k <= m_geodyn_nk; k++)
      for (int j = 1; j <= m_geodyn_nj; j++)
        for (int c = 1; c <= 3; c++)
          geodyndata1[side](c, j, k, 1) = geodyndata2[side](c, j, k, 1);
  }
  for (int side = 2; side <= 3; side++) {
#pragma omp parallel for
    for (int k = 1; k <= m_geodyn_nk; k++)
      for (int i = 1; i <= m_geodyn_ni; i++)
        for (int c = 1; c <= 3; c++)
          geodyndata1[side](c, i, k, 1) = geodyndata2[side](c, i, k, 1);
  }
  int llim = 4;
  if (m_geodyn_faces == 5) llim = 5;
  for (int side = llim; side <= 5; side++) {
#pragma omp parallel for
    for (int j = 1; j <= m_geodyn_nj; j++)
      for (int i = 1; i <= m_geodyn_ni; i++)
        for (int c = 1; c <= 3; c++)
          geodyndata1[side](c, i, j, 1) = geodyndata2[side](c, i, j, 1);
  }
}

//-----------------------------------------------------------------------
void EW::geodyn_second_ghost_point(vector<Sarray>& rho, vector<Sarray>& mu,
                                   vector<Sarray>& lambda,
                                   vector<Sarray>& forcing, float_sw4 t,
                                   vector<Sarray>& U, vector<Sarray>& Um,
                                   int crf) {
  //   int n1 = static_cast<int>(floor(t/m_geodyn_dt));
  //   m_geodyn_step = n1;
  if (m_do_geodynbc) {
    float_sw4 twgh = ((m_geodyn_step + 1) * m_geodyn_dt - t) / m_geodyn_dt;
    float_sw4 d2i = 1 / (mDt * mDt);
    //      cout << "twgh = " << twgh << " geostep = " << m_geodyn_step << endl;
    for (int g = 0; g < mNumberOfCartesianGrids; g++) {
      int i0 = m_geodyn_dims[g][0];
      int i1 = m_geodyn_dims[g][1];
      int j0 = m_geodyn_dims[g][2];
      int j1 = m_geodyn_dims[g][3];
      int k0 = m_geodyn_dims[g][4];
      int k1 = m_geodyn_dims[g][5];
      float_sw4 h = mGridSize[g];
      float_sw4 h2 = h * h;

      bool low_interior, high_interior;
      low_interior = m_iStartInt[g] <= i0 + 1 && i0 + 1 <= m_iEndInt[g];
      high_interior = m_iStartInt[g] <= i1 - 1 && i1 - 1 <= m_iEndInt[g];
      bool surface_correction = k0 <= 1 && g == mNumberOfGrids - 1;

      Sarray Lu0(3, i0, i0, j0, j1, k0, k1);
      if (low_interior)
        evalLu_Dim(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu0.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu0, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i0, j0,
                   j1, k0, k1);
      Sarray Lu1(3, i1, i1, j0, j1, k0, k1);
      if (high_interior)
        evalLu_Dip(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu1.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu1, mu[g].c_ptr(), lambda[g].c_ptr(), h, i1, i1, j0,
                   j1, k0, k1);
      int kstart = k0 + 1;
      if (surface_correction) {
        // Special at corner between free surface and Geodyn cube
        if (low_interior)
          evalLu_DkpDim(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                        m_kStart[g], m_kEnd[g],
                        //			U[g].c_ptr(), Lu0.c_ptr(),
                        //mu[g].c_ptr(), lambda[g].c_ptr(),
                        U[g], Lu0, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i0,
                        j0, j1, k0, k1);
        if (high_interior)
          evalLu_DkpDip(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                        m_kStart[g], m_kEnd[g],
                        //			U[g].c_ptr(), Lu1.c_ptr(),
                        //mu[g].c_ptr(), lambda[g].c_ptr(),
                        U[g], Lu1, mu[g].c_ptr(), lambda[g].c_ptr(), h, i1, i1,
                        j0, j1, k0, k1);
        kstart = k0;
        // k0 should be 1 here.
      }

// Side with i=const.
#pragma omp parallel for
      for (int k = kstart; k <= k1 - 1; k++)
        for (int j = j0 + 1; j <= j1 - 1; j++) {
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(
              floor(((k - 1) * h - m_geodyn_origin[2]) / m_geodyn_h + 1));
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghk =
              ((k - 1) * h - (m_geodyn_origin[2] + (kg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 bnd0[3], bnd1[3];
          for (int c = 1; c <= 3; c++) {
            bnd0[c - 1] =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[0](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[0](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[0](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[0](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[0](c, jg0 + 1, kg0 + 1, 1));
            bnd1[c - 1] =
                twgh *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data1[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data1[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data1[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data1[1](c, jg0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghj) * (1 - wghk) *
                         m_geodyn_data2[1](c, jg0, kg0, 1) +
                     wghj * (1 - wghk) * m_geodyn_data2[1](c, jg0 + 1, kg0, 1) +
                     (1 - wghj) * wghk * m_geodyn_data2[1](c, jg0, kg0 + 1, 1) +
                     wghj * wghk * m_geodyn_data2[1](c, jg0 + 1, kg0 + 1, 1));
          }

          // Lower bndry
          float_sw4 res1, res2, res3;
          if (low_interior) {
            res1 = crf * rho[g](i0, j, k) *
                       (bnd0[0] - 2 * U[g](1, i0, j, k) + Um[g](1, i0, j, k)) *
                       d2i -
                   Lu0(1, i0, j, k) - forcing[g](1, i0, j, k);
            res2 = crf * rho[g](i0, j, k) *
                       (bnd0[1] - 2 * U[g](2, i0, j, k) + Um[g](2, i0, j, k)) *
                       d2i -
                   Lu0(2, i0, j, k) - forcing[g](2, i0, j, k);
            res3 = crf * rho[g](i0, j, k) *
                       (bnd0[2] - 2 * U[g](3, i0, j, k) + Um[g](3, i0, j, k)) *
                       d2i -
                   Lu0(3, i0, j, k) - forcing[g](3, i0, j, k);

            U[g](1, i0 + 1, j, k) =
                U[g](1, i0 + 1, j, k) +
                h2 * res1 /
                    (mu[g](i0 + 1, j, k) + mu[g](i0, j, k) +
                     0.5 * (lambda[g](i0 + 1, j, k) + lambda[g](i0, j, k)));
            U[g](2, i0 + 1, j, k) =
                U[g](2, i0 + 1, j, k) +
                2 * h2 * res2 / (mu[g](i0 + 1, j, k) + mu[g](i0, j, k));
            U[g](3, i0 + 1, j, k) =
                U[g](3, i0 + 1, j, k) +
                2 * h2 * res3 / (mu[g](i0 + 1, j, k) + mu[g](i0, j, k));
          }
          // Upper bndry
          if (high_interior) {
            res1 = crf * rho[g](i1, j, k) *
                       (bnd1[0] - 2 * U[g](1, i1, j, k) + Um[g](1, i1, j, k)) *
                       d2i -
                   Lu1(1, i1, j, k) - forcing[g](1, i1, j, k);
            res2 = crf * rho[g](i1, j, k) *
                       (bnd1[1] - 2 * U[g](2, i1, j, k) + Um[g](2, i1, j, k)) *
                       d2i -
                   Lu1(2, i1, j, k) - forcing[g](2, i1, j, k);
            res3 = crf * rho[g](i1, j, k) *
                       (bnd1[2] - 2 * U[g](3, i1, j, k) + Um[g](3, i1, j, k)) *
                       d2i -
                   Lu1(3, i1, j, k) - forcing[g](3, i1, j, k);

            U[g](1, i1 - 1, j, k) =
                U[g](1, i1 - 1, j, k) +
                h2 * res1 /
                    (mu[g](i1 - 1, j, k) + mu[g](i1, j, k) +
                     0.5 * (lambda[g](i1 - 1, j, k) + lambda[g](i1, j, k)));
            U[g](2, i1 - 1, j, k) =
                U[g](2, i1 - 1, j, k) +
                2 * h2 * res2 / (mu[g](i1 - 1, j, k) + mu[g](i1, j, k));
            U[g](3, i1 - 1, j, k) =
                U[g](3, i1 - 1, j, k) +
                2 * h2 * res3 / (mu[g](i1 - 1, j, k) + mu[g](i1, j, k));
          }
        }

      // Side with j=const
      low_interior = m_jStartInt[g] <= j0 + 1 && j0 + 1 <= m_jEndInt[g];
      high_interior = m_jStartInt[g] <= j1 - 1 && j1 - 1 <= m_jEndInt[g];

      if (low_interior) {
        Lu0.define(3, i0, i1, j0, j0, k0, k1);
        evalLu_Djm(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu0.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu0, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1, j0,
                   j0, k0, k1);
      }
      if (high_interior) {
        Lu1.define(3, i0, i1, j1, j1, k0, k1);
        evalLu_Djp(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu1.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu1, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1, j1,
                   j1, k0, k1);
      }
      if (surface_correction) {
        if (low_interior)
          evalLu_DkpDjm(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                        m_kStart[g], m_kEnd[g],
                        //			U[g].c_ptr(), Lu0.c_ptr(),
                        //mu[g].c_ptr(), lambda[g].c_ptr(),
                        U[g], Lu0, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1,
                        j0, j0, k0, k1);
        if (high_interior)
          evalLu_DkpDjp(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                        m_kStart[g], m_kEnd[g],
                        //			U[g].c_ptr(), Lu1.c_ptr(),
                        //mu[g].c_ptr(), lambda[g].c_ptr(),
                        U[g], Lu1, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1,
                        j1, j1, k0, k1);
      }
#pragma omp parallel for
      for (int k = kstart; k <= k1 - 1; k++)
        for (int i = i0 + 1; i <= i1 - 1; i++) {
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int kg0 = static_cast<int>(
              floor(((k - 1) * h - m_geodyn_origin[2]) / m_geodyn_h + 1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
          if (kg0 <= 0) kg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghk =
              ((k - 1) * h - (m_geodyn_origin[2] + (kg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 bnd0[3], bnd1[3];
          for (int c = 1; c <= 3; c++) {
            bnd0[c - 1] =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[2](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[2](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[2](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[2](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[2](c, ig0 + 1, kg0 + 1, 1));
            bnd1[c - 1] =
                twgh *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data1[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data1[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data1[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data1[3](c, ig0 + 1, kg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghk) *
                         m_geodyn_data2[3](c, ig0, kg0, 1) +
                     wghi * (1 - wghk) * m_geodyn_data2[3](c, ig0 + 1, kg0, 1) +
                     (1 - wghi) * wghk * m_geodyn_data2[3](c, ig0, kg0 + 1, 1) +
                     wghi * wghk * m_geodyn_data2[3](c, ig0 + 1, kg0 + 1, 1));
          }
          float_sw4 res1, res2, res3;
          // Lower bndry
          if (low_interior) {
            res1 = crf * rho[g](i, j0, k) *
                       (bnd0[0] - 2 * U[g](1, i, j0, k) + Um[g](1, i, j0, k)) *
                       d2i -
                   Lu0(1, i, j0, k) - forcing[g](1, i, j0, k);
            res2 = crf * rho[g](i, j0, k) *
                       (bnd0[1] - 2 * U[g](2, i, j0, k) + Um[g](2, i, j0, k)) *
                       d2i -
                   Lu0(2, i, j0, k) - forcing[g](2, i, j0, k);
            res3 = crf * rho[g](i, j0, k) *
                       (bnd0[2] - 2 * U[g](3, i, j0, k) + Um[g](3, i, j0, k)) *
                       d2i -
                   Lu0(3, i, j0, k) - forcing[g](3, i, j0, k);

            U[g](1, i, j0 + 1, k) =
                U[g](1, i, j0 + 1, k) +
                2 * h2 * res1 / (mu[g](i, j0 + 1, k) + mu[g](i, j0, k));
            U[g](2, i, j0 + 1, k) =
                U[g](2, i, j0 + 1, k) +
                h2 * res2 /
                    (mu[g](i, j0 + 1, k) + mu[g](i, j0, k) +
                     0.5 * (lambda[g](i, j0 + 1, k) + lambda[g](i, j0, k)));
            U[g](3, i, j0 + 1, k) =
                U[g](3, i, j0 + 1, k) +
                2 * h2 * res3 / (mu[g](i, j0 + 1, k) + mu[g](i, j0, k));
          }
          // Upper bndry
          if (high_interior) {
            res1 = crf * rho[g](i, j1, k) *
                       (bnd1[0] - 2 * U[g](1, i, j1, k) + Um[g](1, i, j1, k)) *
                       d2i -
                   Lu1(1, i, j1, k) - forcing[g](1, i, j1, k);
            res2 = crf * rho[g](i, j1, k) *
                       (bnd1[1] - 2 * U[g](2, i, j1, k) + Um[g](2, i, j1, k)) *
                       d2i -
                   Lu1(2, i, j1, k) - forcing[g](2, i, j1, k);
            res3 = crf * rho[g](i, j1, k) *
                       (bnd1[2] - 2 * U[g](3, i, j1, k) + Um[g](3, i, j1, k)) *
                       d2i -
                   Lu1(3, i, j1, k) - forcing[g](3, i, j1, k);

            U[g](1, i, j1 - 1, k) =
                U[g](1, i, j1 - 1, k) +
                2 * h2 * res1 / (mu[g](i, j1 - 1, k) + mu[g](i, j1, k));
            U[g](2, i, j1 - 1, k) =
                U[g](2, i, j1 - 1, k) +
                h2 * res2 /
                    (mu[g](i, j1 - 1, k) + mu[g](i, j1, k) +
                     0.5 * (lambda[g](i, j1 - 1, k) + lambda[g](i, j1, k)));
            U[g](3, i, j1 - 1, k) =
                U[g](3, i, j1 - 1, k) +
                2 * h2 * res3 / (mu[g](i, j1 - 1, k) + mu[g](i, j1, k));
          }
        }

      // Side with k=const
      low_interior = (m_kStartInt[g] <= k0 + 1 && k0 + 1 <= m_kEndInt[g]) &&
                     !surface_correction;
      high_interior = m_kStartInt[g] <= k1 - 1 && k1 - 1 <= m_kEndInt[g];
      if (m_geodyn_faces == 6 && low_interior) {
        Lu0.define(3, i0, i1, j0, j1, k0, k0);
        evalLu_Dkm(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu0.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu0, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1, j0,
                   j1, k0, k0);
      }

      if (high_interior) {
        Lu1.define(3, i0, i1, j0, j1, k1, k1);
        evalLu_Dkp(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g], m_kStart[g],
                   m_kEnd[g],
                   //			U[g].c_ptr(), Lu1.c_ptr(), mu[g].c_ptr(),
                   //lambda[g].c_ptr(),
                   U[g], Lu1, mu[g].c_ptr(), lambda[g].c_ptr(), h, i0, i1, j0,
                   j1, k1, k1);
      }
#pragma omp parallel for
      for (int j = j0 + 1; j <= j1 - 1; j++)
        for (int i = i0 + 1; i <= i1 - 1; i++) {
          int ig0 = static_cast<int>(
              floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
          int jg0 = static_cast<int>(
              floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
          if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
          if (ig0 <= 0) ig0 = 1;
          if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
          if (jg0 <= 0) jg0 = 1;
          float_sw4 wghi =
              ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 wghj =
              ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
              m_geodyn_h;
          float_sw4 bnd0[3], bnd1[3];
          for (int c = 1; c <= 3; c++) {
            if (m_geodyn_faces == 6 && low_interior) {
              bnd0[c - 1] =
                  twgh * ((1 - wghi) * (1 - wghj) *
                              m_geodyn_data1[4](c, ig0, jg0, 1) +
                          wghi * (1 - wghj) *
                              m_geodyn_data1[4](c, ig0 + 1, jg0, 1) +
                          (1 - wghi) * wghj *
                              m_geodyn_data1[4](c, ig0, jg0 + 1, 1) +
                          wghi * wghj *
                              m_geodyn_data1[4](c, ig0 + 1, jg0 + 1, 1)) +
                  (1 - twgh) *
                      ((1 - wghi) * (1 - wghj) *
                           m_geodyn_data2[4](c, ig0, jg0, 1) +
                       wghi * (1 - wghj) *
                           m_geodyn_data2[4](c, ig0 + 1, jg0, 1) +
                       (1 - wghi) * wghj *
                           m_geodyn_data2[4](c, ig0, jg0 + 1, 1) +
                       wghi * wghj * m_geodyn_data2[4](c, ig0 + 1, jg0 + 1, 1));
            }
            bnd1[c - 1] =
                twgh *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data1[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data1[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data1[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data1[5](c, ig0 + 1, jg0 + 1, 1)) +
                (1 - twgh) *
                    ((1 - wghi) * (1 - wghj) *
                         m_geodyn_data2[5](c, ig0, jg0, 1) +
                     wghi * (1 - wghj) * m_geodyn_data2[5](c, ig0 + 1, jg0, 1) +
                     (1 - wghi) * wghj * m_geodyn_data2[5](c, ig0, jg0 + 1, 1) +
                     wghi * wghj * m_geodyn_data2[5](c, ig0 + 1, jg0 + 1, 1));
          }
          // Upper bndry
          float_sw4 res1, res2, res3;
          if (high_interior) {
            res1 = crf * rho[g](i, j, k1) *
                       (bnd1[0] - 2 * U[g](1, i, j, k1) + Um[g](1, i, j, k1)) *
                       d2i -
                   Lu1(1, i, j, k1) - forcing[g](1, i, j, k1);
            res2 = crf * rho[g](i, j, k1) *
                       (bnd1[1] - 2 * U[g](2, i, j, k1) + Um[g](2, i, j, k1)) *
                       d2i -
                   Lu1(2, i, j, k1) - forcing[g](2, i, j, k1);
            res3 = crf * rho[g](i, j, k1) *
                       (bnd1[2] - 2 * U[g](3, i, j, k1) + Um[g](3, i, j, k1)) *
                       d2i -
                   Lu1(3, i, j, k1) - forcing[g](3, i, j, k1);

            U[g](1, i, j, k1 - 1) =
                U[g](1, i, j, k1 - 1) +
                2 * h2 * res1 / (mu[g](i, j, k1 - 1) + mu[g](i, j, k1));
            U[g](2, i, j, k1 - 1) =
                U[g](2, i, j, k1 - 1) +
                2 * h2 * res2 / (mu[g](i, j, k1 - 1) + mu[g](i, j, k1));
            U[g](3, i, j, k1 - 1) =
                U[g](3, i, j, k1 - 1) +
                h2 * res3 /
                    (mu[g](i, j, k1 - 1) + mu[g](i, j, k1) +
                     0.5 * (lambda[g](i, j, k1 - 1) + lambda[g](i, j, k1)));
          }
          // Lower bndry
          if (m_geodyn_faces == 6 && low_interior) {
            res1 = crf * rho[g](i, j, k0) *
                       (bnd0[0] - 2 * U[g](1, i, j, k0) + Um[g](1, i, j, k0)) *
                       d2i -
                   Lu0(1, i, j, k0) - forcing[g](1, i, j, k0);
            res2 = crf * rho[g](i, j, k0) *
                       (bnd0[1] - 2 * U[g](2, i, j, k0) + Um[g](2, i, j, k0)) *
                       d2i -
                   Lu0(2, i, j, k0) - forcing[g](2, i, j, k0);
            res3 = crf * rho[g](i, j, k0) *
                       (bnd0[2] - 2 * U[g](3, i, j, k0) + Um[g](3, i, j, k0)) *
                       d2i -
                   Lu0(3, i, j, k0) - forcing[g](3, i, j, k0);

            U[g](1, i, j, k0 + 1) =
                U[g](1, i, j, k0 + 1) +
                2 * h2 * res1 / (mu[g](i, j, k0 + 1) + mu[g](i, j, k0));
            U[g](2, i, j, k0 + 1) =
                U[g](2, i, j, k0 + 1) +
                2 * h2 * res2 / (mu[g](i, j, k0 + 1) + mu[g](i, j, k0));
            U[g](3, i, j, k0 + 1) =
                U[g](3, i, j, k0 + 1) +
                h2 * res3 /
                    (mu[g](i, j, k0 + 1) + mu[g](i, j, k0) +
                     0.5 * (lambda[g](i, j, k0 + 1) + lambda[g](i, j, k0)));
          }
        }
    }

    // Topographic case, grid is curvilinear
    if (topographyExists()) {
      geodyn_second_ghost_point_curvilinear(rho, mu, lambda, forcing, t, U, Um,
                                            crf);
    }
  }
}

//-----------------------------------------------------------------------
extern "C" {
void F77_FUNC(dgesv, DGESV)(int* n, int* nrhs, double* A, int* lda, int* ipiv,
                            double* b, int* ldb, int* info);
}

//-----------------------------------------------------------------------
void EW::geodyn_second_ghost_point_curvilinear(vector<Sarray>& rho,
                                               vector<Sarray>& mu,
                                               vector<Sarray>& lambda,
                                               vector<Sarray>& forcing,
                                               float_sw4 t, vector<Sarray>& U,
                                               vector<Sarray>& Um, int crf) {
  float_sw4 twgh = ((m_geodyn_step + 1) * m_geodyn_dt - t) / m_geodyn_dt;
  float_sw4 d2i = 1 / (mDt * mDt);
  int g = mNumberOfGrids - 1;
  float_sw4 h = mGridSize[g];

  int i0 = m_geodyn_dims[g][0];
  int i1 = m_geodyn_dims[g][1];
  int j0 = m_geodyn_dims[g][2];
  int j1 = m_geodyn_dims[g][3];
  int k0 = m_geodyn_dims[g][4];
  int k1 = m_geodyn_dims[g][5];

  float_sw4 zcubelen = (m_geodyn_nk - 1) * m_geodyn_h;
  //   bool   at_surface = k0==1;

  bool low_interior, high_interior;
  low_interior = m_iStartInt[g] <= i0 + 1 && i0 + 1 <= m_iEndInt[g];
  high_interior = m_iStartInt[g] <= i1 - 1 && i1 - 1 <= m_iEndInt[g];
  bool surface_correction = k0 <= 1 && g == mNumberOfGrids - 1;

  Sarray Lu0(3, i0, i0, j0, j1, k0, k1);
  if (low_interior)
    evalLuCurv<0, 1, 1, 1, 1, 1>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu0,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i0, i0, j0, j1, k0, k1);
  Sarray Lu1(3, i1, i1, j0, j1, k0, k1);
  if (high_interior)
    evalLuCurv<1, 0, 1, 1, 1, 1>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu1,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i1, i1, j0, j1, k0, k1);

  int kstart = k0 + 1;
  if (surface_correction) {
    // Special at corner between free surface and Geodyn cube
    if (low_interior)
      evalLuCurv<0, 1, 1, 1, 1, 0>(m_iStart[g], m_iEnd[g], m_jStart[g],
                                   m_jEnd[g], m_kStart[g], m_kEnd[g], U[g], Lu0,
                                   mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                   mJ[g], i0, i0, j0, j1, k0, k0);
    if (high_interior)
      evalLuCurv<1, 0, 1, 1, 1, 0>(m_iStart[g], m_iEnd[g], m_jStart[g],
                                   m_jEnd[g], m_kStart[g], m_kEnd[g], U[g], Lu1,
                                   mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                   mJ[g], i1, i1, j0, j1, k0, k0);
    kstart = k0;
  }
// Side with i=const.
#pragma omp parallel for
  for (int k = kstart; k <= k1 - 1; k++)
    for (int j = j0 + 1; j <= j1 - 1; j++) {
      float_sw4 strfact = (mZ[g](i0, j, k1) - mZ[g](i0, j, 1)) / zcubelen;
      int jg0 = static_cast<int>(
          floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
      int kg0 = static_cast<int>(floor(
          (mZ[g](i0, j, k) - mZ[g](i0, j, 1)) / (strfact * m_geodyn_h) + 1));
      //               int kg0 = static_cast<int>(floor(((k-1)*h -
      //               m_geodyn_origin[2])/m_geodyn_h+1));
      if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
      if (jg0 <= 0) jg0 = 1;
      if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
      if (kg0 <= 0) kg0 = 1;
      float_sw4 wghj =
          ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
          m_geodyn_h;
      float_sw4 wghk = ((mZ[g](i0, j, k) - mZ[g](i0, j, 1)) / strfact -
                        ((kg0 - 1) * m_geodyn_h)) /
                       m_geodyn_h;
      float_sw4 bnd0[3], bnd1[3];
      //	 if( i0+1==86 && j==102 && k==25)
      //	    cout << "boundary data: wk " << wghk << " sf " << strfact <<
      //" zk1 "
      //		 << mZ[g](i0,j,k1) << " z1 " << mZ[g](i0,j,1) << " k1 "
      //<< k1 <<  endl;
      for (int c = 1; c <= 3; c++) {
        bnd0[c - 1] =
            twgh *
                ((1 - wghj) * (1 - wghk) * m_geodyn_data1[0](c, jg0, kg0, 1) +
                 wghj * (1 - wghk) * m_geodyn_data1[0](c, jg0 + 1, kg0, 1) +
                 (1 - wghj) * wghk * m_geodyn_data1[0](c, jg0, kg0 + 1, 1) +
                 wghj * wghk * m_geodyn_data1[0](c, jg0 + 1, kg0 + 1, 1)) +
            (1 - twgh) *
                ((1 - wghj) * (1 - wghk) * m_geodyn_data2[0](c, jg0, kg0, 1) +
                 wghj * (1 - wghk) * m_geodyn_data2[0](c, jg0 + 1, kg0, 1) +
                 (1 - wghj) * wghk * m_geodyn_data2[0](c, jg0, kg0 + 1, 1) +
                 wghj * wghk * m_geodyn_data2[0](c, jg0 + 1, kg0 + 1, 1));
      }
      strfact = (mZ[g](i1, j, k1) - mZ[g](i1, j, 1)) / zcubelen;
      kg0 = static_cast<int>(floor(
          (mZ[g](i1, j, k) - mZ[g](i1, j, 1)) / (strfact * m_geodyn_h) + 1));
      if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
      if (kg0 <= 0) kg0 = 1;
      wghk = ((mZ[g](i1, j, k) - mZ[g](i1, j, 1)) / strfact -
              ((kg0 - 1) * m_geodyn_h)) /
             m_geodyn_h;
      for (int c = 1; c <= 3; c++) {
        bnd1[c - 1] =
            twgh *
                ((1 - wghj) * (1 - wghk) * m_geodyn_data1[1](c, jg0, kg0, 1) +
                 wghj * (1 - wghk) * m_geodyn_data1[1](c, jg0 + 1, kg0, 1) +
                 (1 - wghj) * wghk * m_geodyn_data1[1](c, jg0, kg0 + 1, 1) +
                 wghj * wghk * m_geodyn_data1[1](c, jg0 + 1, kg0 + 1, 1)) +
            (1 - twgh) *
                ((1 - wghj) * (1 - wghk) * m_geodyn_data2[1](c, jg0, kg0, 1) +
                 wghj * (1 - wghk) * m_geodyn_data2[1](c, jg0 + 1, kg0, 1) +
                 (1 - wghj) * wghk * m_geodyn_data2[1](c, jg0, kg0 + 1, 1) +
                 wghj * wghk * m_geodyn_data2[1](c, jg0 + 1, kg0 + 1, 1));
      }
      // Lower bndry
      float_sw4 res1, res2, res3;
      if (low_interior) {
        res1 = crf * rho[g](i0, j, k) *
                   (bnd0[0] - 2 * U[g](1, i0, j, k) + Um[g](1, i0, j, k)) *
                   d2i -
               Lu0(1, i0, j, k) - forcing[g](1, i0, j, k);
        res2 = crf * rho[g](i0, j, k) *
                   (bnd0[1] - 2 * U[g](2, i0, j, k) + Um[g](2, i0, j, k)) *
                   d2i -
               Lu0(2, i0, j, k) - forcing[g](2, i0, j, k);
        res3 = crf * rho[g](i0, j, k) *
                   (bnd0[2] - 2 * U[g](3, i0, j, k) + Um[g](3, i0, j, k)) *
                   d2i -
               Lu0(3, i0, j, k) - forcing[g](3, i0, j, k);

        //	    U[g](1,i0+1,j,k) = U[g](1,i0+1,j,k) +
        //h2*res1/(mu[g](i0+1,j,k)+mu[g](i0,j,k)+
        //								 0.5*(lambda[g](i0+1,j,k)+lambda[g](i0,j,k)));
        //	    U[g](2,i0+1,j,k) = U[g](2,i0+1,j,k) +
        //2*h2*res2/(mu[g](i0+1,j,k)+mu[g](i0,j,k)); 	    U[g](3,i0+1,j,k) =
        //U[g](3,i0+1,j,k) + 2*h2*res3/(mu[g](i0+1,j,k)+mu[g](i0,j,k));
        U[g](1, i0 + 1, j, k) =
            U[g](1, i0 + 1, j, k) +
            mJ[g](i0, j, k) * res1 /
                ((mu[g](i0 + 1, j, k) + 0.5 * lambda[g](i0 + 1, j, k)) *
                     SQR(mMetric[g](1, i0 + 1, j, k)) +
                 (mu[g](i0, j, k) + 0.5 * lambda[g](i0, j, k)) *
                     SQR(mMetric[g](1, i0, j, k)));
        U[g](2, i0 + 1, j, k) =
            U[g](2, i0 + 1, j, k) +
            2 * mJ[g](i0, j, k) * res2 /
                (mu[g](i0 + 1, j, k) * SQR(mMetric[g](1, i0 + 1, j, k)) +
                 mu[g](i0, j, k) * SQR(mMetric[g](1, i0, j, k)));
        U[g](3, i0 + 1, j, k) =
            U[g](3, i0 + 1, j, k) +
            2 * mJ[g](i0, j, k) * res3 /
                (mu[g](i0 + 1, j, k) * SQR(mMetric[g](1, i0 + 1, j, k)) +
                 mu[g](i0, j, k) * SQR(mMetric[g](1, i0, j, k)));
        //	    if( i0+1==86 && j==102 && k==25)
        //		     cout << "In geodyn bc "  << m_myRank << " " <<
        //U[g](1,i0+1,j,k) << 			" " << Um[g](1,i0,j,k) << " " << bnd0[0] << " " <<
        //Lu0(1,i0,j,k) << endl;
      }
      // Upper bndry
      if (high_interior) {
        res1 = crf * rho[g](i1, j, k) *
                   (bnd1[0] - 2 * U[g](1, i1, j, k) + Um[g](1, i1, j, k)) *
                   d2i -
               Lu1(1, i1, j, k) - forcing[g](1, i1, j, k);
        res2 = crf * rho[g](i1, j, k) *
                   (bnd1[1] - 2 * U[g](2, i1, j, k) + Um[g](2, i1, j, k)) *
                   d2i -
               Lu1(2, i1, j, k) - forcing[g](2, i1, j, k);
        res3 = crf * rho[g](i1, j, k) *
                   (bnd1[2] - 2 * U[g](3, i1, j, k) + Um[g](3, i1, j, k)) *
                   d2i -
               Lu1(3, i1, j, k) - forcing[g](3, i1, j, k);

        //	    U[g](1,i1-1,j,k) = U[g](1,i1-1,j,k) +
        //h2*res1/(mu[g](i1-1,j,k)+mu[g](i1,j,k)+
        //								 0.5*(lambda[g](i1-1,j,k)+lambda[g](i1,j,k)));
        //	    U[g](2,i1-1,j,k) = U[g](2,i1-1,j,k) +
        //2*h2*res2/(mu[g](i1-1,j,k)+mu[g](i1,j,k)); 	    U[g](3,i1-1,j,k) =
        //U[g](3,i1-1,j,k) + 2*h2*res3/(mu[g](i1-1,j,k)+mu[g](i1,j,k));
        U[g](1, i1 - 1, j, k) =
            U[g](1, i1 - 1, j, k) +
            mJ[g](i1, j, k) * res1 /
                ((mu[g](i1 - 1, j, k) + 0.5 * lambda[g](i1 - 1, j, k)) *
                     SQR(mMetric[g](1, i1 - 1, j, k)) +
                 (mu[g](i1, j, k) + 0.5 * lambda[g](i1, j, k)) *
                     SQR(mMetric[g](1, i1, j, k)));
        U[g](2, i1 - 1, j, k) =
            U[g](2, i1 - 1, j, k) +
            2 * mJ[g](i1, j, k) * res2 /
                (mu[g](i1 - 1, j, k) * SQR(mMetric[g](1, i1 - 1, j, k)) +
                 mu[g](i1, j, k) * SQR(mMetric[g](1, i1, j, k)));
        U[g](3, i1 - 1, j, k) =
            U[g](3, i1 - 1, j, k) +
            2 * mJ[g](i1, j, k) * res3 /
                (mu[g](i1 - 1, j, k) * SQR(mMetric[g](1, i1 - 1, j, k)) +
                 mu[g](i1, j, k) * SQR(mMetric[g](1, i1, j, k)));
      }
    }
  // Side with j=const
  low_interior = m_jStartInt[g] <= j0 + 1 && j0 + 1 <= m_jEndInt[g];
  high_interior = m_jStartInt[g] <= j1 - 1 && j1 - 1 <= m_jEndInt[g];
  if (low_interior) {
    Lu0.define(3, i0, i1, j0, j0, k0, k1);
    evalLuCurv<1, 1, 0, 1, 1, 1>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu0,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i0, i1, j0, j0, k0, k1);
  }
  if (high_interior) {
    Lu1.define(3, i0, i1, j1, j1, k0, k1);
    evalLuCurv<1, 1, 1, 0, 1, 1>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu1,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i0, i1, j1, j1, k0, k1);
  }
  if (surface_correction) {
    if (low_interior)
      evalLuCurv<1, 1, 0, 1, 1, 0>(m_iStart[g], m_iEnd[g], m_jStart[g],
                                   m_jEnd[g], m_kStart[g], m_kEnd[g], U[g], Lu0,
                                   mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                   mJ[g], i0, i1, j0, j0, k0, k0);
    if (high_interior)
      evalLuCurv<1, 1, 1, 0, 1, 0>(m_iStart[g], m_iEnd[g], m_jStart[g],
                                   m_jEnd[g], m_kStart[g], m_kEnd[g], U[g], Lu1,
                                   mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                   mJ[g], i0, i1, j1, j1, k0, k0);
  }
#pragma omp parallel for
  for (int k = kstart; k <= k1 - 1; k++)
    for (int i = i0 + 1; i <= i1 - 1; i++) {
      float_sw4 strfact = (mZ[g](i, j0, k1) - mZ[g](i, j0, 1)) / zcubelen;
      int ig0 = static_cast<int>(
          floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
      int kg0 = static_cast<int>(floor(
          (mZ[g](i, j0, k) - mZ[g](i, j0, 1)) / (strfact * m_geodyn_h) + 1));
      //	       int kg0 = static_cast<int>(floor(((k-1)*h -
      //m_geodyn_origin[2])/m_geodyn_h+1));
      if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
      if (ig0 <= 0) ig0 = 1;
      if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
      if (kg0 <= 0) kg0 = 1;
      float_sw4 wghi =
          ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
          m_geodyn_h;
      //	       double wghk = ((k-1)*h -
      //(m_geodyn_origin[2]+(kg0-1)*m_geodyn_h))/m_geodyn_h;
      float_sw4 wghk = ((mZ[g](i, j0, k) - mZ[g](i, j0, 1)) / strfact -
                        ((kg0 - 1) * m_geodyn_h)) /
                       m_geodyn_h;
      float_sw4 bnd0[3], bnd1[3];
      for (int c = 1; c <= 3; c++) {
        bnd0[c - 1] =
            twgh *
                ((1 - wghi) * (1 - wghk) * m_geodyn_data1[2](c, ig0, kg0, 1) +
                 wghi * (1 - wghk) * m_geodyn_data1[2](c, ig0 + 1, kg0, 1) +
                 (1 - wghi) * wghk * m_geodyn_data1[2](c, ig0, kg0 + 1, 1) +
                 wghi * wghk * m_geodyn_data1[2](c, ig0 + 1, kg0 + 1, 1)) +
            (1 - twgh) *
                ((1 - wghi) * (1 - wghk) * m_geodyn_data2[2](c, ig0, kg0, 1) +
                 wghi * (1 - wghk) * m_geodyn_data2[2](c, ig0 + 1, kg0, 1) +
                 (1 - wghi) * wghk * m_geodyn_data2[2](c, ig0, kg0 + 1, 1) +
                 wghi * wghk * m_geodyn_data2[2](c, ig0 + 1, kg0 + 1, 1));
      }
      strfact = (mZ[g](i, j1, k1) - mZ[g](i, j1, 1)) / zcubelen;
      kg0 = static_cast<int>(floor(
          (mZ[g](i, j1, k) - mZ[g](i, j1, 1)) / (strfact * m_geodyn_h) + 1));
      if (kg0 >= m_geodyn_nk) kg0 = m_geodyn_nk - 1;
      if (kg0 <= 0) kg0 = 1;
      wghk = ((mZ[g](i, j1, k) - mZ[g](i, j1, 1)) / strfact -
              ((kg0 - 1) * m_geodyn_h)) /
             m_geodyn_h;
      for (int c = 1; c <= 3; c++) {
        bnd1[c - 1] =
            twgh *
                ((1 - wghi) * (1 - wghk) * m_geodyn_data1[3](c, ig0, kg0, 1) +
                 wghi * (1 - wghk) * m_geodyn_data1[3](c, ig0 + 1, kg0, 1) +
                 (1 - wghi) * wghk * m_geodyn_data1[3](c, ig0, kg0 + 1, 1) +
                 wghi * wghk * m_geodyn_data1[3](c, ig0 + 1, kg0 + 1, 1)) +
            (1 - twgh) *
                ((1 - wghi) * (1 - wghk) * m_geodyn_data2[3](c, ig0, kg0, 1) +
                 wghi * (1 - wghk) * m_geodyn_data2[3](c, ig0 + 1, kg0, 1) +
                 (1 - wghi) * wghk * m_geodyn_data2[3](c, ig0, kg0 + 1, 1) +
                 wghi * wghk * m_geodyn_data2[3](c, ig0 + 1, kg0 + 1, 1));
      }
      float_sw4 res1, res2, res3;
      // Lower bndry
      if (low_interior) {
        res1 = crf * rho[g](i, j0, k) *
                   (bnd0[0] - 2 * U[g](1, i, j0, k) + Um[g](1, i, j0, k)) *
                   d2i -
               Lu0(1, i, j0, k) - forcing[g](1, i, j0, k);
        res2 = crf * rho[g](i, j0, k) *
                   (bnd0[1] - 2 * U[g](2, i, j0, k) + Um[g](2, i, j0, k)) *
                   d2i -
               Lu0(2, i, j0, k) - forcing[g](2, i, j0, k);
        res3 = crf * rho[g](i, j0, k) *
                   (bnd0[2] - 2 * U[g](3, i, j0, k) + Um[g](3, i, j0, k)) *
                   d2i -
               Lu0(3, i, j0, k) - forcing[g](3, i, j0, k);

        //	    U[g](1,i,j0+1,k) = U[g](1,i,j0+1,k) +
        //2*h2*res1/(mu[g](i,j0+1,k)+mu[g](i,j0,k)); 	    U[g](2,i,j0+1,k) =
        //U[g](2,i,j0+1,k) +   h2*res2/(mu[g](i,j0+1,k)+mu[g](i,j0,k)+
        //							     0.5*(lambda[g](i,j0+1,k)+lambda[g](i,j0,k)));
        //	    U[g](3,i,j0+1,k) = U[g](3,i,j0+1,k) +
        //2*h2*res3/(mu[g](i,j0+1,k)+mu[g](i,j0,k));
        U[g](1, i, j0 + 1, k) =
            U[g](1, i, j0 + 1, k) +
            2 * mJ[g](i, j0, k) * res1 /
                (mu[g](i, j0 + 1, k) * SQR(mMetric[g](1, i, j0 + 1, k)) +
                 mu[g](i, j0, k) * SQR(mMetric[g](1, i, j0, k)));
        U[g](2, i, j0 + 1, k) =
            U[g](2, i, j0 + 1, k) +
            mJ[g](i, j0, k) * res2 /
                ((mu[g](i, j0 + 1, k) + 0.5 * lambda[g](i, j0 + 1, k)) *
                     SQR(mMetric[g](1, i, j0 + 1, k)) +
                 (mu[g](i, j0, k) + 0.5 * lambda[g](i, j0, k)) *
                     SQR(mMetric[g](1, i, j0, k)));
        U[g](3, i, j0 + 1, k) =
            U[g](3, i, j0 + 1, k) +
            2 * mJ[g](i, j0, k) * res3 /
                (mu[g](i, j0 + 1, k) * SQR(mMetric[g](1, i, j0 + 1, k)) +
                 mu[g](i, j0, k) * SQR(mMetric[g](1, i, j0, k)));
      }
      // Upper bndry
      if (high_interior) {
        res1 = crf * rho[g](i, j1, k) *
                   (bnd1[0] - 2 * U[g](1, i, j1, k) + Um[g](1, i, j1, k)) *
                   d2i -
               Lu1(1, i, j1, k) - forcing[g](1, i, j1, k);
        res2 = crf * rho[g](i, j1, k) *
                   (bnd1[1] - 2 * U[g](2, i, j1, k) + Um[g](2, i, j1, k)) *
                   d2i -
               Lu1(2, i, j1, k) - forcing[g](2, i, j1, k);
        res3 = crf * rho[g](i, j1, k) *
                   (bnd1[2] - 2 * U[g](3, i, j1, k) + Um[g](3, i, j1, k)) *
                   d2i -
               Lu1(3, i, j1, k) - forcing[g](3, i, j1, k);

        //	    U[g](1,i,j1-1,k) = U[g](1,i,j1-1,k) +
        //2*h2*res1/(mu[g](i,j1-1,k)+mu[g](i,j1,k)); 	    U[g](2,i,j1-1,k) =
        //U[g](2,i,j1-1,k) +   h2*res2/(mu[g](i,j1-1,k)+mu[g](i,j1,k)+
        //							     0.5*(lambda[g](i,j1-1,k)+lambda[g](i,j1,k)));
        //	    U[g](3,i,j1-1,k) = U[g](3,i,j1-1,k) +
        //2*h2*res3/(mu[g](i,j1-1,k)+mu[g](i,j1,k));
        U[g](1, i, j1 - 1, k) =
            U[g](1, i, j1 - 1, k) +
            2 * mJ[g](i, j1, k) * res1 /
                (mu[g](i, j1 - 1, k) * SQR(mMetric[g](1, i, j1 - 1, k)) +
                 mu[g](i, j1, k) * SQR(mMetric[g](1, i, j1, k)));
        U[g](2, i, j1 - 1, k) =
            U[g](2, i, j1 - 1, k) +
            mJ[g](i, j1, k) * res2 /
                ((mu[g](i, j1 - 1, k) + 0.5 * lambda[g](i, j1 - 1, k)) *
                     SQR(mMetric[g](1, i, j1 - 1, k)) +
                 (mu[g](i, j1, k) + 0.5 * lambda[g](i, j1, k)) *
                     SQR(mMetric[g](1, i, j1, k)));
        U[g](3, i, j1 - 1, k) =
            U[g](3, i, j1 - 1, k) +
            2 * mJ[g](i, j1, k) * res3 /
                (mu[g](i, j1 - 1, k) * SQR(mMetric[g](1, i, j1 - 1, k)) +
                 mu[g](i, j1, k) * SQR(mMetric[g](1, i, j1, k)));
      }
    }
  // Side with k=const
  low_interior = (m_kStartInt[g] <= k0 + 1 && k0 + 1 <= m_kEndInt[g]) &&
                 !surface_correction;
  high_interior = m_kStartInt[g] <= k1 - 1 && k1 - 1 <= m_kEndInt[g];
  if (m_geodyn_faces == 6 && low_interior) {
    //      cout << "geodyn low_interior curvi " << endl;
    Lu0.define(3, i0, i1, j0, j1, k0, k0);
    evalLuCurv<1, 1, 1, 1, 0, 1>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu0,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i0, i1, j0, j1, k0, k0);
  }
  if (high_interior) {
    //      cout << "geodyn high_interior curvi " << endl;
    Lu1.define(3, i0, i1, j0, j1, k1, k1);
    evalLuCurv<1, 1, 1, 1, 1, 0>(m_iStart[g], m_iEnd[g], m_jStart[g], m_jEnd[g],
                                 m_kStart[g], m_kEnd[g], U[g], Lu1,
                                 mu[g].c_ptr(), lambda[g].c_ptr(), mMetric[g],
                                 mJ[g], i0, i1, j0, j1, k1, k1);
  }
#pragma omp parallel for
  for (int j = j0 + 1; j <= j1 - 1; j++)
    for (int i = i0 + 1; i <= i1 - 1; i++) {
      int ig0 = static_cast<int>(
          floor(((i - 1) * h - m_geodyn_origin[0]) / m_geodyn_h + 1));
      int jg0 = static_cast<int>(
          floor(((j - 1) * h - m_geodyn_origin[1]) / m_geodyn_h + 1));
      if (ig0 >= m_geodyn_ni) ig0 = m_geodyn_ni - 1;
      if (ig0 <= 0) ig0 = 1;
      if (jg0 >= m_geodyn_nj) jg0 = m_geodyn_nj - 1;
      if (jg0 <= 0) jg0 = 1;
      float_sw4 wghi =
          ((i - 1) * h - (m_geodyn_origin[0] + (ig0 - 1) * m_geodyn_h)) /
          m_geodyn_h;
      float_sw4 wghj =
          ((j - 1) * h - (m_geodyn_origin[1] + (jg0 - 1) * m_geodyn_h)) /
          m_geodyn_h;
      float_sw4 bnd0[3], bnd1[3];
      for (int c = 1; c <= 3; c++) {
        if (m_geodyn_faces == 6 && low_interior) {
          bnd0[c - 1] =
              twgh *
                  ((1 - wghi) * (1 - wghj) * m_geodyn_data1[4](c, ig0, jg0, 1) +
                   wghi * (1 - wghj) * m_geodyn_data1[4](c, ig0 + 1, jg0, 1) +
                   (1 - wghi) * wghj * m_geodyn_data1[4](c, ig0, jg0 + 1, 1) +
                   wghi * wghj * m_geodyn_data1[4](c, ig0 + 1, jg0 + 1, 1)) +
              (1 - twgh) *
                  ((1 - wghi) * (1 - wghj) * m_geodyn_data2[4](c, ig0, jg0, 1) +
                   wghi * (1 - wghj) * m_geodyn_data2[4](c, ig0 + 1, jg0, 1) +
                   (1 - wghi) * wghj * m_geodyn_data2[4](c, ig0, jg0 + 1, 1) +
                   wghi * wghj * m_geodyn_data2[4](c, ig0 + 1, jg0 + 1, 1));
        }
        bnd1[c - 1] =
            twgh *
                ((1 - wghi) * (1 - wghj) * m_geodyn_data1[5](c, ig0, jg0, 1) +
                 wghi * (1 - wghj) * m_geodyn_data1[5](c, ig0 + 1, jg0, 1) +
                 (1 - wghi) * wghj * m_geodyn_data1[5](c, ig0, jg0 + 1, 1) +
                 wghi * wghj * m_geodyn_data1[5](c, ig0 + 1, jg0 + 1, 1)) +
            (1 - twgh) *
                ((1 - wghi) * (1 - wghj) * m_geodyn_data2[5](c, ig0, jg0, 1) +
                 wghi * (1 - wghj) * m_geodyn_data2[5](c, ig0 + 1, jg0, 1) +
                 (1 - wghi) * wghj * m_geodyn_data2[5](c, ig0, jg0 + 1, 1) +
                 wghi * wghj * m_geodyn_data2[5](c, ig0 + 1, jg0 + 1, 1));
      }
      // Upper bndry
      float_sw4 res1, res2, res3;
      if (high_interior) {
        res1 = crf * rho[g](i, j, k1) *
                   (bnd1[0] - 2 * U[g](1, i, j, k1) + Um[g](1, i, j, k1)) *
                   d2i -
               Lu1(1, i, j, k1) - forcing[g](1, i, j, k1);
        res2 = crf * rho[g](i, j, k1) *
                   (bnd1[1] - 2 * U[g](2, i, j, k1) + Um[g](2, i, j, k1)) *
                   d2i -
               Lu1(2, i, j, k1) - forcing[g](2, i, j, k1);
        res3 = crf * rho[g](i, j, k1) *
                   (bnd1[2] - 2 * U[g](3, i, j, k1) + Um[g](3, i, j, k1)) *
                   d2i -
               Lu1(3, i, j, k1) - forcing[g](3, i, j, k1);

        //	    U[g](1,i,j,k1-1) = U[g](1,i,j,k1-1) +
        //2*h2*res1/(mu[g](i,j,k1-1)+mu[g](i,j,k1)); 	    U[g](2,i,j,k1-1) =
        //U[g](2,i,j,k1-1) + 2*h2*res2/(mu[g](i,j,k1-1)+mu[g](i,j,k1));
        //	    U[g](3,i,j,k1-1) = U[g](3,i,j,k1-1) +
        //h2*res3/(mu[g](i,j,k1-1)+mu[g](i,j,k1)+
        //							     0.5*(lambda[g](i,j,k1-1)+lambda[g](i,j,k1)));
        double x[3], amat_[9];
#define a(i, j) amat_[((i)-1) + 3 * ((j)-1)]
        x[0] = 2 * mJ[g](i, j, k1) * res1;
        x[1] = 2 * mJ[g](i, j, k1) * res2;
        x[2] = 2 * mJ[g](i, j, k1) * res3;
        a(1, 1) = (2 * mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      SQR(mMetric[g](2, i, j, k1 - 1)) +
                  (2 * mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      SQR(mMetric[g](2, i, j, k1)) +
                  mu[g](i, j, k1 - 1) * (SQR(mMetric[g](3, i, j, k1 - 1)) +
                                         SQR(mMetric[g](4, i, j, k1 - 1))) +
                  mu[g](i, j, k1) * (SQR(mMetric[g](3, i, j, k1)) +
                                     SQR(mMetric[g](4, i, j, k1)));

        a(1, 2) = (mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      mMetric[g](2, i, j, k1 - 1) *
                      mMetric[g](3, i, j, k1 - 1) +
                  (mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      mMetric[g](2, i, j, k1) * mMetric[g](3, i, j, k1);
        a(1, 3) = (mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      mMetric[g](2, i, j, k1 - 1) *
                      mMetric[g](4, i, j, k1 - 1) +
                  (mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      mMetric[g](2, i, j, k1) * mMetric[g](4, i, j, k1);

        a(2, 1) = a(1, 2);
        a(2, 2) = (2 * mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      SQR(mMetric[g](3, i, j, k1 - 1)) +
                  (2 * mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      SQR(mMetric[g](3, i, j, k1)) +
                  mu[g](i, j, k1 - 1) * (SQR(mMetric[g](2, i, j, k1 - 1)) +
                                         SQR(mMetric[g](4, i, j, k1 - 1))) +
                  mu[g](i, j, k1) * (SQR(mMetric[g](2, i, j, k1)) +
                                     SQR(mMetric[g](4, i, j, k1)));
        a(2, 3) = (mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      mMetric[g](3, i, j, k1 - 1) *
                      mMetric[g](4, i, j, k1 - 1) +
                  (mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      mMetric[g](3, i, j, k1) * mMetric[g](4, i, j, k1);
        a(3, 1) = a(1, 3);
        a(3, 2) = a(2, 3);
        a(3, 3) = (2 * mu[g](i, j, k1 - 1) + lambda[g](i, j, k1 - 1)) *
                      SQR(mMetric[g](4, i, j, k1 - 1)) +
                  (2 * mu[g](i, j, k1) + lambda[g](i, j, k1)) *
                      SQR(mMetric[g](4, i, j, k1)) +
                  mu[g](i, j, k1 - 1) * (SQR(mMetric[g](2, i, j, k1 - 1)) +
                                         SQR(mMetric[g](3, i, j, k1 - 1))) +
                  mu[g](i, j, k1) * (SQR(mMetric[g](2, i, j, k1)) +
                                     SQR(mMetric[g](3, i, j, k1)));
        int ipiv[3], info, three = 3, one = 1;
        F77_FUNC(dgesv, DGESV)
        (&three, &one, amat_, &three, ipiv, x, &three, &info);
        VERIFY2(
            info == 0,
            "ERROR: info = "
                << info
                << " from DGESV in geodyn_second_ghost_point_curvilinear\n");
        //	    if( info != 0 )
        //	       cout << "ERROR: info = " << info << " from DGESV in
        //geodyn_second_ghost_point_curvilinear " << endl;
        U[g](1, i, j, k1 - 1) = U[g](1, i, j, k1 - 1) + x[0];
        U[g](2, i, j, k1 - 1) = U[g](2, i, j, k1 - 1) + x[1];
        U[g](3, i, j, k1 - 1) = U[g](3, i, j, k1 - 1) + x[2];
      }
      // Lower bndry
      if (m_geodyn_faces == 6 && low_interior) {
        res1 = crf * rho[g](i, j, k0) *
                   (bnd0[0] - 2 * U[g](1, i, j, k0) + Um[g](1, i, j, k0)) *
                   d2i -
               Lu0(1, i, j, k0) - forcing[g](1, i, j, k0);
        res2 = crf * rho[g](i, j, k0) *
                   (bnd0[1] - 2 * U[g](2, i, j, k0) + Um[g](2, i, j, k0)) *
                   d2i -
               Lu0(2, i, j, k0) - forcing[g](2, i, j, k0);
        res3 = crf * rho[g](i, j, k0) *
                   (bnd0[2] - 2 * U[g](3, i, j, k0) + Um[g](3, i, j, k0)) *
                   d2i -
               Lu0(3, i, j, k0) - forcing[g](3, i, j, k0);

        //	    U[g](1,i,j,k0+1) = U[g](1,i,j,k0+1) +
        //2*h2*res1/(mu[g](i,j,k0+1)+mu[g](i,j,k0)); 	    U[g](2,i,j,k0+1) =
        //U[g](2,i,j,k0+1) + 2*h2*res2/(mu[g](i,j,k0+1)+mu[g](i,j,k0));
        //	    U[g](3,i,j,k0+1) = U[g](3,i,j,k0+1) +
        //h2*res3/(mu[g](i,j,k0+1)+mu[g](i,j,k0)+
        //							     0.5*(lambda[g](i,j,k0+1)+lambda[g](i,j,k0)));
        double x[3], amat_[9];
        x[0] = 2 * mJ[g](i, j, k0) * res1;
        x[1] = 2 * mJ[g](i, j, k0) * res2;
        x[2] = 2 * mJ[g](i, j, k0) * res3;
        a(1, 1) = (2 * mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      SQR(mMetric[g](2, i, j, k0 + 1)) +
                  (2 * mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      SQR(mMetric[g](2, i, j, k0)) +
                  mu[g](i, j, k0 + 1) * (SQR(mMetric[g](3, i, j, k0 + 1)) +
                                         SQR(mMetric[g](4, i, j, k0 + 1))) +
                  mu[g](i, j, k0) * (SQR(mMetric[g](3, i, j, k0)) +
                                     SQR(mMetric[g](4, i, j, k0)));

        a(1, 2) = (mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      mMetric[g](2, i, j, k0 + 1) *
                      mMetric[g](3, i, j, k0 + 1) +
                  (mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      mMetric[g](2, i, j, k0) * mMetric[g](3, i, j, k0);
        a(1, 3) = (mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      mMetric[g](2, i, j, k0 + 1) *
                      mMetric[g](4, i, j, k0 + 1) +
                  (mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      mMetric[g](2, i, j, k0) * mMetric[g](4, i, j, k0);

        a(2, 1) = a(1, 2);
        a(2, 2) = (2 * mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      SQR(mMetric[g](3, i, j, k0 + 1)) +
                  (2 * mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      SQR(mMetric[g](3, i, j, k0)) +
                  mu[g](i, j, k0 + 1) * (SQR(mMetric[g](2, i, j, k0 + 1)) +
                                         SQR(mMetric[g](4, i, j, k0 + 1))) +
                  mu[g](i, j, k0) * (SQR(mMetric[g](2, i, j, k0)) +
                                     SQR(mMetric[g](4, i, j, k0)));
        a(2, 3) = (mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      mMetric[g](3, i, j, k0 + 1) *
                      mMetric[g](4, i, j, k0 + 1) +
                  (mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      mMetric[g](3, i, j, k0) * mMetric[g](4, i, j, k0);
        a(3, 1) = a(1, 3);
        a(3, 2) = a(2, 3);
        a(3, 3) = (2 * mu[g](i, j, k0 + 1) + lambda[g](i, j, k0 + 1)) *
                      SQR(mMetric[g](4, i, j, k0 + 1)) +
                  (2 * mu[g](i, j, k0) + lambda[g](i, j, k0)) *
                      SQR(mMetric[g](4, i, j, k0)) +
                  mu[g](i, j, k0 + 1) * (SQR(mMetric[g](2, i, j, k0 + 1)) +
                                         SQR(mMetric[g](3, i, j, k0 + 1))) +
                  mu[g](i, j, k0) * (SQR(mMetric[g](2, i, j, k0)) +
                                     SQR(mMetric[g](3, i, j, k0)));
        int ipiv[3], info, three = 3, one = 1;
        F77_FUNC(dgesv, DGESV)
        (&three, &one, amat_, &three, ipiv, x, &three, &info);
        VERIFY2(
            info == 0,
            "ERROR: info = "
                << info
                << " from DGESV in geodyn_second_ghost_point_curvilinear\n");
        //	    if( info != 0 )
        //	       cout << "ERROR: info = " << info << " from DGESV in
        //geodyn_second_ghost_point_curvilinear " << endl;
        U[g](1, i, j, k0 + 1) = U[g](1, i, j, k0 + 1) + x[0];
        U[g](2, i, j, k0 + 1) = U[g](2, i, j, k0 + 1) + x[1];
        U[g](3, i, j, k0 + 1) = U[g](3, i, j, k0 + 1) + x[2];
      }
#undef a
    }
}

//-----------------------------------------------------------------------
void evalLu_Dip(
    int ib, int ie, int jb, int je, int kb, int ke,
    //		 double* a_u, double* a_lu, double* a_mu, double* a_la,
    Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la, float_sw4 h,
    int ilb, int ile, int jlb, int jle, int klb, int kle) {
  // Suggested change: Input Sarray& a_lu instead and define
  //   const long int lb= a_lu.m_base;
  //   const size_t loi=a_lu.m_offi;
  //   const size_t loj=a_lu.m_offj;
  //   const size_t lok=a_lu.m_offk;
  //   const size_t loc=a_lu.m_offc;
  ////   float_sw4* a_lupt = a_lu.c_ptr();
  ////#define lu(c,i,j,k) a_lupt[lb+loc*(c)+loi*(i)+loj*(j)+lok*(k)]
  //   float_sw4* a_lupt = &(a_lu.m_data[lb]);
  //#define lu(c,i,j,k) a_lupt[loc*(c)+loi*(i)+loj*(j)+lok*(k)]
  //   const long int b=a_u.m_base;
  //   const size_t oi=a_u.m_offi;
  //   const size_t oj=a_u.m_offj;
  //   const size_t ok=a_u.m_offk;
  //   const size_t oc=a_u.m_offc;
  //   float_sw4* a_upt=a_u.c_ptr();
  //   float_sw4* a_mupt=a_mu.c_ptr();
  //   float_sw4* a_lapt=a_la.c_ptr();
  //#define u(c,i,j,k) a_upt[b+oc*(c)+oi*(i)+oj*(j)+ok*(k)]
  //#define mu(i,j,k) a_mupt[b+oi*(i)+oj*(j)+ok*(k)]
  //#define la(i,j,k) a_lapt[b+oi*(i)+oj*(j)+ok*(k)]
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  const size_t nijk = nij * (ke - kb + 1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb + 1; k <= kle - 1; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb; i <= ile; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +
             half * (la(i + 1, j, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k) +
                          u(3, i + 1, j, k + 1) - u(3, i + 1, j, k - 1)) -
                     la(i, j, k) * (u(2, i, j + 1, k) - u(2, i, j - 1, k) +
                                    u(3, i, j, k + 1) - u(3, i, j, k - 1))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i + 1, j - 1, k) - u(2, i, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i + 1, j, k + 1) - u(3, i, j, k + 1)) -
                     mu(i, j, k - 1) *
                         (u(3, i + 1, j, k - 1) - u(3, i, j, k - 1))) +
             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +
             fourth * (la(i, j + 1, k) *
                           (2 * (u(1, i + 1, j + 1, k) - u(1, i, j + 1, k)) +
                            u(3, i, j + 1, k + 1) - u(3, i, j + 1, k - 1)) -
                       la(i, j - 1, k) *
                           (2 * (u(1, i + 1, j - 1, k) - u(1, i, j - 1, k)) +
                            u(3, i, j - 1, k + 1) - u(3, i, j - 1, k - 1))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k)) -
                     mu(i, j, k) * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) +
             fourth * (mu(i, j, k + 1) *
                           (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1)) -
                       mu(i, j, k - 1) *
                           (u(3, i, j + 1, k - 1) - u(3, i, j - 1, k - 1))) +
             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +
             fourth * (la(i, j, k + 1) *
                           (2 * (u(1, i + 1, j, k + 1) - u(1, i, j, k + 1)) +
                            u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1)) -
                       la(i, j, k - 1) *
                           (2 * (u(1, i + 1, j, k - 1) - u(1, i, j, k - 1)) +
                            u(2, i, j + 1, k - 1) - u(2, i, j - 1, k - 1))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k - 1)) -
                     mu(i, j, k) * (u(1, i, j, k + 1) - u(1, i, j, k - 1))) +
             fourth * (mu(i, j + 1, k) *
                           (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k - 1)) -
                       mu(i, j - 1, k) *
                           (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k - 1))) +
             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_Dim(int ib, int ie, int jb, int je, int kb, int ke,
                //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                //float_sw4* a_la,
                Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb + 1; k <= kle - 1; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb; i <= ile; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +
             half * (la(i, j, k) * (u(2, i, j + 1, k) - u(2, i, j - 1, k) +
                                    u(3, i, j, k + 1) - u(3, i, j, k - 1)) -
                     la(i - 1, j, k) *
                         (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k) +
                          u(3, i - 1, j, k + 1) - u(3, i - 1, j, k - 1))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k) - u(2, i - 1, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k) - u(2, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j, k + 1) - u(3, i - 1, j, k + 1)) -
                     mu(i, j, k - 1) *
                         (u(3, i, j, k - 1) - u(3, i - 1, j, k - 1))) +
             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +
             fourth * (la(i, j + 1, k) *
                           (2 * (u(1, i, j + 1, k) - u(1, i - 1, j + 1, k)) +
                            u(3, i, j + 1, k + 1) - u(3, i, j + 1, k - 1)) -
                       la(i, j - 1, k) *
                           (2 * (u(1, i, j - 1, k) - u(1, i - 1, j - 1, k)) +
                            u(3, i, j - 1, k + 1) - u(3, i, j - 1, k - 1))) +
             half * (mu(i, j, k) * (u(1, i, j + 1, k) - u(1, i, j - 1, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k))) +
             fourth * (mu(i, j, k + 1) *
                           (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1)) -
                       mu(i, j, k - 1) *
                           (u(3, i, j + 1, k - 1) - u(3, i, j - 1, k - 1))) +
             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +
             fourth * (la(i, j, k + 1) *
                           (2 * (u(1, i, j, k + 1) - u(1, i - 1, j, k + 1)) +
                            u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1)) -
                       la(i, j, k - 1) *
                           (2 * (u(1, i, j, k - 1) - u(1, i - 1, j, k - 1)) +
                            u(2, i, j + 1, k - 1) - u(2, i, j - 1, k - 1))) +
             half * (mu(i, j, k) * (u(1, i, j, k + 1) - u(1, i, j, k - 1)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k - 1))) +
             fourth * (mu(i, j + 1, k) *
                           (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k - 1)) -
                       mu(i, j - 1, k) *
                           (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k - 1))) +
             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_Djp(int ib, int ie, int jb, int je, int kb, int ke,
                //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                //float_sw4* a_la,
                Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb + 1; k <= kle - 1; k++)
    for (int j = jlb; j <= jle; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             fourth * (la(i + 1, j, k) *
                           (2 * (u(2, i + 1, j + 1, k) - u(2, i + 1, j, k)) +
                            u(3, i + 1, j, k + 1) - u(3, i + 1, j, k - 1)) -
                       la(i - 1, j, k) *
                           (2 * (u(2, i - 1, j + 1, k) - u(2, i - 1, j, k)) +
                            u(3, i - 1, j, k + 1) - u(3, i - 1, j, k - 1))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k)) -
                     mu(i, j, k) * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) +
             fourth * (mu(i, j, k + 1) *
                           (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1)) -
                       mu(i, j, k - 1) *
                           (u(3, i + 1, j, k - 1) - u(3, i - 1, j, k - 1))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j + 1, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k) +
                          u(3, i, j + 1, k + 1) - u(3, i, j + 1, k - 1)) -
                     la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    u(3, i, j, k + 1) - u(3, i, j, k - 1))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i + 1, j, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j + 1, k) - u(1, i - 1, j, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j + 1, k + 1) - u(3, i, j, k + 1)) -
                     mu(i, j, k - 1) *
                         (u(3, i, j + 1, k - 1) - u(3, i, j, k - 1))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             fourth * (la(i, j, k + 1) *
                           (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1) +
                            2 * (u(2, i, j + 1, k + 1) - u(2, i, j, k + 1))) -
                       la(i, j, k - 1) *
                           (u(1, i + 1, j, k - 1) - u(1, i - 1, j, k - 1) +
                            2 * (u(2, i, j + 1, k - 1) - u(2, i, j, k - 1)))) +
             fourth * (mu(i + 1, j, k) *
                           (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k - 1)) -
                       mu(i - 1, j, k) *
                           (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k - 1))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k - 1)) -
                     mu(i, j, k) * (u(2, i, j, k + 1) - u(2, i, j, k - 1))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_Djm(int ib, int ie, int jb, int je, int kb, int ke,
                //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                //float_sw4* a_la,
                Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb + 1; k <= kle - 1; k++)
    for (int j = jlb; j <= jle; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             fourth * (la(i + 1, j, k) *
                           (2 * (u(2, i + 1, j, k) - u(2, i + 1, j - 1, k)) +
                            u(3, i + 1, j, k + 1) - u(3, i + 1, j, k - 1)) -
                       la(i - 1, j, k) *
                           (2 * (u(2, i - 1, j, k) - u(2, i - 1, j - 1, k)) +
                            u(3, i - 1, j, k + 1) - u(3, i - 1, j, k - 1))) +
             half * (mu(i, j, k) * (u(2, i + 1, j, k) - u(2, i - 1, j, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k))) +
             fourth * (mu(i, j, k + 1) *
                           (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1)) -
                       mu(i, j, k - 1) *
                           (u(3, i + 1, j, k - 1) - u(3, i - 1, j, k - 1))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    u(3, i, j, k + 1) - u(3, i, j, k - 1)) -
                     la(i, j - 1, k) *
                         (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k) +
                          u(3, i, j - 1, k + 1) - u(3, i, j - 1, k - 1))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k) - u(1, i + 1, j - 1, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k) - u(1, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j, k + 1) - u(3, i, j - 1, k + 1)) -
                     mu(i, j, k - 1) *
                         (u(3, i, j, k - 1) - u(3, i, j - 1, k - 1))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +
             fourth * (la(i, j, k + 1) *
                           (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1) +
                            2 * (u(2, i, j, k + 1) - u(2, i, j - 1, k + 1))) -
                       la(i, j, k - 1) *
                           (u(1, i + 1, j, k - 1) - u(1, i - 1, j, k - 1) +
                            2 * (u(2, i, j, k - 1) - u(2, i, j - 1, k - 1)))) +
             fourth * (mu(i + 1, j, k) *
                           (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k - 1)) -
                       mu(i - 1, j, k) *
                           (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k - 1))) +
             half * (mu(i, j, k) * (u(2, i, j, k + 1) - u(2, i, j, k - 1)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k - 1))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_Dkp(int ib, int ie, int jb, int je, int kb, int ke,
                //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                //float_sw4* a_la,
                Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= kle; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             fourth * (la(i + 1, j, k) *
                           (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k) +
                            2 * (u(3, i + 1, j, k + 1) - u(3, i + 1, j, k))) -
                       la(i - 1, j, k) *
                           (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k) +
                            2 * (u(3, i - 1, j, k + 1) - u(3, i - 1, j, k)))) +
             fourth * (mu(i, j + 1, k) *
                           (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k)) -
                       mu(i, j - 1, k) *
                           (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1)) -
                     mu(i, j, k) * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             fourth * (la(i, j + 1, k) *
                           (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k) +
                            2 * (u(3, i, j + 1, k + 1) - u(3, i, j + 1, k))) -
                       la(i, j - 1, k) *
                           (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k) +
                            2 * (u(3, i, j - 1, k + 1) - u(3, i, j - 1, k)))) +
             fourth * (mu(i + 1, j, k) *
                           (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k)) -
                       mu(i - 1, j, k) *
                           (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1)) -
                     mu(i, j, k) * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k + 1) *
                         (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1) +
                          u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1)) -
                     la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    u(2, i, j + 1, k) - u(2, i, j - 1, k))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_Dkm(int ib, int ie, int jb, int je, int kb, int ke,
                //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                //float_sw4* a_la,
                Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= kle; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             fourth * (la(i + 1, j, k) *
                           (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k) +
                            2 * (u(3, i + 1, j, k) - u(3, i + 1, j, k - 1))) -
                       la(i - 1, j, k) *
                           (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k) +
                            2 * (u(3, i - 1, j, k) - u(3, i - 1, j, k - 1)))) +
             fourth * (mu(i, j + 1, k) *
                           (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k)) -
                       mu(i, j - 1, k) *
                           (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k))) +
             half * (mu(i, j, k) * (u(3, i + 1, j, k) - u(3, i - 1, j, k)) -
                     mu(i, j, k - 1) *
                         (u(3, i + 1, j, k - 1) - u(3, i - 1, j, k - 1))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             fourth * (la(i, j + 1, k) *
                           (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k) +
                            2 * (u(3, i, j + 1, k) - u(3, i, j + 1, k - 1))) -
                       la(i, j - 1, k) *
                           (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k) +
                            2 * (u(3, i, j - 1, k) - u(3, i, j - 1, k - 1)))) +
             fourth * (mu(i + 1, j, k) *
                           (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k)) -
                       mu(i - 1, j, k) *
                           (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k))) +
             half * (mu(i, j, k) * (u(3, i, j + 1, k) - u(3, i, j - 1, k)) -
                     mu(i, j, k - 1) *
                         (u(3, i, j + 1, k - 1) - u(3, i, j - 1, k - 1))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    u(2, i, j + 1, k) - u(2, i, j - 1, k)) -
                     la(i, j, k - 1) *
                         (u(1, i + 1, j, k - 1) - u(1, i - 1, j, k - 1) +
                          u(2, i, j + 1, k - 1) - u(2, i, j - 1, k - 1))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k) - u(1, i + 1, j, k - 1)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k) - u(1, i - 1, j, k - 1))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k) - u(2, i, j + 1, k - 1)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k) - u(2, i, j - 1, k - 1))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_DkpDip(int ib, int ie, int jb, int je, int kb, int ke,
                   //		    float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                   //float_sw4* a_la,
                   Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                   float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                   int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= klb; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb; i <= ile; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             half * (la(i + 1, j, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i + 1, j - 1, k) +
                          2 * (u(3, i + 1, j, k + 1) - u(3, i + 1, j, k))) -
                     la(i, j, k) * (u(2, i, j + 1, k) - u(2, i, j - 1, k) +
                                    2 * (u(3, i, j, k + 1) - u(3, i, j, k)))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i + 1, j - 1, k) - u(2, i, j - 1, k))) +
             (mu(i, j, k + 1) * (u(3, i + 1, j, k + 1) - u(3, i, j, k + 1)) -
              mu(i, j, k) * (u(3, i + 1, j, k) - u(3, i, j, k))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j + 1, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i, j + 1, k) +
                          (u(3, i, j + 1, k + 1) - u(3, i, j + 1, k))) -
                     la(i, j - 1, k) *
                         (u(1, i + 1, j - 1, k) - u(1, i, j - 1, k) +
                          (u(3, i, j - 1, k + 1) - u(3, i, j - 1, k)))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i + 1, j - 1, k)) -
                     mu(i, j, k) * (u(1, i, j + 1, k) - u(1, i, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1)) -
                     mu(i, j, k) * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k + 1) *
                         (2 * (u(1, i + 1, j, k + 1) - u(1, i, j, k + 1)) +
                          u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1)) -
                     la(i, j, k) * (2 * (u(1, i + 1, j, k) - u(1, i, j, k)) +
                                    u(2, i, j + 1, k) - u(2, i, j - 1, k))) +
             (mu(i + 1, j, k) * (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k)) -
              mu(i, j, k) * (u(1, i, j, k + 1) - u(1, i, j, k))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k))) +
             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_DkpDim(int ib, int ie, int jb, int je, int kb, int ke,
                   //		    float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                   //float_sw4* a_la,
                   Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                   float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                   int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= klb; k++)
    for (int j = jlb + 1; j <= jle - 1; j++)
      for (int i = ilb; i <= ile; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             half * (la(i, j, k) * (u(2, i, j + 1, k) - u(2, i, j - 1, k) +
                                    2 * (u(3, i, j, k + 1) - u(3, i, j, k))) -
                     la(i - 1, j, k) *
                         (u(2, i - 1, j + 1, k) - u(2, i - 1, j - 1, k) +
                          2 * (u(3, i - 1, j, k + 1) - u(3, i - 1, j, k)))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k) - u(2, i - 1, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k) - u(2, i - 1, j - 1, k))) +
             (mu(i, j, k + 1) * (u(3, i, j, k + 1) - u(3, i - 1, j, k + 1)) -
              mu(i, j, k) * (u(3, i, j, k) - u(3, i - 1, j, k))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j + 1, k) *
                         (u(1, i, j + 1, k) - u(1, i - 1, j + 1, k) +
                          (u(3, i, j + 1, k + 1) - u(3, i, j + 1, k))) -
                     la(i, j - 1, k) *
                         (u(1, i, j - 1, k) - u(1, i - 1, j - 1, k) +
                          (u(3, i, j - 1, k + 1) - u(3, i, j - 1, k)))) +
             half * (mu(i, j, k) * (u(1, i, j + 1, k) - u(1, i, j - 1, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j + 1, k) - u(1, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i, j + 1, k + 1) - u(3, i, j - 1, k + 1)) -
                     mu(i, j, k) * (u(3, i, j + 1, k) - u(3, i, j - 1, k))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k + 1) *
                         (2 * (u(1, i, j, k + 1) - u(1, i - 1, j, k + 1)) +
                          u(2, i, j + 1, k + 1) - u(2, i, j - 1, k + 1)) -
                     la(i, j, k) * (2 * (u(1, i, j, k) - u(1, i - 1, j, k)) +
                                    u(2, i, j + 1, k) - u(2, i, j - 1, k))) +
             (mu(i, j, k) * (u(1, i, j, k + 1) - u(1, i, j, k)) -
              mu(i - 1, j, k) * (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_DkpDjp(int ib, int ie, int jb, int je, int kb, int ke,
                   //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                   //float_sw4* a_la,
                   Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                   float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                   int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= klb; k++)
    for (int j = jlb; j <= jle; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             half * (la(i + 1, j, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i + 1, j, k) +
                          (u(3, i + 1, j, k + 1) - u(3, i + 1, j, k))) -
                     la(i - 1, j, k) *
                         (u(2, i - 1, j + 1, k) - u(2, i - 1, j, k) +
                          (u(3, i - 1, j, k + 1) - u(3, i - 1, j, k)))) +
             half * (mu(i, j + 1, k) *
                         (u(2, i + 1, j + 1, k) - u(2, i - 1, j + 1, k)) -
                     mu(i, j, k) * (u(2, i + 1, j, k) - u(2, i - 1, j, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1)) -
                     mu(i, j, k) * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j + 1, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i - 1, j + 1, k) +
                          2 * (u(3, i, j + 1, k + 1) - u(3, i, j + 1, k))) -
                     la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    2 * (u(3, i, j, k + 1) - u(3, i, j, k)))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j + 1, k) - u(1, i + 1, j, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j + 1, k) - u(1, i - 1, j, k))) +
             (mu(i, j, k + 1) * (u(3, i, j + 1, k + 1) - u(3, i, j, k + 1)) -
              mu(i, j, k) * (u(3, i, j + 1, k) - u(3, i, j, k))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k + 1) *
                         (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1) +
                          2 * (u(2, i, j + 1, k + 1) - u(2, i, j, k + 1))) -
                     la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    2 * (u(2, i, j + 1, k) - u(2, i, j, k)))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k))) +
             (mu(i, j + 1, k) * (u(2, i, j + 1, k + 1) - u(2, i, j + 1, k)) -
              mu(i, j, k) * (u(2, i, j, k + 1) - u(2, i, j, k))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void evalLu_DkpDjm(int ib, int ie, int jb, int je, int kb, int ke,
                   //		 float_sw4* a_u, float_sw4* a_lu, float_sw4* a_mu,
                   //float_sw4* a_la,
                   Sarray& u, Sarray& lu, float_sw4* a_mu, float_sw4* a_la,
                   float_sw4 h, int ilb, int ile, int jlb, int jle, int klb,
                   int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  const float_sw4 ih2 = 1 / (h * h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;

  for (int k = klb; k <= klb; k++)
    for (int j = jlb; j <= jle; j++)
      for (int i = ilb + 1; i <= ile - 1; i++) {
        float_sw4 mupx = half * (mu(i, j, k) + mu(i + 1, j, k));
        float_sw4 mumx = half * (mu(i, j, k) + mu(i - 1, j, k));
        float_sw4 mupy = half * (mu(i, j + 1, k) + mu(i, j, k));
        float_sw4 mumy = half * (mu(i, j - 1, k) + mu(i, j, k));
        float_sw4 mupz = half * (mu(i, j, k + 1) + mu(i, j, k));
        float_sw4 mumz = half * (mu(i, j, k - 1) + mu(i, j, k));
        lu(1, i, j, k) =
            ih2 *
            ((2 * mupx + half * (la(i, j, k) + la(i + 1, j, k))) *
                 (u(1, i + 1, j, k) - u(1, i, j, k)) -
             (2 * mumx + half * (la(i, j, k) + la(i - 1, j, k))) *
                 (u(1, i, j, k) - u(1, i - 1, j, k)) +

             half * (la(i + 1, j, k) *
                         (u(2, i + 1, j, k) - u(2, i + 1, j - 1, k) +
                          (u(3, i + 1, j, k + 1) - u(3, i + 1, j, k))) -
                     la(i - 1, j, k) *
                         (u(2, i - 1, j, k) - u(2, i - 1, j - 1, k) +
                          (u(3, i - 1, j, k + 1) - u(3, i - 1, j, k)))) +
             half * (mu(i, j, k) * (u(2, i + 1, j, k) - u(2, i - 1, j, k)) -
                     mu(i, j - 1, k) *
                         (u(2, i + 1, j - 1, k) - u(2, i - 1, j - 1, k))) +
             half * (mu(i, j, k + 1) *
                         (u(3, i + 1, j, k + 1) - u(3, i - 1, j, k + 1)) -
                     mu(i, j, k) * (u(3, i + 1, j, k) - u(3, i - 1, j, k))) +

             mupy * (u(1, i, j + 1, k) - u(1, i, j, k)) -
             mumy * (u(1, i, j, k) - u(1, i, j - 1, k)) +
             mupz * (u(1, i, j, k + 1) - u(1, i, j, k)) -
             mumz * (u(1, i, j, k) - u(1, i, j, k - 1)));

        lu(2, i, j, k) =
            ih2 *
            (mupx * (u(2, i + 1, j, k) - u(2, i, j, k)) -
             mumx * (u(2, i, j, k) - u(2, i - 1, j, k)) +

             half * (la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    2 * (u(3, i, j, k + 1) - u(3, i, j, k))) -
                     la(i, j - 1, k) *
                         (u(1, i + 1, j - 1, k) - u(1, i - 1, j - 1, k) +
                          2 * (u(3, i, j - 1, k + 1) - u(3, i, j - 1, k)))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k) - u(1, i + 1, j - 1, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k) - u(1, i - 1, j - 1, k))) +
             (mu(i, j, k + 1) * (u(3, i, j, k + 1) - u(3, i, j - 1, k + 1)) -
              mu(i, j, k) * (u(3, i, j, k) - u(3, i, j - 1, k))) +

             (2 * mupy + half * (la(i, j, k) + la(i, j + 1, k))) *
                 (u(2, i, j + 1, k) - u(2, i, j, k)) -
             (2 * mumy + half * (la(i, j, k) + la(i, j - 1, k))) *
                 (u(2, i, j, k) - u(2, i, j - 1, k)) +
             mupz * (u(2, i, j, k + 1) - u(2, i, j, k)) -
             mumz * (u(2, i, j, k) - u(2, i, j, k - 1)));

        lu(3, i, j, k) =
            ih2 *
            (mupx * (u(3, i + 1, j, k) - u(3, i, j, k)) -
             mumx * (u(3, i, j, k) - u(3, i - 1, j, k)) +

             half * (la(i, j, k + 1) *
                         (u(1, i + 1, j, k + 1) - u(1, i - 1, j, k + 1) +
                          2 * (u(2, i, j, k + 1) - u(2, i, j - 1, k + 1))) -
                     la(i, j, k) * (u(1, i + 1, j, k) - u(1, i - 1, j, k) +
                                    2 * (u(2, i, j, k) - u(2, i, j - 1, k)))) +
             half * (mu(i + 1, j, k) *
                         (u(1, i + 1, j, k + 1) - u(1, i + 1, j, k)) -
                     mu(i - 1, j, k) *
                         (u(1, i - 1, j, k + 1) - u(1, i - 1, j, k))) +
             (mu(i, j, k) * (u(2, i, j, k + 1) - u(2, i, j, k)) -
              mu(i, j - 1, k) * (u(2, i, j - 1, k + 1) - u(2, i, j - 1, k))) +

             mupy * (u(3, i, j + 1, k) - u(3, i, j, k)) -
             mumy * (u(3, i, j, k) - u(3, i, j - 1, k)) +
             (2 * mupz + half * (la(i, j, k + 1) + la(i, j, k))) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             (2 * mumz + half * (la(i, j, k - 1) + la(i, j, k))) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));
      }
#undef mu
#undef la
#undef u
#undef lu
}

//-----------------------------------------------------------------------
void EW::geodyn_up_from_uacc(vector<Sarray>& Up, vector<Sarray>& Uacc,
                             vector<Sarray>& U, vector<Sarray>& Um,
                             float_sw4 dt) {
  if (m_do_geodynbc) {
    float_sw4 dt2 = dt * dt;
    for (int g = 0; g < mNumberOfGrids; g++) {
      int i0 = m_geodyn_dims[g][0];
      int i1 = m_geodyn_dims[g][1];
      int j0 = m_geodyn_dims[g][2];
      int j1 = m_geodyn_dims[g][3];
      int k0 = m_geodyn_dims[g][4];
      int k1 = m_geodyn_dims[g][5];
      bool low_interior, high_interior;
      low_interior = m_iStartInt[g] <= i0 + 1 && i0 + 1 <= m_iEndInt[g];
      high_interior = m_iStartInt[g] <= i1 - 1 && i1 - 1 <= m_iEndInt[g];
      bool surface_correction = k0 <= 1 && g == mNumberOfGrids - 1;
      int kstart = k0 + 1;
      if (surface_correction) kstart = k0;

      // Side with i=const.
      for (int k = kstart; k <= k1 - 1; k++)
        for (int j = j0 + 1; j <= j1 - 1; j++) {
          if (low_interior) {
            Up[g](1, i0 + 1, j, k) = dt2 * Uacc[g](1, i0 + 1, j, k) +
                                     2 * U[g](1, i0 + 1, j, k) -
                                     Um[g](1, i0 + 1, j, k);
            Up[g](2, i0 + 1, j, k) = dt2 * Uacc[g](2, i0 + 1, j, k) +
                                     2 * U[g](2, i0 + 1, j, k) -
                                     Um[g](2, i0 + 1, j, k);
            Up[g](3, i0 + 1, j, k) = dt2 * Uacc[g](3, i0 + 1, j, k) +
                                     2 * U[g](3, i0 + 1, j, k) -
                                     Um[g](3, i0 + 1, j, k);
          }
          if (high_interior) {
            Up[g](1, i1 - 1, j, k) = dt2 * Uacc[g](1, i1 - 1, j, k) +
                                     2 * U[g](1, i1 - 1, j, k) -
                                     Um[g](1, i1 - 1, j, k);
            Up[g](2, i1 - 1, j, k) = dt2 * Uacc[g](2, i1 - 1, j, k) +
                                     2 * U[g](2, i1 - 1, j, k) -
                                     Um[g](2, i1 - 1, j, k);
            Up[g](3, i1 - 1, j, k) = dt2 * Uacc[g](3, i1 - 1, j, k) +
                                     2 * U[g](3, i1 - 1, j, k) -
                                     Um[g](3, i1 - 1, j, k);
          }
        }
      // Side with j=const
      low_interior = m_jStartInt[g] <= j0 + 1 && j0 + 1 <= m_jEndInt[g];
      high_interior = m_jStartInt[g] <= j1 - 1 && j1 - 1 <= m_jEndInt[g];
      for (int k = kstart; k <= k1 - 1; k++)
        for (int i = i0 + 1; i <= i1 - 1; i++) {
          if (low_interior) {
            Up[g](1, i, j0 + 1, k) = dt2 * Uacc[g](1, i, j0 + 1, k) +
                                     2 * U[g](1, i, j0 + 1, k) -
                                     Um[g](1, i, j0 + 1, k);
            Up[g](2, i, j0 + 1, k) = dt2 * Uacc[g](2, i, j0 + 1, k) +
                                     2 * U[g](2, i, j0 + 1, k) -
                                     Um[g](2, i, j0 + 1, k);
            Up[g](3, i, j0 + 1, k) = dt2 * Uacc[g](3, i, j0 + 1, k) +
                                     2 * U[g](3, i, j0 + 1, k) -
                                     Um[g](3, i, j0 + 1, k);
          }
          if (high_interior) {
            Up[g](1, i, j1 - 1, k) = dt2 * Uacc[g](1, i, j1 - 1, k) +
                                     2 * U[g](1, i, j1 - 1, k) -
                                     Um[g](1, i, j1 - 1, k);
            Up[g](2, i, j1 - 1, k) = dt2 * Uacc[g](2, i, j1 - 1, k) +
                                     2 * U[g](2, i, j1 - 1, k) -
                                     Um[g](2, i, j1 - 1, k);
            Up[g](3, i, j1 - 1, k) = dt2 * Uacc[g](3, i, j1 - 1, k) +
                                     2 * U[g](3, i, j1 - 1, k) -
                                     Um[g](3, i, j1 - 1, k);
          }
        }
      // Side with k=const
      low_interior = m_kStartInt[g] <= k0 + 1 && k0 + 1 <= m_kEndInt[g];
      high_interior = m_kStartInt[g] <= k1 - 1 && k1 - 1 <= m_kEndInt[g];
      for (int j = j0 + 1; j <= j1 - 1; j++)
        for (int i = i0 + 1; i <= i1 - 1; i++) {
          if (low_interior && m_geodyn_faces == 6) {
            Up[g](1, i, j, k0 + 1) = dt2 * Uacc[g](1, i, j, k0 + 1) +
                                     2 * U[g](1, i, j, k0 + 1) -
                                     Um[g](1, i, j, k0 + 1);
            Up[g](2, i, j, k0 + 1) = dt2 * Uacc[g](2, i, j, k0 + 1) +
                                     2 * U[g](2, i, j, k0 + 1) -
                                     Um[g](2, i, j, k0 + 1);
            Up[g](3, i, j, k0 + 1) = dt2 * Uacc[g](3, i, j, k0 + 1) +
                                     2 * U[g](3, i, j, k0 + 1) -
                                     Um[g](3, i, j, k0 + 1);
          }
          if (high_interior) {
            Up[g](1, i, j, k1 - 1) = dt2 * Uacc[g](1, i, j, k1 - 1) +
                                     2 * U[g](1, i, j, k1 - 1) -
                                     Um[g](1, i, j, k1 - 1);
            Up[g](2, i, j, k1 - 1) = dt2 * Uacc[g](2, i, j, k1 - 1) +
                                     2 * U[g](2, i, j, k1 - 1) -
                                     Um[g](2, i, j, k1 - 1);
            Up[g](3, i, j, k1 - 1) = dt2 * Uacc[g](3, i, j, k1 - 1) +
                                     2 * U[g](3, i, j, k1 - 1) -
                                     Um[g](3, i, j, k1 - 1);
          }
        }
    }
  }
}

//-----------------------------------------------------------------------
void EW::save_geoghost(vector<Sarray>& U) {
  //   if( m_do_geodynbc && m_geodyn_faces==5 )
  if (m_do_geodynbc) {
    int g = mNumberOfGrids - 1;
    int i0 = m_geodyn_dims[g][0];
    int i1 = m_geodyn_dims[g][1];
    int j0 = m_geodyn_dims[g][2];
    int j1 = m_geodyn_dims[g][3];
    int k0 = m_geodyn_dims[g][4];
    int k1 = m_geodyn_dims[g][5];
    if (k0 <= 1 && k1 - k0 + 1 > 0) {
      for (int j = j0; j <= j1; j++)
        for (int c = 1; c <= 3; c++) {
          m_geo_usgh[0][c - 1 + 3 * (j - j0)] = U[g](c, i0, j, 0);
          m_geo_usgh[1][c - 1 + 3 * (j - j0)] = U[g](c, i1, j, 0);
        }
      for (int i = i0; i <= i1; i++)
        for (int c = 1; c <= 3; c++) {
          m_geo_usgh[2][c - 1 + 3 * (i - i0)] = U[g](c, i, j0, 0);
          m_geo_usgh[3][c - 1 + 3 * (i - i0)] = U[g](c, i, j1, 0);
        }
    }
  }
}

//-----------------------------------------------------------------------
void EW::restore_geoghost(vector<Sarray>& U) {
  //   if( m_do_geodynbc && m_geodyn_faces==5 )
  {
    int g = mNumberOfGrids - 1;
    int i0 = m_geodyn_dims[g][0];
    int i1 = m_geodyn_dims[g][1];
    int j0 = m_geodyn_dims[g][2];
    int j1 = m_geodyn_dims[g][3];
    int k0 = m_geodyn_dims[g][4];
    int k1 = m_geodyn_dims[g][5];
    if (k0 <= 1 && k1 - k0 + 1 > 0) {
      for (int j = j0; j <= j1; j++)
        for (int c = 1; c <= 3; c++) {
          U[g](c, i0, j, 0) = m_geo_usgh[0][c - 1 + 3 * (j - j0)];
          U[g](c, i1, j, 0) = m_geo_usgh[1][c - 1 + 3 * (j - j0)];
        }
      for (int i = i0; i <= i1; i++)
        for (int c = 1; c <= 3; c++) {
          U[g](c, i, j0, 0) = m_geo_usgh[2][c - 1 + 3 * (i - i0)];
          U[g](c, i, j1, 0) = m_geo_usgh[3][c - 1 + 3 * (i - i0)];
        }
    }
  }
}

//-----------------------------------------------------------------------
void EW::bcsurf_curvilinear_2nd_order(int side, int i0, int i1, int j0, int j1,
                                      int k0, int g, Sarray& u,
                                      float_sw4* bforcing) {
  if (side == 0) {
    //  side  i = i0;
    i1 = i0;
  } else if (side == 1) {
    // side i=i1;
    i0 = i1;
  } else if (side == 2) {
    // side j=j0
    j1 = j0;
  } else if (side == 3) {
    // side j=j1
    j0 = j1;
  }
  int ib = m_iStart[g];
  int jb = m_jStart[g];
  int ni = m_iEnd[g] - m_iStart[g] + 1;
  int nj = m_jEnd[g] - m_jStart[g] + 1;
  Sarray& met = mMetric[g];  // Rename mMetric, makes formulas shorter.
  for (int j = j0; j <= j1; j++)
    for (int i = i0; i <= i1; i++) {
      size_t qq = i - ib + ni * (j - jb);

      // One sided x-derivatives
      float_sw4 ux, vx, wx, uy, vy, wy;
      if (side == 0) {
        ux = u(1, i, j, k0) - u(1, i - 1, j, k0);
        vx = u(2, i, j, k0) - u(2, i - 1, j, k0);
        wx = u(3, i, j, k0) - u(3, i - 1, j, k0);
      } else if (side == 1) {
        ux = u(1, i + 1, j, k0) - u(1, i, j, k0);
        vx = u(2, i + 1, j, k0) - u(2, i, j, k0);
        wx = u(3, i + 1, j, k0) - u(3, i, j, k0);
      } else {
        // Centered x-derivatives
        ux = 0.5 * (u(1, i + 1, j, k0) - u(1, i - 1, j, k0));
        vx = 0.5 * (u(2, i + 1, j, k0) - u(2, i - 1, j, k0));
        wx = 0.5 * (u(3, i + 1, j, k0) - u(3, i - 1, j, k0));
      }
      // One sided y-derivatives
      if (side == 2) {
        uy = u(1, i, j, k0) - u(1, i, j - 1, k0);
        vy = u(2, i, j, k0) - u(2, i, j - 1, k0);
        wy = u(3, i, j, k0) - u(3, i, j - 1, k0);

      } else if (side == 3) {
        uy = u(1, i, j + 1, k0) - u(1, i, j, k0);
        vy = u(2, i, j + 1, k0) - u(2, i, j, k0);
        wy = u(3, i, j + 1, k0) - u(3, i, j, k0);
      } else {
        // Centered y-derivatives
        uy = 0.5 * (u(1, i, j + 1, k0) - u(1, i, j - 1, k0));
        vy = 0.5 * (u(2, i, j + 1, k0) - u(2, i, j - 1, k0));
        wy = 0.5 * (u(3, i, j + 1, k0) - u(3, i, j - 1, k0));
      }

      // Boundary stresses tangential derivatives
      float_sw4 rhs1 =
          (2 * mMu[g](i, j, k0) + mLambda[g](i, j, k0)) * met(1, i, j, k0) *
              met(2, i, j, k0) * ux +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(3, i, j, k0) * vx +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(4, i, j, k0) * wx +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(3, i, j, k0) * uy +
          mLambda[g](i, j, k0) * met(1, i, j, k0) * met(2, i, j, k0) * vy;
      float_sw4 rhs2 =
          mLambda[g](i, j, k0) * met(1, i, j, k0) * met(3, i, j, k0) * ux +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(2, i, j, k0) * vx +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(2, i, j, k0) * uy +
          (2 * mMu[g](i, j, k0) + mLambda[g](i, j, k0)) * met(1, i, j, k0) *
              met(3, i, j, k0) * vy +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(4, i, j, k0) * wy;
      float_sw4 rhs3 =
          mLambda[g](i, j, k0) * met(1, i, j, k0) * met(4, i, j, k0) * ux +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(2, i, j, k0) * wx +
          mMu[g](i, j, k0) * met(1, i, j, k0) * met(3, i, j, k0) * wy +
          mLambda[g](i, j, k0) * met(1, i, j, k0) * met(4, i, j, k0) * vy;

      float_sw4 acp = met(2, i, j, k0 + 1) * met(2, i, j, k0 + 1) +
                      met(3, i, j, k0 + 1) * met(3, i, j, k0 + 1) +
                      met(4, i, j, k0 + 1) * met(4, i, j, k0 + 1);
      float_sw4 ac = met(2, i, j, k0) * met(2, i, j, k0) +
                     met(3, i, j, k0) * met(3, i, j, k0) +
                     met(4, i, j, k0) * met(4, i, j, k0);
      float_sw4 acm = met(2, i, j, k0 - 1) * met(2, i, j, k0 - 1) +
                      met(3, i, j, k0 - 1) * met(3, i, j, k0 - 1) +
                      met(4, i, j, k0 - 1) * met(4, i, j, k0 - 1);

      float_sw4 muplap = mMu[g](i, j, k0 + 1) + mLambda[g](i, j, k0 + 1);
      float_sw4 mupla = mMu[g](i, j, k0) + mLambda[g](i, j, k0);
      float_sw4 muplam = mMu[g](i, j, k0 - 1) + mLambda[g](i, j, k0 - 1);
      double b[3], amat_[9];
#define amat(i, j) amat_[(i)-1 + 3 * ((j)-1)]

      b[0] = 2 * rhs1 - 2 * bforcing[3 * qq] +
             0.5 *
                 (mMu[g](i, j, k0 + 1) * acp + mMu[g](i, j, k0) * ac +
                  muplap * met(2, i, j, k0 + 1) * met(2, i, j, k0 + 1) +
                  mupla * met(2, i, j, k0) * met(2, i, j, k0)) *
                 (u(1, i, j, k0 + 1) - u(1, i, j, k0)) +
             0.5 *
                 (muplap * met(2, i, j, k0 + 1) * met(3, i, j, k0 + 1) +
                  mupla * met(2, i, j, k0) * met(3, i, j, k0)) *
                 (u(2, i, j, k0 + 1) - u(2, i, j, k0)) +
             0.5 *
                 (muplap * met(2, i, j, k0 + 1) * met(4, i, j, k0 + 1) +
                  mupla * met(2, i, j, k0) * met(4, i, j, k0)) *
                 (u(3, i, j, k0 + 1) - u(3, i, j, k0));
      amat(1, 1) = 0.5 * (mMu[g](i, j, k0 - 1) * acm + mMu[g](i, j, k0) * ac +
                          muplam * met(2, i, j, k0 - 1) * met(2, i, j, k0 - 1) +
                          mupla * met(2, i, j, k0) * met(2, i, j, k0));
      amat(1, 2) = 0.5 * (muplam * met(2, i, j, k0 - 1) * met(3, i, j, k0 - 1) +
                          mupla * met(2, i, j, k0) * met(3, i, j, k0));
      amat(1, 3) = 0.5 * (muplam * met(2, i, j, k0 - 1) * met(4, i, j, k0 - 1) +
                          mupla * met(2, i, j, k0) * met(4, i, j, k0));

      b[1] = 2 * rhs2 - 2 * bforcing[1 + 3 * qq] +
             0.5 *
                 (muplap * met(2, i, j, k0 + 1) * met(3, i, j, k0 + 1) +
                  mupla * met(2, i, j, k0) * met(3, i, j, k0)) *
                 (u(1, i, j, k0 + 1) - u(1, i, j, k0)) +
             0.5 *
                 (mMu[g](i, j, k0 + 1) * acp + mMu[g](i, j, k0) * ac +
                  muplap * met(3, i, j, k0 + 1) * met(3, i, j, k0 + 1) +
                  mupla * met(3, i, j, k0) * met(3, i, j, k0)) *
                 (u(2, i, j, k0 + 1) - u(2, i, j, k0)) +
             0.5 *
                 (muplap * met(4, i, j, k0 + 1) * met(3, i, j, k0 + 1) +
                  mupla * met(4, i, j, k0) * met(3, i, j, k0)) *
                 (u(3, i, j, k0 + 1) - u(3, i, j, k0));

      amat(2, 1) = 0.5 * (muplam * met(2, i, j, k0 - 1) * met(3, i, j, k0 - 1) +
                          mupla * met(2, i, j, k0) * met(3, i, j, k0));
      amat(2, 2) = 0.5 * (mMu[g](i, j, k0 - 1) * acm + mMu[g](i, j, k0) * ac +
                          muplam * met(3, i, j, k0 - 1) * met(3, i, j, k0 - 1) +
                          mupla * met(3, i, j, k0) * met(3, i, j, k0));
      amat(2, 3) = 0.5 * (muplam * met(4, i, j, k0 - 1) * met(3, i, j, k0 - 1) +
                          mupla * met(4, i, j, k0) * met(3, i, j, k0));

      b[2] = 2 * rhs3 - 2 * bforcing[2 + 3 * qq] +
             0.5 *
                 (muplap * met(2, i, j, k0 + 1) * met(4, i, j, k0 + 1) +
                  mupla * met(2, i, j, k0) * met(4, i, j, k0)) *
                 (u(1, i, j, k0 + 1) - u(1, i, j, k0)) +
             0.5 *
                 (muplap * met(3, i, j, k0 + 1) * met(4, i, j, k0 + 1) +
                  mupla * met(3, i, j, k0) * met(4, i, j, k0)) *
                 (u(2, i, j, k0 + 1) - u(2, i, j, k0)) +
             0.5 *
                 (mMu[g](i, j, k0 + 1) * acp + mMu[g](i, j, k0) * ac +
                  muplap * met(4, i, j, k0 + 1) * met(4, i, j, k0 + 1) +
                  mupla * met(4, i, j, k0) * met(4, i, j, k0)) *
                 (u(3, i, j, k0 + 1) - u(3, i, j, k0));

      amat(3, 1) = 0.5 * (muplam * met(2, i, j, k0 - 1) * met(4, i, j, k0 - 1) +
                          mupla * met(2, i, j, k0) * met(4, i, j, k0));
      amat(3, 2) = 0.5 * (muplam * met(3, i, j, k0 - 1) * met(4, i, j, k0 - 1) +
                          mupla * met(3, i, j, k0) * met(4, i, j, k0));
      amat(3, 3) = 0.5 * (mMu[g](i, j, k0 - 1) * acm + mMu[g](i, j, k0) * ac +
                          muplam * met(4, i, j, k0 - 1) * met(4, i, j, k0 - 1) +
                          mupla * met(4, i, j, k0) * met(4, i, j, k0));
#undef amat
      // solve linear 3x3 system:
      int ipiv[3], info, three = 3, one = 1;
      F77_FUNC(dgesv, DGESV)
      (&three, &one, amat_, &three, ipiv, b, &three, &info);
      VERIFY2(info == 0,
              "ERROR: info = "
                  << info << " from DGESV in bcsurf_curvilinear_2nd_order\n");
      //	 if( info != 0 )
      //	    cout << "ERROR: info = " << info << " from DGESV in
      //bcsurf_curvilinear_2nd_order " << endl;

      u(1, i, j, k0 - 1) = u(1, i, j, k0) + b[0];
      u(2, i, j, k0 - 1) = u(2, i, j, k0) + b[1];
      u(3, i, j, k0 - 1) = u(3, i, j, k0) + b[2];
    }
}

//-----------------------------------------------------------------------
template <int iu, int il, int ju, int jl, int ku, int kl>
void evalLuCurv(int ib, int ie, int jb, int je, int kb, int ke, Sarray& u,
                Sarray& lu, float_sw4* a_mu, float_sw4* a_la, Sarray& met,
                Sarray& jac, int ilb, int ile, int jlb, int jle, int klb,
                int kle) {
#define mu(i, j, k) a_mu[i - ib + ni * (j - jb) + nij * (k - kb)]
#define la(i, j, k) a_la[i - ib + ni * (j - jb) + nij * (k - kb)]
  //#define u(c,i,j,k)  a_u[i-ib+ni*(j-jb)+nij*(k-kb)+(c-1)*nijk]
  //#define lu(c,i,j,k) a_lu[i-ilb+nli*(j-jlb)+nlij*(k-klb)+(c-1)*nlijk]
  const size_t ni = ie - ib + 1;
  const size_t nij = ni * (je - jb + 1);
  //   const size_t nijk=nij*(ke-kb+1);
  //   const size_t nli=ile-ilb+1;
  //   const size_t nlij=nli*(jle-jlb+1);
  //   const size_t nlijk=nlij*(kle-klb+1);
  //   const float_sw4 ih2 = 1/(h*h);
  const float_sw4 half = 0.5;
  const float_sw4 fourth = 0.25;
  // differences in mixed terms are computed as (u(i+iu)-u(i-il))/(iu+il), so
  // that iu=1,il=0 gives forward diff, iu=0,il=1 gives backward diff, iu=1,il=1
  // gives centered diff.

  const int ok = ku == kl ? 1 : 0;
  const int oj = ju == jl ? 1 : 0;
  const int oi = iu == il ? 1 : 0;
  for (int k = klb + ok; k <= kle - ok; k++)
    for (int j = jlb + oj; j <= jle - oj; j++)
      for (int i = ilb + oi; i <= ile - oi; i++)
      //   for( int k=klb+1 ; k <= kle-1; k++ )
      //      for( int j=jlb ; j <= jle; j++ )
      //	 for( int i=ilb+1 ; i <= ile-1; i++ )
      {
        float_sw4 r1 = 0, r2 = 0, r3 = 0;
        float_sw4 ijac = 1.0 / jac(i, j, k);
        // U-equation D+D-
        float_sw4 mup = 0.5 * ((2 * mu(i + 1, j, k) + la(i + 1, j, k)) *
                                   met(1, i + 1, j, k) * met(1, i + 1, j, k) +
                               (2 * mu(i, j, k) + la(i, j, k)) *
                                   met(1, i, j, k) * met(1, i, j, k));

        float_sw4 mum = 0.5 * ((2 * mu(i, j, k) + la(i, j, k)) *
                                   met(1, i, j, k) * met(1, i, j, k) +
                               (2 * mu(i - 1, j, k) + la(i - 1, j, k)) *
                                   met(1, i - 1, j, k) * met(1, i - 1, j, k));

        r1 +=
            mup * (u(1, i + 1, j, k) - u(1, i, j, k)) -
            mum * (u(1, i, j, k) - u(1, i - 1, j, k)) +
            0.5 *
                ((mu(i, j + 1, k) * met(1, i, j + 1, k) * met(1, i, j + 1, k) +
                  mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k)) *
                     (u(1, i, j + 1, k) - u(1, i, j, k)) -
                 (mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k) +
                  mu(i, j - 1, k) * met(1, i, j - 1, k) * met(1, i, j - 1, k)) *
                     (u(1, i, j, k) - u(1, i, j - 1, k)));

        mup = 0.5 *
              ((2 * mu(i, j, k + 1) + la(i, j, k + 1)) * met(2, i, j, k + 1) *
                   met(2, i, j, k + 1) +
               mu(i, j, k + 1) * (met(3, i, j, k + 1) * met(3, i, j, k + 1) +
                                  met(4, i, j, k + 1) * met(4, i, j, k + 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) *
                   met(2, i, j, k) +
               mu(i, j, k) * (met(3, i, j, k) * met(3, i, j, k) +
                              met(4, i, j, k) * met(4, i, j, k)));
        mum = 0.5 *
              ((2 * mu(i, j, k - 1) + la(i, j, k - 1)) * met(2, i, j, k - 1) *
                   met(2, i, j, k - 1) +
               mu(i, j, k - 1) * (met(3, i, j, k - 1) * met(3, i, j, k - 1) +
                                  met(4, i, j, k - 1) * met(4, i, j, k - 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) *
                   met(2, i, j, k) +
               mu(i, j, k) * (met(3, i, j, k) * met(3, i, j, k) +
                              met(4, i, j, k) * met(4, i, j, k)));

        r1 += mup * (u(1, i, j, k + 1) - u(1, i, j, k)) -
              mum * (u(1, i, j, k) - u(1, i, j, k - 1));

        r1 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(2, i, j, k + 1) *
                  met(3, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(3, i, j, k)) *
                 (u(2, i, j, k + 1) - u(2, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(3, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(2, i, j, k - 1) *
                  met(3, i, j, k - 1)) *
                 (u(2, i, j, k) - u(2, i, j, k - 1)));

        r1 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(2, i, j, k + 1) *
                  met(4, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(4, i, j, k)) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(4, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(2, i, j, k - 1) *
                  met(4, i, j, k - 1)) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));

        // V-equation D+D-
        mup = 0.5 * ((2 * mu(i, j + 1, k) + la(i, j + 1, k)) *
                         met(1, i, j + 1, k) * met(1, i, j + 1, k) +
                     (2 * mu(i, j, k) + la(i, j, k)) * met(1, i, j, k) *
                         met(1, i, j, k));

        mum = 0.5 * ((2 * mu(i, j, k) + la(i, j, k)) * met(1, i, j, k) *
                         met(1, i, j, k) +
                     (2 * mu(i, j - 1, k) + la(i, j - 1, k)) *
                         met(1, i, j - 1, k) * met(1, i, j - 1, k));

        r2 +=
            mup * (u(2, i, j + 1, k) - u(2, i, j, k)) -
            mum * (u(2, i, j, k) - u(2, i, j - 1, k)) +
            0.5 *
                ((mu(i + 1, j, k) * met(1, i + 1, j, k) * met(1, i + 1, j, k) +
                  mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k)) *
                     (u(2, i + 1, j, k) - u(2, i, j, k)) -
                 (mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k) +
                  mu(i - 1, j, k) * met(1, i - 1, j, k) * met(1, i - 1, j, k)) *
                     (u(2, i, j, k) - u(2, i - 1, j, k)));

        r2 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(2, i, j, k + 1) *
                  met(3, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(3, i, j, k)) *
                 (u(1, i, j, k + 1) - u(1, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(3, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(2, i, j, k - 1) *
                  met(3, i, j, k - 1)) *
                 (u(1, i, j, k) - u(1, i, j, k - 1)));

        mup = 0.5 *
              ((2 * mu(i, j, k + 1) + la(i, j, k + 1)) * met(3, i, j, k + 1) *
                   met(3, i, j, k + 1) +
               mu(i, j, k + 1) * (met(2, i, j, k + 1) * met(2, i, j, k + 1) +
                                  met(4, i, j, k + 1) * met(4, i, j, k + 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) *
                   met(3, i, j, k) +
               mu(i, j, k) * (met(2, i, j, k) * met(2, i, j, k) +
                              met(4, i, j, k) * met(4, i, j, k)));
        mum = 0.5 *
              ((2 * mu(i, j, k - 1) + la(i, j, k - 1)) * met(3, i, j, k - 1) *
                   met(3, i, j, k - 1) +
               mu(i, j, k - 1) * (met(2, i, j, k - 1) * met(2, i, j, k - 1) +
                                  met(4, i, j, k - 1) * met(4, i, j, k - 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) *
                   met(3, i, j, k) +
               mu(i, j, k) * (met(2, i, j, k) * met(2, i, j, k) +
                              met(4, i, j, k) * met(4, i, j, k)));

        r2 += mup * (u(2, i, j, k + 1) - u(2, i, j, k)) -
              mum * (u(2, i, j, k) - u(2, i, j, k - 1));

        r2 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(3, i, j, k + 1) *
                  met(4, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) * met(4, i, j, k)) *
                 (u(3, i, j, k + 1) - u(3, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) * met(4, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(3, i, j, k - 1) *
                  met(4, i, j, k - 1)) *
                 (u(3, i, j, k) - u(3, i, j, k - 1)));

        // W-equation D+D-
        r3 += 0.5 *
              ((mu(i, j + 1, k) * met(1, i, j + 1, k) * met(1, i, j + 1, k) +
                mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k)) *
                   (u(3, i, j + 1, k) - u(3, i, j, k)) -
               (mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k) +
                mu(i, j - 1, k) * met(1, i, j - 1, k) * met(1, i, j - 1, k)) *
                   (u(3, i, j, k) - u(3, i, j - 1, k)));
        r3 += 0.5 *
              ((mu(i + 1, j, k) * met(1, i + 1, j, k) * met(1, i + 1, j, k) +
                mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k)) *
                   (u(3, i + 1, j, k) - u(3, i, j, k)) -
               (mu(i, j, k) * met(1, i, j, k) * met(1, i, j, k) +
                mu(i - 1, j, k) * met(1, i - 1, j, k) * met(1, i - 1, j, k)) *
                   (u(3, i, j, k) - u(3, i - 1, j, k)));

        r3 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(2, i, j, k + 1) *
                  met(4, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(4, i, j, k)) *
                 (u(1, i, j, k + 1) - u(1, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(2, i, j, k) * met(4, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(2, i, j, k - 1) *
                  met(4, i, j, k - 1)) *
                 (u(1, i, j, k) - u(1, i, j, k - 1)));

        r3 +=
            0.5 *
            (((mu(i, j, k + 1) + la(i, j, k + 1)) * met(3, i, j, k + 1) *
                  met(4, i, j, k + 1) +
              (mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) * met(4, i, j, k)) *
                 (u(2, i, j, k + 1) - u(2, i, j, k)) -
             ((mu(i, j, k) + la(i, j, k)) * met(3, i, j, k) * met(4, i, j, k) +
              (mu(i, j, k - 1) + la(i, j, k - 1)) * met(3, i, j, k - 1) *
                  met(4, i, j, k - 1)) *
                 (u(2, i, j, k) - u(2, i, j, k - 1)));

        mup = 0.5 *
              ((2 * mu(i, j, k + 1) + la(i, j, k + 1)) * met(4, i, j, k + 1) *
                   met(4, i, j, k + 1) +
               mu(i, j, k + 1) * (met(2, i, j, k + 1) * met(2, i, j, k + 1) +
                                  met(3, i, j, k + 1) * met(3, i, j, k + 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(4, i, j, k) *
                   met(4, i, j, k) +
               mu(i, j, k) * (met(2, i, j, k) * met(2, i, j, k) +
                              met(3, i, j, k) * met(3, i, j, k)));
        mum = 0.5 *
              ((2 * mu(i, j, k - 1) + la(i, j, k - 1)) * met(4, i, j, k - 1) *
                   met(4, i, j, k - 1) +
               mu(i, j, k - 1) * (met(2, i, j, k - 1) * met(2, i, j, k - 1) +
                                  met(3, i, j, k - 1) * met(3, i, j, k - 1)) +
               (2 * mu(i, j, k) + la(i, j, k)) * met(4, i, j, k) *
                   met(4, i, j, k) +
               mu(i, j, k) * (met(2, i, j, k) * met(2, i, j, k) +
                              met(3, i, j, k) * met(3, i, j, k)));

        r3 += mup * (u(3, i, j, k + 1) - u(3, i, j, k)) -
              mum * (u(3, i, j, k) - u(3, i, j, k - 1));

        // U-equation, mixed derivatives
        // pq-derivatives
        r1 += (mu(i, j + ju, k) * met(1, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(2, i + iu, j + ju, k) - u(2, i - il, j + ju, k)) -
               mu(i, j - jl, k) * met(1, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(2, i + iu, j - jl, k) - u(2, i - il, j - jl, k))) /
              ((iu + il) * (ju + jl));

        // qp-derivatives
        r1 += (la(i + iu, j, k) * met(1, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(2, i + iu, j + ju, k) - u(2, i + iu, j - jl, k)) -
               la(i - il, j, k) * met(1, i - il, j, k) * met(1, i - il, j, k) *
                   (u(2, i - il, j + ju, k) - u(2, i - il, j - jl, k))) /
              ((iu + il) * (ju + jl));

        // pr-derivatives
        r1 += ((2 * mu(i, j, k + ku) + la(i, j, k + ku)) *
                   met(2, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(1, i + iu, j, k + ku) - u(1, i - il, j, k + ku)) -
               (2 * mu(i, j, k - kl) + la(i, j, k - kl)) *
                   met(2, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(1, i + iu, j, k - kl) - u(1, i - il, j, k - kl)) +
               mu(i, j, k + ku) * met(3, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(2, i + iu, j, k + ku) - u(2, i - il, j, k + ku)) -
               mu(i, j, k - kl) * met(3, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(2, i + iu, j, k - kl) - u(2, i - il, j, k - kl)) +
               mu(i, j, k + ku) * met(4, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(3, i + iu, j, k + ku) - u(3, i - il, j, k + ku)) -
               mu(i, j, k - kl) * met(4, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(3, i + iu, j, k - kl) - u(3, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));

        // rp derivatives
        r1 += ((2 * mu(i + iu, j, k) + la(i + iu, j, k)) *
                   met(2, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(1, i + iu, j, k + ku) - u(1, i + iu, j, k - kl)) -
               (2 * mu(i - il, j, k) + la(i - il, j, k)) *
                   met(2, i - il, j, k) * met(1, i - il, j, k) *
                   (u(1, i - il, j, k + ku) - u(1, i - il, j, k - kl)) +
               la(i + iu, j, k) * met(3, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(2, i + iu, j, k + ku) - u(2, i + iu, j, k - kl)) -
               la(i - il, j, k) * met(3, i - il, j, k) * met(1, i - il, j, k) *
                   (u(2, i - il, j, k + ku) - u(2, i - il, j, k - kl)) +
               la(i + iu, j, k) * met(4, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(3, i + iu, j, k + ku) - u(3, i + iu, j, k - kl)) -
               la(i - il, j, k) * met(4, i - il, j, k) * met(1, i - il, j, k) *
                   (u(3, i - il, j, k + ku) - u(3, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));

        // qr derivatives
        r1 += (mu(i, j, k + ku) * met(3, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(1, i, j + ju, k + ku) - u(1, i, j - jl, k + ku)) +
               la(i, j, k + ku) * met(2, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j - jl, k + ku)) -
               mu(i, j, k - kl) * met(3, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(1, i, j + ju, k - kl) - u(1, i, j - jl, k - kl)) -
               la(i, j, k - kl) * met(2, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(2, i, j + ju, k - kl) - u(2, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));

        // rq derivatives
        r1 += (mu(i, j + ju, k) * met(3, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(1, i, j + ju, k + ku) - u(1, i, j + ju, k - kl)) +
               mu(i, j + ju, k) * met(2, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j + ju, k - kl)) -
               mu(i, j - jl, k) * met(3, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(1, i, j - jl, k + ku) - u(1, i, j - jl, k - kl)) -
               mu(i, j - jl, k) * met(2, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(2, i, j - jl, k + ku) - u(2, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));

        lu(1, i, j, k) = r1 * ijac;

        // V-equation, mixed derivatives

        // pq-derivatives
        r2 += (la(i, j + ju, k) * met(1, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(1, i + iu, j + ju, k) - u(1, i - il, j + ju, k)) -
               la(i, j - jl, k) * met(1, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(1, i + iu, j - jl, k) - u(1, i - il, j - jl, k))) /
              ((iu + il) * (ju + jl));

        // qp-derivatives
        r2 += (mu(i + iu, j, k) * met(1, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(1, i + iu, j + ju, k) - u(1, i + iu, j - jl, k)) -
               mu(i - il, j, k) * met(1, i - il, j, k) * met(1, i - il, j, k) *
                   (u(1, i - il, j + ju, k) - u(1, i - il, j - jl, k))) /
              ((iu + il) * (ju + jl));

        // pr-derivatives
        r2 += (la(i, j, k + ku) * met(3, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(1, i + iu, j, k + ku) - u(1, i - il, j, k + ku)) +
               mu(i, j, k + ku) * met(2, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(2, i + iu, j, k + ku) - u(2, i - il, j, k + ku)) -
               la(i, j, k - kl) * met(3, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(1, i + iu, j, k - kl) - u(1, i - il, j, k - kl)) -
               mu(i, j, k - kl) * met(2, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(2, i + iu, j, k - kl) - u(2, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));
        // rp derivatives
        r2 += (mu(i + iu, j, k) * met(3, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(1, i + iu, j, k + ku) - u(1, i + iu, j, k - kl)) +
               mu(i + iu, j, k) * met(2, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(2, i + iu, j, k + ku) - u(2, i + iu, j, k - kl)) -
               mu(i - il, j, k) * met(3, i - il, j, k) * met(1, i - il, j, k) *
                   (u(1, i - il, j, k + ku) - u(1, i - il, j, k - kl)) -
               mu(i - il, j, k) * met(2, i - il, j, k) * met(1, i - il, j, k) *
                   (u(2, i - il, j, k + ku) - u(2, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));

        // qr derivatives
        r2 += (mu(i, j, k + ku) * met(2, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(1, i, j + ju, k + ku) - u(1, i, j - jl, k + ku)) +
               (2 * mu(i, j, k + ku) + la(i, j, k + ku)) *
                   met(3, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j - jl, k + ku)) +
               mu(i, j, k + ku) * met(4, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(3, i, j + ju, k + ku) - u(3, i, j - jl, k + ku)) -
               mu(i, j, k - kl) * met(2, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(1, i, j + ju, k - kl) - u(1, i, j - jl, k - kl)) -
               (2 * mu(i, j, k - kl) + la(i, j, k - kl)) *
                   met(3, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(2, i, j + ju, k - kl) - u(2, i, j - jl, k - kl)) -
               mu(i, j, k - kl) * met(4, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(3, i, j + ju, k - kl) - u(3, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));

        // rq derivatives
        r2 += (la(i, j + ju, k) * met(2, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(1, i, j + ju, k + ku) - u(1, i, j + ju, k - kl)) +
               (2 * mu(i, j + ju, k) + la(i, j + ju, k)) *
                   met(3, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j + ju, k - kl)) +
               la(i, j + ju, k) * met(4, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(3, i, j + ju, k + ku) - u(3, i, j + ju, k - kl)) -
               la(i, j - jl, k) * met(2, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(1, i, j - jl, k + ku) - u(1, i, j - jl, k - kl)) -
               (2 * mu(i, j - jl, k) + la(i, j - jl, k)) *
                   met(3, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(2, i, j - jl, k + ku) - u(2, i, j - jl, k - kl)) -
               la(i, j - jl, k) * met(4, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(3, i, j - jl, k + ku) - u(3, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));

        lu(2, i, j, k) = r2 * ijac;

        // W-equation, mixed derivatives

        // pr-derivatives
        r3 += (la(i, j, k + ku) * met(4, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(1, i + iu, j, k + ku) - u(1, i - il, j, k + ku)) +
               mu(i, j, k + ku) * met(2, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(3, i + iu, j, k + ku) - u(3, i - il, j, k + ku)) -
               la(i, j, k - kl) * met(4, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(1, i + iu, j, k - kl) - u(1, i - il, j, k - kl)) -
               mu(i, j, k - kl) * met(2, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(3, i + iu, j, k - kl) - u(3, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));
        // rp derivatives
        r3 += (mu(i + iu, j, k) * met(4, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(1, i + iu, j, k + ku) - u(1, i + iu, j, k - kl)) +
               mu(i + iu, j, k) * met(2, i + iu, j, k) * met(1, i + iu, j, k) *
                   (u(3, i + iu, j, k + ku) - u(3, i + iu, j, k - kl)) -
               mu(i - il, j, k) * met(4, i - il, j, k) * met(1, i - il, j, k) *
                   (u(1, i - il, j, k + ku) - u(1, i - il, j, k - kl)) -
               mu(i - il, j, k) * met(2, i - il, j, k) * met(1, i - il, j, k) *
                   (u(3, i - il, j, k + ku) - u(3, i - il, j, k - kl))) /
              ((iu + il) * (ku + kl));
        // qr derivatives
        r3 += (mu(i, j, k + ku) * met(3, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(3, i, j + ju, k + ku) - u(3, i, j - jl, k + ku)) +
               la(i, j, k + ku) * met(4, i, j, k + ku) * met(1, i, j, k + ku) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j - jl, k + ku)) -
               mu(i, j, k - kl) * met(3, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(3, i, j + ju, k - kl) - u(3, i, j - jl, k - kl)) -
               la(i, j, k - kl) * met(4, i, j, k - kl) * met(1, i, j, k - kl) *
                   (u(2, i, j + ju, k - kl) - u(2, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));
        // rq derivatives
        r3 += (mu(i, j + ju, k) * met(3, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(3, i, j + ju, k + ku) - u(3, i, j + ju, k - kl)) +
               mu(i, j + ju, k) * met(4, i, j + ju, k) * met(1, i, j + ju, k) *
                   (u(2, i, j + ju, k + ku) - u(2, i, j + ju, k - kl)) -
               mu(i, j - jl, k) * met(3, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(3, i, j - jl, k + ku) - u(3, i, j - jl, k - kl)) -
               mu(i, j - jl, k) * met(4, i, j - jl, k) * met(1, i, j - jl, k) *
                   (u(2, i, j - jl, k + ku) - u(2, i, j - jl, k - kl))) /
              ((ju + jl) * (ku + kl));

        lu(3, i, j, k) = r3 * ijac;
      }
#undef mu
#undef la
}
