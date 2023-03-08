#include "caliper.h"
#include <fcntl.h>
#include <unistd.h>
#include <cmath>

#include <iostream>

#include "MParGridFile.h"

//-----------------------------------------------------------------------
MParGridFile::MParGridFile(std::string fname) { read_mparcartfile(fname); }

//-----------------------------------------------------------------------
MParGridFile::MParGridFile(float_sw4* xms, int vars, int nx, int ny, int nz,
                           float_sw4 hx, float_sw4 hy, float_sw4 hz,
                           float_sw4 xmin, float_sw4 ymin, float_sw4 zmin,
                           int update) {

  SW4_MARK_FUNCTION;
  
  m_update = update;
  m_vars = vars;
  m_nxr = nx;
  m_nyr = ny;
  m_nzr = nz;
  m_hxr = hx;
  m_hyr = hy;
  m_hzr = hz;
  m_xrmin = xmin;
  m_yrmin = ymin;
  m_zrmin = zmin;
  int nc = 3;
  if (m_vars == 3)
    nc = 1;
  else if (m_vars == 4)
    nc = 2;
  m_xms = new double[nc * m_nxr * m_nyr * m_nzr];
  for (int i = 0; i < nc * m_nxr * m_nyr * m_nzr; i++) m_xms[i] = xms[i];
}

//-----------------------------------------------------------------------
MParGridFile::~MParGridFile() { delete[] m_xms; }

//-----------------------------------------------------------------------
void MParGridFile::read_mparcartfile(std::string fname) {

  SW4_MARK_FUNCTION;
  
  int fd = open(fname.c_str(), O_RDONLY);
  int mag;
  size_t nr = read(fd, &mag, sizeof(int));
  nr = read(fd, &m_update, sizeof(int));
  nr = read(fd, &m_vars, sizeof(int));
  nr = read(fd, &m_nxr, sizeof(int));
  nr = read(fd, &m_nyr, sizeof(int));
  nr = read(fd, &m_nzr, sizeof(int));
  nr = read(fd, &m_hxr, sizeof(double));
  nr = read(fd, &m_hyr, sizeof(double));
  nr = read(fd, &m_hzr, sizeof(double));
  nr = read(fd, &m_xrmin, sizeof(double));
  nr = read(fd, &m_yrmin, sizeof(double));
  nr = read(fd, &m_zrmin, sizeof(double));
  int nc = 3;
  if (m_vars == 3)
    nc = 1;
  else if (m_vars == 4)
    nc = 2;
  m_xms = new double[nc * m_nxr * m_nyr * m_nzr];
  nr = read(fd, m_xms, nc * m_nxr * m_nyr * m_nzr * sizeof(double));
  close(fd);
}

//-----------------------------------------------------------------------
void MParGridFile::write_mparcartfile(std::string fname) {

  SW4_MARK_FUNCTION;
  
  int fd = open(fname.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0660);
  int mag = 1;
  size_t nr = write(fd, &mag, sizeof(int));
  nr = write(fd, &m_update, sizeof(int));
  nr = write(fd, &m_vars, sizeof(int));
  nr = write(fd, &m_nxr, sizeof(int));
  nr = write(fd, &m_nyr, sizeof(int));
  nr = write(fd, &m_nzr, sizeof(int));
  nr = write(fd, &m_hxr, sizeof(double));
  nr = write(fd, &m_hyr, sizeof(double));
  nr = write(fd, &m_hzr, sizeof(double));
  nr = write(fd, &m_xrmin, sizeof(double));
  nr = write(fd, &m_yrmin, sizeof(double));
  nr = write(fd, &m_zrmin, sizeof(double));
  int nc = 3;
  if (m_vars == 3)
    nc = 1;
  else if (m_vars == 4)
    nc = 2;
  nr = write(fd, m_xms, nc * m_nxr * m_nyr * m_nzr * sizeof(double));
  close(fd);
}

//-----------------------------------------------------------------------
void MParGridFile::interpolate_to_other(float_sw4* xms, int vars, int nx,
                                        int ny, int nz, float_sw4 hx,
                                        float_sw4 hy, float_sw4 hz,
                                        float_sw4 xmin, float_sw4 ymin,
                                        float_sw4 zmin) {

  SW4_MARK_FUNCTION;

  
  int nc;
  if (vars == 1 || vars == 2)
    nc = 3;
  else if (vars == 3)
    nc = 1;
  else if (vars == 4)
    nc = 2;

  bool fail = false;
  float_sw4* xmsrc;
  if (vars == m_vars)
    xmsrc = m_xms;
  else {
    xmsrc = new float_sw4[nc * m_nxr * m_nyr * m_nzr];
    for (int ind = 0; ind < m_nxr * m_nyr * m_nzr; ind++) {
      if (vars == 1) {
        if (m_vars == 2) {
          xmsrc[3 * ind] = m_xms[3 * ind];  // rho=rho
          xmsrc[3 * ind + 1] = m_xms[3 * ind] * m_xms[3 * ind + 1] *
                               m_xms[3 * ind + 1];  // mu=rho*cs*cs
          xmsrc[3 * ind + 1] =
              m_xms[3 * ind] * m_xms[3 * ind + 2] * m_xms[3 * ind + 2] -
              2 * xmsrc[3 * ind + 1];  // la=rho*cp*cp-2*mu
        } else
          fail = true;
      } else if (vars == 2) {
        if (m_vars == 1) {
          xmsrc[3 * ind] = m_xms[3 * ind];  // rho=rho
          xmsrc[3 * ind + 1] =
              sqrt(m_xms[3 * ind + 1] / m_xms[3 * ind]);  // cs = sqrt(mu/rho)
          xmsrc[3 * ind + 2] =
              sqrt((2 * m_xms[3 * ind + 1] + m_xms[3 * ind + 2]) /
                   m_xms[3 * ind]);  // cp=sqrt((2*mu+la)/rho)
        }

        else
          fail = true;
      } else if (vars == 3) {
        if (m_vars == 1)
          xmsrc[ind] = sqrt((2 * m_xms[3 * ind + 1] + m_xms[3 * ind + 2]) /
                            m_xms[3 * ind]);
        else if (m_vars == 2)
          xmsrc[ind] = m_xms[3 * ind + 2];
        else if (m_vars == 3)
          xmsrc[ind] = m_xms[ind];
        else if (m_vars == 4)
          xmsrc[ind] = m_xms[2 * ind + 1];
      } else if (vars == 4) {
        if (m_vars == 1) {
          if (std::isnan(m_xms[3 * ind + 1]) || std::isnan(m_xms[3 * ind]) ||
              m_xms[3 * ind + 1] < 0 || m_xms[3 * ind] < 0)
            std::cout << "error for xmsrc (rho,mu,lambda)= " << m_xms[3 * ind]
                      << " " << m_xms[3 * ind + 1] << " " << m_xms[3 * ind + 2]
                      << std::endl;
          xmsrc[2 * ind] =
              sqrt(m_xms[3 * ind + 1] / m_xms[3 * ind]);  // cs = sqrt(mu/rho)
          xmsrc[2 * ind + 1] =
              sqrt((2 * m_xms[3 * ind + 1] + m_xms[3 * ind + 2]) /
                   m_xms[3 * ind]);  // cp=sqrt((2*mu+la)/rho)
        } else if (m_vars == 2) {
          xmsrc[2 * ind] = m_xms[3 * ind + 1];
          xmsrc[2 * ind + 1] = m_xms[3 * ind + 2];
        } else
          fail = true;
      }
    }
  }
  if (fail) {
    std::cout << "ERROR: MParGridFile, attempt to extract too many components"
              << std::endl;
    if (vars != m_vars) delete[] xmsrc;
    return;
  }
  int nxyr = m_nxr * m_nyr;
  size_t indpar = 0;
  //   std::cout << " nx,ny,nz = " << nx << " " << ny << " " <<  nz << " m_ " <<
  //     m_nxr << " " << m_nyr << " " << m_nzr << std::endl;
  for (int k = 1; k <= nz; k++) {
    double z = zmin + (k - 1) * hz;
    int kr = static_cast<int>((z - m_zrmin) / m_hzr + 1);
    if (kr < 1) kr = 1;
    if (kr > m_nzr - 1) kr = m_nzr - 1;
    double zr = m_zrmin + (kr - 1) * m_hzr;
    float_sw4 wghz = (z - zr) / m_hzr;
    for (int j = 1; j <= ny; j++) {
      double y = ymin + (j - 1) * hy;
      int jr = static_cast<int>((y - m_yrmin) / m_hyr + 1);
      if (jr < 1) jr = 1;
      if (jr > m_nyr - 1) jr = m_nyr - 1;
      double yr = m_yrmin + (jr - 1) * m_hyr;
      float_sw4 wghy = (y - yr) / m_hyr;
      for (int i = 1; i <= nx; i++) {
        double x = xmin + (i - 1) * hx;
        int ir = static_cast<int>((x - m_xrmin) / m_hxr + 1);
        if (ir < 1) ir = 1;
        if (ir > m_nxr - 1) ir = m_nxr - 1;
        double xr = m_xrmin + (ir - 1) * m_hxr;
        float_sw4 wghx = (x - xr) / m_hxr;
        size_t ind = ir - 1 + m_nxr * (jr - 1) + nxyr * (kr - 1);
        //	    std::cout << "ind= " << ind << " ind+ " << ind+1+m_nxr+nxyr
        //<< " indpar = " << indpar << std::endl; 	    std::cout << "ind " << ind <<
        //" wgh " << wghx << " " << wghy << " " << wghz << std::endl;
        for (int c = 0; c < nc; c++) {
          xms[c + nc * indpar] =
              (1 - wghz) * ((1 - wghy) * ((1 - wghx) * xmsrc[c + nc * ind] +
                                          wghx * xmsrc[c + nc * (ind + 1)]) +
                            wghy * ((1 - wghx) * xmsrc[c + nc * (ind + m_nxr)] +
                                    wghx * xmsrc[c + nc * (ind + 1 + m_nxr)])) +
              wghz *
                  ((1 - wghy) * ((1 - wghx) * xmsrc[c + nc * (ind + nxyr)] +
                                 wghx * xmsrc[c + nc * (ind + 1 + nxyr)]) +
                   wghy * ((1 - wghx) * xmsrc[c + nc * (ind + m_nxr + nxyr)] +
                           wghx * xmsrc[c + nc * (ind + 1 + m_nxr + nxyr)]));
          if (std::isnan(xms[c + nc * indpar]))
            std::cout << "nan in interp " << c << " " << indpar
                      << " ind= " << ind << " wgh= " << wghx << " " << wghy
                      << " " << wghz << " " << xmsrc[c + nc * ind] << " "
                      << xmsrc[c + nc * (ind + 1)] << " "
                      << xmsrc[c + nc * (ind + m_nxr)] << " "
                      << xmsrc[c + nc * (ind + 1 + m_nxr)] << std::endl;
        }
        indpar++;
      }
    }
  }
  if (vars != m_vars) delete[] xmsrc;
}
