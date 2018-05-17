#ifndef __OPTIMIZEDCUDA_RHS__
#define __OPTIMIZEDCUDA_RHS__
#include <cuda_runtime.h>
#define BX 16
#define BY 16
#define u(c,idx) (c_order ? a_u[idx+(c)*nijk] : a_u[3*(idx)+c])
#define up(c,idx) (c_order ? a_up[idx+(c)*nijk] : a_up[3*(idx)+c])


template <bool c_order, bool pred>
__launch_bounds__(BX*BY)
  __global__ void rhs4_v2 (int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			   int ni, int nj, int nk,
			   float_sw4* __restrict__ a_up,
			   const float_sw4* __restrict__ a_u,
			   const float_sw4* __restrict__ a_mu,
			   const float_sw4* __restrict__ a_lambda,
			   const float_sw4* __restrict__ a_strx,
			   const float_sw4* __restrict__ a_stry,
			   const float_sw4* __restrict__ a_strz,
			   const float_sw4 h,const int a1fac,const float_sw4 cof) {
  int index;

  const float_sw4 i6   = 1.0/6;
  const float_sw4 i144 = 1.0/144;
  const float_sw4 tf   = 0.75;
  
  float_sw4 strx_im2, strx_im1, strx_i, strx_ip1, strx_ip2;
  float_sw4 stry_jm2, stry_jm1, stry_j, stry_jp1, stry_jp2;
  float_sw4 strz_km2, strz_km1, strz_k, strz_kp1, strz_kp2;

  float_sw4 la_km2, la_km1, la_ijk, la_kp1, la_kp2;
  float_sw4 la_im2, la_im1, la_ip1, la_ip2;
  float_sw4 la_jm2, la_jm1, la_jp1, la_jp2;

  float_sw4 mu_km2, mu_km1, mu_ijk, mu_kp1, mu_kp2;
  float_sw4 mu_im2, mu_im1, mu_ip1, mu_ip2;
  float_sw4 mu_jm2, mu_jm1, mu_jp1, mu_jp2;

  // float_sw4 fo0, fo1, fo2, rho;
  // float_sw4 up0, up1, up2;
  // float_sw4 um0, um1, um2;

  // Prefetch variables
  float_sw4 a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2;
  
  float_sw4 r1, r2, r3;
  int active=0, loader=0, loady2=0, loadx2=0;

  __shared__ float_sw4 shu[3][5][BY+4][BX+4];

  // i starts at 0, while j starts at jfirst
  const int i = threadIdx.x + blockIdx.x * BX;
  const int j = jfirst + threadIdx.y + blockIdx.y * BY;

  const int ti = threadIdx.x + 2;
  const int tj = threadIdx.y + 2;
  //const int tk = 2;

  const int nij = ni * nj;
  const int nijk = nij * nk;

  int kthm0=3, kthm1=2, kthm2=1, kthm3=0, kthm4=4, kthtmp;

  // Active threads doing the computation at (i,j,k)
  if (i >= ifirst && i <= ilast && j <= jlast)
    active = 1;

  // Threads loading data with halos (4 regions)
  if (j-2 <= jlast + 2) {
    if (i-2 >= 0 && i-2 < ni)
      loader = 1;
    if (threadIdx.x < 4 && i+BX-2 < ni)
      loadx2 = 1;
    if (threadIdx.y < 4 && j+BY-2 <= jlast+2)
      loady2 = 1;
  }

  // Loading index for data with halos, starting at (i-2,j-2,kfirst-2)
  index = (kfirst - 2) * nij + (j - 2) * ni + i - 2;

  // Load the first 4 plans of U in the shared memory
  // and prefetch the 5th plan in a,b,c,d
  if (loader) {
    int idx = index;
    shu[0][0][threadIdx.y][threadIdx.x] = u(0,idx);
    shu[1][0][threadIdx.y][threadIdx.x] = u(1,idx);
    shu[2][0][threadIdx.y][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][1][threadIdx.y][threadIdx.x] = u(0,idx);
    shu[1][1][threadIdx.y][threadIdx.x] = u(1,idx);
    shu[2][1][threadIdx.y][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][2][threadIdx.y][threadIdx.x] = u(0,idx);
    shu[1][2][threadIdx.y][threadIdx.x] = u(1,idx);
    shu[2][2][threadIdx.y][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][3][threadIdx.y][threadIdx.x] = u(0,idx);
    shu[1][3][threadIdx.y][threadIdx.x] = u(1,idx);
    shu[2][3][threadIdx.y][threadIdx.x] = u(2,idx);
    idx += nij;
    a0 = u(0,idx);
    a1 = u(1,idx);
    a2 = u(2,idx);
  }
  if (loadx2) {
    int idx = index + BX;
    shu[0][0][threadIdx.y][threadIdx.x+BX] = u(0,idx);
    shu[1][0][threadIdx.y][threadIdx.x+BX] = u(1,idx);
    shu[2][0][threadIdx.y][threadIdx.x+BX] = u(2,idx);
    idx += nij;		         
    shu[0][1][threadIdx.y][threadIdx.x+BX] = u(0,idx);
    shu[1][1][threadIdx.y][threadIdx.x+BX] = u(1,idx);
    shu[2][1][threadIdx.y][threadIdx.x+BX] = u(2,idx);
    idx += nij;		         
    shu[0][2][threadIdx.y][threadIdx.x+BX] = u(0,idx);
    shu[1][2][threadIdx.y][threadIdx.x+BX] = u(1,idx);
    shu[2][2][threadIdx.y][threadIdx.x+BX] = u(2,idx);
    idx += nij;		         
    shu[0][3][threadIdx.y][threadIdx.x+BX] = u(0,idx);
    shu[1][3][threadIdx.y][threadIdx.x+BX] = u(1,idx);
    shu[2][3][threadIdx.y][threadIdx.x+BX] = u(2,idx);
    idx += nij;
    b0 = u(0,idx);
    b1 = u(1,idx);
    b2 = u(2,idx);
  }
  if (loader && loady2) {
    int idx = index + BY * ni;
    shu[0][0][threadIdx.y+BY][threadIdx.x] = u(0,idx);
    shu[1][0][threadIdx.y+BY][threadIdx.x] = u(1,idx);
    shu[2][0][threadIdx.y+BY][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][1][threadIdx.y+BY][threadIdx.x] = u(0,idx);
    shu[1][1][threadIdx.y+BY][threadIdx.x] = u(1,idx);
    shu[2][1][threadIdx.y+BY][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][2][threadIdx.y+BY][threadIdx.x] = u(0,idx);
    shu[1][2][threadIdx.y+BY][threadIdx.x] = u(1,idx);
    shu[2][2][threadIdx.y+BY][threadIdx.x] = u(2,idx);
    idx += nij;
    shu[0][3][threadIdx.y+BY][threadIdx.x] = u(0,idx);
    shu[1][3][threadIdx.y+BY][threadIdx.x] = u(1,idx);
    shu[2][3][threadIdx.y+BY][threadIdx.x] = u(2,idx);
    idx += nij;
    c0 = u(0,idx);
    c1 = u(1,idx);
    c2 = u(2,idx);
  }
  if (loadx2 && loady2) {
    int idx = index + BY * ni + BX;
    shu[0][0][threadIdx.y+BY][threadIdx.x+BX] = u(0,idx);
    shu[1][0][threadIdx.y+BY][threadIdx.x+BX] = u(1,idx);
    shu[2][0][threadIdx.y+BY][threadIdx.x+BX] = u(2,idx);
    idx += nij;
    shu[0][1][threadIdx.y+BY][threadIdx.x+BX] = u(0,idx);
    shu[1][1][threadIdx.y+BY][threadIdx.x+BX] = u(1,idx);
    shu[2][1][threadIdx.y+BY][threadIdx.x+BX] = u(2,idx);
    idx += nij;
    shu[0][2][threadIdx.y+BY][threadIdx.x+BX] = u(0,idx);
    shu[1][2][threadIdx.y+BY][threadIdx.x+BX] = u(1,idx);
    shu[2][2][threadIdx.y+BY][threadIdx.x+BX] = u(2,idx);
    idx += nij;
    shu[0][3][threadIdx.y+BY][threadIdx.x+BX] = u(0,idx);
    shu[1][3][threadIdx.y+BY][threadIdx.x+BX] = u(1,idx);
    shu[2][3][threadIdx.y+BY][threadIdx.x+BX] = u(2,idx);
    idx += nij;
    d0 = u(0,idx);
    d1 = u(1,idx);
    d2 = u(2,idx);
  }

  if (active) {
    // Load strx & stry, centered at i and j
    strx_im2 = a_strx[i-2];
    strx_im1 = a_strx[i-1];
    strx_i   = a_strx[i];
    strx_ip1 = a_strx[i+1];
    strx_ip2 = a_strx[i+2];
    
    stry_jm2 = a_stry[j-2];
    stry_jm1 = a_stry[j-1];
    stry_j   = a_stry[j];
    stry_jp1 = a_stry[j+1];
    stry_jp2 = a_stry[j+2];

    // Load first 4 lambda and mu values 
    int idx = index + 2 * ni + 2;
    la_km2 = a_lambda[idx];
    mu_km2 = a_mu[idx];
    idx += nij;
    la_km1 = a_lambda[idx];
    mu_km1 = a_mu[idx];
    idx += nij;
    la_ijk = a_lambda[idx];
    mu_ijk = a_mu[idx];
    idx += nij;
    la_kp1 = a_lambda[idx];
    mu_kp1 = a_mu[idx];
  }

  // Move the loading index at kfirst + 2
  index += 4* nij; 

  // Load the first 4 values of strz
  strz_km2 = a_strz[kfirst-2];
  strz_km1 = a_strz[kfirst-1];
  strz_k   = a_strz[kfirst];
  strz_kp1 = a_strz[kfirst+1];

  //cof = 1.0/(h*h);

  // Main loop on K dimension
  for (int k = kfirst; k <= klast; k++) {

    // Rotate the shared memory indices 
    kthtmp = kthm4;
    kthm4  = kthm3;
    kthm3  = kthm2;
    kthm2  = kthm1;
    kthm1  = kthm0;
    kthm0  = kthtmp;

    // Rotate the shared memory and store the prefetched variables
    __syncthreads();

    shu[0][kthm0][threadIdx.y][threadIdx.x] = a0;
    shu[1][kthm0][threadIdx.y][threadIdx.x] = a1;
    shu[2][kthm0][threadIdx.y][threadIdx.x] = a2;

    if (threadIdx.x < 4) {
      shu[0][kthm0][threadIdx.y][threadIdx.x+BX] = b0;
      shu[1][kthm0][threadIdx.y][threadIdx.x+BX] = b1;
      shu[2][kthm0][threadIdx.y][threadIdx.x+BX] = b2;
    }
    if (threadIdx.y < 4) {
      shu[0][kthm0][threadIdx.y+BY][threadIdx.x] = c0;
      shu[1][kthm0][threadIdx.y+BY][threadIdx.x] = c1;
      shu[2][kthm0][threadIdx.y+BY][threadIdx.x] = c2;
      if (threadIdx.x < 4) {
        shu[0][kthm0][threadIdx.y+BY][threadIdx.x+BX] = d0;
        shu[1][kthm0][threadIdx.y+BY][threadIdx.x+BX] = d1;
        shu[2][kthm0][threadIdx.y+BY][threadIdx.x+BX] = d2;
      }
    }
    __syncthreads();
    
    // Prefetch the next iteration variables for U to better hide the latency
    if (k < klast) {
      int idx = index + nij;
      if (loader) {
	a0 = u(0,idx);
	a1 = u(1,idx);
	a2 = u(2,idx);
      }
      idx += BX;
      if (loadx2) {
	b0 = u(0,idx);
	b1 = u(1,idx);
	b2 = u(2,idx);
      }
      idx = index + nij + BY * ni;
      if (loader & loady2) {
	c0 = u(0,idx);
	c1 = u(1,idx);
	c2 = u(2,idx);
      }
      idx += BX;
      if (loadx2 && loady2) {
	d0 = u(0,idx);
	d1 = u(1,idx);
	d2 = u(2,idx);
      }
    }

    // Active threads load center values and do the computation
    if (active) {
      // Load a new value of strz
      strz_kp2 = a_strz[k+2];

      // Load new values of lambda and mu, and the (i,j) halos
      {
	int idx = index + 2 * ni + 2;
	la_kp2 = a_lambda[idx];
	mu_kp2 = a_mu[idx];
	// Halos at k
	idx -= 2 * nij;
	la_im2 = a_lambda[idx-2];
	la_im1 = a_lambda[idx-1];
	la_ip1 = a_lambda[idx+1];
	la_ip2 = a_lambda[idx+2];
	la_jm2 = a_lambda[idx-2*ni];
	la_jm1 = a_lambda[idx-1*ni];
	la_jp1 = a_lambda[idx+ni];
	la_jp2 = a_lambda[idx+2*ni];
	mu_im2 = a_mu[idx-2];
	mu_im1 = a_mu[idx-1];
	mu_ip1 = a_mu[idx+1];
	mu_ip2 = a_mu[idx+2];
	mu_jm2 = a_mu[idx-2*ni];
	mu_jm1 = a_mu[idx-ni];
	mu_jp1 = a_mu[idx+ni];
	mu_jp2 = a_mu[idx+2*ni];
	// Read rho, fo, um, up in advance
	// fo0 = fo(0,idx);
	// fo1 = fo(1,idx);
	// fo2 = fo(2,idx);
	//rho = a_rho[idx];
	// if (pred) {
	//   // um0 = um(0,idx);
	//   // um1 = um(1,idx);
	//   // um2 = um(2,idx);
	// }
	// else {
	//   up0 = up(0,idx);
	//   up1 = up(1,idx);
	//   up2 = up(2,idx);
	// }
      }

      // X contribution
      {
	float_sw4 mux1, mux2, mux3, mux4;
	mux1 = mu_im1 * strx_im1 - tf * (mu_ijk * strx_i + mu_im2 * strx_im2);
	mux2 = mu_im2 * strx_im2 + mu_ip1 * strx_ip1 + 3 * (mu_ijk * strx_i   + mu_im1 * strx_im1);
	mux3 = mu_im1 * strx_im1 + mu_ip2 * strx_ip2 + 3 * (mu_ip1 * strx_ip1 + mu_ijk * strx_i  );
	mux4 = mu_ip1 * strx_ip1 - tf * (mu_ijk * strx_i + mu_ip2 * strx_ip2);

	r1 = strx_i * ((2 * mux1 + la_im1 * strx_im1 - tf * (la_ijk * strx_i   + la_im2 * strx_im2)) *
		       (shu[0][kthm2][tj+0][ti-2] - shu[0][kthm2][tj+0][ti+0]) +
		       (2 * mux2 + la_im2 * strx_im2 + la_ip1 * strx_ip1 +
			3 * (la_ijk * strx_i   + la_im1 * strx_im1)) *
		       (shu[0][kthm2][tj+0][ti-1] - shu[0][kthm2][tj+0][ti+0]) +
		       (2 * mux3 + la_im1 * strx_im1 + la_ip2 * strx_ip2 +
			3 * (la_ip1 * strx_ip1 + la_ijk * strx_i  )) *
		       (shu[0][kthm2][tj+0][ti+1] - shu[0][kthm2][tj+0][ti+0]) +
		       (2 * mux4 + la_ip1 * strx_ip1 - tf * (la_ijk * strx_i   + la_ip2 * strx_ip2)) *
		       (shu[0][kthm2][tj+0][ti+2] - shu[0][kthm2][tj+0][ti+0]));
	r2 = strx_i * (mux1 * (shu[1][kthm2][tj+0][ti-2] - shu[1][kthm2][tj+0][ti+0]) +
		       mux2 * (shu[1][kthm2][tj+0][ti-1] - shu[1][kthm2][tj+0][ti+0]) +
		       mux3 * (shu[1][kthm2][tj+0][ti+1] - shu[1][kthm2][tj+0][ti+0]) +
		       mux4 * (shu[1][kthm2][tj+0][ti+2] - shu[1][kthm2][tj+0][ti+0]));
	r3 = strx_i * (mux1 * (shu[2][kthm2][tj+0][ti-2] - shu[2][kthm2][tj+0][ti+0]) +
		       mux2 * (shu[2][kthm2][tj+0][ti-1] - shu[2][kthm2][tj+0][ti+0]) +
		       mux3 * (shu[2][kthm2][tj+0][ti+1] - shu[2][kthm2][tj+0][ti+0]) +
		       mux4 * (shu[2][kthm2][tj+0][ti+2] - shu[2][kthm2][tj+0][ti+0]));
      }
      // Y contribution
      {
	float_sw4 muy1, muy2, muy3, muy4;
	muy1 = mu_jm1 * stry_jm1 - tf * (mu_ijk * stry_j + mu_jm2 * stry_jm2);
	muy2 = mu_jm2 * stry_jm2 + mu_jp1 * stry_jp1 + 3 * (mu_ijk * stry_j   + mu_jm1 * stry_jm1);
	muy3 = mu_jm1 * stry_jm1 + mu_jp2 * stry_jp2 + 3 * (mu_jp1 * stry_jp1 + mu_ijk * stry_j  );
	muy4 = mu_jp1 * stry_jp1 - tf * (mu_ijk * stry_j + mu_jp2 * stry_jp2);

	r1 += stry_j * (muy1 * (shu[0][kthm2][tj-2][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muy2 * (shu[0][kthm2][tj-1][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muy3 * (shu[0][kthm2][tj+1][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muy4 * (shu[0][kthm2][tj+2][ti+0] - shu[0][kthm2][tj+0][ti+0]));
	r2 += stry_j * ((2 * muy1 + la_jm1 * stry_jm1 - tf * (la_ijk * stry_j + la_jm2 * stry_jm2)) *
			(shu[1][kthm2][tj-2][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			(2 * muy2 + la_jm2 * stry_jm2 + la_jp1 * stry_jp1 +
			 3 * (la_ijk * stry_j   + la_jm1 * stry_jm1)) *
			(shu[1][kthm2][tj-1][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			(2 * muy3 + la_jm1 * stry_jm1 + la_jp2 * stry_jp2 +
			 3 * (la_jp1 * stry_jp1 + la_ijk * stry_j  )) *
			(shu[1][kthm2][tj+1][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			(2 * muy4 + la_jp1 * stry_jp1 - tf * (la_ijk * stry_j  + la_jp2 * stry_jp2)) *
			(shu[1][kthm2][tj+2][ti+0] - shu[1][kthm2][tj+0][ti+0]));
	r3 += stry_j  *(muy1 * (shu[2][kthm2][tj-2][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			muy2 * (shu[2][kthm2][tj-1][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			muy3 * (shu[2][kthm2][tj+1][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			muy4 * (shu[2][kthm2][tj+2][ti+0] - shu[2][kthm2][tj+0][ti+0]));
      }
      // Z contribution
      {
	float_sw4 muz1, muz2, muz3, muz4;
	muz1 = mu_km1 * strz_km1 - tf * (mu_ijk * strz_k + mu_km2 * strz_km2);
	muz2 = mu_km2 * strz_km2 + mu_kp1 * strz_kp1 + 3 * (mu_ijk * strz_k   + mu_km1 * strz_km1);
	muz3 = mu_km1 * strz_km1 + mu_kp2 * strz_kp2 + 3 * (mu_kp1 * strz_kp1 + mu_ijk * strz_k);
	muz4 = mu_kp1 * strz_kp1 - tf * (mu_ijk * strz_k + mu_kp2 * strz_kp2);
	r1 += strz_k * (muz1 * (shu[0][kthm4][tj+0][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muz2 * (shu[0][kthm3][tj+0][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muz3 * (shu[0][kthm1][tj+0][ti+0] - shu[0][kthm2][tj+0][ti+0]) +
			muz4 * (shu[0][kthm0][tj+0][ti+0] - shu[0][kthm2][tj+0][ti+0]));
	r2 += strz_k * (muz1 * (shu[1][kthm4][tj+0][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			muz2 * (shu[1][kthm3][tj+0][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			muz3 * (shu[1][kthm1][tj+0][ti+0] - shu[1][kthm2][tj+0][ti+0]) +
			muz4 * (shu[1][kthm0][tj+0][ti+0] - shu[1][kthm2][tj+0][ti+0]));
	r3 += strz_k * ((2 * muz1 + la_km1 * strz_km1 - tf * (la_ijk * strz_k + la_km2 * strz_km2)) *
			(shu[2][kthm4][tj+0][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			(2 * muz2 + la_km2 * strz_km2 + la_kp1 * strz_kp1 +
			 3 * (la_ijk * strz_k   + la_km1 * strz_km1)) *
			(shu[2][kthm3][tj+0][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			(2 * muz3 + la_km1 * strz_km1 + la_kp2 * strz_kp2 +
			 3 * (la_kp1 * strz_kp1 + la_ijk * strz_k  )) *
			(shu[2][kthm1][tj+0][ti+0] - shu[2][kthm2][tj+0][ti+0]) +
			(2 * muz4 + la_kp1 * strz_kp1 - tf * (la_ijk * strz_k + la_kp2 * strz_kp2)) *
			(shu[2][kthm0][tj+0][ti+0] - shu[2][kthm2][tj+0][ti+0]));
      }
      r1 *= i6;
      r2 *= i6;
      r3 *= i6;


      /* Mixed derivatives: */
      /*   (la*v_y)_x */
      r1 += strx_i  * stry_j * i144 * 
	(la_im2 * (shu[1][kthm2][tj-2][ti-2] - shu[1][kthm2][tj+2][ti-2] + 8 *
		   (-shu[1][kthm2][tj-1][ti-2] + shu[1][kthm2][tj+1][ti-2])) - 8 *
	 (la_im1 * (shu[1][kthm2][tj-2][ti-1] - shu[1][kthm2][tj+2][ti-1] + 8 *
		    (-shu[1][kthm2][tj-1][ti-1] + shu[1][kthm2][tj+1][ti-1]))) + 8 *
	 (la_ip1 * (shu[1][kthm2][tj-2][ti+1] - shu[1][kthm2][tj+2][ti+1] + 8 *
		    (-shu[1][kthm2][tj-1][ti+1] + shu[1][kthm2][tj+1][ti+1]))) - 
	 (la_ip2 * (shu[1][kthm2][tj-2][ti+2] - shu[1][kthm2][tj+2][ti+2] + 8 *
		    (-shu[1][kthm2][tj-1][ti+2] + shu[1][kthm2][tj+1][ti+2]))))
	/*   (la*w_z)_x */
        + strx_i * strz_k * i144 *
	(la_im2 * (shu[2][kthm4][tj+0][ti-2] - shu[2][kthm0][tj+0][ti-2] + 8 *
		   (-shu[2][kthm3][tj+0][ti-2] + shu[2][kthm1][tj+0][ti-2])) - 8 *
	 (la_im1 * (shu[2][kthm4][tj+0][ti-1] - shu[2][kthm0][tj+0][ti-1] + 8 *
		    (-shu[2][kthm3][tj+0][ti-1] + shu[2][kthm1][tj+0][ti-1]))) + 8 *
	 (la_ip1 * (shu[2][kthm4][tj+0][ti+1] - shu[2][kthm0][tj+0][ti+1] + 8 *
		    (-shu[2][kthm3][tj+0][ti+1] + shu[2][kthm1][tj+0][ti+1]))) - 
	 (la_ip2 * (shu[2][kthm4][tj+0][ti+2] - shu[2][kthm0][tj+0][ti+2] + 8 *
		    (-shu[2][kthm3][tj+0][ti+2] + shu[2][kthm1][tj+0][ti+2]))))
	/*   (mu*v_x)_y */
        + strx_i * stry_j * i144 *
	(mu_jm2 * (shu[1][kthm2][tj-2][ti-2] - shu[1][kthm2][tj-2][ti+2] + 8 *
		   (-shu[1][kthm2][tj-2][ti-1] + shu[1][kthm2][tj-2][ti+1])) - 8 *
	 (mu_jm1 * (shu[1][kthm2][tj-1][ti-2] - shu[1][kthm2][tj-1][ti+2] + 8 *
		    (-shu[1][kthm2][tj-1][ti-1] + shu[1][kthm2][tj-1][ti+1]))) + 8 *
	 (mu_jp1 * (shu[1][kthm2][tj+1][ti-2] - shu[1][kthm2][tj+1][ti+2] + 8 *
		    (-shu[1][kthm2][tj+1][ti-1] + shu[1][kthm2][tj+1][ti+1]))) - 
	 (mu_jp2 * (shu[1][kthm2][tj+2][ti-2] - shu[1][kthm2][tj+2][ti+2] + 8 *
		    (-shu[1][kthm2][tj+2][ti-1] + shu[1][kthm2][tj+2][ti+1]))))
	/*   (mu*w_x)_z */
        + strx_i * strz_k * i144 *
	(mu_km2 * (shu[2][kthm4][tj+0][ti-2] - shu[2][kthm4][tj+0][ti+2] + 8 *
		   (-shu[2][kthm4][tj+0][ti-1] + shu[2][kthm4][tj+0][ti+1])) - 8 *
 	 (mu_km1 * (shu[2][kthm3][tj+0][ti-2] - shu[2][kthm3][tj+0][ti+2] + 8 *
		    (-shu[2][kthm3][tj+0][ti-1] + shu[2][kthm3][tj+0][ti+1]))) + 8 *
 	 (mu_kp1 * (shu[2][kthm1][tj+0][ti-2] - shu[2][kthm1][tj+0][ti+2] + 8 *
		    (-shu[2][kthm1][tj+0][ti-1] + shu[2][kthm1][tj+0][ti+1]))) -
	 (mu_kp2 * (shu[2][kthm0][tj+0][ti-2] - shu[2][kthm0][tj+0][ti+2] + 8 *
		    (-shu[2][kthm0][tj+0][ti-1] + shu[2][kthm0][tj+0][ti+1]))));

      /*   (mu*u_y)_x */
      r2 += strx_i  *stry_j  *i144*
	(mu_im2 * (shu[0][kthm2][tj-2][ti-2] - shu[0][kthm2][tj+2][ti-2] + 8 *
		   (-shu[0][kthm2][tj-1][ti-2] + shu[0][kthm2][tj+1][ti-2])) - 8 *
	 (mu_im1 * (shu[0][kthm2][tj-2][ti-1] - shu[0][kthm2][tj+2][ti-1] + 8 *
		    (-shu[0][kthm2][tj-1][ti-1] + shu[0][kthm2][tj+1][ti-1]))) + 8 *
	 (mu_ip1 * (shu[0][kthm2][tj-2][ti+1] - shu[0][kthm2][tj+2][ti+1] + 8 *
		    (-shu[0][kthm2][tj-1][ti+1] + shu[0][kthm2][tj+1][ti+1]))) -
	 (mu_ip2 * (shu[0][kthm2][tj-2][ti+2] - shu[0][kthm2][tj+2][ti+2] + 8 *
		    (-shu[0][kthm2][tj-1][ti+2] + shu[0][kthm2][tj+1][ti+2]))))
	/* (la*u_x)_y */
        + strx_i * stry_j * i144 *
	(la_jm2 * (shu[0][kthm2][tj-2][ti-2] - shu[0][kthm2][tj-2][ti+2] + 8 *
		   (-shu[0][kthm2][tj-2][ti-1] + shu[0][kthm2][tj-2][ti+1])) - 8 *
	 (la_jm1 * (shu[0][kthm2][tj-1][ti-2] - shu[0][kthm2][tj-1][ti+2] + 8 *
		    (-shu[0][kthm2][tj-1][ti-1] + shu[0][kthm2][tj-1][ti+1]))) + 8 *
	 (la_jp1 * (shu[0][kthm2][tj+1][ti-2] - shu[0][kthm2][tj+1][ti+2] + 8 *
		    (-shu[0][kthm2][tj+1][ti-1] + shu[0][kthm2][tj+1][ti+1]))) -
	 (la_jp2 * (shu[0][kthm2][tj+2][ti-2] - shu[0][kthm2][tj+2][ti+2] + 8 *
		    (-shu[0][kthm2][tj+2][ti-1] + shu[0][kthm2][tj+2][ti+1]))))
	/* (la*w_z)_y */
        + stry_j * strz_k * i144 *
	(la_jm2 * (shu[2][kthm4][tj-2][ti+0] - shu[2][kthm0][tj-2][ti+0] + 8 *
		   (-shu[2][kthm3][tj-2][ti+0] + shu[2][kthm1][tj-2][ti+0])) - 8 *
	 (la_jm1 * (shu[2][kthm4][tj-1][ti+0] - shu[2][kthm0][tj-1][ti+0] + 8 *
		    (-shu[2][kthm3][tj-1][ti+0] + shu[2][kthm1][tj-1][ti+0]))) + 8 *
	 (la_jp1 * (shu[2][kthm4][tj+1][ti+0] - shu[2][kthm0][tj+1][ti+0] + 8 *
		    (-shu[2][kthm3][tj+1][ti+0] + shu[2][kthm1][tj+1][ti+0]))) -
	 (la_jp2 * (shu[2][kthm4][tj+2][ti+0] - shu[2][kthm0][tj+2][ti+0] + 8 *
		    (-shu[2][kthm3][tj+2][ti+0] + shu[2][kthm1][tj+2][ti+0]))))
	/* (mu*w_y)_z */
        + stry_j * strz_k * i144 *
	(mu_km2 * (shu[2][kthm4][tj-2][ti+0] - shu[2][kthm4][tj+2][ti+0] + 8 *
		   (-shu[2][kthm4][tj-1][ti+0] + shu[2][kthm4][tj+1][ti+0])) - 8 *
	 (mu_km1 * (shu[2][kthm3][tj-2][ti+0] - shu[2][kthm3][tj+2][ti+0] + 8 *
		    (-shu[2][kthm3][tj-1][ti+0] + shu[2][kthm3][tj+1][ti+0]))) + 8 *
	 (mu_kp1 * (shu[2][kthm1][tj-2][ti+0] - shu[2][kthm1][tj+2][ti+0] + 8 *
		    (-shu[2][kthm1][tj-1][ti+0] + shu[2][kthm1][tj+1][ti+0]))) -
	 (mu_kp2 * (shu[2][kthm0][tj-2][ti+0] - shu[2][kthm0][tj+2][ti+0] + 8 *
		    (-shu[2][kthm0][tj-1][ti+0] + shu[2][kthm0][tj+1][ti+0]))));

      /*  (mu*u_z)_x */
      r3 += strx_i * strz_k * i144 *
	(mu_im2 * (shu[0][kthm4][tj+0][ti-2] - shu[0][kthm0][tj+0][ti-2] + 8 *
		   (-shu[0][kthm3][tj+0][ti-2] + shu[0][kthm1][tj+0][ti-2])) - 8 *
	 (mu_im1 * (shu[0][kthm4][tj+0][ti-1] - shu[0][kthm0][tj+0][ti-1] + 8 *
		    (-shu[0][kthm3][tj+0][ti-1] + shu[0][kthm1][tj+0][ti-1]))) + 8 *
	 (mu_ip1 * (shu[0][kthm4][tj+0][ti+1] - shu[0][kthm0][tj+0][ti+1] + 8 *
		    (-shu[0][kthm3][tj+0][ti+1] + shu[0][kthm1][tj+0][ti+1]))) -
	 (mu_ip2 * (shu[0][kthm4][tj+0][ti+2] - shu[0][kthm0][tj+0][ti+2] + 8 *
		    (-shu[0][kthm3][tj+0][ti+2] + shu[0][kthm1][tj+0][ti+2]))))
	/* (mu*v_z)_y */
        + stry_j  *strz_k*i144*
	(mu_jm2 * (shu[1][kthm4][tj-2][ti+0] - shu[1][kthm0][tj-2][ti+0] + 8 *
		   (- shu[1][kthm3][tj-2][ti+0]+shu[1][kthm1][tj-2][ti+0])) - 8 * 
	 (mu_jm1 * (shu[1][kthm4][tj-1][ti+0] - shu[1][kthm0][tj-1][ti+0] + 8 *
		    (- shu[1][kthm3][tj-1][ti+0]+shu[1][kthm1][tj-1][ti+0]))) + 8 * 
	 (mu_jp1 * (shu[1][kthm4][tj+1][ti+0] - shu[1][kthm0][tj+1][ti+0] + 8 *
		    (- shu[1][kthm3][tj+1][ti+0]+shu[1][kthm1][tj+1][ti+0]))) -
	 (mu_jp2 * (shu[1][kthm4][tj+2][ti+0] - shu[1][kthm0][tj+2][ti+0] + 8 *
		    (- shu[1][kthm3][tj+2][ti+0]+shu[1][kthm1][tj+2][ti+0]))))
	/*   (la*u_x)_z */
        + strx_i  *strz_k*i144*
	(la_km2*(shu[0][kthm4][tj+0][ti-2] - shu[0][kthm4][tj+0][ti+2] + 8 *
		 (- shu[0][kthm4][tj+0][ti-1]+shu[0][kthm4][tj+0][ti+1])) - 8*
	 (la_km1*(shu[0][kthm3][tj+0][ti-2] - shu[0][kthm3][tj+0][ti+2] + 8 *
		  (- shu[0][kthm3][tj+0][ti-1]+shu[0][kthm3][tj+0][ti+1]))) + 8*
	 (la_kp1*(shu[0][kthm1][tj+0][ti-2] - shu[0][kthm1][tj+0][ti+2] + 8 *
		  (- shu[0][kthm1][tj+0][ti-1]+shu[0][kthm1][tj+0][ti+1]))) -
	 (la_kp2*(shu[0][kthm0][tj+0][ti-2] - shu[0][kthm0][tj+0][ti+2] + 8 *
		  (- shu[0][kthm0][tj+0][ti-1]+shu[0][kthm0][tj+0][ti+1]))))
	/* (la*v_y)_z */
        + stry_j  *strz_k*i144*
	(la_km2*(shu[1][kthm4][tj-2][ti+0] - shu[1][kthm4][tj+2][ti+0] + 8 *
		 (- shu[1][kthm4][tj-1][ti+0]+shu[1][kthm4][tj+1][ti+0])) - 8*
	 (la_km1*(shu[1][kthm3][tj-2][ti+0] - shu[1][kthm3][tj+2][ti+0] + 8 *
		  (- shu[1][kthm3][tj-1][ti+0]+shu[1][kthm3][tj+1][ti+0]))) + 8*
	 (la_kp1*(shu[1][kthm1][tj-2][ti+0] - shu[1][kthm1][tj+2][ti+0] + 8 *
		  (- shu[1][kthm1][tj-1][ti+0]+shu[1][kthm1][tj+1][ti+0]))) -
	 (la_kp2*(shu[1][kthm0][tj-2][ti+0] - shu[1][kthm0][tj+2][ti+0] + 8 *
		  (- shu[1][kthm0][tj-1][ti+0]+shu[1][kthm0][tj+1][ti+0]))));

      
      // Use the result at point (i,j,k)
      if (pred) {
	// Apply predictor
	// int idx = k * nij + j * ni + i;
	// float_sw4 fact = dt*dt / rho;
	// up(0,idx) = 2 * shu[0][kthm2][tj][ti] - um0 + fact * (cof * r1 + fo0);
	// up(1,idx) = 2 * shu[1][kthm2][tj][ti] - um1 + fact * (cof * r2 + fo1);
	// up(2,idx) = 2 * shu[2][kthm2][tj][ti] - um2 + fact * (cof * r3 + fo2);
      }
      else {
	// Apply 4th order correction
	int idx = k * nij + j * ni + i;
	//float_sw4 fact = dt*dt*dt*dt / (12 * rho);
	up(0,idx) = a1fac*up(0,idx)+ cof * r1;
	up(1,idx) = a1fac*up(1,idx) + cof * r2;
	up(2,idx) = a1fac*up(2,idx) + cof * r3;
      }
      
      // Rotate the register queues
      strz_km2 = strz_km1;
      strz_km1 = strz_k;
      strz_k = strz_kp1;
      strz_kp1 = strz_kp2;

      la_km2 = la_km1;
      la_km1 = la_ijk;
      la_ijk = la_kp1;
      la_kp1 = la_kp2;

      mu_km2 = mu_km1;
      mu_km1 = mu_ijk;
      mu_ijk = mu_kp1;
      mu_kp1 = mu_kp2;

    }  // End active

    // Move the index to the next K
    index += nij;

  } // End K loop
}
#undef u
#undef up
void rhs4th3fortsgstr_ciopt( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast,
			     int ni,int nj,int nk,
			     float_sw4* a_lu, float_sw4* a_u,
			     float_sw4* a_mu, float_sw4* a_lambda, 
			     float_sw4 h, float_sw4* a_strx, float_sw4* a_stry, 
			     float_sw4* a_strz, char op );
 
#endif
