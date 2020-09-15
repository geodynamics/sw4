/* Fast marching main interface. */
/*
  Copyright (C) 2004 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// modified to use in FWI
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <float.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>



#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))

#define SF_ABS(a)   ((a) >= 0  ? (a) : (-(a)))
#define SF_SIG(a)   ((a) >= 0  ?  1  :  -1 )

#define SF_NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#define SF_ODD(n)  ((n) & 1)
#define SF_EVEN(n) (!(SF_ODD(n)))

#define SF_PI (3.14159265358979323846264338328)

#define SF_EPS FLT_EPSILON
#define SF_HUGE FLT_MAX

float **x, **xn, **x1;


enum {SF_IN, SF_FRONT, SF_OUT};

/*------------------------------------------------------------*/
/*@out@*/ void *sf_alloc (size_t n    /* number of elements */,
              size_t size /* size of one element */)
      /*< output-checking allocation >*/
{
    void *ptr;

    size *= n;

    if (0>=size) fprintf(stderr, "%s: illegal allocation (%ld bytes)",__FILE__, size);

    ptr = malloc (size);

    if (NULL == ptr)
    fprintf(stderr, "%s: cannot allocate %lu bytes:", __FILE__,size);

    return ptr;
}


void sf_pqueue_init (int n)
/*< Initialize heap with the maximum size >*/
{
    //fprintf(stderr, "pqueue_init sizeof(float*)=%ld\n", sizeof(float*));

    x = (float **) sf_alloc ((n+1),sizeof (float *));
    if(x==NULL) {
        fprintf(stderr, "unable to sf_alloc **x\n");
    }
}

void sf_pqueue_start (void)
/*< Set starting values >*/
{
    xn = x;
    x1 = x+1;
}

void sf_pqueue_close (void)
/*< Free the allocated storage >*/
{
    free (x);
}

void sf_pqueue_insert (float* v)
/*< Insert an element (smallest first) >*/
{
    float **xi, **xq;
    unsigned int q;

    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
    xq = x + q;
    if (*v > **xq) break;
    *xi = *xq; xi = xq;
    }
    *xi = v;
}

void sf_pqueue_insert2 (float* v)
/*< Insert an element (largest first) >*/
{
    float **xi, **xq;
    unsigned int q;

    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
    xq = x + q;
    if (*v < **xq) break;
    *xi = *xq; xi = xq;
    }
    *xi = v;
}

float* sf_pqueue_extract (void)
/*< Extract the smallest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;

    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
    xc = x + c;
    if (c < (unsigned int) n && **xc > **(xc+1)) {
        c++; xc++;
    }
    if (*t <= **xc) break;
    *xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

float* sf_pqueue_extract2 (void)
/*< Extract the largest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;

    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
    xc = x + c;
    if (c < (unsigned int) n && **xc < **(xc+1)) {
        c++; xc++;
    }
    if (*t >= **xc) break;
    *xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

void sf_pqueue_update (float **v)
/*< restore the heap: the value has been altered >*/
{
  unsigned int c;
  int n;
  float **xc, **xi;

  xi = v;
  n = (int) (xn-x); c = (unsigned int) (xi-x);
  for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
      xc = x + c;
      if (c < (unsigned int) n && **xc > **(xc+1)) {
      c++; xc++;
      }
      if (**v <= **xc) break;
      *xi = *xc; xi = xc;
  }
  xi = v; c = (unsigned int) (xi-x);
  for (c >>= 1; c > 0; c >>= 1) {
      xc = x + c;
      if (**v > **xc) break;
      *xi = *xc; xi = xc;
  }
  *xi = *v;
}

struct Upd {
    double stencil, value;
    double delta;
};

static int update (float value, int i);
static int update2 (float value, int i);
static float qsolve(int i);
static float qsolve2(int i);
static void stencil (float t, struct Upd *x);
static bool updaten (int m, float* res, struct Upd *v[]);
static bool updaten2 (int m, float* res, struct Upd *v[]);
static void grid (int *i, const int *n);

static int *in, *n, s[3], order;
static float *ttime, *vv, rdx[3];
static double v1;

void sf_neighbors_init (int *in1     /* status flag [n[0]*n[1]*n[2]] */,
            float *rdx1  /* grid sampling [3] */,
            int *n1      /* grid samples [3] */,
            int order1   /* accuracy order */,
            float *time1 /* traveltime [n[0]*n[1]*n[2]] */)
/*< Initialize >*/
{
    in = in1; ttime = time1;
    n = n1; order = order1;
    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];
    rdx[0] = 1./(rdx1[0]*rdx1[0]);
    rdx[1] = 1./(rdx1[1]*rdx1[1]);
    rdx[2] = 1./(rdx1[2]*rdx1[2]);
}

int  sf_neighbours(int i)
/*< Update neighbors of gridpoint i, return number of updated points >*/
{
    int j, k, ix, npoints;

    npoints = 0;
    for (j=0; j < 3; j++) {
    ix = (i/s[j])%n[j];
    if (ix+1 <= n[j]-1) {
        k = i+s[j];
        if (in[k] != SF_IN) npoints += update(qsolve(k),k);
    }
    if (ix-1 >= 0  ) {
        k = i-s[j];
        if (in[k] != SF_IN) npoints += update(qsolve(k),k);
    }
    }
    return npoints;
}

int  sf_neighbours2(int i)
/*< Update neighbors of gridpoint i, return number of updated points >*/
{
    int j, k, ix, npoints;

    npoints = 0;
    for (j=0; j < 3; j++) {
    ix = (i/s[j])%n[j];
    if (ix+1 <= n[j]-1) {
        k = i+s[j];
        if (in[k] != SF_IN) npoints += update2(qsolve2(k),k);
    }
    if (ix-1 >= 0  ) {
        k = i-s[j];
        if (in[k] != SF_IN) npoints += update2(qsolve2(k),k);
    }
    }
    return npoints;
}

static int update (float value, int i)
/* update gridpoint i with new value */
{
    if (value < ttime[i]) {
    ttime[i]   = value;
    if (in[i] == SF_OUT) {
        in[i] = SF_FRONT;
        sf_pqueue_insert (ttime+i);
        return 1;
    }
/*	sf_pqueue_update (&(ttime+i)); */
    }

    return 0;
}

static int update2 (float value, int i)
/* update gridpoint i with new value */
{
    if (value > ttime[i]) {
    ttime[i]   = value;
    if (in[i] == SF_OUT) {
        in[i] = SF_FRONT;
        sf_pqueue_insert2 (ttime+i);
        return 1;
    }
/*	sf_pqueue_update (&(ttime+i)); */
    }

    return 0;
}

static float qsolve(int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix;
    float a, b, t, res;
    struct Upd *v[3], x[3], *xj;

    for (j=0; j<3; j++) {
    ix = (i/s[j])%n[j];

    if (ix > 0) {
        k = i-s[j];
        a = ttime[k];
    } else {
        a = SF_HUGE;
    }

    if (ix < n[j]-1) {
        k = i+s[j];
        b = ttime[k];
    } else {
        b = SF_HUGE;
    }

    xj = x+j;
    xj->delta = rdx[j];

    if (a < b) {
        xj->stencil = xj->value = a;
    } else {
        xj->stencil = xj->value = b;
    }

    if (order > 1) {
        if (a < b  && ix-2 >= 0) {
        k = i-2*s[j];
        if (in[k] != SF_OUT && a >= (t=ttime[k]))
            stencil(t,xj);
        }
        if (a > b && ix+2 <= n[j]-1) {
        k = i+2*s[j];
        if (in[k] != SF_OUT && b >= (t=ttime[k]))
            stencil(t,xj);
        }
    }
    }

    if (x[0].value <= x[1].value) {
    if (x[1].value <= x[2].value) {
        v[0] = x; v[1] = x+1; v[2] = x+2;
    } else if (x[2].value <= x[0].value) {
        v[0] = x+2; v[1] = x; v[2] = x+1;
    } else {
        v[0] = x; v[1] = x+2; v[2] = x+1;
    }
    } else {
    if (x[0].value <= x[2].value) {
        v[0] = x+1; v[1] = x; v[2] = x+2;
    } else if (x[2].value <= x[1].value) {
        v[0] = x+2; v[1] = x+1; v[2] = x;
    } else {
        v[0] = x+1; v[1] = x+2; v[2] = x;
    }
    }

    v1=vv[i];

    if(v[2]->value < SF_HUGE) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
    if (updaten(3, &res, v) ||
        updaten(2, &res, v) ||
        updaten(1, &res, v)) return res;

    } else if(v[1]->value < SF_HUGE) { /* TWO DIRECTIONS CONTRIBUTE */
    if (updaten(2, &res, v) ||
        updaten(1, &res, v)) return res;

    } else if(v[0]->value < SF_HUGE) { /* ONE DIRECTION CONTRIBUTES */
    if (updaten(1, &res, v)) return res;

    }

    return SF_HUGE;
}

static float qsolve2(int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix;
    float a, b, t, res;
    struct Upd *v[3], x[3], *xj;

    for (j=0; j<3; j++) {
    ix = (i/s[j])%n[j];

    if (ix > 0) {
        k = i-s[j];
        a = ttime[k];
    } else {
        a = 0.;
    }

    if (ix < n[j]-1) {
        k = i+s[j];
        b = ttime[k];
    } else {
        b = 0.;
    }

    xj = x+j;
    xj->delta = rdx[j];

    if (a > b) {
        xj->stencil = xj->value = a;
    } else {
        xj->stencil = xj->value = b;
    }

    if (order > 1) {
        if (a > b  && ix-2 >= 0) {
        k = i-2*s[j];
        if (in[k] != SF_OUT && a <= (t=ttime[k]))
            stencil(t,xj);
        }
        if (a < b && ix+2 <= n[j]-1) {
        k = i+2*s[j];
        if (in[k] != SF_OUT && b <= (t=ttime[k]))
            stencil(t,xj);
        }
    }
    }

    if (x[0].value >= x[1].value) {
    if (x[1].value >= x[2].value) {
        v[0] = x; v[1] = x+1; v[2] = x+2;
    } else if (x[2].value >= x[0].value) {
        v[0] = x+2; v[1] = x; v[2] = x+1;
    } else {
        v[0] = x; v[1] = x+2; v[2] = x+1;
    }
    } else {
    if (x[0].value >= x[2].value) {
        v[0] = x+1; v[1] = x; v[2] = x+2;
    } else if (x[2].value >= x[1].value) {
        v[0] = x+2; v[1] = x+1; v[2] = x;
    } else {
        v[0] = x+1; v[1] = x+2; v[2] = x;
    }
    }

    v1=vv[i];

    if(v[2]->value > 0) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
    if (updaten2(3, &res, v) ||
        updaten2(2, &res, v) ||
        updaten2(1, &res, v)) return res;
    } else if(v[1]->value > 0) { /* TWO DIRECTIONS CONTRIBUTE */
    if (updaten2(2, &res, v) ||
        updaten2(1, &res, v)) return res;
    } else if(v[0]->value > 0) { /* ONE DIRECTION CONTRIBUTES */
    if (updaten2(1, &res, v)) return res;
    }

    return 0.;
}

static void stencil (float t, struct Upd *x)
/* second-order stencil */
{
    x->delta *= 2.25;
    x->stencil = (4.0*x->value - t)/3.0;
}

static bool updaten (int m, float* res, struct Upd *v[])
/* updating */
{
    double a, b, c, discr, t;
    int j;

    a = b = c = 0.;

    for (j=0; j<m; j++) {
    a += v[j]->delta;
    b += v[j]->stencil*v[j]->delta;
    c += v[j]->stencil*v[j]->stencil*v[j]->delta;
    }
    b /= a;

    discr=b*b+(v1-c)/a;

    if (discr < 0.) return false;

    t = b + sqrt(discr);
    if (t <= v[m-1]->value) return false;

    *res = t;
    return true;
}

static bool updaten2 (int m, float* res, struct Upd *v[])
/* updating */
{
    double a, b, c, discr, t;
    int j;

    a = b = c = 0.;

    for (j=0; j<m; j++) {
    a += v[j]->delta;
    b += v[j]->stencil*v[j]->delta;
    c += v[j]->stencil*v[j]->stencil*v[j]->delta;
    }
    b /= a;

    discr=b*b+(v1-c)/a;

    if (discr < 0.) return false;

    t = b - sqrt(discr);
    if (t >= v[m-1]->value) return false;

    *res = t;
    return true;
}

static void grid (int *i, const int *n)
/* restrict i[3] to the grid n[3] */
{
    int j;

    for (j=0; j < 3; j++) {
    if (i[j] < 0) {
        i[j]=0;
    } else if (i[j] >= n[j]) {
        i[j]=n[j]-1;
    }
    }
}

static int dist(int k, float x1, float x2, float x3)
/* assign distance to a neighboring grid point */
{
    float ti;

    ti = sqrtf(vv[k])*hypotf(x1,hypotf(x2,x3));
    if (SF_OUT == in[k]) {
    in[k] = SF_IN;
    ttime[k] = ti;
    sf_pqueue_insert (ttime+k);
    return 1;
    } else if (ti < ttime[k]) {
    ttime[k] = ti;
    }

    return 0;
}

int sf_neighbors_distance(int np         /* number of points */,
              float *vv1     /* slowness squared */,
              float **points /* point coordinates[np][3] */,
              float *d       /* grid sampling [3] */,
              float *o       /* grid origin [3] */)
/*< initialize distance computation >*/
{
    int ip, i, j, n123, ix[3], k;
    float x[3];

    n123 = n[0]*n[1]*n[2];

    vv = vv1;

    /* initialize everywhere */
    for (i=0; i < n123; i++) {
    in[i] = SF_OUT;
    ttime[i] = SF_HUGE;
    }

    for (ip=0; ip < np; ip++) {
    for (j=0; j < 3; j++) {
        x[j] = (points[ip][j]-o[j])/d[j];
        ix[j] = floorf(x[j]);
    }
    if (x[0] < 0. || ix[0] >= n[0] ||
        x[1] < 0. || ix[1] >= n[1] ||
        x[2] < 0. || ix[2] >= n[2]) continue;
    k = 0;
    for (j=0; j < 3; j++) {
        x[j] = (x[j]-ix[j])*d[j];
        k += ix[j]*s[j];
    }
    n123 -= dist(k,x[0],x[1],x[2]);
    if (ix[0] != n[0]-1) {
        n123 -= dist(k+s[0],d[0]-x[0],x[1],x[2]);
        if (ix[1] != n[1]-1) {
        n123 -= dist(k+s[0]+s[1],d[0]-x[0],d[1]-x[1],x[2]);
        if (ix[2] != n[2]-1)
            n123 -=
            dist(k+s[0]+s[1]+s[2],d[0]-x[0],d[1]-x[1],d[2]-x[2]);
        }
        if (ix[2] != n[2]-1)
        n123 -= dist(k+s[0]+s[2],d[0]-x[0],x[1],d[2]-x[2]);
    }
    if (ix[1] != n[1]-1) {
        n123 -= dist(k+s[1],x[0],d[1]-x[1],x[2]);
        if (ix[2] != n[2]-1)
        n123 -= dist(k+s[1]+s[2],x[0],d[1]-x[1],d[2]-x[2]);
    }
    if (ix[2] != n[2]-1) n123 -= dist(k+s[2],x[0],x[1],d[2]-x[2]);
    }

    return n123;
}

int sf_neighbors_nearsource(float* xs   /* source location [3] */,
                int* b      /* constant-velocity box around it [3] */,
                float* d    /* grid sampling [3] */,
                float* vv1  /* slowness [n[0]*n[1]*n[2]] */,
                bool *plane /* if plane-wave source */)
/*< initialize the source >*/
{
    int npoints, ic, i, j, is, start[3], endx[3], ix, iy, iz;
    double delta[3], delta2;


    /* initialize everywhere */
    for (i=0; i < n[0]*n[1]*n[2]; i++) {
    in[i] = SF_OUT;
    ttime[i] = SF_HUGE;
    }

    vv = vv1;

    /* Find index of the source location and project it to the grid */
    for (j=0; j < 3; j++) {
    is = xs[j]/d[j]+0.5;
    start[j] = is-b[j];
    endx[j]  = is+b[j];
    }

    grid(start, n);
    grid(endx, n);

    ic = (start[0]+endx[0])/2 +
    n[0]*((start[1]+endx[1])/2 +
          n[1]*(start[2]+endx[2])/2);

    v1 = vv[ic];

    /* loop in a small box around the source */
    npoints = n[0]*n[1]*n[2];
    for (ix=start[2]; ix <= endx[2]; ix++) {
    for (iy=start[1]; iy <= endx[1]; iy++) {
        for (iz=start[0]; iz <= endx[0]; iz++) {
        npoints--;
        i = iz + n[0]*(iy + n[1]*ix);

        delta[0] = xs[0]-iz*d[0];
        delta[1] = xs[1]-iy*d[1];
        delta[2] = xs[2]-ix*d[2];

        delta2 = 0.;
        for (j=0; j < 3; j++) {
            if (!plane[2-j]) delta2 += delta[j]*delta[j];
        }

        /* analytical formula (Euclid) */
        ttime[i] = sqrtf(v1*delta2);
        in[i] = SF_IN;

        if ((n[0] > 1 && (iz == start[0] || iz == endx[0])) ||
            (n[1] > 1 && (iy == start[1] || iy == endx[1])) ||
            (n[2] > 1 && (ix == start[2] || ix == endx[2]))) {
            sf_pqueue_insert (ttime+i);
        }
        }
    }
    }

    return npoints;
}

int sf_neighbors_surface(float* vv1  /* slowness [n[0]*n[1]*n[2]] */,
             float* tt0  /* surface traveltime [n[1]*n[2]] */,
             bool forw /* forward or backward continuation */)
/*< initialize the source at the surface >*/
{
    int npoints, i, j, ix, iy;

    /* initialize everywhere */
    for (i=0; i < n[0]*n[1]*n[2]; i++) {
    in[i] = SF_OUT;
    if (forw) {
        ttime[i] = SF_HUGE;
    } else {
        ttime[i] = 0.;
    }
    }

    vv = vv1;

    npoints = (n[0]-1)*n[1]*n[2];

    for (ix=0; ix < n[2]; ix++) {
    for (iy=0; iy < n[1]; iy++) {
        j = iy + n[1]*ix;
        i = j*n[0];

        ttime[i] = tt0[j];
        in[i] = SF_IN;

        if (forw) {
        sf_pqueue_insert (ttime+i);
        } else {
        sf_pqueue_insert2 (ttime+i);
        }
    }
    }

    return npoints;
}

int sf_neighbors_mask(float* vv1  /* slowness [n[0]*n[1]*n[2]] */,
              float* tref /* reference traveltime [n[0]*n[1]*n[2]] */,
              bool* known /* where known [n[0]*n[1]*n[2]] */,
              bool forw   /* forward or backward continuation */)
/*< initialize the source using a mask >*/
{
    int npoints, i, nxy;

    /* save velocity */
    vv = vv1;

    /* total number of points */
    nxy = n[0]*n[1]*n[2];
    npoints = nxy;

    for (i=0; i < nxy; i++) {
    if (known[i]) {
        in[i] = SF_IN;
        ttime[i] = tref[i];

        if (forw) {
        sf_pqueue_insert (ttime+i);
        } else {
        sf_pqueue_insert2 (ttime+i);
        }

        npoints--;
    } else {
        in[i] = SF_OUT;
        if (forw) {
        ttime[i] = SF_HUGE;
        } else {
        ttime[i] = 0.;
        }
    }
    }

    return npoints;
}



void fastmarch_init (int n3,int n2,int n1)
/*< Initialize data dimensions >*/
{
    int maxband;

    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;
    //fprintf(stderr, "fastmarch_init: n1=%d n2=%d n3=%d maxband=%d\n", n1, n2, n3, maxband);
    sf_pqueue_init (100*maxband);
}

void fastmarch (float* time                /* time */,
        float* v                   /* slowness squared */,
        int* in                    /* in/front/out flag */,
        bool* plane                /* if plane source */,
        int   n3,  int n2,  int n1 /* dimensions */,
        float o3,float o2,float o1 /* origin */,
        float d3,float d2,float d1 /* sampling */,
        float s3,float s2,float s1 /* source */,
        int   b3,  int b2,  int b1 /* box around the source */,
        int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;

    n[0] = n1; xs[0] = s1-o1; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = s2-o2; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = s3-o3; b[2] = b3; d[2] = d3;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, time);

    for (npoints =  sf_neighbors_nearsource (xs, b, d, v, plane);
     npoints > 0;
     npoints -= sf_neighbours(i)) {
    /* Pick smallest value in the NarrowBand
       mark as good, decrease points_left */

    /* sf_warning("npoints=%d",npoints); */

        //fprintf(stderr, "npoints=%d\n", npoints);

    p = sf_pqueue_extract();

    if (p == NULL) {
        fprintf(stderr, "%s: heap exausted!\n",__FILE__);
        break;
    }

    i = p - time;

    in[i] = SF_IN;
    }
}

void fastmarch_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}

/* 	$Id$	 */
