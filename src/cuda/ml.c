#include <math.h>
#include <stdlib.h>
#include <functions.h>

/* ml.c for use in the ode file mlnet.ode
 *
 * without the -O3 -m64 it is slower than XPP!!
 * gcc -fPIC -dynamiclib -O3 -o ML.SO ml.c -m64
 *
 * for linux
 *
 * gcc -fPIC -shared  -O3 -o ML.SO ml.c -m64
 */

double *p;
#define iapp p[0]
#define phi p[1]
#define va p[2]
#define vb p[3]
#define vc p[4]
#define vd p[5]
#define gca p[6]
#define vk p[7]
#define vl p[8]
#define gk p[9]
#define gl p[10]
#define tsyn p[11]
#define gsyn p[12]
#define vsyn p[13]
#define vth p[14]
#define vshp p[15]

/*
double iapp=0.1,phi=.333;
double va=-.01,vb=0.15,vc=0.1,vd=0.145,gca=1;
double vk=-.7,vl=-.5,gk=2.0,gl=.5;
double tsyn=1,gsyn=0.,vsyn=.5;
double vth=0,vshp=.05;

*/

double *sum;

int allocflag = 0;
double
cuda_minf(double v) {
    return .5*(1 + tanh((v - va) / vb));
}
double
cuda_ninf(double v) {
    return .5*(1 + tanh((v - vc) / vd));
}
double
cuda_lamn(double v) {
    return phi*cosh((v - vc) / (2*vd));
}
double
cuda_s_inf(double v) {
    return 1.0 / (1.0 + exp(-(v - vth) / vshp));
}

void
cuda_update_sums(double *s, double *wgt, int n) {
    for (int32 i = 0; i < n; i++) {
        sum[i] = 0.0;
        for (int32 j = 0; j < n; j++) {
            sum[i] += (s[j]*wgt[j + i*n]);
        }
    }
}

void
cuda_update_rhs(double *vp, double *wp, double *sp, double *v, double *w, double *s, int n) {
    for (int32 i = 0; i < n; i++) {
        vp[i] = iapp - gl*(v[i] - vl) - gk*w[i]*(v[i] - vk) -
                gca*cuda_minf(v[i])*(v[i] - 1.0) - gsyn*sum[i]*(v[i] - vsyn);
        wp[i] = cuda_lamn(v[i])*(cuda_ninf(v[i]) - w[i]);
        sp[i] = (cuda_s_inf(v[i]) - s[i]) / tsyn;
    }
}

void
alloc_sum(int n) {
    if (allocflag == 1) {
        return;
    }
    sum = malloc(n*sizeof(*sum));
    allocflag = 1;
}
void
ML(int nn, int ivar, double *par, double *var, double *z[50], double *ydot) {
    double *s, *w, *v;
    double *sp, *wp, *vp;
    double *wgt;
    int n = nn / 3;

    double t = var[0];
    v = var + ivar;
    w = var + n + ivar;
    s = var + 2*n + ivar;
    vp = ydot;
    wp = ydot + n;
    sp = ydot + 2*n;
    wgt = z[0];
    p = par;
    alloc_sum(n);
    cuda_update_sums(s, wgt, n);
    cuda_update_rhs(vp, wp, sp, v, w, s, n);
}
