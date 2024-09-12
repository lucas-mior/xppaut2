#ifndef _volterra2_h_
#define _volterra2_h_
#include "integers.h"

/* volterra2.c */
double ker_val(int32 in);
void alloc_v_memory(void);
void allocate_volterra(int32 npts, int32 flag);
void re_evaluate_kernels(void);
void alloc_kernels(int32 flag);
void init_sums(double t0, int32 n, double dt, int32 i0, int32 iend,
               int32 ishift);
double alpha1n(double mu, double dt, double t, double t0);
double alpbetjn(double mu, double dt, int32 l);
double betnn(double mu, double dt, double t0, double t);
void get_kn(double *y, double t);
int32 volterra(double *y, double *t, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 volt_step(double *y, double t, double dt, int32 neq, double *yg,
                double *yp, double *yp2, double *ytemp, double *errvec,
                double *jac);

#endif
