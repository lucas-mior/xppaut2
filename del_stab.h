#ifndef _del_stab_h_
#define _del_stab_h_

#include "integers.h"
typedef struct {
    double r, i;
} COMPLEX;

/* del_stab.c */
void do_delay_sing(double *x, double eps, double err, double big, int32 maxit,
                   int32 n, int32 *ierr, float *stabinfo);
COMPLEX csum(COMPLEX z, COMPLEX w);
COMPLEX cdif(COMPLEX z, COMPLEX w);
COMPLEX cmlt(COMPLEX z, COMPLEX w);
COMPLEX cdivv(COMPLEX z, COMPLEX w);
COMPLEX cexp2(COMPLEX z);
void switch_rows(COMPLEX *z, int32 i1, int32 i2, int32 n);
COMPLEX rtoc(double x, double y);
void cprintn(COMPLEX z);
void cprint(COMPLEX z);
void cprintarr(COMPLEX *z, int32 n, int32 m);
double c_abs(COMPLEX z);
COMPLEX cdeterm(COMPLEX *z, int32 n);
COMPLEX cxdeterm(COMPLEX *z, int32 n);
void make_z(COMPLEX *z, double *delay, int32 n, int32 m, double *coef,
            COMPLEX lambda);
int32 find_positive_root(double *coef, double *delay, int32 n, int32 m, double rad,
                       double err, double eps, double big, int32 maxit,
                       double *rr);
void process_root(double real, double im);
double get_arg(double *delay, double *coef, int32 m, int32 n, COMPLEX lambda);
int32 test_sign(double old, double new);
int32 plot_args(double *coef, double *delay, int32 n, int32 m, int32 npts, double almax,
              double wmax);

#endif
