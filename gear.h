#ifndef _gear_h_
#define _gear_h_
#include "integers.h"

void silent_fixpt(double *x, double eps, double err, double big, int32 maxit,
                  int32 n, double *er, double *em, int32 *ierr);
void do_sing(double *x, double eps, double err, double big, int32 maxit,
             int32 n, int32 *ierr, float *stabinfo);
void do_sing_info(double *x, double eps, double err, double big, int32 maxit,
                  int32 n, double *er, double *em, int32 *ierr);

void shoot_this_now(void);
void pr_evec(double *x, double *ev, int32 n, int32 pr, double eval, int32 type);
void get_complex_evec(double *m, double evr, double evm, double *br, double *bm,
                      int32 n, int32 maxit, double err, int32 *ierr);
void get_evec(double *a, double *anew, double *b, double *bp, int32 n,
              int32 maxit, double err, int32 *ipivot, double eval, int32 *ierr);
void eigen(int32 n, double *a, double *ev, double *work, int32 *ierr);
void hqrx(int32 n, int32 low, int32 igh, double *h, double *ev, int32 *ierr);
void orthesx(int32 n, int32 low, int32 igh, double *a, double *ort);
double sign(double x, double y);
int32 imin(int32 x, int32 y);
double amax(double u, double v);
void getjactrans(double *x, double *y, double *yp, double *xp, double eps,
                 double *d, int32 n);
void getjac(double *x, double *y, double *yp, double *xp, double eps,
            double *dermat, int32 n);
void rooter(double *x, double err, double eps, double big, double *work,
            int32 *ierr, int32 maxit, int32 n);
double sqr2(double z);
int32 gear(int32 n, double *t, double tout, double *y, double hmin, double hmax,
           double eps, int32 mf, double *error, int32 *kflag, int32 *jstart,
           double *work, int32 *iwork);
int32 ggear(int32 n, double *t, double tout, double *y, double hmin,
            double hmax, double eps, int32 mf, double *error, int32 *kflag,
            int32 *jstart, double *work, int32 *iwork);
double sgnum(double x, double y);
double Max(double x, double y);
double Min(double x, double y);
void sgefa(double *a, int32 lda, int32 n, int32 *ipvt, int32 *info);
void sgesl(double *a, int32 lda, int32 n, int32 *ipvt, double *b, int32 job);
void saxpy(int32 n, double sa, double *sx, int32 incx, double *sy, int32 incy);
int32 isamax(int32 n, double *sx, int32 incx);
double sdot(int32 n, double *sx, int32 incx, double *sy, int32 incy);
void sscal(int32 n, double sa, double *sx, int32 incx);

#endif
