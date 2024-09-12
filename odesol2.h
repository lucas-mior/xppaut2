#ifndef _odesol2_h_
#define _odesol2_h_
#include "integers.h"

int32 symplect3(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 discrete(double *y, double *tim, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 bak_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 one_bak_step(double *y, double *t, double dt, int32 neq, double *yg,
                   double *yp, double *yp2, double *ytemp, double *errvec,
                   double *jac, int32 *istart);
void one_step_discrete(double *y, double dt, double *yp, int32 neq, double *t);
void one_step_symp(double *y, double h, double *f, int32 n, double *t);
void one_step_euler(double *y, double dt, double *yp, int32 neq, double *t);
void one_step_rk4(double *y, double dt, double *yval[3], int32 neq,
                  double *tim);
void one_step_heun(double *y, double dt, double *yval[2], int32 neq,
                   double *tim);
int32 euler(double *y, double *tim, double dt, int32 nt, int32 neq,
            int32 *istart, double *work);
int32 mod_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 rung_kut(double *y, double *tim, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 adams(double *y, double *tim, double dt, int32 nstep, int32 neq,
            int32 *ist, double *work);
int32 abmpc(double *y, double *t, double dt, int32 neq);
int32 rb23(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
           double *work, int32 *ierr);
int32 rosen(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
            double *work, int32 *ierr);
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy,
                 int32 neq, double eps, double scal);
void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold,
                  int32 n, double eps, double scal);
int32 bandfac(double *a, int32 ml, int32 mr, int32 n);
void bandsol(double *a, double *b, int32 ml, int32 mr, int32 n);

#endif
