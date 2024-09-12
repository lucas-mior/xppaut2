#ifndef _flags_h_
#define _flags_h_
#include "integers.h"

/* flags.c */
int32 add_global(char *cond, int32 sign, char *rest);
void show_flags(void);
int32 compile_flags(void);
int32 one_flag_step(double *yold, double *ynew, int32 *istart, double told,
                    double *tnew, int32 neq, double *s);
int32 one_flag_step_symp(double *y, double dt, double *work, int32 neq,
                         double *tim, int32 *istart);
int32 one_flag_step_euler(double *y, double dt, double *work, int32 neq,
                          double *tim, int32 *istart);
int32 one_flag_step_discrete(double *y, double dt, double *work, int32 neq,
                             double *tim, int32 *istart);
int32 one_flag_step_heun(double *y, double dt, double *yval[2], int32 neq,
                         double *tim, int32 *istart);
int32 one_flag_step_rk4(double *y, double dt, double *yval[3], int32 neq,
                        double *tim, int32 *istart);
void printflaginfo(void);
int32 one_flag_step_gear(int32 neq, double *t, double tout, double *y,
                         double hmin, double hmax, double eps, int32 mf,
                         double *error, int32 *kflag, int32 *jstart,
                         double *work, int32 *iwork);
int32 one_flag_step_rosen(double *y, double *tstart, double tfinal,
                          int32 *istart, int32 n, double *work, int32 *ierr);
int32 one_flag_step_dp(int32 *istart, double *y, double *t, int32 n,
                       double tout, double *tol, double *atol, int32 flag,
                       int32 *kflag);
int32 one_flag_step_adap(double *y, int32 neq, double *t, double tout,
                         double eps, double *hguess, double hmin, double *work,
                         int32 *ier, double epjac, int32 iflag, int32 *jstart);
int32 one_flag_step_backeul(double *y, double *t, double dt, int32 neq,
                            double *yg, double *yp, double *yp2, double *ytemp,
                            double *errvec, double *jac, int32 *istart);
int32 one_flag_step_cvode(int32 *command, double *y, double *t, int32 n,
                          double tout, int32 *kflag, double *atol,
                          double *rtol);

#endif
