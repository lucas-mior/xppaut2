#ifndef _do_fit_h_
#define _do_fit_h_
#include "integers.h"

typedef struct {
    char file[25];
    char varlist[25], collist[25];
    char parlist1[25], parlist2[25];
    int32 dim, npars, nvars, npts, maxiter;
    int32 icols[50], ipar[50], ivar[50];
    double tol, eps;
} FITINFO;

void init_fit_info(void);
void get_fit_info(double *y, double *a, double *t0, int32 *flag, double eps,
                  double *yfit, double **yderv, int32 npts, int32 npars, int32 nvars,
                  int32 *ivar, int32 *ipar);
void printem(double **yderv, double *yfit, double *t0, int32 npars, int32 nvars,
             int32 npts);
int32 one_step_int(double *y, double t0, double t1, int32 *istart);
void print_fit_info(void);
void test_fit(void);
int32 run_fit(char *filename, int32 npts, int32 npars, int32 nvars, int32 maxiter,
            int32 ndim, double eps, double tol, int32 *ipar, int32 *ivar, int32 *icols,
            double *y0, double *a, double *yfit);
int32 marlevstep(double *t0, double *y0, double *y, double *sig, double *a,
               int32 npts, int32 nvars, int32 npars, int32 *ivar, int32 *ipar,
               double *covar, double *alpha, double *chisq, double *alambda,
               double *work, double **yderv, double *yfit, double *ochisq,
               int32 ictrl, double eps);
int32 mrqcof(double *t0, double *y0, double *y, double *sig, double *a, int32 npts,
           int32 nvars, int32 npars, int32 *ivar, int32 *ipar, double *alpha,
           double *chisq, double *beta, double **yderv, double *yfit,
           double eps);
int32 get_fit_params(void);
void parse_collist(char *collist, int32 *icols, int32 *n);
void parse_varlist(char *varlist, int32 *ivars, int32 *n);
void parse_parlist(char *parlist, int32 *ipars, int32 *n);

#endif
