#ifndef _stiff_h_
#define _stiff_h_
#include "integers.h"

/* stiff.c */
void jacobn(double x, double *y, double *dfdx, double *dermat, double eps,
            double *work, int32 n);
int32 adaptive(double *ystart, int32 nvar, double *xs, double x2, double eps,
               double *hguess, double hmin, double *work, int32 *ier,
               double epjac, int32 iflag, int32 *jstart);
int32 gadaptive(double *ystart, int32 nvar, double *xs, double x2, double eps,
                double *hguess, double hmin, double *work, int32 *ier,
                double epjac, int32 iflag, int32 *jstart);
int32 stiff(double y[], double dydx[], int32 n, double *x, double htry,
            double eps, double yscal[], double *hdid, double *hnext,
            double *work, double epjac, int32 *ier);
int32 rkqs(double *y, double *dydx, int32 n, double *x, double htry, double eps,
           double *yscal, double *hdid, double *hnext, double *work,
           int32 *ier);
void rkck(double *y, double *dydx, int32 n, double x, double h, double *yout,
          double *yerr, double *work);

#endif
