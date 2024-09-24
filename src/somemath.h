#ifndef SOMEMATH_H
#define SOMEMATH_H

#include "complex_math.h"

double d_imag(doublecomplex *z);
double d_lg10(double *x);
double d_sign(double a, double b);
int64 i_dnnt(double *x);
int64 i_nint(double *x);
double pow_dd(double *ap, double *bp);
double pow_di(double *ap, int64 *bp);
int64 pow_ii(int64 ap, int64 bp);
double r_lg10(double x);

#endif
