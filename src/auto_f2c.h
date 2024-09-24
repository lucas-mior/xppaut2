/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

        - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "integers.h"

typedef struct {
    double r;
    double i;
} floatcomplex;
typedef struct {
    double r;
    double i;
} doublecomplex;

/* I/O stuff */

#define VOID void

#ifndef abs
#define abs(x) ((x) >= 0 ? (x) : -(x))
#endif
#ifndef fabs
#define fabs(x) ((x) >= 0 ? (x) : -(x))
#endif

#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))
#define dmin(a, b) (double)min(a, b)
#define dmax(a, b) (double)max(a, b)

#define ARRAY2D(array, i, j) array[(i) + (j)*array##_dim1]
#define ARRAY3D(array, i, j, k)                                                \
    array[(i) + ((j) + (k)*array##_dim2)*array##_dim1]

double f__cabs(double, double imag);
double d_imag(doublecomplex *z);
double d_lg10(double *x);
double d_sign(double a, double b);
int64 i_dnnt(double *x);
int64 i_nint(double *x);
double pow_dd(double *ap, double *bp);
double pow_di(double *ap, int64 *bp);
int64 pow_ii(int64 ap, int64 bp);
double r_lg10(double x);
double z_abs(doublecomplex *z);
void z_exp(doublecomplex *r, doublecomplex *z);
void z_log(doublecomplex *r, doublecomplex *z);

#endif
