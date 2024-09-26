/******************************************************************
 *                                                                *
 * File          : llnlmath.c                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for a C math library.          *
 *                                                                *
 ******************************************************************/

#include <math.h>
#include "functions.h"
#include "integers.h"

#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0

double
llnlmath_unit_roundoff(void) {
    double u;
    volatile double one_plus_u;

    u = ONE;
    one_plus_u = ONE + u;
    while (one_plus_u != ONE) {
        u /= TWO;
        one_plus_u = ONE + u;
    }
    u *= TWO;

    return u;
}

double
llnlmath_rpower_i(double base, int32 exponent) {
    int32 expt;
    double prod;

    prod = ONE;
    expt = ABS(exponent);
    for (int32 i = 1; i <= expt; i++)
        prod *= base;
    if (exponent < 0)
        prod = ONE / prod;
    return prod;
}

double
llnlmath_rpower_r(double base, double exponent) {
    if (base <= ZERO)
        return ZERO;

    return (double)pow((double)base, (double)exponent);
}

double
llnlmath_rsqrt(double x) {
    if (x <= ZERO)
        return ZERO;

    return (double)sqrt((double)x);
}
