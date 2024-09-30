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
#include "llnlmath.h"

double
llnlmath_unit_roundoff(void) {
    double u;
    volatile double one_plus_u;

    u = 1.0;
    one_plus_u = 1.0 + u;
    while (one_plus_u != 1.0) {
        u /= 2.0;
        one_plus_u = 1.0 + u;
    }
    u *= 2.0;

    return u;
}

double
llnlmath_rpower_i(double base, int32 exponent) {
    int32 expt;
    double prod;

    prod = 1.0;
    expt = ABS(exponent);
    for (int32 i = 1; i <= expt; i++) {
        prod *= base;
    }
    if (exponent < 0) {
        prod = 1.0 / prod;
    }
    return prod;
}

double
llnlmath_rpower_r(double base, double exponent) {
    if (base <= 0.0) {
        return 0.0;
    }

    return (double)pow((double)base, (double)exponent);
}

double
llnlmath_rsqrt(double x) {
    if (x <= 0.0) {
        return 0.0;
    }

    return (double)sqrt((double)x);
}
