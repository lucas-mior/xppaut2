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
#include "llnltyps.h"
#include "integers.h"

#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)

double
UnitRoundoff(void) {
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
RPowerI(double base, int32 exponent) {
    int32 i;
    int32 expt;
    double prod;

    prod = ONE;
    expt = ABS(exponent);
    for (i = 1; i <= expt; i++)
        prod *= base;
    if (exponent < 0)
        prod = ONE / prod;
    return prod;
}

double
RPowerR(double base, double exponent) {
    if (base <= ZERO)
        return ZERO;

    return (double)pow((double)base, (double)exponent);
}

double
RSqrt(double x) {
    if (x <= ZERO)
        return ZERO;

    return (double)sqrt((double)x);
}
