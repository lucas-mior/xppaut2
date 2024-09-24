#include <math.h>
#include "auto_f2c.h"
#include "somemath.h"
#include "integers.h"
#include "complex_math.h"

#define log10e 0.43429448190325182765

double
d_imag(doublecomplex *z) {
    return z->i;
}

double
d_lg10(double *x) {
    return log10e*log(*x);
}

double
pow_dd(double *ap, double *bp) {
    return pow(*ap, *bp);
}

double
d_sign(double a, double b) {
    double x;
    x = (a >= 0 ? a : -a);
    return b >= 0 ? x : -x;
}

int64
i_nint(double *x) {
    return (int64)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

int64
i_dnnt(double *x) {
    return (int64)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}

double
r_lg10(double x) {
    return log10e*log(x);
}

double
pow_di(double *ap, int64 *bp) {
    double pow;
    double x;
    int64 n;
    ulong u;

    pow = 1;
    x = *ap;
    n = *bp;

    if (n != 0) {
        if (n < 0) {
            n = -n;
            x = 1 / x;
        }
        for (u = (ulong)n;;) {
            if (u & 01)
                pow *= x;
            if (u >>= 1)
                x *= x;
            else
                break;
        }
    }
    return pow;
}

int64
pow_ii(int64 ap, int64 bp) {
    int64 pow, x, n;
    ulong u;

    x = ap;
    n = bp;

    if (n <= 0) {
        if (n == 0 || x == 1)
            return 1;
        if (x != -1)
            return x == 0 ? 1 / x : 0;
        n = -n;
    }
    u = (ulong)n;
    for (pow = 1;;) {
        if (u & 01)
            pow *= x;
        if (u >>= 1)
            x *= x;
        else
            break;
    }
    return pow;
}
