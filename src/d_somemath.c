#include "auto_f2c.h"
#include "math.h"

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
d_sign(double a, double b) {
    double x;
    x = (a >= 0 ? a : -a);
    return b >= 0 ? x : -x;
}
