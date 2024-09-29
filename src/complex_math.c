#include "math.h"
#include "complex_math.h"

double
z_abs(doublecomplex *z) {
    double temp;
    double real = z->r;
    double imag = z->i;

    if (real < 0) {
        real = -real;
    }
    if (imag < 0) {
        imag = -imag;
    }
    if (imag > real) {
        temp = real;
        real = imag;
        imag = temp;
    }
    if ((real + imag) == real) {
        return real;
    }

    temp = imag / real;
    temp = real*sqrt(1.0 + temp*temp);  // overflow!!
    return temp;
}

void
z_exp(doublecomplex *r, doublecomplex *z) {
    double expx;
    double zi = z->i;

    expx = exp(z->r);
    r->r = expx*cos(zi);
    r->i = expx*sin(zi);
    return;
}

void
z_log(doublecomplex *r, doublecomplex *z) {
    double zi = z->i;
    double zr = z->r;

    r->i = atan2(zi, zr);
    r->r = log(z_abs(z));
    return;
}
