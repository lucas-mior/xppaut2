#include "auto_f2c.h"
#include "math.h"
#include "complex_math.h"

double
f__cabs(double real, double imag) {
    double temp;

    if (real < 0)
        real = -real;
    if (imag < 0)
        imag = -imag;
    if (imag > real) {
        temp = real;
        real = imag;
        imag = temp;
    }
    if ((real + imag) == real)
        return real;

    temp = imag / real;
    temp = real*sqrt(1.0 + temp*temp); /*overflow!!*/
    return temp;
}

double
z_abs(doublecomplex *z) {
    return f__cabs(z->r, z->i);
}

void
z_exp(doublecomplex *r, doublecomplex *z) {
    double expx, zi = z->i;

    expx = exp(z->r);
    r->r = expx*cos(zi);
    r->i = expx*sin(zi);
    return;
}

void
z_log(doublecomplex *r, doublecomplex *z) {
    double zi = z->i, zr = z->r;
    r->i = atan2(zi, zr);
    r->r = log(f__cabs(zr, zi));
    return;
}
