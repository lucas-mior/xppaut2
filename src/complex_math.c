#include "auto_f2c.h"
#include "math.h"

double f__cabs(double, double);

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
