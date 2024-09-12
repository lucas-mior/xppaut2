#include "auto_f2c.h"

#include "math.h"
extern double f__cabs(double, double);
void
z_log(doublecomplex *r, doublecomplex *z) {
    double zi = z->i, zr = z->r;
    r->i = atan2(zi, zr);
    r->r = log(f__cabs(zr, zi));
    return;
}
