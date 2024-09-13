#include "auto_f2c.h"

#include "math.h"

void
z_exp(doublecomplex *r, doublecomplex *z) {
    double expx, zi = z->i;

    expx = exp(z->r);
    r->r = expx*cos(zi);
    r->i = expx*sin(zi);
    return;
}
