#include "auto_f2c.h"

double
pow_di(double *ap, int64 *bp) {
    double pow, x;
    int64 n;
    unsigned long u;

    pow = 1;
    x = *ap;
    n = *bp;

    if (n != 0) {
        if (n < 0) {
            n = -n;
            x = 1 / x;
        }
        for (u = n;;) {
            if (u & 01)
                pow *= x;
            if (u >>= 1)
                x *= x;
            else
                break;
        }
    }
    return (pow);
}
