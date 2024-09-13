#include <math.h>

#include "integers.h"

int64
i_dnnt(double *x) {
    return (int64)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}
