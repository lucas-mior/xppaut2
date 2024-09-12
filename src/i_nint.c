#include "auto_f2c.h"

#include "math.h"
int64
i_nint(float *x) {
    return (int64)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}
