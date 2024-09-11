#include "auto_f2c.h"

#ifdef KR_headers
double floor();
int64
i_dnnt(x)
double *x;
#else
#undef abs
#include "math.h"
int64
i_dnnt(double *x)
#endif
{ return (int64)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x)); }
