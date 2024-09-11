#include "auto_f2c.h"

#ifdef KR_headers
double floor();
integer i_dnnt(x)
double *x;
#else
#undef abs
#include "math.h"
integer i_dnnt(double *x)
#endif
{ return (integer)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x)); }
