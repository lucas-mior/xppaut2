#include "auto_f2c.h"

#ifdef KR_headers
double floor();
integer i_nint(x)
float *x;
#else
#undef abs
#include "math.h"
integer i_nint(float *x)
#endif
{ return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x)); }
