#include "auto_f2c.h"

#ifdef KR_headers
double pow();
double
pow_dd(ap, bp)
double *ap, *bp;
#else
#undef abs
#include "math.h"
double
pow_dd(double *ap, double *bp)
#endif
{ return (pow(*ap, *bp)); }
