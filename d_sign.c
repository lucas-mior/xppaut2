#include "auto_f2c.h"

#ifdef KR_headers
double
d_sign(a, b)
double a, b;
#else
double
d_sign(double a, double b)
#endif
{
    double x;
    x = (a >= 0 ? a : -a);
    return (b >= 0 ? x : -x);
}
