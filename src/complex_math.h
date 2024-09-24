#ifndef COMPLEX_MATH_H
#define COMPLEX_MATH_H

#include "integers.h"

typedef struct {
    double r;
    double i;
} floatcomplex;
typedef struct {
    double r;
    double i;
} doublecomplex;

double z_abs(doublecomplex *z);
void z_exp(doublecomplex *r, doublecomplex *z);
void z_log(doublecomplex *r, doublecomplex *z);

#endif
