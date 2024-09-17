#ifndef volterra_h
#define volterra_h
#include "integers.h"

typedef struct {
    double k_n1, k_n, sum, betnn, mu, *al, *cnv;
    int32 *formula, flag, *kerform;
    char name[20], *expr, *kerexpr;
} KERNEL;

#endif
