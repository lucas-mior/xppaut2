#ifndef volterra_h_
#define volterra_h_
#include "integers.h"

#define KN_OLD 1
#define KN 0
typedef struct {
    double k_n1, k_n, sum, betnn, mu, *al, *cnv;
    int32 *formula, flag, *kerform;
    char name[20], *expr, *kerexpr;
} KERNEL;

#endif
