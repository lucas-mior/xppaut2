#include "integers.h"
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include "xpplim.h"

double **storage;
double *WORK;
int32 IWORK[10000];

#define BACKEUL 7
#define VOLTERRA 6
#define STIFF 9
#define GEAR 5
#define RB23 13

void
storage_init_alloc_info(void) {
    xpv.node = NODE + NMarkov;
    xpv.nvec = 0; /* this is just for now */
    xpv.x = xmalloc((usize)(xpv.nvec + xpv.node)*sizeof(*(xpv.x)));
    for (int32 i = xpv.node; i < (xpv.nvec + xpv.node); i++) {
        xpv.x[i] = 0.0;
    }
    return;
}

void
storage_alloc_meth(void) {
    int32 nn = xpv.node + xpv.nvec;
    int32 sz = 30*nn;
    switch (METHOD) {
    case STIFF:
        sz = 2*nn*nn + 13*nn + 100;

        break;
    case GEAR:
        sz = 30*nn + nn*nn + 100;
        break;
    case BACKEUL:
    case VOLTERRA:
        sz = 10*nn + nn*nn + 100;
        break;
    case RB23:
        sz = 12*nn + 100 + nn*nn;
        break;
    default:
        break;
    }
    if (WORK) {
        free(WORK);
    }
    WORK = xmalloc((usize)sz*sizeof(*WORK));
    return;
}

int32
storage_realloc(int32 ncol, int32 nrow) {
    int32 i = 0;
    while ((storage[i] = realloc(storage[i], (usize)nrow*sizeof(*storage))) !=
           NULL) {
        i++;
        if (i == ncol) {
            return 1;
        }
    }
    ggets_err_msg("Cannot allocate sufficient storage");
    return 0;
}

void
storage_init_stor(int32 nrow, int32 ncol) {
    int32 i;
    WORK = NULL;
    storage = xmalloc((MAX_ODE + 1)*sizeof(double *));
    MAXSTOR = nrow;
    storind = 0;
    if (storage != NULL) {
        i = 0;
        while ((storage[i] = xmalloc((usize)nrow*sizeof(*storage))) != NULL) {
            i++;
            if (i == ncol) {
                return;
            }
        }
    }
    ggets_err_msg("Cannot allocate sufficient storage");
    exit(0);
}
