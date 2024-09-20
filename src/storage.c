#include "integers.h"
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>
#include "xpplim.h"
double **storage;
double *WORK;
extern int32 MAXSTOR, storind;
int32 IWORK[10000];
extern int32 NODE, NMarkov;
extern int32 METHOD;

#define BACKEUL 7
#define VOLTERRA 6
#define STIFF 9
#define GEAR 5
#define RB23 13
typedef struct {
    int32 nvec, node;
    double *x;
} XPPVEC;

extern XPPVEC xpv;

void
init_alloc_info(void) {
    int32 i;
    xpv.node = NODE + NMarkov;
    xpv.nvec = 0; /* this is just for now */
    xpv.x = malloc((usize)(xpv.nvec + xpv.node)*sizeof(*(xpv.x)));
    for (i = xpv.node; i < (xpv.nvec + xpv.node); i++)
        xpv.x[i] = 0.0;
    return;
}

void
alloc_meth(void) {
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
        fprintf(stderr, "%s: METHOD not in switch case.\n", __func__);
        break;
    }
    if (WORK)
        free(WORK);
    WORK = malloc((usize)sz*sizeof(*WORK));
    return;
}

int32
reallocstor(int32 ncol, int32 nrow) {
    int32 i = 0;
    while ((storage[i] = realloc(storage[i], (usize)nrow*sizeof(*storage))) != NULL) {
        i++;
        if (i == ncol)
            return 1;
    }
    err_msg("Cannot allocate sufficient storage");
    return 0;
}

void
init_stor(int32 nrow, int32 ncol) {
    int32 i;
    WORK = NULL;
    storage = malloc((MAX_ODE + 1)*sizeof(double *));
    MAXSTOR = nrow;
    storind = 0;
    if (storage != NULL) {
        i = 0;
        while ((storage[i] = malloc((usize)nrow*sizeof(*storage))) != NULL) {
            i++;
            if (i == ncol)
                return;
        }
    }
    err_msg("Cannot allocate sufficient storage");
    exit(0);
}

void
free_storage(int32 ncol) {
    int32 i;
    for (i = 0; i < ncol; i++)
        free(storage[i]);
    free(storage);
    if (WORK)
        free(WORK);
    return;
}
