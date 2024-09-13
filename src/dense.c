/******************************************************************
 *                                                                *
 * File          : dense.c                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for a generic DENSE linear     *
 * solver package.                                                *
 *                                                                *
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "llnltyps.h"
#include "vector.h"
#include "integers.h"

#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)

/* Implementation */

DenseMat
DenseAllocMat(int64 N) {
    DenseMat A;

    if (N <= 0)
        return NULL;

    A = malloc(sizeof *A);
    if (A == NULL)
        return NULL;

    A->data = denalloc(N);
    if (A->data == NULL) {
        free(A);
        return NULL;
    }

    A->size = N;

    return A;
}

int64 *
DenseAllocPiv(int64 N) {
    if (N <= 0)
        return NULL;

    return malloc(N * sizeof(int64));
}

int64
DenseFactor(DenseMat A, int64 *p) {
    return gefa(A->data, A->size, p);
}

void
DenseBacksolve(DenseMat A, int64 *p, N_Vector b) {
    gesl(A->data, A->size, p, N_VDATA(b));
    return;
}

void
DenseZero(DenseMat A) {
    denzero(A->data, A->size);
    return;
}

void
DenseCopy(DenseMat A, DenseMat B) {
    dencopy(A->data, B->data, A->size);
    return;
}

void
DenseScale(double c, DenseMat A) {
    denscale(c, A->data, A->size);
    return;
}

void
DenseAddI(DenseMat A) {
    denaddI(A->data, A->size);
    return;
}

void
DenseFreeMat(DenseMat A) {
    denfree(A->data);
    free(A);
    return;
}

void
DenseFreePiv(int64 *p) {
    free(p);
    return;
}

void
DensePrint(DenseMat A) {
    denprint(A->data, A->size);
    return;
}

double **
denalloc(int64 n) {
    int64 j;
    double **a;

    if (n <= 0)
        return NULL;

    a = malloc(n*sizeof(double *));
    if (a == NULL)
        return NULL;

    a[0] = malloc(n*n*sizeof(double));
    if (a[0] == NULL) {
        free(a);
        return NULL;
    }

    for (j = 1; j < n; j++)
        a[j] = a[0] + j*n;

    return a;
}

int64 *
denallocpiv(int64 n) {
    if (n <= 0)
        return NULL;

    return malloc(n*sizeof(int64));
}

int64
gefa(double **a, int64 n, int64 *p) {
    int64 i, j, k, l;
    double *col_j, *col_k, *diag_k;
    double temp, mult, a_kj;
    bool swap;

    /* k = elimination step number */

    for (k = 0; k < n - 1; k++, p++) {

        col_k = a[k];
        diag_k = col_k + k;

        /* find l = pivot row number */

        l = k;
        for (i = k + 1; i < n; i++)
            if (ABS(col_k[i]) > ABS(col_k[l]))
                l = i;
        *p = l;

        /* check for zero pivot element */

        if (col_k[l] == ZERO)
            return k + 1;

        /* swap a(l,k) and a(k,k) if necessary */

        if ((swap = (l != k))) {
            temp = col_k[l];
            col_k[l] = *diag_k;
            *diag_k = temp;
        }

        /* Scale the elements below the diagonal in         */
        /* column k by -1.0 / a(k,k). After the above swap, */
        /* a(k,k) holds the pivot element. This scaling     */
        /* stores the pivot row multipliers -a(i,k)/a(k,k)  */
        /* in a(i,k), i=k+1, ..., n-1.                      */

        mult = -ONE / (*diag_k);
        for (i = k + 1; i < n; i++)
            col_k[i] *= mult;

        /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., n-1 */
        /* row k is the pivot row after swapping with row l.      */
        /* The computation is done one column at a time,          */
        /* column j=k+1, ..., n-1.                                */

        for (j = k + 1; j < n; j++) {

            col_j = a[j];
            a_kj = col_j[l];

            /* Swap the elements a(k,j) and a(k,l) if l!=k. */

            if (swap) {
                col_j[l] = col_j[k];
                col_j[k] = a_kj;
            }

            /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
            /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */

            if (a_kj != ZERO) {
                for (i = k + 1; i < n; i++)
                    col_j[i] += a_kj*col_k[i];
            }
        }
    }

    /* set the last pivot row to be n-1 and check for a zero pivot */

    *p = n - 1;
    if (a[n - 1][n - 1] == ZERO)
        return n;

    /* return 0 to indicate success */

    return 0;
}

void
gesl(double **a, int64 n, int64 *p, double *b) {
    int64 k, l, i;
    double mult, *col_k;

    /* Solve Ly = Pb, store solution y in b */

    for (k = 0; k < n - 1; k++) {
        l = p[k];
        mult = b[l];
        if (l != k) {
            b[l] = b[k];
            b[k] = mult;
        }
        col_k = a[k];
        for (i = k + 1; i < n; i++)
            b[i] += mult*col_k[i];
    }

    /* Solve Ux = y, store solution x in b */

    for (k = n - 1; k >= 0; k--) {
        col_k = a[k];
        b[k] /= col_k[k];
        mult = -b[k];
        for (i = 0; i < k; i++)
            b[i] += mult*col_k[i];
    }
    return;
}

void
denzero(double **a, int64 n) {
    int64 i, j;
    double *col_j;

    for (j = 0; j < n; j++) {
        col_j = a[j];
        for (i = 0; i < n; i++)
            col_j[i] = ZERO;
    }
    return;
}

void
dencopy(double **a, double **b, int64 n) {
    int64 i, j;
    double *a_col_j, *b_col_j;

    for (j = 0; j < n; j++) {
        a_col_j = a[j];
        b_col_j = b[j];
        for (i = 0; i < n; i++)
            b_col_j[i] = a_col_j[i];
    }
    return;
}

void
denscale(double c, double **a, int64 n) {
    int64 i, j;
    double *col_j;

    for (j = 0; j < n; j++) {
        col_j = a[j];
        for (i = 0; i < n; i++)
            col_j[i] *= c;
    }
    return;
}

void
denaddI(double **a, int64 n) {
    int64 i;

    for (i = 0; i < n; i++)
        a[i][i] += ONE;
    return;
}

void
denfreepiv(int64 *p) {
    free(p);
    return;
}

void
denfree(double **a) {
    free(a[0]);
    free(a);
    return;
}

void
denprint(double **a, int64 n) {
    int64 i, j;

    plintf("\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            plintf("%10g", a[j][i]);
        }
        plintf("\n");
    }
    plintf("\n");
}
