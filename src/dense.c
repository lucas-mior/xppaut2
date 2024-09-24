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
dense_alloc_mat(int64 N) {
    DenseMat A;

    if (N <= 0)
        return NULL;

    A = xmalloc(sizeof *A);
    if (A == NULL)
        return NULL;

    A->data = dense_alloc2(N);
    if (A->data == NULL) {
        free(A);
        return NULL;
    }

    A->size = N;

    return A;
}

int64 *
dense_alloc_piv(int64 N) {
    if (N <= 0)
        return NULL;

    return xmalloc((usize)N*sizeof(int64));
}

int64
dense_factor(DenseMat A, int64 *p) {
    return dense_gefa(A->data, A->size, p);
}

void
dense_back_solve(DenseMat A, int64 *p, N_Vector b) {
    dense_gesl(A->data, A->size, p, N_VDATA(b));
    return;
}

void
dense_zero(DenseMat A) {
    den_zero(A->data, A->size);
    return;
}

void
dense_copy(DenseMat A, DenseMat B) {
    den_copy(A->data, B->data, A->size);
    return;
}

void
dense_scal(double c, DenseMat A) {
    den_scale(c, A->data, A->size);
    return;
}

void
dense_add_i(DenseMat A) {
    den_add_i(A->data, A->size);
    return;
}

void
dense_free_mat(DenseMat A) {
    den_free(A->data);
    free(A);
    return;
}

void
dense_free_piv(int64 *p) {
    free(p);
    return;
}

void
dense_print(DenseMat A) {
    den_print(A->data, A->size);
    return;
}

double **
dense_alloc2(int64 n) {
    int64 j;
    double **a;

    if (n <= 0)
        return NULL;

    a = xmalloc((usize)n*sizeof(double *));
    if (a == NULL)
        return NULL;

    a[0] = xmalloc((usize)(n*n)*sizeof(*(a[0])));
    if (a[0] == NULL) {
        free(a);
        return NULL;
    }

    for (j = 1; j < n; j++)
        a[j] = a[0] + j*n;

    return a;
}

int64 *
dense_alloc_piv2(int64 n) {
    if (n <= 0)
        return NULL;

    return xmalloc((usize)n*sizeof(int64));
}

int64
dense_gefa(double **a, int64 n, int64 *p) {
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
dense_gesl(double **a, int64 n, int64 *p, double *b) {
    int64 k, l, i;
    double mult;
    double *col_k;

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
den_zero(double **a, int64 n) {
    int64 i;
    int64 j;
    double *col_j;

    for (j = 0; j < n; j++) {
        col_j = a[j];
        for (i = 0; i < n; i++)
            col_j[i] = ZERO;
    }
    return;
}

void
den_copy(double **a, double **b, int64 n) {
    int64 i;
    int64 j;
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
den_scale(double c, double **a, int64 n) {
    int64 i;
    int64 j;
    double *col_j;

    for (j = 0; j < n; j++) {
        col_j = a[j];
        for (i = 0; i < n; i++)
            col_j[i] *= c;
    }
    return;
}

void
den_add_i(double **a, int64 n) {
    int64 i;

    for (i = 0; i < n; i++)
        a[i][i] += ONE;
    return;
}

void
den_free_piv(int64 *p) {
    free(p);
    return;
}

void
den_free(double **a) {
    free(a[0]);
    free(a);
    return;
}

void
den_print(double **a, int64 n) {
    int64 i;
    int64 j;

    ggets_plintf("\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            ggets_plintf("%10g", a[j][i]);
        }
        ggets_plintf("\n");
    }
    ggets_plintf("\n");
}
