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
#include "dense.h"
#include "vector.h"
#include "integers.h"
#include "xmalloc.h"

/* Implementation */

DenseMat
dense_alloc_mat(int64 N) {
    DenseMat A;

    if (N <= 0) {
        return NULL;
    }

    A = xmalloc(sizeof *A);
    if (A == NULL) {
        return NULL;
    }

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
    if (N <= 0) {
        return NULL;
    }

    return xmalloc((usize)N*sizeof(int64));
}

int64
dense_factor(DenseMat A, int64 *p) {
    return dense_gefa(A->data, A->size, p);
}

void
dense_back_solve(DenseMat A, int64 *p, Vector b) {
    dense_gesl(A->data, A->size, p, N_VDATA(b));
    return;
}

void
dense_zero(DenseMat A) {
    dense_zero2(A->data, A->size);
    return;
}

void
dense_copy(DenseMat A, DenseMat B) {
    dense_copy2(A->data, B->data, A->size);
    return;
}

void
dense_scal(double c, DenseMat A) {
    dense_scale2(c, A->data, A->size);
    return;
}

void
dense_add_i(DenseMat A) {
    dense_add_i2(A->data, A->size);
    return;
}

void
dense_free_mat(DenseMat A) {
    dense_free2(A->data);
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
    dense_print2(A->data, A->size);
    return;
}

double **
dense_alloc2(int64 n) {
    double **a;

    if (n <= 0) {
        return NULL;
    }

    a = xmalloc((usize)n*sizeof(double *));
    if (a == NULL) {
        return NULL;
    }

    a[0] = xmalloc((usize)(n*n)*sizeof(*(a[0])));
    if (a[0] == NULL) {
        free(a);
        return NULL;
    }

    for (int32 j = 1; j < n; j++) {
        a[j] = a[0] + j*n;
    }

    return a;
}

int64 *
dense_alloc_piv2(int64 n) {
    if (n <= 0) {
        return NULL;
    }

    return xmalloc((usize)n*sizeof(int64));
}

int64
dense_gefa(double **a, int64 n, int64 *p) {
    int64 l;
    double *col_j, *col_k, *diag_k;
    double temp;
    double mult;
    double a_kj;
    bool swap;

    // k = elimination step number

    for (int32 k = 0; k < n - 1; k++, p++) {
        col_k = a[k];
        diag_k = col_k + k;

        // find l = pivot row number

        l = k;
        for (int64 i = k + 1; i < n; i++) {
            if (ABS(col_k[i]) > ABS(col_k[l])) {
                l = i;
            }
        }
        *p = l;

        // check for zero pivot element

        if (col_k[l] == 0.0) {
            return k + 1;
        }

        // swap a(l,k) and a(k,k) if necessary

        if ((swap = (l != k))) {
            temp = col_k[l];
            col_k[l] = *diag_k;
            *diag_k = temp;
        }

        // Scale the elements below the diagonal in
        // column k by -1.0 / a(k,k). After the above swap,
        // a(k,k) holds the pivot element. This scaling
        // stores the pivot row multipliers -a(i,k)/a(k,k)
        // in a(i,k), i=k+1, ..., n-1.

        mult = -1.0 / (*diag_k);
        for (int64 i = k + 1; i < n; i++) {
            col_k[i] *= mult;
        }

        // row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., n-1
        // row k is the pivot row after swapping with row l.
        // The computation is done one column at a time,
        // column j=k+1, ..., n-1.

        for (int32 j = k + 1; j < n; j++) {
            col_j = a[j];
            a_kj = col_j[l];

            // Swap the elements a(k,j) and a(k,l) if l!=k.

            if (swap) {
                col_j[l] = col_j[k];
                col_j[k] = a_kj;
            }

            // a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)
            // a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k)

            if (a_kj != 0.0) {
                for (int64 i = k + 1; i < n; i++) {
                    col_j[i] += a_kj*col_k[i];
                }
            }
        }
    }

    // set the last pivot row to be n-1 and check for a zero pivot

    *p = n - 1;
    if (a[n - 1][n - 1] == 0.0) {
        return n;
    }

    // return 0 to indicate success

    return 0;
}

void
dense_gesl(double **a, int64 n, int64 *p, double *b) {
    int64 l;
    double mult;
    double *col_k;

    // Solve Ly = Pb, store solution y in b

    for (int32 k = 0; k < n - 1; k++) {
        l = p[k];
        mult = b[l];
        if (l != k) {
            b[l] = b[k];
            b[k] = mult;
        }
        col_k = a[k];
        for (int64 i = k + 1; i < n; i++) {
            b[i] += mult*col_k[i];
        }
    }

    // Solve Ux = y, store solution x in b

    for (int64 k = n - 1; k >= 0; k--) {
        col_k = a[k];
        b[k] /= col_k[k];
        mult = -b[k];
        for (int32 i = 0; i < k; i++) {
            b[i] += mult*col_k[i];
        }
    }
    return;
}

void
dense_zero2(double **a, int64 n) {
    double *col_j;

    for (int32 j = 0; j < n; j++) {
        col_j = a[j];
        for (int64 i = 0; i < n; i++) {
            col_j[i] = 0.0;
        }
    }
    return;
}

void
dense_copy2(double **a, double **b, int64 n) {
    double *a_col_j, *b_col_j;

    for (int32 j = 0; j < n; j++) {
        a_col_j = a[j];
        b_col_j = b[j];
        for (int64 i = 0; i < n; i++) {
            b_col_j[i] = a_col_j[i];
        }
    }
    return;
}

void
dense_scale2(double c, double **a, int64 n) {
    double *col_j;

    for (int32 j = 0; j < n; j++) {
        col_j = a[j];
        for (int64 i = 0; i < n; i++) {
            col_j[i] *= c;
        }
    }
    return;
}

void
dense_add_i2(double **a, int64 n) {
    for (int64 i = 0; i < n; i++) {
        a[i][i] += 1.0;
    }
    return;
}

void
dense_free_piv2(int64 *p) {
    free(p);
    return;
}

void
dense_free2(double **a) {
    free(a[0]);
    free(a);
    return;
}

void
dense_print2(double **a, int64 n) {
    ggets_plintf("\n");
    for (int64 i = 0; i < n; i++) {
        for (int32 j = 0; j < n; j++) {
            ggets_plintf("%10g", a[j][i]);
        }
        ggets_plintf("\n");
    }
    ggets_plintf("\n");
}
