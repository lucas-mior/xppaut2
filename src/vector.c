/*******************************************************************************
 *
 * File          : vector.c
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL
 * Last Modified : 1 September 1994
 * -----------------------------------------------------------------------------
 *
 * This is the implementation file for a generic VECTOR
 * package. It contains the implementation of the Vector
 * kernels listed in vector.h.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vector.h"
#include "xmalloc.h"
#include "functions.h"
#include "llnlmath.h"

static void vector_copy(Vector x, Vector z);
static void vector_sum(Vector x, Vector y, Vector z);
static void vector_diff(Vector x, Vector y, Vector z);
static void vector_neg(Vector x, Vector z);
static void vector_scaleSum(double c, Vector x, Vector y, Vector z);
static void vector_scaleDiff(double c, Vector x, Vector y, Vector z);
static void vector_lin1(double a, Vector x, Vector y, Vector z);
static void vector_lin2(double a, Vector x, Vector y, Vector z);
static void vector_axpy(double a, Vector x, Vector y);
static void vector_scaleBy(double a, Vector x);

Vector
vector_new(int64 N) {
    Vector v;

    if (N <= 0) {
        return NULL;
    }

    v = xmalloc(sizeof *v);
    if (v == NULL) {
        return NULL;
    }

    v->data = xmalloc((usize)N*sizeof(*(v->data)));
    if (v->data == NULL) {
        free(v);
        return NULL;
    }

    v->length = N;

    return v;
}

void
vector_free(Vector x) {
    free(x->data);
    free(x);
    return;
}

void
vector_linear_sum(double a, Vector x, double b, Vector y, Vector z) {
    int64 N;
    double c, *xd, *yd, *zd;
    Vector v1;
    Vector v2;
    bool test;

    if ((b == 1.0) && (z == y)) {  // BLAS usage: axpy y <- ax+y
        vector_axpy(a, x, y);
        return;
    }

    if ((a == 1.0) && (z == x)) {  // BLAS usage: axpy x <- by+x
        vector_axpy(b, y, x);
        return;
    }

    // Case: a == b == 1.0

    if ((a == 1.0) && (b == 1.0)) {
        vector_sum(x, y, z);
        return;
    }

    // Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0

    if ((test = ((a == 1.0) && (b == -1.0))) || ((a == -1.0) && (b == 1.0))) {
        v1 = test ? y : x;
        v2 = test ? x : y;
        vector_diff(v2, v1, z);
        return;
    }

    /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0
     */
    // if a or b is 0.0, then user should have called N_VScale

    if ((test = (a == 1.0)) || (b == 1.0)) {
        c = test ? b : a;
        v1 = test ? y : x;
        v2 = test ? x : y;
        vector_lin1(c, v1, v2, z);
        return;
    }

    // Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0

    if ((test = (a == -1.0)) || (b == -1.0)) {
        c = test ? b : a;
        v1 = test ? y : x;
        v2 = test ? x : y;
        vector_lin2(c, v1, v2, z);
        return;
    }

    // Case: a == b
    // catches case both a and b are 0.0 - user should have called N_VConst

    if (a == b) {
        vector_scaleSum(a, x, y, z);
        return;
    }

    // Case: a == -b

    if (a == -b) {
        vector_scaleDiff(a, x, y, z);
        return;
    }

    /* Do all cases not handled above:
       (1) a == other, b == 0.0 - user should have called N_VScale
       (2) a == 0.0, b == other - user should have called N_VScale
       (3) a,b == other, a !=b, a != -b */

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int64 i = 0; i < N; i++) {
        *zd++ = a*(*xd++) + b*(*yd++);
    }
    return;
}

void
vector_const(double c, Vector z) {
    int64 N;
    double *zd;

    N = z->length;
    zd = z->data;

    for (int64 i = 0; i < N; i++) {
        *zd++ = c;
    }
    return;
}

void
vector_prod(Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int64 i = 0; i < N; i++) {
        *zd++ = (*xd++)*(*yd++);
    }
    return;
}

void
vector_div(Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int64 i = 0; i < N; i++) {
        *zd++ = (*xd++) / (*yd++);
    }
    return;
}

void
vector_scale(double c, Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    if (z == x) {  // BLAS usage: scale x <- cx
        vector_scaleBy(c, x);
        return;
    }

    if (c == 1.0) {
        vector_copy(x, z);
    } else if (c == -1.0) {
        vector_neg(x, z);
    } else {
        N = x->length;
        xd = x->data;
        zd = z->data;
        for (int32 i = 0; i < N; i++) {
            *zd++ = c*(*xd++);
        }
    }
    return;
}

void
vector_abs(Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int64 i = 0; i < N; i++, xd++, zd++) {
        *zd = ABS(*xd);
    }
    return;
}

void
vector_inv(Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = 1.0 / (*xd++);
    }
    return;
}

void
vector_add_const(Vector x, double b, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = (*xd++) + b;
    }
    return;
}

double
vector_dot_prod(Vector x, Vector y) {
    int64 N;
    double sum = 0.0, *xd, *yd;

    N = x->length;
    xd = x->data;
    yd = y->data;

    for (int32 i = 0; i < N; i++) {
        sum += (*xd++)*(*yd++);
    }

    return sum;
}

double
vector_max_norm(Vector x) {
    int64 N;
    double max = 0.0, *xd;

    N = x->length;
    xd = x->data;

    for (int32 i = 0; i < N; i++, xd++) {
        if (ABS(*xd) > max) {
            max = ABS(*xd);
        }
    }

    return max;
}

double
vector_wrms_norm(Vector x, Vector w) {
    int64 N;
    double sum = 0.0, prodi, *xd, *wd;

    N = x->length;
    xd = x->data;
    wd = w->data;

    for (int32 i = 0; i < N; i++) {
        prodi = (*xd++)*(*wd++);
        sum += prodi*prodi;
    }

    return llnlmath_rsqrt(sum / (double)N);
}

double
vector_min(Vector x) {
    int64 N;
    double min, *xd;

    N = x->length;
    xd = x->data;
    min = xd[0];

    for (int32 i = 0; i < N; i++, xd++) {
        if ((*xd) < min) {
            min = *xd;
        }
    }

    return min;
}

void
vector_compare(double c, Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++, xd++, zd++) {
        *zd = (ABS(*xd) >= c) ? 1.0 : 0.0;
    }
    return;
}

bool
vector_inv_test(Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        if (*xd == 0.0) {
            return false;
        }
        *zd++ = 1.0 / (*xd++);
    }

    return true;
}

void
vector_print(Vector x) {
    int64 N;
    double *xd;

    N = x->length;
    xd = x->data;

    for (int32 i = 0; i < N; i++) {
        ggets_plintf("%g\n", *xd++);
    }

    ggets_plintf("\n");
    return;
}

void
vector_copy(Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = *xd++;
    }
    return;
}

void
vector_sum(Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = (*xd++) + (*yd++);
    }
    return;
}

void
vector_diff(Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = (*xd++) - (*yd++);
    }
    return;
}

void
vector_neg(Vector x, Vector z) {
    int64 N;
    double *xd, *zd;

    N = x->length;
    xd = x->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = -(*xd++);
    }
    return;
}

void
vector_scaleSum(double c, Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = c*((*xd++) + (*yd++));
    }
    return;
}

void
vector_scaleDiff(double c, Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = c*((*xd++) - (*yd++));
    }
    return;
}

void
vector_lin1(double a, Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = a*(*xd++) + (*yd++);
    }
    return;
}

void
vector_lin2(double a, Vector x, Vector y, Vector z) {
    int64 N;
    double *xd, *yd, *zd;

    N = x->length;
    xd = x->data;
    yd = y->data;
    zd = z->data;

    for (int32 i = 0; i < N; i++) {
        *zd++ = a*(*xd++) - (*yd++);
    }
    return;
}

void
vector_axpy(double a, Vector x, Vector y) {
    int64 N;
    double *xd, *yd;

    N = x->length;
    xd = x->data;
    yd = y->data;

    if (a == 1.0) {
        for (int32 i = 0; i < N; i++) {
            *yd++ += (*xd++);
        }
        return;
    }

    if (a == -1.0) {
        for (int32 i = 0; i < N; i++) {
            *yd++ -= (*xd++);
        }
        return;
    }

    for (int32 i = 0; i < N; i++) {
        *yd++ += a*(*xd++);
    }
    return;
}

void
vector_scaleBy(double a, Vector x) {
    int64 N;
    double *xd;

    N = x->length;
    xd = x->data;

    for (int32 i = 0; i < N; i++) {
        *xd++ *= a;
    }
}
