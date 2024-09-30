/******************************************************************
 *                                                                *
 * File          : iterativ.c                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the iterativ.h header      *
 * file. It contains the implementation of functions that may be  *
 * useful for many different iterative solvers of A x = b.        *
 *                                                                *
 ******************************************************************/

#include "iterativ.h"
#include "llnlmath.h"
#include "vector.h"
#include "integers.h"

#define FACTOR 1000.0

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define ABS(A) ((A > 0) ? (A) : -(A))
#define SQR(A) ((A)*(A))

/************************* modified_gs ***********************************
 This implementation of modified_gs is a slight modification of a previous
 modified Gram-Schmidt routine (called mgs) written by Milo Dorr.
*************************************************************************/

int32
iterativ_modified_gs(Vector *v, double **h, int32 k, int32 p,
                     double *new_vk_norm) {
    int32 k_minus_1;
    int32 i0;
    double new_norm_2;
    double new_product;
    double vk_norm;
    double temp;

    vk_norm = llnlmath_rsqrt(vector_dot_prod(v[k], v[k]));
    k_minus_1 = k - 1;
    i0 = MAX(k - p, 0);

    // Perform modified Gram-Schmidt

    for (int32 i = i0; i < k; i++) {
        h[i][k_minus_1] = vector_dot_prod(v[i], v[k]);
        vector_linear_sum(1.0, v[k], -h[i][k_minus_1], v[i], v[k]);
    }

    // Compute the norm of the new vector at v[k].

    *new_vk_norm = llnlmath_rsqrt(vector_dot_prod(v[k], v[k]));

    /* If the norm of the new vector at v[k] is less than
       FACTOR (== 1000) times unit roundoff times the norm of the
       input vector v[k], then the vector will be reorthogonalized
       in order to ensure that nonorthogonality is not being masked
       by a very small vector length.           */

    temp = FACTOR*vk_norm;
    if ((temp + (*new_vk_norm)) != temp) {
        return 0;
    }

    new_norm_2 = 0.0;

    for (int32 i = i0; i < k; i++) {
        new_product = vector_dot_prod(v[i], v[k]);
        temp = FACTOR*h[i][k_minus_1];
        if ((temp + new_product) == temp) {
            continue;
        }
        h[i][k_minus_1] += new_product;
        vector_linear_sum(1.0, v[k], -new_product, v[i], v[k]);
        new_norm_2 += SQR(new_product);
    }

    if (new_norm_2 != 0.0) {
        new_product = SQR(*new_vk_norm) - new_norm_2;
        *new_vk_norm =
            (new_product > 0.0) ? llnlmath_rsqrt(new_product) : 0.0;
    }

    return 0;
}

/************************ classical_gs ********************************
 This implementation of classical_gs was contributed to by Homer Walker
 and Peter Brown.
**********************************************************************/

int32
iterativ_classical_gs(Vector *v, double **h, int32 k, int32 p,
                      double *new_vk_norm, Vector temp, double *s) {
    int32 k_minus_1;
    int32 i0;
    double vk_norm;

    k_minus_1 = k - 1;

    // Perform Classical Gram-Schmidt

    vk_norm = llnlmath_rsqrt(vector_dot_prod(v[k], v[k]));

    i0 = MAX(k - p, 0);
    for (int32 i = i0; i < k; i++) {
        h[i][k_minus_1] = vector_dot_prod(v[i], v[k]);
    }

    for (int32 i = i0; i < k; i++) {
        vector_linear_sum(1.0, v[k], -h[i][k_minus_1], v[i], v[k]);
    }

    // Compute the norm of the new vector at v[k].

    *new_vk_norm = llnlmath_rsqrt(vector_dot_prod(v[k], v[k]));

    // Reorthogonalize if necessary

    if ((FACTOR*(*new_vk_norm)) < vk_norm) {
        for (int32 i = i0; i < k; i++) {
            s[i] = vector_dot_prod(v[i], v[k]);
        }

        if (i0 < k) {
            vector_scale(s[i0], v[i0], temp);
            h[i0][k_minus_1] += s[i0];
        }
        for (int32 i = i0 + 1; i < k; i++) {
            vector_linear_sum(s[i], v[i], 1.0, temp, temp);
            h[i][k_minus_1] += s[i];
        }
        vector_linear_sum(1.0, v[k], -1.0, temp, v[k]);

        *new_vk_norm = llnlmath_rsqrt(vector_dot_prod(v[k], v[k]));
    }

    return 0;
}

/*************** QRfact **********************************************
 This implementation of QRfact is a slight modification of a previous
 routine (called qrfact) written by Milo Dorr.
**********************************************************************/

int32
iterativ_qr_fact(int32 n, double **h, double *q, int32 job) {
    double c;
    double s;
    double temp1;
    double temp2;
    double temp3;
    int32 i;
    int32 q_ptr;
    int32 n_minus_1;
    int32 code = 0;

    switch (job) {
    case 0:
        // Compute a new factorization of H.
        code = 0;
        for (int32 k = 0; k < n; k++) {
            // Multiply column k by the previous k-1 Givens rotations.
            for (int32 j = 0; j < k - 1; j++) {
                i = 2*j;
                temp1 = h[j][k];
                temp2 = h[j + 1][k];
                c = q[i];
                s = q[i + 1];
                h[j][k] = c*temp1 - s*temp2;
                h[j + 1][k] = s*temp1 + c*temp2;
            }

            // Compute the Givens rotation components c and s
            q_ptr = 2*k;
            temp1 = h[k][k];
            temp2 = h[k + 1][k];
            if (temp2 == 0.0) {
                c = 1.0;
                s = 0.0;
            } else if (ABS(temp2) >= ABS(temp1)) {
                temp3 = temp1 / temp2;
                s = -1.0 / llnlmath_rsqrt(1.0 + SQR(temp3));
                c = -s*temp3;
            } else {
                temp3 = temp2 / temp1;
                c = 1.0 / llnlmath_rsqrt(1.0 + SQR(temp3));
                s = -c*temp3;
            }
            q[q_ptr] = c;
            q[q_ptr + 1] = s;
            if ((h[k][k] = c*temp1 - s*temp2) == 0.0) {
                code = k + 1;
            }
        }
        break;

    default:
        // Update the factored H to which a new column has been added.
        n_minus_1 = n - 1;
        code = 0;

        // Multiply the new column by the previous n-1 Givens rotations.
        for (int32 k = 0; k < n_minus_1; k++) {
            i = 2*k;
            temp1 = h[k][n_minus_1];
            temp2 = h[k + 1][n_minus_1];
            c = q[i];
            s = q[i + 1];
            h[k][n_minus_1] = c*temp1 - s*temp2;
            h[k + 1][n_minus_1] = s*temp1 + c*temp2;
        }

        /* Compute new Givens rotation and multiply it times the last two
           entries in the new column of H.  Note that the second entry of
           this product will be 0, so it is not necessary to compute it. */
        temp1 = h[n_minus_1][n_minus_1];
        temp2 = h[n][n_minus_1];
        if (temp2 == 0.0) {
            c = 1.0;
            s = 0.0;
        } else if (ABS(temp2) >= ABS(temp1)) {
            temp3 = temp1 / temp2;
            s = -1.0 / llnlmath_rsqrt(1.0 + SQR(temp3));
            c = -s*temp3;
        } else {
            temp3 = temp2 / temp1;
            c = 1.0 / llnlmath_rsqrt(1.0 + SQR(temp3));
            s = -c*temp3;
        }
        q_ptr = 2*n_minus_1;
        q[q_ptr] = c;
        q[q_ptr + 1] = s;
        if ((h[n_minus_1][n_minus_1] = c*temp1 - s*temp2) == 0.0) {
            code = n;
        }
    }

    return code;
}

/*************** QRsol ************************************************
 This implementation of QRsol is a slight modification of a previous
 routine (called qrsol) written by Milo Dorr.
**********************************************************************/

int32
iterativ_qr_sol(int32 n, double **h, double *q, double *b) {
    double c;
    double s;
    double temp1;
    double temp2;
    int32 q_ptr;
    int32 code = 0;

    // Compute Q*b.

    for (int32 k = 0; k < n; k++) {
        q_ptr = 2*k;
        c = q[q_ptr];
        s = q[q_ptr + 1];
        temp1 = b[k];
        temp2 = b[k + 1];
        b[k] = c*temp1 - s*temp2;
        b[k + 1] = s*temp1 + c*temp2;
    }

    // Solve  R*x = Q*b.

    for (int32 k = n - 1; k >= 0; k--) {
        if (h[k][k] == 0.0) {
            code = k + 1;
            break;
        }
        b[k] /= h[k][k];
        for (int32 i = 0; i < k; i++) {
            b[i] -= b[k]*h[i][k];
        }
    }

    return code;
}
