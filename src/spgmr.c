/******************************************************************
 *                                                                *
 * File          : spgmr.c                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the scaled preconditioned  *
 * GMRES (SPGMR) iterative linear solver.                         *
 *                                                                *
 ******************************************************************/

#include "integers.h"
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include "vector.h"
#include <stdbool.h>

#define ZERO 0.0
#define ONE 1.0

/*************** Private Helper Function Prototype *******************/

static void spgmr_free_vector_array(Vector *A, int32 indMax);

/* Implementation of Spgmr algorithm */

/*************** SpgmrMalloc *****************************************/

SpgmrMem
spgmr_malloc(int64 N, int32 l_max) {
    SpgmrMem mem;
    Vector *V, xcor, vtemp;
    double **Hes, *givens, *yg;

    /* Check the input parameters */

    if ((N <= 0) || (l_max <= 0))
        return NULL;

    /* Get memory for the Krylov basis vectors V[0], ..., V[l_max] */

    V = xmalloc((usize)(l_max + 1)*sizeof(*V));
    if (V == NULL)
        return NULL;

    for (int32 k = 0; k <= l_max; k++) {
        V[k] = vector_new(N);
        if (V[k] == NULL) {
            spgmr_free_vector_array(V, k - 1);
            return NULL;
        }
    }

    /* Get memory for the Hessenberg matrix Hes */

    Hes = xmalloc((usize)(l_max + 1)*sizeof(*Hes));
    if (Hes == NULL) {
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    for (int32 k = 0; k <= l_max; k++) {
        Hes[k] = xmalloc((usize)l_max*sizeof(*(Hes[k])));
        if (Hes[k] == NULL) {
            for (int32 i = 0; i < k; i++)
                free(Hes[i]);
            spgmr_free_vector_array(V, l_max);
            return NULL;
        }
    }

    /* Get memory for Givens rotation components */

    givens = xmalloc(2*(usize)l_max*sizeof(*givens));
    if (givens == NULL) {
        for (int32 i = 0; i <= l_max; i++)
            free(Hes[i]);
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    /* Get memory to hold the correction to z_tilde */

    xcor = vector_new(N);
    if (xcor == NULL) {
        free(givens);
        for (int32 i = 0; i <= l_max; i++)
            free(Hes[i]);
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    /* Get memory to hold SPGMR y and g vectors */

    yg = xmalloc(((usize)l_max + 1)*sizeof(*yg));
    if (yg == NULL) {
        vector_free(xcor);
        free(givens);
        for (int32 i = 0; i <= l_max; i++)
            free(Hes[i]);
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    /* Get an array to hold a temporary vector */

    vtemp = vector_new(N);
    if (vtemp == NULL) {
        free(yg);
        vector_free(xcor);
        free(givens);
        for (int32 i = 0; i <= l_max; i++)
            free(Hes[i]);
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    /* Get memory for an SpgmrMemRec containing SPGMR matrices and vectors */

    mem = xmalloc(sizeof(*mem));
    if (mem == NULL) {
        vector_free(vtemp);
        free(yg);
        vector_free(xcor);
        free(givens);
        for (int32 i = 0; i <= l_max; i++)
            free(Hes[i]);
        spgmr_free_vector_array(V, l_max);
        return NULL;
    }

    /* Set the fields of mem */

    mem->N = N;
    mem->l_max = l_max;
    mem->V = V;
    mem->Hes = Hes;
    mem->givens = givens;
    mem->xcor = xcor;
    mem->yg = yg;
    mem->vtemp = vtemp;

    /* Return the pointer to SPGMR memory */

    return mem;
}

/*************** spgmr_solve ******************************************/

int32
spgmr_solve(SpgmrMem mem, void *A_data, Vector x, Vector b, int32 pretype,
            int32 gstype, double delta, int32 max_restarts, void *P_data,
            Vector sx, Vector sb, ATimesFn atimes, PSolveFn psolve,
            double *res_norm, int32 *nli, int32 *nps) {
    Vector *V, xcor, vtemp;
    double **Hes, *givens, *yg;
    /*double s_r0_norm, beta, rotation_product, r_norm, s_product, rho;*/
    double beta, rotation_product, r_norm, s_product, rho = 0.0;
    bool preOnLeft, preOnRight, scale_x, scale_b, converged;
    int32 j, l, l_plus_1, l_max, krydim = 0, ier, ntries;

    if (mem == NULL)
        return SPGMR_MEM_NULL;

    /* Make local copies of mem variables */

    l_max = mem->l_max;
    V = mem->V;
    Hes = mem->Hes;
    givens = mem->givens;
    xcor = mem->xcor;
    yg = mem->yg;
    vtemp = mem->vtemp;

    *nli = *nps = 0;   /* Initialize counters */
    converged = false; /* Initialize converged flag */

    if (max_restarts < 0)
        max_restarts = 0;

    if ((pretype != PRE_LEFT) && (pretype != PRE_RIGHT) &&
        (pretype != PRE_BOTH))
        pretype = PRE_NONE;

    preOnLeft = ((pretype == PRE_LEFT) || (pretype == PRE_BOTH));
    preOnRight = ((pretype == PRE_RIGHT) || (pretype == PRE_BOTH));
    scale_x = (sx != NULL);
    scale_b = (sb != NULL);

    /* Set vtemp and V[0] to initial (unscaled) residual r_0 = b - A*x_0  */

    if (vector_dot_prod(x, x) == ZERO) {
        vector_scale(ONE, b, vtemp);
    } else {
        if (atimes(A_data, x, vtemp) != 0)
            return SPGMR_ATIMES_FAIL;
        vector_linear_sum(ONE, b, -ONE, vtemp, vtemp);
    }
    vector_scale(ONE, vtemp, V[0]);

    /* Apply b-scaling to vtemp, get L2 norm of sb r_0, and return if small */
    /*
      if (scale_b) vector_prod(sb, vtemp, vtemp);
      s_r0_norm = llnlmath_rsqrt(vector_dot_prod(vtemp, vtemp));
      if (s_r0_norm <= delta) return SPGMR_SUCCESS;
    */
    /* Apply left preconditioner and b-scaling to V[0] = r_0 */

    if (preOnLeft) {
        ier = psolve(P_data, V[0], vtemp, PRE_LEFT);
        (*nps)++;
        if (ier != 0)
            return ((ier < 0) ? SPGMR_PSOLVE_FAIL_UNREC
                              : SPGMR_PSOLVE_FAIL_REC);
    } else {
        vector_scale(ONE, V[0], vtemp);
    }

    if (scale_b) {
        vector_prod(sb, vtemp, V[0]);
    } else {
        vector_scale(ONE, vtemp, V[0]);
    }

    /* Set r_norm = beta to L2 norm of V[0] = sb P1_inv r_0, and
       return if small  */

    *res_norm = r_norm = beta = llnlmath_rsqrt(vector_dot_prod(V[0], V[0]));
    if (r_norm <= delta)
        return SPGMR_SUCCESS;

    /* Set xcor = 0 */

    vector_const(ZERO, xcor);

    /* Begin outer iterations: up to (max_restarts + 1) attempts */

    for (ntries = 0; ntries <= max_restarts; ntries++) {
        /* Initialize the Hessenberg matrix Hes and Givens rotation
           product.  Normalize the initial vector V[0].             */

        for (int32 i = 0; i <= l_max; i++)
            for (j = 0; j < l_max; j++)
                Hes[i][j] = ZERO;

        rotation_product = ONE;

        vector_scale(ONE / r_norm, V[0], V[0]);

        /* Inner loop: generate Krylov sequence and Arnoldi basis */

        for (l = 0; l < l_max; l++) {
            (*nli)++;

            krydim = l_plus_1 = l + 1;

            /* Generate A-tilde V[l], where A-tilde = sb P1_inv A P2_inv sx_inv
             */

            /* Apply x-scaling: vtemp = sx_inv V[l] */
            if (scale_x) {
                vector_div(V[l], sx, vtemp);
            } else {
                vector_scale(ONE, V[l], vtemp);
            }

            /* Apply right precoditioner: vtemp = P2_inv sx_inv V[l] */
            vector_scale(ONE, vtemp, V[l_plus_1]);
            if (preOnRight) {
                ier = psolve(P_data, V[l_plus_1], vtemp, PRE_RIGHT);
                (*nps)++;
                if (ier != 0)
                    return ((ier < 0) ? SPGMR_PSOLVE_FAIL_UNREC
                                      : SPGMR_PSOLVE_FAIL_REC);
            }

            /* Apply A: V[l+1] = A P2_inv sx_inv V[l] */
            if (atimes(A_data, vtemp, V[l_plus_1]) != 0)
                return SPGMR_ATIMES_FAIL;

            /* Apply left preconditioning: vtemp = P1_inv A P2_inv sx_inv V[l]
             */
            if (preOnLeft) {
                ier = psolve(P_data, V[l_plus_1], vtemp, PRE_LEFT);
                (*nps)++;
                if (ier != 0)
                    return ((ier < 0) ? SPGMR_PSOLVE_FAIL_UNREC
                                      : SPGMR_PSOLVE_FAIL_REC);
            } else {
                vector_scale(ONE, V[l_plus_1], vtemp);
            }

            /* Apply b-scaling: V[l+1] = sb P1_inv A P2_inv sx_inv V[l] */
            if (scale_b) {
                vector_prod(sb, vtemp, V[l_plus_1]);
            } else {
                vector_scale(ONE, vtemp, V[l_plus_1]);
            }

            /*  Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde. */

            if (gstype == CLASSICAL_GS) {
                if (iterativ_classical_gs(V, Hes, l_plus_1, l_max,
                                          &(Hes[l_plus_1][l]), vtemp, yg) != 0)
                    return SPGMR_GS_FAIL;
            } else {
                if (iterativ_modified_gs(V, Hes, l_plus_1, l_max,
                                         &(Hes[l_plus_1][l])) != 0)
                    return SPGMR_GS_FAIL;
            }

            /*  Update the QR factorization of Hes  */

            if (iterativ_qr_fact(krydim, Hes, givens, l) != 0)
                return SPGMR_QRFACT_FAIL;

            /*  Update residual norm estimate; break if convergence test passes
             */

            rotation_product *= givens[2*l + 1];

            if ((*res_norm = rho = ABS(rotation_product*r_norm)) <= delta) {
                converged = true;
                break;
            }

            /* Normalize V[l+1] with norm value from the Gram-Schmidt */
            vector_scale(ONE / Hes[l_plus_1][l], V[l_plus_1], V[l_plus_1]);
        }

        /* Inner loop is done.  Compute the new correction vector xcor */

        /* Construct g, then solve for y */
        yg[0] = r_norm;
        for (int32 i = 1; i <= krydim; i++)
            yg[i] = ZERO;
        if (iterativ_qr_sol(krydim, Hes, givens, yg) != 0)
            return SPGMR_QRSOL_FAIL;

        /* Add correction vector V_l y to xcor */
        for (int32 k = 0; k < krydim; k++)
            vector_linear_sum(yg[k], V[k], ONE, xcor, xcor);

        /* If converged, construct the final solution vector x */
        if (converged) {
            /* Apply x-scaling and right precond.: vtemp = P2_inv sx_inv xcor */

            if (scale_x)
                vector_div(xcor, sx, xcor);
            if (preOnRight) {
                ier = psolve(P_data, xcor, vtemp, PRE_RIGHT);
                (*nps)++;
                if (ier != 0)
                    return ((ier < 0) ? SPGMR_PSOLVE_FAIL_UNREC
                                      : SPGMR_PSOLVE_FAIL_REC);
            } else {
                vector_scale(ONE, xcor, vtemp);
            }

            /* Add correction to initial x to get final solution x, and return
             */
            vector_linear_sum(ONE, x, ONE, vtemp, x);

            return SPGMR_SUCCESS;
        }

        /* Not yet converged; if allowed, prepare for restart */

        if (ntries == max_restarts)
            break;

        /* Construct last column of Q in yg */
        s_product = ONE;
        for (int32 i = krydim; i > 0; i--) {
            yg[i] = s_product*givens[2*i - 2];
            s_product *= givens[2*i - 1];
        }
        yg[0] = s_product;

        /* Scale r_norm and yg */
        r_norm *= s_product;
        for (int32 i = 0; i <= krydim; i++)
            yg[i] *= r_norm;
        r_norm = ABS(r_norm);

        /* Multiply yg by V_(krydim+1) to get last residual vector; restart */
        vector_scale(yg[0], V[0], V[0]);
        for (int32 k = 1; k <= krydim; k++)
            vector_linear_sum(yg[k], V[k], ONE, V[0], V[0]);
    }

    /* Failed to converge, even after allowed restarts.
       If the residual norm was reduced below its initial value, compute
       and return x anyway.  Otherwise return failure flag.              */

    if (rho < beta) {
        /* Apply the x-scaling and right precond.: vtemp = P2_inv sx_inv xcor */

        if (scale_x)
            vector_div(xcor, sx, xcor);
        if (preOnRight) {
            ier = psolve(P_data, xcor, vtemp, PRE_RIGHT);
            (*nps)++;
            if (ier != 0)
                return ((ier < 0) ? SPGMR_PSOLVE_FAIL_UNREC
                                  : SPGMR_PSOLVE_FAIL_REC);
        } else {
            vector_scale(ONE, xcor, vtemp);
        }

        /* Add vtemp to initial x to get final solution x, and return */
        vector_linear_sum(ONE, x, ONE, vtemp, x);

        return SPGMR_RES_REDUCED;
    }

    return SPGMR_CONV_FAIL;
}

/*************** SpgmrFree *******************************************/

void
spgmr_free(SpgmrMem mem) {
    int32 l_max;
    double **Hes;

    if (mem == NULL)
        return;

    l_max = mem->l_max;
    Hes = mem->Hes;

    spgmr_free_vector_array(mem->V, l_max);
    for (int32 i = 0; i <= l_max; i++)
        free(Hes[i]);
    free(Hes);
    free(mem->givens);
    vector_free(mem->xcor);
    free(mem->yg);
    vector_free(mem->vtemp);

    free(mem);
    return;
}

/*************** Private Helper Function: FreeVectorArray ************/

static void
spgmr_free_vector_array(Vector *A, int32 indMax) {
    int32 j;

    for (j = 0; j <= indMax; j++)
        vector_free(A[j]);

    free(A);
}
