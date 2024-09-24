/******************************************************************
 *                                                                *
 * File          : cvdense.c                                      *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the CVODE dense linear     *
 * solver, CVDENSE.                                               *
 *                                                                *
 ******************************************************************/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "cvdense.h"
#include "cvode.h"
#include "functions.h"
#include "llnltyps.h"
#include "vector.h"
#include "integers.h"

/* Error Messages */

#define CVDENSE_INIT "cv_dense_init-- "

#define MSG_MEM_FAIL CVDENSE_INIT "A memory request failed.\n\n"

/* Other Constants */

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)

/******************************************************************
 *                                                                *
 * Types : CVDenseMemRec, CVDenseMem                              *
 *----------------------------------------------------------------*
 * The type CVDenseMem is pointer to a CVDenseMemRec. This        *
 * structure contains CVDense solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {
    CVDenseJacFn d_jac; /* jac = Jacobian routine to be called    */

    DenseMat d_M; /* M = I - gamma J, gamma = h / l1        */

    int64 *d_pivots; /* pivots = pivot array for PM = LU       */

    DenseMat d_savedJ; /* savedJ = old Jacobian                  */

    int32 d_nstlj; /* nstlj = nst at last Jacobian eval.     */

    int32 d_nje; /* nje = no. of calls to jac              */

    void *d_J_data; /* J_data is passed to jac                */

} CVDenseMemRec, *CVDenseMem;

/* CVDENSE linit, lsetup, lsolve, and lfree routines */

static int32 cv_dense_init(CVodeMem cv_mem, bool *setupNonNull);

static int32 cv_dense_setup(CVodeMem cv_mem, int32 convfail, Vector ypred,
                          Vector fpred, bool *jcurPtr, Vector vtemp1,
                          Vector vtemp2, Vector vtemp3);

static int32 cv_dense_solve(CVodeMem cv_mem, Vector b, Vector ycur,
                          Vector fcur);

static void cv_dense_free(CVodeMem cv_mem);

/*************** cv_dense_dq_jac ****************************************

 This routine generates a dense difference quotient approximation to
 the Jacobian of f(t,y). It assumes that a dense matrix of type
 DenseMat is stored column-wise, and that elements within each column
 are contiguous. The address of the jth column of J is obtained via
 the macro DENSE_COL and an Vector with the jth column as the
 component array is created using N_VMAKE and N_VDATA. Finally, the
 actual computation of the jth column of the Jacobian is done with a
 call to N_VLinearSum.

**********************************************************************/

void
cv_dense_dq_jac(int64 N, DenseMat J, RhsFn f, void *f_data, double tn, Vector y,
             Vector fy, Vector ewt, double h, double uround, void *jac_data,
             int32 *nfePtr, Vector vtemp1, Vector vtemp2, Vector vtemp3) {
    double fnorm, minInc, inc, inc_inv, yjsaved, srur;
    double *y_data;
    double *ewt_data;
    Vector ftemp, jthCol;
    int64 j;

    (void)jac_data;
    (void)vtemp2;
    (void)vtemp3;

    ftemp = vtemp1; /* Rename work vector for use as f vector value */

    /* Obtain pointers to the data for ewt, y */
    ewt_data = N_VDATA(ewt);
    y_data = N_VDATA(y);

    /* Set minimum increment based on uround and norm of f */
    srur = llnlmath_rsqrt(uround);
    fnorm = vector_wrms_norm(fy, ewt);
    minInc = (fnorm != ZERO)
                 ? (MIN_INC_MULT*ABS(h)*uround*(double)N*fnorm)
                 : ONE;

    N_VMAKE(jthCol, NULL, N);

    /* this is the only for loop for 0..N-1 in CVODE */
    for (j = 0; j < N; j++) {

        /* Generate the jth col of J(tn,y) */

        N_VDATA(jthCol) = DENSE_COL(J, j);
        yjsaved = y_data[j];
        inc = MAX(srur*ABS(yjsaved), minInc / ewt_data[j]);
        y_data[j] += inc;
        f(N, tn, y, ftemp, f_data);
        inc_inv = ONE / inc;
        vector_linear_sum(inc_inv, ftemp, -inc_inv, fy, jthCol);
        y_data[j] = yjsaved;
    }

    N_VDISPOSE(jthCol);

    /* Increment counter nfe = *nfePtr */
    *nfePtr += N;
    return;
}

/* Readability Replacements */

#define N (cv_mem->cv_N)
#define lmm (cv_mem->cv_lmm)
#define f (cv_mem->cv_f)
#define f_data (cv_mem->cv_f_data)
#define uround (cv_mem->cv_uround)
#define nst (cv_mem->cv_nst)
#define tn (cv_mem->cv_tn)
#define h (cv_mem->cv_h)
#define gamma (cv_mem->cv_gamma)
#define gammap (cv_mem->cv_gammap)
#define gamrat (cv_mem->cv_gamrat)
#define ewt (cv_mem->cv_ewt)
#define nfe (cv_mem->cv_nfe)
#define errfp (cv_mem->cv_errfp)
#define iopt (cv_mem->cv_iopt)
#define linit (cv_mem->cv_linit)
#define lsetup (cv_mem->cv_lsetup)
#define lsolve (cv_mem->cv_lsolve)
#define lfree (cv_mem->cv_lfree)
#define lmem (cv_mem->cv_lmem)

#define jac (cvdense_mem->d_jac)
#define M (cvdense_mem->d_M)
#define pivots (cvdense_mem->d_pivots)
#define savedJ (cvdense_mem->d_savedJ)
#define nstlj (cvdense_mem->d_nstlj)
#define nje (cvdense_mem->d_nje)
#define J_data (cvdense_mem->d_J_data)

/*************** CVDense *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the dense linear solver module. CVDense sets the
 cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 to be cv_dense_init, cv_dense_setup, cv_dense_solve, and cv_dense_free,
 respectively. It allocates memory for a structure of type
 CVDenseMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure. Finally, it sets d_J_data field in the
 CVDenseMemRec structure to be the input parameter jac_data and the
 d_jac field to be:

 (1) the input parameter djac if djac != NULL or

 (2) cv_dense_dq_jac if djac == NULL.

**********************************************************************/

void
cv_dense(void *cvode_mem, CVDenseJacFn djac, void *jac_data) {
    CVodeMem cv_mem;
    CVDenseMem cvdense_mem;

    /* Return immediately if cvode_mem is NULL */
    cv_mem = (CVodeMem)cvode_mem;
    if (cv_mem == NULL)
        return; /* CVode reports this error */

    /* Set four main function fields in cv_mem */
    linit = cv_dense_init;
    lsetup = cv_dense_setup;
    lsolve = cv_dense_solve;
    lfree = cv_dense_free;

    /* Get memory for CVDenseMemRec */
    lmem = cvdense_mem = xmalloc(sizeof(CVDenseMemRec));
    if (cvdense_mem == NULL)
        return; /* cv_dense_init reports this error */

    /* Set Jacobian routine field to user's djac or cv_dense_dq_jac */
    if (djac == NULL) {
        jac = cv_dense_dq_jac;
    } else {
        jac = djac;
    }
    J_data = jac_data;
    return;
}

/*************** cv_dense_init *****************************************

 This routine initializes remaining memory specific to the dense
 linear solver.  If any memory request fails, all memory previously
 allocated is freed, and an error message printed, before returning.

**********************************************************************/

static int32
cv_dense_init(CVodeMem cv_mem, bool *setupNonNull) {
    CVDenseMem cvdense_mem;

    cvdense_mem = (CVDenseMem)lmem;

    /* Print error message and return if cvdense_mem is NULL */
    if (cvdense_mem == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }

    /* Set flag setupNonNull = TRUE */
    *setupNonNull = TRUE;

    /* Allocate memory for M, savedJ, and pivot array */

    M = dense_alloc_mat(N);
    if (M == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }
    savedJ = dense_alloc_mat(N);
    if (savedJ == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        dense_free_mat(M);
        return LINIT_ERR;
    }
    pivots = dense_alloc_piv(N);
    if (pivots == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        dense_free_mat(M);
        dense_free_mat(savedJ);
        return LINIT_ERR;
    }

    /* Initialize nje and nstlj, and set workspace lengths */

    nje = 0;
    if (iopt != NULL) {
        iopt[DENSE_NJE] = nje;
        iopt[DENSE_LRW] = (int32)(2*N*N);
        iopt[DENSE_LIW] = (int32)N;
    }
    nstlj = 0;

    return LINIT_OK;
}

/*************** cv_dense_setup ****************************************

 This routine does the setup operations for the dense linear solver.
 It makes a decision whether or not to call the Jacobian evaluation
 routine based on various state variables, and if not it uses the
 saved copy.  In any case, it constructs the Newton matrix
 M = I - gamma*J, updates counters, and calls the dense LU
 factorization routine.

**********************************************************************/

static int32
cv_dense_setup(CVodeMem cv_mem, int32 convfail, Vector ypred, Vector fpred,
             bool *jcurPtr, Vector vtemp1, Vector vtemp2, Vector vtemp3) {
    bool jbad;
    bool jok;
    double dgamma;
    int64 ier;
    CVDenseMem cvdense_mem;

    cvdense_mem = (CVDenseMem)lmem;

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */

    dgamma = ABS((gamma / gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
           ((convfail == FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
           (convfail == FAIL_OTHER);
    jok = !jbad;

    if (jok) {
        /* If jok = TRUE, use saved copy of J */
        *jcurPtr = FALSE;
        dense_copy(savedJ, M);
    } else {
        /* If jok = FALSE, call jac routine for new J value */
        nje++;
        if (iopt != NULL)
            iopt[DENSE_NJE] = nje;
        nstlj = nst;
        *jcurPtr = TRUE;
        dense_zero(M);
        jac(N, M, f, f_data, tn, ypred, fpred, ewt, h, uround, J_data, &nfe,
            vtemp1, vtemp2, vtemp3);
        dense_copy(M, savedJ);
    }

    /* Scale and add I to get M = I - gamma*J */
    dense_scal(-gamma, M);
    dense_add_i(M);

    /* Do LU factorization of M */
    ier = dense_factor(M, pivots);

    /* Return 0 if the LU was complete; otherwise return 1 */
    if (ier > 0)
        return 1;
    return 0;
}

/*************** cv_dense_solve ****************************************

 This routine handles the solve operation for the dense linear solver
 by calling the dense backsolve routine.  The returned value is 0.

**********************************************************************/

static int32
cv_dense_solve(CVodeMem cv_mem, Vector b, Vector ycur, Vector fcur) {
    CVDenseMem cvdense_mem;
    (void)ycur;
    (void)fcur;

    cvdense_mem = (CVDenseMem)lmem;

    dense_back_solve(M, pivots, b);

    /* If BDF, scale the correction to account for change in gamma */
    if ((lmm == BDF) && (gamrat != ONE)) {
        vector_scale(TWO / (ONE + gamrat), b, b);
    }

    return 0;
}

/*************** cv_dense_free *****************************************

 This routine frees memory specific to the dense linear solver.

**********************************************************************/

static void
cv_dense_free(CVodeMem cv_mem) {
    CVDenseMem cvdense_mem;

    cvdense_mem = (CVDenseMem)lmem;

    dense_free_mat(M);
    dense_free_mat(savedJ);
    dense_free_piv(pivots);
    free(lmem);
}
