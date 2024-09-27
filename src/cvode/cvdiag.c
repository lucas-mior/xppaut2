/******************************************************************
 *                                                                *
 * File          : cvdiag.c                                       *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the CVODE diagonal linear  *
 * solver, CVDIAG.                                                *
 *                                                                *
 ******************************************************************/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "cvdiag.h"
#include "cvode.h"
#include "vector.h"
#include "integers.h"
#include "functions.h"

/* Error Messages */

#define CVDIAG_INIT "cv_diag_init-- "

#define MSG_MEM_FAIL CVDIAG_INIT "A memory request failed.\n\n"

/* Other Constants */

#define FRACT 0.1
#define ONE 1.0

/******************************************************************
 *                                                                *
 * Types : CVDiagMemRec, CVDiagMem                                *
 *----------------------------------------------------------------*
 * The type CVDiagMem is pointer to a CVDiagMemRec. This          *
 * structure contains CVDiag solver-specific data.                *
 *                                                                *
 ******************************************************************/

typedef struct {
    double di_gammasv; /* gammasv = gamma at the last call to setup */
                       /* or solve                                  */

    Vector di_M; /* M = (I - gamma J)^{-1} , gamma = h / l1   */

    Vector di_bit; /* temporary storage vector                  */

    Vector di_bitcomp; /* temporary storage vector                  */

} CVDiagMemRec, *CVDiagMem;

/* CVDIAG linit, lsetup, lsolve, and lfree routines */

static int32 cv_diag_init(CVodeMem cv_mem, bool *setupNonNull);

static int32 cv_diag_setup(CVodeMem cv_mem, int32 convfail, Vector ypred,
                           Vector fpred, bool *jcurPtr, Vector vtemp1,
                           Vector vtemp2, Vector vtemp3);

static int32 cv_diag_solve(CVodeMem cv_mem, Vector b, Vector ycur, Vector fcur);

static void cv_diag_free(CVodeMem cv_mem);

/* Readability Replacements */

#define N (cv_mem->cv_N)
#define f (cv_mem->cv_f)
#define f_data (cv_mem->cv_f_data)
#define uround (cv_mem->cv_uround)
#define tn (cv_mem->cv_tn)
#define h (cv_mem->cv_h)
#define rl1 (cv_mem->cv_rl1)
#define gamma (cv_mem->cv_gamma)
#define ewt (cv_mem->cv_ewt)
#define nfe (cv_mem->cv_nfe)
#define errfp (cv_mem->cv_errfp)
#define iopt (cv_mem->cv_iopt)
#define zn (cv_mem->cv_zn)
#define linit (cv_mem->cv_linit)
#define lsetup (cv_mem->cv_lsetup)
#define lsolve (cv_mem->cv_lsolve)
#define lfree (cv_mem->cv_lfree)
#define lmem (cv_mem->cv_lmem)

#define gammasv (cvdiag_mem->di_gammasv)
#define M (cvdiag_mem->di_M)
#define bit (cvdiag_mem->di_bit)
#define bitcomp (cvdiag_mem->di_bitcomp)

/*************** CVDiag **********************************************

 This routine initializes the memory record and sets various function
 fields specific to the diagonal linear solver module. CVDiag sets the
 cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 to be cv_diag_init, cv_diag_setup, cv_diag_solve, and cv_diag_free,
 respectively. It allocates memory for a structure of type
 CVDiagMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure.

**********************************************************************/

void
cv_diag(void *cvode_mem) {
    CVodeMem cv_mem;
    CVDiagMem cvdiag_mem;

    /* Return immediately if cvode_mem is NULL */
    cv_mem = (CVodeMem)cvode_mem;
    if (cv_mem == NULL) {
        return; /* CVode reports this error */
    }

    /* Set four main function fields in cv_mem */
    linit = cv_diag_init;
    lsetup = cv_diag_setup;
    lsolve = cv_diag_solve;
    lfree = cv_diag_free;

    /* Get memory for CVDiagMemRec */
    lmem = cvdiag_mem = xmalloc(sizeof(CVDiagMemRec));
    if (cvdiag_mem == NULL) {
        return; /* cv_diag_init reports this error */
    }
    return;
}

/*************** cv_diag_init ******************************************

 This routine initializes remaining memory specific to the diagonal
 linear solver.  If any memory request fails, all memory previously
 allocated is freed, and an error message printed, before returning.

**********************************************************************/

static int32
cv_diag_init(CVodeMem cv_mem, bool *setupNonNull) {
    CVDiagMem cvdiag_mem;

    cvdiag_mem = (CVDiagMem)lmem;

    /* Print error message and return if cvdiag_mem is NULL */
    if (cvdiag_mem == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }

    /* Set flag setupNonNull = true */
    *setupNonNull = true;

    /* Allocate memory for M, bit, and bitcomp */

    M = vector_new(N);
    if (M == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }
    bit = vector_new(N);
    if (bit == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        vector_free(M);
        return LINIT_ERR;
    }
    bitcomp = vector_new(N);
    if (bitcomp == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        vector_free(M);
        vector_free(bit);
        return LINIT_ERR;
    }

    /* Set workspace lengths */
    if (iopt != NULL) {
        iopt[DIAG_LRW] = (int32)N*3;
        iopt[DIAG_LIW] = 0;
    }

    return LINIT_OK;
}

/*************** cv_diag_setup *****************************************

 This routine does the setup operations for the diagonal linear
 solver.  It constructs a diagonal approximation to the Newton matrix
 M = I - gamma*J, updates counters, and inverts M.

**********************************************************************/

static int32
cv_diag_setup(CVodeMem cv_mem, int32 convfail, Vector ypred, Vector fpred,
              bool *jcurPtr, Vector vtemp1, Vector vtemp2, Vector vtemp3) {
    double r;
    Vector ftemp;
    Vector y;
    bool invOK;
    CVDiagMem cvdiag_mem;
    (void)convfail;
    (void)vtemp3;

    cvdiag_mem = (CVDiagMem)lmem;

    /* Rename work vectors for use as temporary values of y and f */
    ftemp = vtemp1;
    y = vtemp2;

    r = FRACT*rl1;
    vector_linear_sum(h, fpred, -ONE, zn[1], ftemp);
    vector_linear_sum(r, ftemp, ONE, ypred, y);

    /* Evaluate f at perturbed y */
    f(N, tn, y, M, f_data);
    nfe++;

    /* Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i) */
    vector_linear_sum(ONE, M, -ONE, fpred, M);
    vector_linear_sum(FRACT, ftemp, -h, M, M);
    vector_prod(ftemp, ewt, y);
    /* Protect against deltay_i being at roundoff level */
    vector_compare(uround, y, bit);
    vector_add_const(bit, -ONE, bitcomp);
    vector_prod(ftemp, bit, y);
    vector_linear_sum(FRACT, y, -ONE, bitcomp, y);
    vector_div(M, y, M);
    vector_prod(M, bit, M);
    vector_linear_sum(ONE, M, -ONE, bitcomp, M);

    /* Invert M with test for zero components */
    invOK = vector_inv_test(M, M);
    if (!invOK) {
        return 1;
    }

    /* Set jcur = true, save gamma in gammasv, and return */
    *jcurPtr = true;
    gammasv = gamma;
    return 0;
}

/*************** cv_diag_solve *****************************************

 This routine performs the solve operation for the diagonal linear
 solver.  If necessary it first updates gamma in M = I - gamma*J.

**********************************************************************/

static int32
cv_diag_solve(CVodeMem cv_mem, Vector b, Vector ycur, Vector fcur) {
    bool invOK;
    double r;
    CVDiagMem cvdiag_mem;
    (void)ycur;
    (void)fcur;

    cvdiag_mem = (CVDiagMem)lmem;

    /* If gamma has changed, update factor in M, and save gamma value */

    if (gammasv != gamma) {
        r = gamma / gammasv;
        vector_inv(M, M);
        vector_add_const(M, -ONE, M);
        vector_scale(r, M, M);
        vector_add_const(M, ONE, M);
        invOK = vector_inv_test(M, M);
        if (!invOK) {
            return 1;
        }

        gammasv = gamma;
    }

    /* Apply M-inverse to b */
    vector_prod(b, M, b);
    return 0;
}

/*************** cv_diag_free ******************************************

 This routine frees memory specific to the diagonal linear solver.

**********************************************************************/

static void
cv_diag_free(CVodeMem cv_mem) {
    CVDiagMem cvdiag_mem;

    cvdiag_mem = (CVDiagMem)lmem;

    vector_free(M);
    vector_free(bit);
    vector_free(bitcomp);
    free(lmem);
}
