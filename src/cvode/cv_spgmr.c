/******************************************************************
 *                                                                *
 * File          : cvspgmr.c                                      *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the CVODE scaled,          *
 * preconditioned GMRES linear solver, CVSPGMR.                   *
 *                                                                *
 ******************************************************************/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "cv_spgmr.h"
#include "cvode.h"
#include "vector.h"
#include "functions.h"
#include "integers.h"

/* Error Messages */

#define CVSPGMR_INIT "cv_spgmr_init-- "

#define MSG_MEM_FAIL CVSPGMR_INIT "A memory request failed.\n\n"

#define MSG_BAD_PRETYPE_1 CVSPGMR_INIT "pretype=%d illegal.\n"
#define MSG_BAD_PRETYPE_2 "The legal values are PRE_NONE=%d, PRE_LEFT=%d, "
#define MSG_BAD_PRETYPE_3 "PRE_RIGHT=%d, and PRE_BOTH=%d.\n\n"
#define MSG_BAD_PRETYPE MSG_BAD_PRETYPE_1 MSG_BAD_PRETYPE_2 MSG_BAD_PRETYPE_3

#define MSG_PSOLVE_REQ_1 CVSPGMR_INIT "pretype!=NONE, but PSOLVE=NULL is "
#define MSG_PSOLVE_REQ_2 "illegal.\n\n"
#define MSG_PSOLVE_REQ MSG_PSOLVE_REQ_1 MSG_PSOLVE_REQ_2

#define MSG_BAD_GSTYPE_1 CVSPGMR_INIT "gstype=%d illegal.\n"
#define MSG_BAD_GSTYPE_2 "The legal values are MODIFIED_GS=%d and "
#define MSG_BAD_GSTYPE_3 "CLASSICAL_GS=%d.\n\n"
#define MSG_BAD_GSTYPE MSG_BAD_GSTYPE_1 MSG_BAD_GSTYPE_2 MSG_BAD_GSTYPE_3

/* Other Constants */

#define ZERO 0.0
#define ONE 1.0

/******************************************************************
 *                                                                *
 * Types : CVSpgmrMemRec, CVSpgmrMem                              *
 *----------------------------------------------------------------*
 * The type CVSpgmrMem is pointer to a CVSpgmrMemRec. This        *
 * structure contains cv_spgmr solver-specific data.               *
 *                                                                *
 ******************************************************************/

typedef struct {
    int32 g_pretype; /* type of preconditioning                      */
    int32 g_gstype;  /* type of Gram-Schmidt orthogonalization       */
    double g_srqtN;  /* sqrt(N)                                      */
    double g_delt;   /* delt = user specified or DELT_DEFAULT        */
    double g_deltar; /* deltar = delt*tq4                          */
    double g_delta;  /* delta = deltar*sqrtN                       */
    int32 g_maxl;    /* maxl = maximum dimension of the Krylov space */

    int32 g_nstlpre; /* value of nst at the last precond call       */
    int32 g_npe;     /* npe = total number of precond calls         */
    int32 g_nli;     /* nli = total number of linear iterations     */
    int32 g_nps;     /* nps = total number of psolve calls          */
    int32 g_ncfl;    /* ncfl = total number of convergence failures */

    Vector g_ytemp; /* temp vector used by CVAtimesDQ              */
    Vector g_x;     /* temp vector used by cv_spgmr_solve            */
    Vector g_ycur;  /* CVODE current y vector in Newton Iteration  */
    Vector g_fcur;  /* fcur = f(tn, ycur)                          */

    CVSpgmrPrecondFn g_precond; /* precond = user-supplied routine to   */
                                /* compute a preconditioner             */

    CVSpgmrPSolveFn g_psolve; /* psolve = user-supplied routine to    */
                              /* solve preconditioner linear system   */

    void *g_P_data;       /* P_data passed to psolve and precond   */
    SpgmrMem g_spgmr_mem; /* spgmr_mem is memory used by the       */
                          /* generic Spgmr solver                  */

} CVSpgmrMemRec, *CVSpgmrMem;

/* CVSPGMR linit, lsetup, lsolve, and lfree routines */

static int32 cv_spgmr_init(CVodeMem cv_mem, bool *setupNonNull);

static int32 cv_spgmr_setup(CVodeMem cv_mem, int32 convfail, Vector ypred,
                            Vector fpred, bool *jcurPtr, Vector vtemp1,
                            Vector vtemp2, Vector vtemp3);

static int32 cv_spgmr_solve(CVodeMem cv_mem, Vector b, Vector ycur,
                            Vector fcur);

static void cv_spgmr_free(CVodeMem cv_mem);

/* CVSPGMR Atimes and PSolve routines called by generic SPGMR solver */

static int32 cv_spgmr_atimes_dq(void *lin_mem, Vector v, Vector z);

static int32 cv_spgmr_psolve(void *lin_mem, Vector r, Vector z, int32 lr);

/* Readability Replacements */

#define N (cv_mem->cv_N)
#define uround (cv_mem->cv_uround)
#define tq (cv_mem->cv_tq)
#define nst (cv_mem->cv_nst)
#define tn (cv_mem->cv_tn)
#define h (cv_mem->cv_h)
#define gamma (cv_mem->cv_gamma)
#define gammap (cv_mem->cv_gammap)
#define nfe (cv_mem->cv_nfe)
#define f (cv_mem->cv_f)
#define f_data (cv_mem->cv_f_data)
#define ewt (cv_mem->cv_ewt)
#define errfp (cv_mem->cv_errfp)
#define mnewt (cv_mem->cv_mnewt)
#define iopt (cv_mem->cv_iopt)
#define linit (cv_mem->cv_linit)
#define lsetup (cv_mem->cv_lsetup)
#define lsolve (cv_mem->cv_lsolve)
#define lfree (cv_mem->cv_lfree)
#define lmem (cv_mem->cv_lmem)

#define sqrtN (cvspgmr_mem->g_srqtN)
#define ytemp (cvspgmr_mem->g_ytemp)
#define x (cvspgmr_mem->g_x)
#define ycur (cvspgmr_mem->g_ycur)
#define fcur (cvspgmr_mem->g_fcur)
#define delta (cvspgmr_mem->g_delta)
#define deltar (cvspgmr_mem->g_deltar)
#define npe (cvspgmr_mem->g_npe)
#define nli (cvspgmr_mem->g_nli)
#define nps (cvspgmr_mem->g_nps)
#define ncfl (cvspgmr_mem->g_ncfl)
#define nstlpre (cvspgmr_mem->g_nstlpre)
#define spgmr_mem (cvspgmr_mem->g_spgmr_mem)

/*************** cv_spgmr *********************************************

 This routine initializes the memory record and sets various function
 fields specific to the Spgmr linear solver module. cv_spgmr sets the
 cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 to be cv_spgmr_init, cv_spgmr_setup, cv_spgmr_solve, and cv_spgmr_free,
 respectively. It allocates memory for a structure of type
 CVSpgmrMemRec and sets the cv_lmem field in (*cvode_mem) to the
 address of this structure. cv_spgmr sets the following fields in the
 CVSpgmrMemRec structure:

   g_pretype = pretype
   g_maxl    = MIN(N,CVSPGMR_MAXL)  if maxl <= 0
             = maxl                 if maxl > 0
   g_delt    = CVSPGMR_DELT if delt == 0.0
             = delt         if delt != 0.0
   g_P_data  = P_data
   g_precond = precond
   g_psolve  = psolve

**********************************************************************/

void
cv_spgmr(void *cvode_mem, int32 pretype, int32 gstype, int32 maxl, double delt,
         CVSpgmrPrecondFn precond, CVSpgmrPSolveFn psolve, void *P_data) {
    CVodeMem cv_mem;
    CVSpgmrMem cvspgmr_mem;

    /* Return immediately if cvode_mem is NULL */
    cv_mem = (CVodeMem)cvode_mem;
    if (cv_mem == NULL)
        return; /* CVode reports this error */

    /* Set four main function fields in cv_mem */
    linit = cv_spgmr_init;
    lsetup = cv_spgmr_setup;
    lsolve = cv_spgmr_solve;
    lfree = cv_spgmr_free;

    /* Get memory for CVSpgmrMemRec */
    lmem = cvspgmr_mem = xmalloc(sizeof(CVSpgmrMemRec));
    if (cvspgmr_mem == NULL)
        return; /* cv_spgmr_init reports this error */

    /* Set Spgmr parameters that have been passed in call sequence */
    cvspgmr_mem->g_pretype = pretype;
    cvspgmr_mem->g_gstype = gstype;
    cvspgmr_mem->g_maxl = (int32)((maxl <= 0) ? MIN(CVSPGMR_MAXL, N) : maxl);
    cvspgmr_mem->g_delt = (delt == ZERO) ? CVSPGMR_DELT : delt;
    cvspgmr_mem->g_P_data = P_data;
    cvspgmr_mem->g_precond = precond;
    cvspgmr_mem->g_psolve = psolve;
    return;
}

/* Additional readability Replacements */

#define pretype (cvspgmr_mem->g_pretype)
#define gstype (cvspgmr_mem->g_gstype)
#define delt (cvspgmr_mem->g_delt)
#define maxl (cvspgmr_mem->g_maxl)
#define psolve (cvspgmr_mem->g_psolve)
#define precond (cvspgmr_mem->g_precond)
#define P_data (cvspgmr_mem->g_P_data)

/*************** cv_spgmr_init *****************************************

 This routine initializes remaining memory specific to the Spgmr
 linear solver.  If any memory request fails, all memory previously
 allocated is freed, and an error message printed, before returning.

**********************************************************************/

static int32
cv_spgmr_init(CVodeMem cv_mem, bool *setupNonNull) {
    CVSpgmrMem cvspgmr_mem;

    cvspgmr_mem = (CVSpgmrMem)lmem;

    /* Print error message and return if cvspgmr_mem is NULL */
    if (cvspgmr_mem == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }

    /* Check for legal pretype, precond, and psolve */
    if ((pretype != PRE_NONE) && (pretype != PRE_LEFT) &&
        (pretype != PRE_RIGHT) && (pretype != PRE_BOTH)) {
        fprintf(errfp, MSG_BAD_PRETYPE, pretype, PRE_NONE, PRE_LEFT, PRE_RIGHT,
                PRE_BOTH);
        return LINIT_ERR;
    }
    if ((pretype != PRE_NONE) && (psolve == NULL)) {
        fprintf(errfp, MSG_PSOLVE_REQ);
        return LINIT_ERR;
    }

    /* Check for legal gstype */
    if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
        fprintf(errfp, MSG_BAD_GSTYPE, gstype, MODIFIED_GS, CLASSICAL_GS);
        return LINIT_ERR;
    }

    /* Allocate memory for ytemp and x */
    ytemp = vector_new(N);
    if (ytemp == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        return LINIT_ERR;
    }
    x = vector_new(N);
    if (x == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        vector_free(ytemp);
        return LINIT_ERR;
    }

    /* Call SpgmrMalloc to allocate workspace for Spgmr */
    spgmr_mem = spgmr_malloc(N, maxl);
    if (spgmr_mem == NULL) {
        fprintf(errfp, MSG_MEM_FAIL);
        vector_free(ytemp);
        vector_free(x);
        return LINIT_ERR;
    }

    /* Initialize sqrtN and counters, and set workspace lengths */

    sqrtN = llnlmath_rsqrt((double)N);
    npe = nli = nps = ncfl = nstlpre = 0;

    if (iopt != NULL) {
        iopt[SPGMR_NPE] = npe;
        iopt[SPGMR_NLI] = nli;
        iopt[SPGMR_NPS] = nps;
        iopt[SPGMR_NCFL] = ncfl;
        iopt[SPGMR_LRW] = (int32)(N*(maxl + 5) + maxl*(maxl + 4) + 1);
        iopt[SPGMR_LIW] = 0;
    }

    /* Set setupNonNull to true iff there is preconditioning        */
    /* (pretype != PRE_NONE) and there is a preconditioning setup phase */
    /* (precond != NULL)                                            */
    *setupNonNull = (pretype != PRE_NONE) && (precond != NULL);

    return LINIT_OK;
}

/*************** cv_spgmr_setup ****************************************

 This routine does the setup operations for the Spgmr linear solver.
 It makes a decision as to whether or not to signal for re-evaluation
 of Jacobian data in the precond routine, based on various state
 variables, then it calls precond.  If we signal for re-evaluation,
 then we reset jcur = *jcurPtr to true, regardless of the precond output.
 In any case, if jcur == true, we increment npe and save nst in nstlpre.

**********************************************************************/

static int32
cv_spgmr_setup(CVodeMem cv_mem, int32 convfail, Vector ypred, Vector fpred,
               bool *jcurPtr, Vector vtemp1, Vector vtemp2, Vector vtemp3) {
    bool jbad;
    bool jok;
    double dgamma;
    int32 ier;
    CVSpgmrMem cvspgmr_mem;

    cvspgmr_mem = (CVSpgmrMem)lmem;

    /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
    dgamma = ABS((gamma / gammap) - ONE);
    jbad = (nst == 0) || (nst > nstlpre + CVSPGMR_MSBPRE) ||
           ((convfail == FAIL_BAD_J) && (dgamma < CVSPGMR_DGMAX)) ||
           (convfail == FAIL_OTHER);
    *jcurPtr = jbad;
    jok = !jbad;

    /* Call precond routine and possibly reset jcur */
    ier = precond(N, tn, ypred, fpred, jok, jcurPtr, gamma, ewt, h, uround,
                  &nfe, P_data, vtemp1, vtemp2, vtemp3);
    if (jbad)
        *jcurPtr = true;

    /* If jcur = true, increment npe and save nst value */
    if (*jcurPtr) {
        npe++;
        nstlpre = nst;
    }

    /* Set npe, and return the same value ier that precond returned */
    if (iopt != NULL)
        iopt[SPGMR_NPE] = npe;
    return ier;
}

/*************** cv_spgmr_solve ****************************************

 This routine handles the call to the generic SPGMR solver spgmr_solve
 for the solution of the linear system Ax = b.

 If the WRMS norm of b is small, we return x = b (if this is the first
 Newton iteration) or x = 0 (if a later Newton iteration).

 Otherwise, we set the tolerance parameter and initial guess (x = 0),
 call spgmr_solve, and copy the solution x into b.  The x-scaling and
 b-scaling arrays are both equal to ewt, and no restarts are allowed.

 The counters nli, nps, and ncfl are incremented, and the return value
 is set according to the success of spgmr_solve.  The success flag is
 returned if spgmr_solve converged, or if this is the first Newton
 iteration and the residual norm was reduced below its initial value.

**********************************************************************/

static int32
cv_spgmr_solve(CVodeMem cv_mem, Vector b, Vector ynow, Vector fnow) {
    double bnorm;
    double res_norm;
    CVSpgmrMem cvspgmr_mem;
    int32 nli_inc;
    int32 nps_inc;
    int32 ier;

    cvspgmr_mem = (CVSpgmrMem)lmem;

    /* Test norm(b); if small, return x = 0 or x = b */
    deltar = delt*tq[4];
    bnorm = vector_wrms_norm(b, ewt);
    if (bnorm <= deltar) {
        if (mnewt > 0)
            vector_const(ZERO, b);
        return 0;
    }

    /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
    ycur = ynow;
    fcur = fnow;

    /* Set inputs delta and initial guess x = 0 to spgmr_solve */
    delta = deltar*sqrtN;
    vector_const(ZERO, x);

    /* Call spgmr_solve and copy x to b */
    ier = spgmr_solve(spgmr_mem, cv_mem, x, b, pretype, gstype, delta, 0,
                      cv_mem, ewt, ewt, cv_spgmr_atimes_dq, cv_spgmr_psolve,
                      &res_norm, &nli_inc, &nps_inc);
    vector_scale(ONE, x, b);

    /* Increment counters nli, nps, and ncfl */
    nli += nli_inc;
    nps += nps_inc;
    if (iopt != NULL) {
        iopt[SPGMR_NLI] = nli;
        iopt[SPGMR_NPS] = nps;
    }
    if (ier != 0) {
        ncfl++;
        if (iopt != NULL)
            iopt[SPGMR_NCFL] = ncfl;
    }

    /* Set return value to -1, 0, or 1 */
    if (ier < 0)
        return -1;
    if ((ier == SPGMR_SUCCESS) || ((ier == SPGMR_RES_REDUCED) && (mnewt == 0)))
        return 0;
    return 1;
}

/*************** cv_spgmr_free *****************************************

 This routine frees memory specific to the Spgmr linear solver.

**********************************************************************/

static void
cv_spgmr_free(CVodeMem cv_mem) {
    CVSpgmrMem cvspgmr_mem;

    cvspgmr_mem = (CVSpgmrMem)lmem;

    vector_free(ytemp);
    vector_free(x);
    spgmr_free(spgmr_mem);
    free(lmem);
    return;
}

/*************** cv_spgmr_atimes_dq *************************************

 This routine generates the matrix-vector product z = Mv, where
 M = I - gamma*J, by using a difference quotient approximation to
 the product Jv.  The approximation is Jv = rho[f(y + v/rho) - f(y)],
 where rho = (WRMS norm of v), i.e. the WRMS norm of v/rho is 1.

**********************************************************************/

static int32
cv_spgmr_atimes_dq(void *cvode_mem, Vector v, Vector z) {
    double rho;
    CVodeMem cv_mem;
    CVSpgmrMem cvspgmr_mem;

    cv_mem = (CVodeMem)cvode_mem;
    cvspgmr_mem = (CVSpgmrMem)lmem;

    /* If rho = norm(v) is 0, return z = 0 */
    rho = vector_wrms_norm(v, ewt);
    if (rho == ZERO) {
        vector_const(ZERO, z);
        return 0;
    }

    /* Set ytemp = ycur + (1/rho) v */
    vector_linear_sum(ONE / rho, v, ONE, ycur, ytemp);

    /* Set z = f(tn, ytemp) */
    f(N, tn, ytemp, z, f_data);
    nfe++;

    /* Replace z by v - (gamma*rho)(z - fcur) */
    vector_linear_sum(ONE, z, -ONE, fcur, z);
    vector_linear_sum(-gamma*rho, z, ONE, v, z);

    return 0;
}

/*************** cv_spgmr_psolve ***************************************

 This routine interfaces between the generic spgmr_solve routine and
 the user's psolve routine.  It passes to psolve all required state
 information from cvode_mem.  Its return value is the same as that
 returned by psolve. Note that the generic SPGMR solver guarantees
 that cv_spgmr_psolve will not be called in the case in which
 preconditioning is not done. This is the only case in which the
 user's psolve routine is allowed to be NULL.

**********************************************************************/

static int32
cv_spgmr_psolve(void *cvode_mem, Vector r, Vector z, int32 lr) {
    CVodeMem cv_mem;
    CVSpgmrMem cvspgmr_mem;
    int32 ier;

    cv_mem = (CVodeMem)cvode_mem;
    cvspgmr_mem = (CVSpgmrMem)lmem;

    ier = psolve(N, tn, ycur, fcur, ytemp, gamma, ewt, delta, &nfe, r, lr,
                 P_data, z);
    /* This call is counted in nps within the cv_spgmr_solve routine */

    return ier;
}
