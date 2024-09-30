/******************************************************************
 *                                                                *
 * File          : cvode.c                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for the main CVODE integrator. *
 * It is independent of the CVODE linear solver in use.           *
 *                                                                *
 ******************************************************************/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "cvode.h"
#include "xmalloc.h"
#include "llnlmath.h"
#include "vector.h"
#include "integers.h"

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define ABS(A) ((A > 0) ? (A) : -(A))

/************************************************************/
/************** BEGIN CVODE Private Constants ***************/
/************************************************************/

#define HALF 0.5     // double 0.5
#define TWELVE 12.0  // double 12.0

/***************************************************************/
/************** BEGIN Default Constants ************************/
/***************************************************************/

#define HMIN_DEFAULT 0.0      // hmin default value
#define HMAX_INV_DEFAULT 0.0  // hmax_inv default value
#define MXHNIL_DEFAULT 10      // mxhnil default value
#define MXSTEP_DEFAULT 2000    // mxstep default value

/***************************************************************/
/*************** END Default Constants *************************/
/***************************************************************/

/***************************************************************/
/************ BEGIN Routine-Specific Constants *****************/
/***************************************************************/

/* CVodeDky */

#define FUZZ_FACTOR 100.0

/* CVHin */

#define HLB_FACTOR 100.0
#define HUB_FACTOR 0.1
#define H_BIAS HALF
#define MAX_ITERS 4

/* CVSet */

#define CORTES 0.1

/* CVStep return values */

#define SUCCESS_STEP 0
#define REP_ERR_FAIL -1
#define REP_CONV_FAIL -2
#define SETUP_FAILED -3
#define SOLVE_FAILED -4

/* CVStep control constants */

#define PREDICT_AGAIN -5
#define DO_ERROR_TEST 1

/* CVStep */

#define THRESH 1.5
#define ETAMX1 10000.0
#define ETAMX2 10.0
#define ETAMX3 10.0
#define ETAMXF 0.2
#define ETAMIN 0.1
#define ETACF 0.25
#define ADDON 0.000001
#define BIAS1 6.0
#define BIAS2 6.0
#define BIAS3 10.0
#define ONEPSM 1.000001

#define SMALL_NST 10  // nst > SMALL_NST => use ETAMX3
#define MXNCF                                                                  \
    10  // max no. of convergence failures during
        // one step try
#define MXNEF                                                                  \
    7  // max no. of error test failures during
       // one step try
#define MXNEF1                                                                 \
    3  // max no. of error test failures before
       // forcing a reduction of order
#define SMALL_NEF                                                              \
    2  // if an error failure occurs and
       // SMALL_NEF <= nef <= MXNEF1, then
       // reset eta =  MIN(eta, ETAMXF)
#define LONG_WAIT                                                              \
    10  // number of steps to wait before
        // considering an order change when
        // q==1 and MXNEF1 error test failures
        // have occurred

/* CVnls return values */

#define SOLVED 0
#define CONV_FAIL -1
#define SETUP_FAIL_UNREC -2
#define SOLVE_FAIL_UNREC -3

/* CVnls input flags */

#define FIRST_CALL 0
#define PREV_CONV_FAIL -1
#define PREV_ERR_FAIL -2

/* CVnls other constants */

#define FUNC_MAXCOR                                                            \
    3  // maximum no. of corrector iterations
       // for iter == FUNCTIONAL
#define NEWT_MAXCOR                                                            \
    3  // maximum no. of corrector iterations
       // for iter == NEWTON

#define CRDOWN                                                                 \
    0.3  // constant used in the estimation of the
         // convergence rate (crate) of the
         // iterates for the nonlinear equation
#define DGMAX                                                                  \
    0.3  // iter == NEWTON, |gamma/gammap-1| > DGMAX
         // => call lsetup

#define RDIV 2.0  // declare divergence if ratio del/delp > RDIV
#define MSBP 20   // max no. of steps between lsetup calls

#define TRY_AGAIN                                                              \
    99  // control constant for CVnlsNewton - should be
        // distinct from CVnls return values

/***************************************************************/
/*************** END Routine-Specific Constants  ***************/
/***************************************************************/

/***************************************************************/
/***************** BEGIN Error Messages ************************/
/***************************************************************/

/* cvode_malloc Error Messages */

#define CVM "cvode_malloc-- "

#define MSG_Y0_NULL CVM "y0=NULL illegal.\n\n"

#define MSG_BAD_N CVM "N=%ld < 1 illegal.\n\n"

#define MSG_BAD_LMM_1 CVM "lmm=%d illegal.\n"
#define MSG_BAD_LMM_2 "The legal values are ADAMS=%d and BDF=%d.\n\n"
#define MSG_BAD_LMM MSG_BAD_LMM_1 MSG_BAD_LMM_2

#define MSG_BAD_ITER_1 CVM "iter=%d illegal.\n"
#define MSG_BAD_ITER_2 "The legal values are FUNCTIONAL=%d "
#define MSG_BAD_ITER_3 "and NEWTON=%d.\n\n"
#define MSG_BAD_ITER MSG_BAD_ITER_1 MSG_BAD_ITER_2 MSG_BAD_ITER_3

#define MSG_BAD_ITOL_1 CVM "itol=%d illegal.\n"
#define MSG_BAD_ITOL_2 "The legal values are SS=%d and SV=%d.\n\n"
#define MSG_BAD_ITOL MSG_BAD_ITOL_1 MSG_BAD_ITOL_2

#define MSG_F_NULL CVM "f=NULL illegal.\n\n"

#define MSG_RELTOL_NULL CVM "reltol=NULL illegal.\n\n"

#define MSG_BAD_RELTOL CVM "*reltol=%g < 0 illegal.\n\n"

#define MSG_ABSTOL_NULL CVM "abstol=NULL illegal.\n\n"

#define MSG_BAD_ABSTOL CVM "Some abstol component < 0.0 illegal.\n\n"

#define MSG_BAD_OPTIN_1 CVM "optIn=%d illegal.\n"
#define MSG_BAD_OPTIN_2 "The legal values are false=%d and true=%d.\n\n"
#define MSG_BAD_OPTIN MSG_BAD_OPTIN_1 MSG_BAD_OPTIN_2

#define MSG_BAD_OPT CVM "optIn=true, but iopt=ropt=NULL.\n\n"

#define MSG_BAD_HMIN_HMAX_1 CVM "Inconsistent step size limits:\n"
#define MSG_BAD_HMIN_HMAX_2 "ropt[CV_HMIN]=%g > ropt[CV_HMAX]=%g.\n\n"
#define MSG_BAD_HMIN_HMAX MSG_BAD_HMIN_HMAX_1 MSG_BAD_HMIN_HMAX_2

#define MSG_MEM_FAIL CVM "A memory request failed.\n\n"

#define MSG_BAD_EWT CVM "Some initial ewt component = 0.0 illegal.\n\n"

/* CVode error messages */

#define CVODE "CVode-- "

#define NO_MEM "cvode_mem=NULL illegal.\n\n"

#define MSG_CVODE_NO_MEM CVODE NO_MEM

#define MSG_LINIT_NULL CVODE "The linear solver's init routine is NULL.\n\n"

#define MSG_LSETUP_NULL CVODE "The linear solver's setup routine is NULL.\n\n"

#define MSG_LSOLVE_NULL CVODE "The linear solver's solve routine is NULL.\n\n"

#define MSG_LFREE_NULL CVODE "The linear solver's free routine is NULL.\n\n"

#define MSG_LINIT_FAIL CVODE "The linear solver's init routine failed.\n\n"

#define MSG_YOUT_NULL CVODE "yout=NULL illegal.\n\n"

#define MSG_T_NULL CVODE "t=NULL illegal.\n\n"

#define MSG_BAD_ITASK_1 CVODE "itask=%d illegal.\nThe legal values are"
#define MSG_BAD_ITASK_2 " NORMAL=%d and ONE_STEP=%d.\n\n"
#define MSG_BAD_ITASK MSG_BAD_ITASK_1 MSG_BAD_ITASK_2

#define MSG_BAD_H0 CVODE "h0=%g and tout-t0=%g inconsistent.\n\n"

#define MSG_BAD_TOUT_1 CVODE "Trouble interpolating at tout = %g.\n"
#define MSG_BAD_TOUT_2 "tout too far back in direction of integration.\n\n"
#define MSG_BAD_TOUT MSG_BAD_TOUT_1 MSG_BAD_TOUT_2

#define MSG_MAX_STEPS_1 CVODE "At t=%g, mxstep=%d steps taken on "
#define MSG_MAX_STEPS_2 "this call before\nreaching tout=%g.\n\n"
#define MSG_MAX_STEPS MSG_MAX_STEPS_1 MSG_MAX_STEPS_2

#define MSG_EWT_NOW_BAD_1 CVODE "At t=%g, "
#define MSG_EWT_NOW_BAD_2 "some ewt component has become <= 0.0.\n\n"
#define MSG_EWT_NOW_BAD MSG_EWT_NOW_BAD_1 MSG_EWT_NOW_BAD_2

#define MSG_TOO_MUCH_ACC CVODE "At t=%g, too much accuracy requested.\n\n"

#define MSG_HNIL_1 CVODE "Warning.. internal t=%g and step size h=%g\n"
#define MSG_HNIL_2 "are such that t + h == t on the next step.\n"
#define MSG_HNIL_3 "The solver will continue anyway.\n\n"
#define MSG_HNIL MSG_HNIL_1 MSG_HNIL_2 MSG_HNIL_3

#define MSG_HNIL_DONE_1 CVODE "The above warning has been issued %d times "
#define MSG_HNIL_DONE_2 "and will not be\nissued again for this problem.\n\n"
#define MSG_HNIL_DONE MSG_HNIL_DONE_1 MSG_HNIL_DONE_2

#define MSG_ERR_FAILS_1 CVODE "At t=%g and step size h=%g, the error test\n"
#define MSG_ERR_FAILS_2 "failed repeatedly or with |h| = hmin.\n\n"
#define MSG_ERR_FAILS MSG_ERR_FAILS_1 MSG_ERR_FAILS_2

#define MSG_CONV_FAILS_1 CVODE "At t=%g and step size h=%g, the corrector\n"
#define MSG_CONV_FAILS_2 "convergence failed repeatedly or "
#define MSG_CONV_FAILS_3 "with |h| = hmin.\n\n"
#define MSG_CONV_FAILS MSG_CONV_FAILS_1 MSG_CONV_FAILS_2 MSG_CONV_FAILS_3

#define MSG_SETUP_FAILED_1 CVODE "At t=%g, the setup routine failed in an "
#define MSG_SETUP_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SETUP_FAILED MSG_SETUP_FAILED_1 MSG_SETUP_FAILED_2

#define MSG_SOLVE_FAILED_1 CVODE "At t=%g, the solve routine failed in an "
#define MSG_SOLVE_FAILED_2 "unrecoverable manner.\n\n"
#define MSG_SOLVE_FAILED MSG_SOLVE_FAILED_1 MSG_SOLVE_FAILED_2

#define MSG_TOO_CLOSE_1 CVODE "tout=%g too close to t0=%g to start"
#define MSG_TOO_CLOSE_2 " integration.\n\n"
#define MSG_TOO_CLOSE MSG_TOO_CLOSE_1 MSG_TOO_CLOSE_2

/* CVodeDky Error Messages */

#define DKY "CVodeDky-- "

#define MSG_DKY_NO_MEM DKY NO_MEM

#define MSG_BAD_K DKY "k=%d illegal.\n\n"

#define MSG_BAD_T_1 DKY "t=%g illegal.\n"
#define MSG_BAD_T_2 "t not in interval tcur-hu=%g to tcur=%g.\n\n"
#define MSG_BAD_T MSG_BAD_T_1 MSG_BAD_T_2

#define MSG_BAD_DKY DKY "dky=NULL illegal.\n\n"

/***************************************************************/
/****************** END Error Messages *************************/
/***************************************************************/

/************************************************************/
/*************** END CVODE Private Constants ****************/
/************************************************************/

/**************************************************************/
/********* BEGIN Private Helper Functions Prototypes **********/
/**************************************************************/

static bool cv_alloc_vectors(CVodeMem cv_mem, int64 neq, int32 maxord);
static void cv_free_vectors(CVodeMem cv_mem, int32 maxord);

static bool cv_ewt_set(CVodeMem cv_mem, double *rtol, void *atol,
                       int32 tol_type, Vector ycur, Vector ewtvec, int64 neq);
static bool cv_ewt_set_ss(CVodeMem cv_mem, double *rtol, double *atol,
                          Vector ycur, Vector ewtvec, int64 neq);
static bool cv_ewt_set_sv(CVodeMem cv_mem, double *rtol, Vector atol,
                          Vector ycur, Vector ewtvec, int64 neq);

static bool cv_hin(CVodeMem cv_mem, double tout);
static double cv_upper_bound_h0(CVodeMem cv_mem, double tdist);
static double cv_ydd_norm(CVodeMem cv_mem, double hg);

static int32 cv_step(CVodeMem cv_mem);

static void cv_adjust_params(CVodeMem cv_mem);
static void cv_adjust_order(CVodeMem cv_mem, int32 deltaq);
static void cv_adjust_adams(CVodeMem cv_mem, int32 deltaq);
static void cv_adjust_bdf(CVodeMem cv_mem, int32 deltaq);
static void cv_increase_bdf(CVodeMem cv_mem);
static void cv_decrese_bdf(CVodeMem cv_mem);

static void cv_rescale(CVodeMem cv_mem);

static void cv_predict(CVodeMem cv_mem);

static void cv_set(CVodeMem cv_mem);
static void cv_set_adams(CVodeMem cv_mem);
static double cv_adams_start(CVodeMem cv_mem, double m[]);
static void cv_adams_finish(CVodeMem cv_mem, double m[], double M[],
                            double hsum);
static double cv_alt_sum(int32 iend, double a[], int32 k);
static void cv_set_bdf(CVodeMem cv_mem);
static void cv_set_tq_bdf(CVodeMem cv_mem, double hsum, double alpha0,
                          double alpha0_hat, double xi_inv, double xistar_inv);

static int32 cv_nls(CVodeMem cv_mem, int32 nflag);
static int32 cv_nls_functional(CVodeMem cv_mem);
static int32 cv_nls_newton(CVodeMem cv_mem, int32 nflag);
static int32 cv_newton_iteration(CVodeMem cv_mem);

static int32 cv_handle_n_flag(CVodeMem cv_mem, int32 *nflagPtr, double saved_t,
                              int32 *ncfPtr);

static void cv_restore(CVodeMem cv_mem, double saved_t);

static bool cv_do_error_test(CVodeMem cv_mem, int32 *nflagPtr, int32 *kflagPtr,
                             double saved_t, int32 *nefPtr, double *dsmPtr);

static void cv_complete_step(CVodeMem cv_mem);

static void cv_prepare_next_step(CVodeMem cv_mem, double dsm);
static void cv_set_eta(CVodeMem cv_mem);
static double cv_compute_etaqm1(CVodeMem cv_mem);
static double cv_compute_etaqp1(CVodeMem cv_mem);
static void cv_choose_eta(CVodeMem cv_mem, double etaqm1, double etaq,
                          double etaqp1);

static int32 cv_handle_failure(CVodeMem cv_mem, int32 kflag);

/**************************************************************/
/********** END Private Helper Functions Prototypes ***********/
/**************************************************************/

/**************************************************************/
/**************** BEGIN Readability Constants *****************/
/**************************************************************/

#define uround (cv_mem->cv_uround)
#define zn (cv_mem->cv_zn)
#define ewt (cv_mem->cv_ewt)
#define y (cv_mem->cv_y)
#define acor (cv_mem->cv_acor)
#define tempv (cv_mem->cv_tempv)
#define ftemp (cv_mem->cv_ftemp)
#define q (cv_mem->cv_q)
#define qprime (cv_mem->cv_qprime)
#define qwait (cv_mem->cv_qwait)
#define L (cv_mem->cv_L)
#define h (cv_mem->cv_h)
#define hprime (cv_mem->cv_hprime)
#define eta (cv_mem->cv_eta)
#define hscale (cv_mem->cv_hscale)
#define tn (cv_mem->cv_tn)
#define tau (cv_mem->cv_tau)
#define tq (cv_mem->cv_tq)
#define l (cv_mem->cv_l)
#define rl1 (cv_mem->cv_rl1)
#define gamma (cv_mem->cv_gamma)
#define gammap (cv_mem->cv_gammap)
#define gamrat (cv_mem->cv_gamrat)
#define crate (cv_mem->cv_crate)
#define acnrm (cv_mem->cv_acnrm)
#define mnewt (cv_mem->cv_mnewt)
#define qmax (cv_mem->cv_qmax)
#define mxstep (cv_mem->cv_mxstep)
#define maxcor (cv_mem->cv_maxcor)
#define mxhnil (cv_mem->cv_mxhnil)
#define hmin (cv_mem->cv_hmin)
#define hmax_inv (cv_mem->cv_hmax_inv)
#define etamax (cv_mem->cv_etamax)
#define nst (cv_mem->cv_nst)
#define nfe (cv_mem->cv_nfe)
#define ncfn (cv_mem->cv_ncfn)
#define netf (cv_mem->cv_netf)
#define nni (cv_mem->cv_nni)
#define nsetups (cv_mem->cv_nsetups)
#define nhnil (cv_mem->cv_nhnil)
#define lrw (cv_mem->cv_lrw)
#define liw (cv_mem->cv_liw)
#define linit (cv_mem->cv_linit)
#define lsetup (cv_mem->cv_lsetup)
#define lsolve (cv_mem->cv_lsolve)
#define lfree (cv_mem->cv_lfree)
#define lmem (cv_mem->cv_lmem)
#define linitOK (cv_mem->cv_linitOK)
#define qu (cv_mem->cv_qu)
#define nstlp (cv_mem->cv_nstlp)
#define hu (cv_mem->cv_hu)
#define saved_tq5 (cv_mem->cv_saved_tq5)
#define jcur (cv_mem->cv_jcur)
#define tolsf (cv_mem->cv_tolsf)
#define setupNonNull (cv_mem->cv_setupNonNull)

/**************************************************************/
/***************** END Readability Constants ******************/
/**************************************************************/

/***************************************************************/
/************* BEGIN CVODE Implementation **********************/
/***************************************************************/

/***************************************************************/
/********* BEGIN Exported Functions Implementation *************/
/***************************************************************/

/******************** cvode_malloc *******************************

 CVode Malloc allocates and initializes memory for a problem. All
 problem specification inputs are checked for errors. If any
 error occurs during initialization, it is reported to the file
 whose file pointer is errfp and NULL is returned. Otherwise, the
 pointer to successfully initialized problem memory is returned.

*****************************************************************/

void *
cvode_malloc(int64 N, RhsFn f, double t0, Vector y0, int32 lmm, int32 iter,
             int32 itol, double *reltol, void *abstol, void *f_data,
             FILE *errfp, bool optIn, int32 iopt[], double ropt[]) {
    bool allocOK;
    bool ioptExists;
    bool roptExists;
    bool neg_abstol;
    bool ewtsetOK;
    int32 maxord;
    CVodeMem cv_mem;
    FILE *fp;

    // Check for legal input parameters

    fp = (errfp == NULL) ? stdout : errfp;

    if (y0 == NULL) {
        fprintf(fp, MSG_Y0_NULL);
        return NULL;
    }

    if (N <= 0) {
        fprintf(fp, MSG_BAD_N, (int64)N);
        return NULL;
    }

    if ((lmm != ADAMS) && (lmm != BDF)) {
        fprintf(fp, MSG_BAD_LMM, lmm, ADAMS, BDF);
        return NULL;
    }

    if ((iter != FUNCTIONAL) && (iter != NEWTON)) {
        fprintf(fp, MSG_BAD_ITER, iter, FUNCTIONAL, NEWTON);
        return NULL;
    }

    if ((itol != SS) && (itol != SV)) {
        fprintf(fp, MSG_BAD_ITOL, itol, SS, SV);
        return NULL;
    }

    if (f == NULL) {
        fprintf(fp, MSG_F_NULL);
        return NULL;
    }

    if (reltol == NULL) {
        fprintf(fp, MSG_RELTOL_NULL);
        return NULL;
    }

    if (*reltol < 0.0) {
        fprintf(fp, MSG_BAD_RELTOL, *reltol);
        return NULL;
    }

    if (abstol == NULL) {
        fprintf(fp, MSG_ABSTOL_NULL);
        return NULL;
    }

    if (itol == SS) {
        neg_abstol = (*((double *)abstol) < 0.0);
    } else {
        neg_abstol = (vector_min((Vector)abstol) < 0.0);
    }
    if (neg_abstol) {
        fprintf(fp, MSG_BAD_ABSTOL);
        return NULL;
    }

    if ((optIn != false) && (optIn != true)) {
        fprintf(fp, MSG_BAD_OPTIN, optIn, false, true);
        return NULL;
    }

    if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
        fprintf(fp, MSG_BAD_OPT);
        return NULL;
    }

    ioptExists = (iopt != NULL);
    roptExists = (ropt != NULL);

    if (optIn && roptExists) {
        if ((ropt[CV_HMAX] > 0.0) && (ropt[CV_HMIN] > ropt[CV_HMAX])) {
            fprintf(fp, MSG_BAD_HMIN_HMAX, ropt[CV_HMIN], ropt[CV_HMAX]);
            return NULL;
        }
    }

    // compute maxord

    maxord = (lmm == ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

    if (optIn && ioptExists) {
        if (iopt[MAXORD] > 0) {
            maxord = MIN(maxord, iopt[MAXORD]);
        }
    }

    cv_mem = xmalloc(sizeof(*cv_mem));
    if (cv_mem == NULL) {
        fprintf(fp, MSG_MEM_FAIL);
        return NULL;
    }

    // Allocate the vectors

    allocOK = cv_alloc_vectors(cv_mem, N, maxord);
    if (!allocOK) {
        fprintf(fp, MSG_MEM_FAIL);
        free(cv_mem);
        return NULL;
    }

    // Set the ewt vector

    ewtsetOK = cv_ewt_set(cv_mem, reltol, abstol, itol, y0, ewt, N);
    if (!ewtsetOK) {
        fprintf(fp, MSG_BAD_EWT);
        cv_free_vectors(cv_mem, maxord);
        free(cv_mem);
        return NULL;
    }

    // All error checking is complete at this point

    // Copy the input parameters into CVODE state

    cv_mem->cv_N = N;  // readability constants defined below cvode_malloc
    cv_mem->cv_f = f;
    cv_mem->cv_f_data = f_data;
    cv_mem->cv_lmm = lmm;
    cv_mem->cv_iter = iter;
    cv_mem->cv_itol = itol;
    cv_mem->cv_reltol = reltol;
    cv_mem->cv_abstol = abstol;
    cv_mem->cv_iopt = iopt;
    cv_mem->cv_ropt = ropt;
    cv_mem->cv_errfp = fp;
    tn = t0;

    // Set step parameters

    q = 1;
    L = 2;
    qwait = L;
    qmax = maxord;
    etamax = ETAMX1;

    // Set uround

    uround = llnlmath_unit_roundoff();

    // Set the linear solver addresses to NULL, linitOK to false

    linit = NULL;
    lsetup = NULL;
    lsolve = NULL;
    lfree = NULL;
    lmem = NULL;
    // We check != NULL later, in CVode and linit, if using NEWTON
    linitOK = false;

    // Initialize the history array zn

    vector_scale(1.0, y0, zn[0]);
    f(N, t0, y0, zn[1], f_data);
    nfe = 1;

    // Handle the remaining optional inputs

    hmin = HMIN_DEFAULT;
    hmax_inv = HMAX_INV_DEFAULT;
    if (optIn && roptExists) {
        if (ropt[CV_HMIN] > 0.0) {
            hmin = ropt[CV_HMIN];
        }
        if (ropt[CV_HMAX] > 0.0) {
            hmax_inv = 1.0 / ropt[CV_HMAX];
        }
    }

    mxhnil = MXHNIL_DEFAULT;
    mxstep = MXSTEP_DEFAULT;
    if (optIn && ioptExists) {
        if (iopt[MXHNIL] > 0) {
            mxhnil = iopt[MXHNIL];
        }
        if (iopt[MXSTEP] > 0) {
            mxstep = iopt[MXSTEP];
        }
    }

    if ((!optIn) && roptExists) {
        ropt[H0] = 0.0;
    }

    // Set maxcor

    maxcor = (iter == NEWTON) ? NEWT_MAXCOR : FUNC_MAXCOR;

    // Initialize all the counters

    nst = ncfn = netf = nni = nsetups = nhnil = nstlp = 0;

    // Initialize all other vars corresponding to optional outputs

    qu = 0;
    hu = 0.0;
    tolsf = 1.0;

    // Initialize optional output locations in iopt, ropt

    if (ioptExists) {
        iopt[NST] = iopt[NFE] = iopt[NSETUPS] = iopt[NNI] = 0;
        iopt[NCFN] = iopt[NETF] = 0;
        iopt[QU] = qu;
        iopt[QCUR] = 0;
        iopt[LENRW] = lrw;
        iopt[LENIW] = liw;
    }

    if (roptExists) {
        ropt[HU] = hu;
        ropt[HCUR] = 0.0;
        ropt[TCUR] = t0;
        ropt[TOLSF] = tolsf;
    }

    // Problem has been successfully initialized

    return (void *)cv_mem;
}

/**************************************************************/
/************** BEGIN More Readability Constants **************/
/**************************************************************/

#define N (cv_mem->cv_N)
#define f (cv_mem->cv_f)
#define f_data (cv_mem->cv_f_data)
#define lmm (cv_mem->cv_lmm)
#define iter (cv_mem->cv_iter)
#define itol (cv_mem->cv_itol)
#define reltol (cv_mem->cv_reltol)
#define abstol (cv_mem->cv_abstol)
#define iopt (cv_mem->cv_iopt)
#define ropt (cv_mem->cv_ropt)
#define errfp (cv_mem->cv_errfp)

/**************************************************************/
/*************** END More Readability Constants ***************/
/**************************************************************/

/********************* CVode ****************************************

 This routine is the main driver of the CVODE package.

 It integrates over a time interval defined by the user, by calling
 CVStep to do internal time steps.

 The first time that CVode is called for a successfully initialized
 problem, it computes a tentative initial step size h.

 CVode supports two modes, specified by itask: NORMAL and ONE_STEP.
 In the NORMAL mode, the solver steps until it reaches or passes tout
 and then interpolates to obtain y(tout).
 In the ONE_STEP mode, it takes one internal step and returns.

********************************************************************/

int32
CVode(void *cvode_mem, double tout, Vector yout, double *t, int32 itask) {
    int32 nstloc;
    int32 kflag;
    int32 istate;
    int32 next_q;
    int32 ier;
    double rh;
    double next_h;
    bool hOK;
    bool ewtsetOK;
    CVodeMem cv_mem;

    // Check for legal inputs in all cases

    cv_mem = (CVodeMem)cvode_mem;
    if (cvode_mem == NULL) {
        fprintf(stdout, MSG_CVODE_NO_MEM);
        return CVODE_NO_MEM;
    }

    if ((y = yout) == NULL) {
        fprintf(errfp, MSG_YOUT_NULL);
        return ILL_INPUT;
    }

    if (t == NULL) {
        fprintf(errfp, MSG_T_NULL);
        return ILL_INPUT;
    }
    *t = tn;

    if ((itask != NORMAL) && (itask != ONE_STEP)) {
        fprintf(errfp, MSG_BAD_ITASK, itask, NORMAL, ONE_STEP);
        return ILL_INPUT;
    }

    // On first call, check solver functions and call linit function

    if (nst == 0) {
        if (iter == NEWTON) {
            if (linit == NULL) {
                fprintf(errfp, MSG_LINIT_NULL);
                return ILL_INPUT;
            }
            if (lsetup == NULL) {
                fprintf(errfp, MSG_LSETUP_NULL);
                return ILL_INPUT;
            }
            if (lsolve == NULL) {
                fprintf(errfp, MSG_LSOLVE_NULL);
                return ILL_INPUT;
            }
            if (lfree == NULL) {
                fprintf(errfp, MSG_LFREE_NULL);
                return ILL_INPUT;
            }
            linitOK = (linit(cv_mem, &(setupNonNull)) == LINIT_OK);
            if (!linitOK) {
                fprintf(errfp, MSG_LINIT_FAIL);
                return ILL_INPUT;
            }
        }

        // On first call, set initial h (from H0 or CVHin) and scale zn[1]

        h = 0.0;
        if (ropt != NULL) {
            h = ropt[H0];
        }
        if ((h != 0.0) && ((tout - tn)*h < 0.0)) {
            fprintf(errfp, MSG_BAD_H0, h, tout - tn);
            return ILL_INPUT;
        }
        if (h == 0.0) {
            hOK = cv_hin(cv_mem, tout);
            if (!hOK) {
                fprintf(errfp, MSG_TOO_CLOSE, tout, tn);
                return ILL_INPUT;
            }
        }
        rh = ABS(h)*hmax_inv;
        if (rh > 1.0) {
            h /= rh;
        }
        if (ABS(h) < hmin) {
            h *= hmin / ABS(h);
        }
        hscale = h;
        vector_scale(h, zn[1], zn[1]);
    }

    // If not the first call, check if tout already reached

    if ((itask == NORMAL) && (nst > 0) && ((tn - tout)*h >= 0.0)) {
        *t = tout;
        ier = cvode_dky(cv_mem, tout, 0, yout);
        if (ier != OKAY) {  // ier must be == BAD_T
            fprintf(errfp, MSG_BAD_TOUT, tout);
            return ILL_INPUT;
        }
        return SUCCESS;
    }

    // Looping point for internal steps

    nstloc = 0;
    while (true) {
        next_h = h;
        next_q = q;

        // Reset and check ewt

        if (nst > 0) {
            ewtsetOK = cv_ewt_set(cv_mem, reltol, abstol, itol, zn[0], ewt, N);
            if (!ewtsetOK) {
                fprintf(errfp, MSG_EWT_NOW_BAD, tn);
                istate = ILL_INPUT;
                *t = tn;
                vector_scale(1.0, zn[0], yout);
                break;
            }
        }

        // Check for too many steps

        if (nstloc >= mxstep) {
            fprintf(errfp, MSG_MAX_STEPS, tn, mxstep, tout);
            istate = TOO_MUCH_WORK;
            *t = tn;
            vector_scale(1.0, zn[0], yout);
            break;
        }

        // Check for too much accuracy requested

        if ((tolsf = uround*vector_wrms_norm(zn[0], ewt)) > 1.0) {
            fprintf(errfp, MSG_TOO_MUCH_ACC, tn);
            istate = TOO_MUCH_ACC;
            *t = tn;
            vector_scale(1.0, zn[0], yout);
            tolsf *= 2.0;
            break;
        }

        // Check for h below roundoff level in tn

        if (tn + h == tn) {
            nhnil++;
            if (nhnil <= mxhnil) {
                fprintf(errfp, MSG_HNIL, tn, h);
            }
            if (nhnil == mxhnil) {
                fprintf(errfp, MSG_HNIL_DONE, mxhnil);
            }
        }

        // Call CVStep to take a step

        kflag = cv_step(cv_mem);

        // Process failed step cases, and exit loop

        if (kflag != SUCCESS_STEP) {
            istate = cv_handle_failure(cv_mem, kflag);
            *t = tn;
            vector_scale(1.0, zn[0], yout);
            break;
        }

        nstloc++;

        // Check if in one-step mode, and if so copy y and exit loop

        if (itask == ONE_STEP) {
            istate = SUCCESS;
            *t = tn;
            vector_scale(1.0, zn[0], yout);
            next_q = qprime;
            next_h = hprime;
            break;
        }

        // Check if tout reached, and if so interpolate and exit loop

        if ((tn - tout)*h >= 0.0) {
            istate = SUCCESS;
            *t = tout;
            (void)cvode_dky(cv_mem, tout, 0, yout);
            next_q = qprime;
            next_h = hprime;
            break;
        }
    }

    // End of step loop; load optional outputs and return

    if (iopt != NULL) {
        iopt[NST] = nst;
        iopt[NFE] = nfe;
        iopt[NSETUPS] = nsetups;
        iopt[NNI] = nni;
        iopt[NCFN] = ncfn;
        iopt[NETF] = netf;
        iopt[QU] = q;
        iopt[QCUR] = next_q;
    }

    if (ropt != NULL) {
        ropt[HU] = h;
        ropt[HCUR] = next_h;
        ropt[TCUR] = tn;
        ropt[TOLSF] = tolsf;
    }

    return istate;
}

/*************** CVodeDky ********************************************

 This routine computes the k-th derivative of the interpolating
 polynomial at the time t and stores the result in the vector dky.
 The formula is:
          q
   dky = SUM c(j,k)*(t - tn)^(j-k)*h^(-j)*zn[j] ,
         j=k
 where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 zn[j] is the j-th column of the Nordsieck history array.

 This function is called by CVode with k = 0 and t = tout, but
 may also be called directly by the user.

**********************************************************************/

int32
cvode_dky(void *cvode_mem, double t, int32 k, Vector dky) {
    double s;
    double c;
    double r;
    double tfuzz;
    double tp;
    double tn1;
    int32 i;
    int32 j;
    CVodeMem cv_mem;

    cv_mem = (CVodeMem)cvode_mem;

    // Check all inputs for legality

    if (cvode_mem == NULL) {
        fprintf(stdout, MSG_DKY_NO_MEM);
        return DKY_NO_MEM;
    }

    if (dky == NULL) {
        fprintf(stdout, MSG_BAD_DKY);
        return BAD_DKY;
    }

    if ((k < 0) || (k > q)) {
        fprintf(errfp, MSG_BAD_K, k);
        return BAD_K;
    }

    tfuzz = FUZZ_FACTOR*uround*(tn + hu);
    tp = tn - hu - tfuzz;
    tn1 = tn + tfuzz;
    if ((t - tp)*(t - tn1) > 0.0) {
        fprintf(errfp, MSG_BAD_T, t, tn - hu, tn);
        return BAD_T;
    }

    // Sum the differentiated interpolating polynomial

    s = (t - tn) / h;
    for (j = q; j >= k; j--) {
        c = 1.0;
        for (i = j; i >= j - k + 1; i--) {
            c *= i;
        }
        if (j == q) {
            vector_scale(c, zn[q], dky);
        } else {
            vector_linear_sum(c, zn[j], s, dky, dky);
        }
    }
    if (k == 0) {
        return OKAY;
    }
    r = llnlmath_rpower_i(h, -k);
    vector_scale(r, dky, dky);
    return OKAY;
}

/********************* CVodeFree **********************************

 This routine frees the problem memory allocated by cvode_malloc.
 Such memory includes all the vectors allocated by CVAllocVectors,
 and the memory lmem for the linear solver (deallocated by a call
 to lfree).

*******************************************************************/

void
cvode_free(void *cvode_mem) {
    CVodeMem cv_mem;

    cv_mem = (CVodeMem)cvode_mem;

    if (cvode_mem == NULL) {
        return;
    }

    cv_free_vectors(cv_mem, qmax);
    if ((iter == NEWTON) && linitOK) {
        lfree(cv_mem);
    }
    free(cv_mem);
    return;
}

/***************************************************************/
/********** END Exported Functions Implementation **************/
/***************************************************************/

/*******************************************************************/
/******** BEGIN Private Helper Functions Implementation ************/
/*******************************************************************/

/****************** CVAllocVectors ***********************************

 This routine allocates the CVODE vectors ewt, acor, tempv, ftemp, and
 zn[0], ..., zn[maxord]. The length of the vectors is the input
 parameter neq and the maximum order (needed to allocate zn) is the
 input parameter maxord. If all memory allocations are successful,
 CVAllocVectors returns true. Otherwise all allocated memory is freed
 and CVAllocVectors returns false.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the double and int64 work spaces
 allocated here.

**********************************************************************/

static bool
cv_alloc_vectors(CVodeMem cv_mem, int64 neq, int32 maxord) {
    int32 i;
    int32 j;

    // Allocate ewt, acor, tempv, ftemp

    ewt = vector_new(neq);
    if (ewt == NULL) {
        return false;
    }
    acor = vector_new(neq);
    if (acor == NULL) {
        vector_free(ewt);
        return false;
    }
    tempv = vector_new(neq);
    if (tempv == NULL) {
        vector_free(ewt);
        vector_free(acor);
        return false;
    }
    ftemp = vector_new(neq);
    if (ftemp == NULL) {
        vector_free(tempv);
        vector_free(ewt);
        vector_free(acor);
        return false;
    }

    // Allocate zn[0] ... zn[maxord]

    for (j = 0; j <= maxord; j++) {
        zn[j] = vector_new(neq);
        if (zn[j] == NULL) {
            vector_free(ewt);
            vector_free(acor);
            vector_free(tempv);
            vector_free(ftemp);
            for (i = 0; i < j; i++) {
                vector_free(zn[i]);
            }
            return false;
        }
    }

    // Set solver workspace lengths

    lrw = (int32)((maxord + 5)*neq);
    liw = 0;

    return true;
}

/***************** CVFreeVectors *********************************

 This routine frees the CVODE vectors allocated in CVAllocVectors.

******************************************************************/

static void
cv_free_vectors(CVodeMem cv_mem, int32 maxord) {
    int32 j;

    vector_free(ewt);
    vector_free(acor);
    vector_free(tempv);
    vector_free(ftemp);
    for (j = 0; j <= maxord; j++) {
        vector_free(zn[j]);
    }
    return;
}

/*********************** cv_ewt_set **************************************

 This routine is responsible for setting the error weight vector
 ewtvec, according to tol_type, as follows:

 (1) ewtvec[i] = 1 / (*rtol*ABS(ycur[i]) + *atol), i=0,...,neq-1
     if tol_type = SS
 (2) ewtvec[i] = 1 / (*rtol*ABS(ycur[i]) + atol[i]), i=0,...,neq-1
     if tol_type = SV

  cv_ewt_set returns true if ewtvec is successfully set as above to a
  positive vector and false otherwise. In the latter case, ewtvec is
  considered undefined after the false return from cv_ewt_set.

  All the double work is done in the routines cv_ewt_set_ss, cv_ewt_set_sv.

***********************************************************************/

static bool
cv_ewt_set(CVodeMem cv_mem, double *rtol, void *atol, int32 tol_type,
           Vector ycur, Vector ewtvec, int64 neq) {
    switch (tol_type) {
    case SS:
        return cv_ewt_set_ss(cv_mem, rtol, (double *)atol, ycur, ewtvec, neq);
    case SV:
        return cv_ewt_set_sv(cv_mem, rtol, (Vector)atol, ycur, ewtvec, neq);
    default:
        break;
    }
    return 0;
}

/*********************** cv_ewt_set_ss *********************************

 This routine sets ewtvec as decribed above in the case tol_type=SS.
 It tests for non-positive components before inverting. cv_ewt_set_ss
 returns true if ewtvec is successfully set to a positive vector
 and false otherwise. In the latter case, ewtvec is considered
 undefined after the false return from cv_ewt_set_ss.

********************************************************************/

static bool
cv_ewt_set_ss(CVodeMem cv_mem, double *rtol, double *atol, Vector ycur,
              Vector ewtvec, int64 neq) {
    double rtoli;
    double atoli;
    (void)neq;

    rtoli = *rtol;
    atoli = *atol;
    vector_abs(ycur, tempv);
    vector_scale(rtoli, tempv, tempv);
    vector_add_const(tempv, atoli, tempv);
    if (vector_min(tempv) <= 0.0) {
        return false;
    }
    vector_inv(tempv, ewtvec);
    return true;
}

/*********************** cv_ewt_set_sv *********************************

 This routine sets ewtvec as decribed above in the case tol_type=SV.
 It tests for non-positive components before inverting. cv_ewt_set_sv
 returns true if ewtvec is successfully set to a positive vector
 and false otherwise. In the latter case, ewtvec is considered
 undefined after the false return from cv_ewt_set_sv.

********************************************************************/

static bool
cv_ewt_set_sv(CVodeMem cv_mem, double *rtol, Vector atol, Vector ycur,
              Vector ewtvec, int64 neq) {
    double rtoli;
    (void)neq;

    rtoli = *rtol;
    vector_abs(ycur, tempv);
    vector_linear_sum(rtoli, tempv, 1.0, atol, tempv);
    if (vector_min(tempv) <= 0.0) {
        return false;
    }
    vector_inv(tempv, ewtvec);
    return true;
}

/******************* CVHin ***************************************

 This routine computes a tentative initial step size h0.
 If tout is too close to tn (= t0), then CVHin returns false and
 h remains uninitialized. Otherwise, CVHin sets h to the chosen
 value h0 and returns true.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1,
 where ydd = estimated second derivative of y.

*****************************************************************/

static bool
cv_hin(CVodeMem cv_mem, double tout) {
    int32 sign;
    int32 count;
    double tdiff;
    double tdist;
    double tround;
    double hlb;
    double hub;
    double hg;
    double hgs;
    double hnew;
    double hrat;
    double h0;
    double yddnrm;

    // Test for tout too close to tn

    if ((tdiff = tout - tn) == 0.0) {
        return false;
    }

    sign = (tdiff > 0.0) ? 1 : -1;
    tdist = ABS(tdiff);
    tround = uround*MAX(ABS(tn), ABS(tout));
    if (tdist < 2.0*tround) {
        return false;
    }

    /* Set lower and upper bounds on h0, and take geometric mean
       Exit with this value if the bounds cross each other       */

    hlb = HLB_FACTOR*tround;
    hub = cv_upper_bound_h0(cv_mem, tdist);
    hg = llnlmath_rsqrt(hlb*hub);
    if (hub < hlb) {
        if (sign == -1) {
            hg = -hg;
        }
        h = hg;
        return true;
    }

    /* Loop up to MAX_ITERS times to find h0.
       Stop if new and previous values differ by a factor < 2.
       Stop if hnew/hg > 2 after one iteration, as this probably means
       that the ydd value is bad because of cancellation error.        */

    count = 0;
    while (true) {
        hgs = hg*sign;
        yddnrm = cv_ydd_norm(cv_mem, hgs);
        hnew = (yddnrm*hub*hub > 2.0) ? llnlmath_rsqrt(2.0 / yddnrm)
                                          : llnlmath_rsqrt(hg*hub);
        count++;
        if (count >= MAX_ITERS) {
            break;
        }
        hrat = hnew / hg;
        if ((hrat > HALF) && (hrat < 2.0)) {
            break;
        }
        if ((count >= 2) && (hrat > 2.0)) {
            hnew = hg;
            break;
        }
        hg = hnew;
    }

    // Apply bounds, bias factor, and attach sign

    h0 = H_BIAS*hnew;
    if (h0 < hlb) {
        h0 = hlb;
    }
    if (h0 > hub) {
        h0 = hub;
    }
    if (sign == -1) {
        h0 = -h0;
    }
    h = h0;
    return true;
}

/******************** CVUpperBoundH0 ******************************

 This routine sets an upper bound on ABS(h0) based on
 tdist = tn - t0 and the values of y[i]/ydot[i].

******************************************************************/

static double
cv_upper_bound_h0(CVodeMem cv_mem, double tdist) {
    double atoli = 0.0;
    double hub_inv;
    double hub;
    bool vectorAtol;
    Vector temp1;
    Vector temp2;

    vectorAtol = (itol == SV);
    if (!vectorAtol) {
        atoli = *((double *)abstol);
    }
    temp1 = tempv;
    temp2 = acor;
    vector_abs(zn[0], temp1);
    vector_abs(zn[1], temp2);
    if (vectorAtol) {
        vector_linear_sum(HUB_FACTOR, temp1, 1.0, (Vector)abstol, temp1);
    } else {
        vector_scale(HUB_FACTOR, temp1, temp1);
        vector_add_const(temp1, atoli, temp1);
    }
    vector_div(temp2, temp1, temp1);
    hub_inv = vector_max_norm(temp1);
    hub = HUB_FACTOR*tdist;
    if (hub*hub_inv > 1.0) {
        hub = 1.0 / hub_inv;
    }
    return hub;
}

/****************** CVYddNorm *************************************

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.

******************************************************************/

static double
cv_ydd_norm(CVodeMem cv_mem, double hg) {
    double yddnrm;

    vector_linear_sum(hg, zn[1], 1.0, zn[0], y);
    f(N, tn + hg, y, tempv, f_data);
    nfe++;
    vector_linear_sum(1.0, tempv, -1.0, zn[1], tempv);
    vector_scale(1.0 / hg, tempv, tempv);

    yddnrm = vector_wrms_norm(tempv, ewt);
    return yddnrm;
}

/********************* CVStep **************************************

 This routine performs one internal cvode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
  * preliminary adjustments if a new step size was chosen;
  * prediction of the Nordsieck history array zn at tn + h;
  * setting of multistep method coefficients and test quantities;
  * solution of the nonlinear system;
  * testing the local error;
  * updating zn and other state data if successful;
  * resetting stepsize and order for the next step.

 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.

********************************************************************/

int32
cv_step(CVodeMem cv_mem) {
    double saved_t;
    double dsm;
    int32 ncf;
    int32 nef;
    int32 nflag;
    int32 kflag;
    bool passed;

    saved_t = tn;
    ncf = nef = 0;
    nflag = FIRST_CALL;

    if ((nst > 0) && (hprime != h)) {
        cv_adjust_params(cv_mem);
    }

    // Looping point for attempts to take a step
    while (true) {
        cv_predict(cv_mem);
        cv_set(cv_mem);

        nflag = cv_nls(cv_mem, nflag);
        kflag = cv_handle_n_flag(cv_mem, &nflag, saved_t, &ncf);
        if (kflag == PREDICT_AGAIN) {
            continue;
        }
        if (kflag != DO_ERROR_TEST) {
            return kflag;
        }
        // Return if nonlinear solve failed and recovery not possible.

        passed = cv_do_error_test(cv_mem, &nflag, &kflag, saved_t, &nef, &dsm);
        if ((!passed) && (kflag == REP_ERR_FAIL)) {
            return kflag;
        }
        // Return if error test failed and recovery not possible.
        if (passed) {
            break;
        }
        // Retry step if error test failed, nflag == PREV_ERR_FAIL
    }

    /* Nonlinear system solve and error test were both successful;
       update data, and consider change of step and/or order       */

    cv_complete_step(cv_mem);
    cv_prepare_next_step(cv_mem, dsm);

    return SUCCESS_STEP;
}

/********************* CVAdjustParams ********************************

 This routine is called when a change in step size was decided upon,
 and it handles the required adjustments to the history array zn.
 If there is to be a change in order, we call CVAdjustOrder and reset
 q, L = q+1, and qwait.  Then in any case, we call CVRescale, which
 resets h and rescales the Nordsieck array.

**********************************************************************/

static void
cv_adjust_params(CVodeMem cv_mem) {
    if (qprime != q) {
        cv_adjust_order(cv_mem, qprime - q);
        q = qprime;
        L = q + 1;
        qwait = L;
    }
    cv_rescale(cv_mem);
    return;
}

/********************* CVAdjustOrder *****************************

  This routine is a high level routine which handles an order
  change by an amount deltaq (= +1 or -1). If a decrease in order
  is requested and q==2, then the routine returns immediately.
  Otherwise CVAdjustAdams or CVAdjustBDF is called to handle the
  order change (depending on the value of lmm).

******************************************************************/

static void
cv_adjust_order(CVodeMem cv_mem, int32 deltaq) {
    if ((q == 2) && (deltaq != 1)) {
        return;
    }

    switch (lmm) {
    case ADAMS:
        cv_adjust_adams(cv_mem, deltaq);
        break;
    case BDF:
        cv_adjust_bdf(cv_mem, deltaq);
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

/*************** CVAdjustAdams ***********************************

 This routine adjusts the history array on a change of order q by
 deltaq, in the case that lmm == ADAMS.

*****************************************************************/

static void
cv_adjust_adams(CVodeMem cv_mem, int32 deltaq) {
    int32 i;
    int32 j;
    double xi;
    double hsum;

    // On an order increase, set new column of zn to zero and return

    if (deltaq == 1) {
        vector_const(0.0, zn[L]);
        return;
    }

    /* On an order decrease, each zn[j] is adjusted by a multiple
       of zn[q].  The coefficients in the adjustment are the
       coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
       integrated, where xi_j = [t_n - t_(n-j)]/h.               */

    for (i = 0; i <= qmax; i++) {
        l[i] = 0.0;
    }
    l[1] = 1.0;
    hsum = 0.0;
    for (j = 1; j <= q - 2; j++) {
        hsum += tau[j];
        xi = hsum / hscale;
        for (i = j + 1; i >= 1; i--) {
            l[i] = l[i]*xi + l[i - 1];
        }
    }

    for (j = 1; j <= q - 2; j++) {
        l[j + 1] = q*(l[j] / (j + 1));
    }

    for (j = 2; j < q; j++) {
        vector_linear_sum(-l[j], zn[q], 1.0, zn[j], zn[j]);
    }
    return;
}

/********************** CVAdjustBDF *******************************

 This is a high level routine which handles adjustments to the
 history array on a change of order by deltaq in the case that
 lmm == BDF.  CVAdjustBDF calls CVIncreaseBDF if deltaq = +1 and
 CVDecreaseBDF if deltaq = -1 to do the actual work.

******************************************************************/

static void
cv_adjust_bdf(CVodeMem cv_mem, int32 deltaq) {
    switch (deltaq) {
    case 1:
        cv_increase_bdf(cv_mem);
        return;
    case -1:
        cv_decrese_bdf(cv_mem);
        return;
    default:
        break;
    }
    return;
}

/******************** CVIncreaseBDF **********************************

 This routine adjusts the history array on an increase in the
 order q in the case that lmm == BDF.
 A new column zn[q+1] is set equal to a multiple of the saved
 vector (= acor) in zn[qmax].  Then each zn[j] is adjusted by
 a multiple of zn[q+1].  The coefficients in the adjustment are the
 coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 where xi_j = [t_n - t_(n-j)]/h.

*********************************************************************/

static void
cv_increase_bdf(CVodeMem cv_mem) {
    double alpha0;
    double alpha1;
    double prod;
    double xi;
    double xiold;
    double hsum;
    double A1;
    int32 i;
    int32 j;

    for (i = 0; i <= qmax; i++) {
        l[i] = 0.0;
    }
    l[2] = alpha1 = prod = xiold = 1.0;
    alpha0 = -1.0;
    hsum = hscale;
    if (q > 1) {
        for (j = 1; j < q; j++) {
            hsum += tau[j + 1];
            xi = hsum / hscale;
            prod *= xi;
            alpha0 -= 1.0 / (j + 1);
            alpha1 += 1.0 / xi;
            for (i = j + 2; i >= 2; i--) {
                l[i] = l[i]*xiold + l[i - 1];
            }
            xiold = xi;
        }
    }
    A1 = (-alpha0 - alpha1) / prod;
    vector_scale(A1, zn[qmax], zn[L]);
    for (j = 2; j <= q; j++) {
        vector_linear_sum(l[j], zn[L], 1.0, zn[j], zn[j]);
    }
    return;
}

/********************* CVDecreaseBDF ******************************

 This routine adjusts the history array on a decrease in the
 order q in the case that lmm == BDF.
 Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 in the adjustment are the coefficients of the polynomial
 x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.

******************************************************************/

static void
cv_decrese_bdf(CVodeMem cv_mem) {
    double hsum;
    double xi;
    int32 i;
    int32 j;

    for (i = 0; i <= qmax; i++) {
        l[i] = 0.0;
    }
    l[2] = 1.0;
    hsum = 0.0;
    for (j = 1; j <= q - 2; j++) {
        hsum += tau[j];
        xi = hsum / hscale;
        for (i = j + 2; i >= 2; i--) {
            l[i] = l[i]*xi + l[i - 1];
        }
    }

    for (j = 2; j < q; j++) {
        vector_linear_sum(-l[j], zn[q], 1.0, zn[j], zn[j]);
    }
    return;
}

/**************** CVRescale ***********************************

  This routine rescales the Nordsieck array by multiplying the
  jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
  h is rescaled by eta, and hscale is reset to h.

***************************************************************/

static void
cv_rescale(CVodeMem cv_mem) {
    int32 j;
    double factor;

    factor = eta;
    for (j = 1; j <= q; j++) {
        vector_scale(factor, zn[j], zn[j]);
        factor *= eta;
    }
    h = hscale*eta;
    hscale = h;
    return;
}

/********************* CVPredict *************************************

 This routine advances tn by the tentative step size h, and computes
 the predicted array z_n(0), which is overwritten on zn.  The
 prediction of zn is done by repeated additions.

*********************************************************************/

static void
cv_predict(CVodeMem cv_mem) {
    int32 j;

    tn += h;
    for (int32 k = 1; k <= q; k++) {
        for (j = q; j >= k; j--) {
            vector_linear_sum(1.0, zn[j - 1], 1.0, zn[j], zn[j - 1]);
        }
    }
    return;
}

/************************** CVSet *********************************

 This routine is a high level routine which calls CVSetAdams or
 CVSetBDF to set the polynomial l, the test quantity array tq,
 and the related variables  rl1, gamma, and gamrat.

******************************************************************/

static void
cv_set(CVodeMem cv_mem) {
    switch (lmm) {
    case ADAMS:
        cv_set_adams(cv_mem);
        break;
    case BDF:
        cv_set_bdf(cv_mem);
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    rl1 = 1.0 / l[1];
    gamma = h*rl1;
    if (nst == 0) {
        gammap = gamma;
    }
    gamrat = (nst > 0) ? gamma / gammap : 1.0;  // protect x / x != 1.0
    return;
}

/******************** CVSetAdams *********************************

 This routine handles the computation of l and tq for the
 case lmm == ADAMS.

 The components of the vector l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                          q-1
 (d/dx) Lambda(x) = c*PRODUCT (1 + x / xi_i) , where
                          i=1
 Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
 Here xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.

*****************************************************************/

static void
cv_set_adams(CVodeMem cv_mem) {
    double m[L_MAX];
    double M[3];
    double hsum;

    if (q == 1) {
        l[0] = l[1] = tq[1] = tq[5] = 1.0;
        tq[2] = 2.0;
        tq[3] = TWELVE;
        tq[4] = CORTES*tq[2];  // = 0.1*tq[2]
        return;
    }

    hsum = cv_adams_start(cv_mem, m);

    M[0] = cv_alt_sum(q - 1, m, 1);
    M[1] = cv_alt_sum(q - 1, m, 2);

    cv_adams_finish(cv_mem, m, M, hsum);
    return;
}

/****************** CVAdamsStart ********************************

 This routine generates in m[] the coefficients of the product
 polynomial needed for the Adams l and tq coefficients for q > 1.

******************************************************************/

static double
cv_adams_start(CVodeMem cv_mem, double m[]) {
    double hsum;
    double xi_inv;
    double sum;
    int32 i;
    int32 j;

    hsum = h;
    m[0] = 1.0;
    for (i = 1; i <= q; i++) {
        m[i] = 0.0;
    }
    for (j = 1; j < q; j++) {
        if ((j == q - 1) && (qwait == 1)) {
            sum = cv_alt_sum(q - 2, m, 2);
            tq[1] = m[q - 2] / (q*sum);
        }
        xi_inv = h / hsum;
        for (i = j; i >= 1; i--) {
            m[i] += m[i - 1]*xi_inv;
        }
        hsum += tau[j];
        // The m[i] are coefficients of product(1 to j) (1 + x/xi_i)
    }
    return hsum;
}

/****************** CVAdamsFinish  *******************************

 This routine completes the calculation of the Adams l and tq.

******************************************************************/

static void
cv_adams_finish(CVodeMem cv_mem, double m[], double M[], double hsum) {
    int32 i;
    double M0_inv;
    double xi;
    double xi_inv;

    M0_inv = 1.0 / M[0];

    l[0] = 1.0;
    for (i = 1; i <= q; i++) {
        l[i] = M0_inv*(m[i - 1] / i);
    }
    xi = hsum / h;
    xi_inv = 1.0 / xi;

    tq[2] = xi*M[0] / M[1];
    tq[5] = xi / l[q];

    if (qwait == 1) {
        for (i = q; i >= 1; i--) {
            m[i] += m[i - 1]*xi_inv;
        }
        M[2] = cv_alt_sum(q, m, 2);
        tq[3] = L*M[0] / M[2];
    }

    tq[4] = CORTES*tq[2];
    return;
}

/****************** CVAltSum **************************************

 CVAltSum returns the value of the alternating sum
   sum (i= 0 ... iend) [ (-1)^i*(a[i] / (i + k)) ].
 If iend < 0 then CVAltSum returns 0.
 This operation is needed to compute the integral, from -1 to 0,
 of a polynomial x^(k-1) M(x) given the coefficients of M(x).

******************************************************************/

static double
cv_alt_sum(int32 iend, double a[], int32 k) {
    int32 i;
    int32 sign;
    double sum;

    if (iend < 0) {
        return 0.0;
    }

    sum = 0.0;
    sign = 1;
    for (i = 0; i <= iend; i++) {
        sum += sign*(a[i] / (i + k));
        sign = -sign;
    }
    return sum;
}

/***************** CVSetBDF **************************************

 This routine computes the coefficients l and tq in the case
 lmm == BDF.  CVSetBDF calls cv_set_tq_bdf to set the test
 quantity vector tq.

 The components of the vector l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                                 q-1
 Lambda(x) = (1 + x / xi*_q)*PRODUCT (1 + x / xi_i) , where
                                 i=1
 xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order.

*****************************************************************/

static void
cv_set_bdf(CVodeMem cv_mem) {
    double alpha0;
    double alpha0_hat;
    double xi_inv;
    double xistar_inv;
    double hsum;
    int32 i;
    int32 j;

    l[0] = l[1] = xi_inv = xistar_inv = 1.0;
    for (i = 2; i <= q; i++) {
        l[i] = 0.0;
    }
    alpha0 = alpha0_hat = -1.0;
    hsum = h;
    if (q > 1) {
        for (j = 2; j < q; j++) {
            hsum += tau[j - 1];
            xi_inv = h / hsum;
            alpha0 -= 1.0 / j;
            for (i = j; i >= 1; i--) {
                l[i] += l[i - 1]*xi_inv;
            }
            // The l[i] are coefficients of product(1 to j) (1 + x/xi_i)
        }

        // j = q
        alpha0 -= 1.0 / q;
        xistar_inv = -l[1] - alpha0;
        hsum += tau[q - 1];
        xi_inv = h / hsum;
        alpha0_hat = -l[1] - xi_inv;
        for (i = q; i >= 1; i--) {
            l[i] += l[i - 1]*xistar_inv;
        }
    }

    cv_set_tq_bdf(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
    return;
}

/****************** cv_set_tq_bdf ************************************

 This routine sets the test quantity vector tq in the case
 lmm == BDF.

******************************************************************/

static void
cv_set_tq_bdf(CVodeMem cv_mem, double hsum, double alpha0, double alpha0_hat,
              double xi_inv, double xistar_inv) {
    double A1;
    double A2;
    double A3;
    double A4;
    double A5;
    double A6;
    double C;
    double CPrime;
    double CPrimePrime;

    A1 = 1.0 - alpha0_hat + alpha0;
    A2 = 1.0 + q*A1;
    tq[2] = ABS(alpha0*(A2 / A1));
    tq[5] = ABS((A2) / (l[q]*xi_inv / xistar_inv));
    if (qwait == 1) {
        C = xistar_inv / l[q];
        A3 = alpha0 + 1.0 / q;
        A4 = alpha0_hat + xi_inv;
        CPrime = A3 / (1.0 - A4 + A3);
        tq[1] = ABS(CPrime / C);
        hsum += tau[q];
        xi_inv = h / hsum;
        A5 = alpha0 - (1.0 / (q + 1));
        A6 = alpha0_hat - xi_inv;
        CPrimePrime = A2 / (1.0 - A6 + A5);
        tq[3] = ABS(CPrimePrime*xi_inv*(q + 2)*A5);
    }
    tq[4] = CORTES*tq[2];
    return;
}

/****************** CVnls *****************************************

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 Depending on iter, it calls CVnlsFunctional or CVnlsNewton
 to do the work.

******************************************************************/

static int32
cv_nls(CVodeMem cv_mem, int32 nflag) {
    switch (iter) {
    case FUNCTIONAL:
        return cv_nls_functional(cv_mem);
    case NEWTON:
        return cv_nls_newton(cv_mem, nflag);
    default:
        break;
    }
    return 0;
}

/***************** CVnlsFunctional ********************************

 This routine attempts to solve the nonlinear system using
 functional iteration (no matrices involved).

******************************************************************/

static int32
cv_nls_functional(CVodeMem cv_mem) {
    int32 m;
    double del;
    double delp = 0;
    double dcon;

    // Initialize counter and evaluate f at predicted y

    crate = 1.0;
    m = 0;
    f(N, tn, zn[0], tempv, f_data);
    nfe++;
    vector_const(0.0, acor);

    // Loop until convergence; accumulate corrections in acor

    while (true) {
        // Correct y directly from the last f value
        vector_linear_sum(h, tempv, -1.0, zn[1], tempv);
        vector_scale(rl1, tempv, tempv);
        vector_linear_sum(1.0, zn[0], 1.0, tempv, y);
        // Get WRMS norm of current correction to use in convergence test
        vector_linear_sum(1.0, tempv, -1.0, acor, acor);
        del = vector_wrms_norm(acor, ewt);
        vector_scale(1.0, tempv, acor);

        /* Test for convergence.  If m > 0, an estimate of the convergence
           rate constant is stored in crate, and used in the test.        */
        if (m > 0) {
            crate = MAX(CRDOWN*crate, del / delp);
        }
        dcon = del*MIN(1.0, crate) / tq[4];
        if (dcon <= 1.0) {
            acnrm = (m == 0) ? del : vector_wrms_norm(acor, ewt);
            return SOLVED;  // Convergence achieved
        }

        // Stop at maxcor iterations or if iter. seems to be diverging
        m++;
        if ((m == maxcor) || ((m >= 2) && (del > RDIV*delp))) {
            return CONV_FAIL;
        }
        // Save norm of correction, evaluate f, and loop again
        delp = del;
        f(N, tn, y, tempv, f_data);
        nfe++;
    }
}

/*********************** CVnlsNewton **********************************

 This routine handles the Newton iteration. It calls lsetup if
 indicated, calls CVNewtonIteration to perform the iteration, and
 retries a failed attempt at Newton iteration if that is indicated.
 See return values at top of this file.

**********************************************************************/

static int32
cv_nls_newton(CVodeMem cv_mem, int32 nflag) {
    Vector vtemp1;
    Vector vtemp2;
    Vector vtemp3;
    int32 convfail;
    int32 ier;
    bool callSetup;

    vtemp1 = acor;   // rename acor as vtemp1 for readability
    vtemp2 = y;      // rename y as vtemp2 for readability
    vtemp3 = tempv;  // rename tempv as vtemp3 for readability

    // Set flag convfail, input to lsetup for its evaluation decision
    convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ? NO_FAILURES
                                                                   : FAIL_OTHER;

    // Decide whether or not to call setup routine (if one exists)
    if (setupNonNull) {
        callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
                    (nst == 0) || (nst >= nstlp + MSBP) ||
                    (ABS(gamrat - 1.0) > DGMAX);
    } else {
        crate = 1.0;
        callSetup = false;
    }

    /* Looping point for the solution of the nonlinear system.
       Evaluate f at the predicted y, call lsetup if indicated, and
       call CVNewtonIteration for the Newton iteration itself.      */

    while (true) {

        f(N, tn, zn[0], ftemp, f_data);
        nfe++;

        if (callSetup) {
            ier = lsetup(cv_mem, convfail, zn[0], ftemp, &jcur, vtemp1, vtemp2,
                         vtemp3);
            nsetups++;
            callSetup = false;
            gamrat = crate = 1.0;
            gammap = gamma;
            nstlp = nst;
            // Return if lsetup failed
            if (ier < 0) {
                return SETUP_FAIL_UNREC;
            }
            if (ier > 0) {
                return CONV_FAIL;
            }
        }

        // Set acor to zero and load prediction into y vector
        vector_const(0.0, acor);
        vector_scale(1.0, zn[0], y);

        // Do the Newton iteration
        ier = cv_newton_iteration(cv_mem);

        /* If there is a convergence failure and the Jacobian-related
           data appears not to be current, loop again with a call to lsetup
           in which convfail=FAIL_BAD_J.  Otherwise return.                 */
        if (ier != TRY_AGAIN) {
            return ier;
        }

        callSetup = true;
        convfail = FAIL_BAD_J;
    }
}

/********************** CVNewtonIteration ****************************

 This routine performs the Newton iteration. If the iteration succeeds,
 it returns the value SOLVED. If not, it may signal the CVnlsNewton
 routine to call lsetup again and reattempt the iteration, by
 returning the value TRY_AGAIN. (In this case, CVnlsNewton must set
 convfail to FAIL_BAD_J before calling setup again).
 Otherwise, this routine returns one of the appropriate values
 SOLVE_FAIL_UNREC or CONV_FAIL back to CVnlsNewton.

*********************************************************************/

static int32
cv_newton_iteration(CVodeMem cv_mem) {
    int32 m;
    int32 ret;
    double del;
    double delp = 0.0;
    double dcon;
    Vector b;

    mnewt = m = 0;

    // Looping point for Newton iteration
    while (true) {

        // Evaluate the residual of the nonlinear system
        vector_linear_sum(rl1, zn[1], 1.0, acor, tempv);
        vector_linear_sum(gamma, ftemp, -1.0, tempv, tempv);

        // Call the lsolve function
        b = tempv;
        ret = lsolve(cv_mem, b, y, ftemp);
        nni++;

        if (ret < 0) {
            return SOLVE_FAIL_UNREC;
        }

        /* If lsolve had a recoverable failure and Jacobian data is
           not current, signal to try the solution again            */
        if (ret > 0) {
            if ((!jcur) && (setupNonNull)) {
                return TRY_AGAIN;
            }
            return CONV_FAIL;
        }
        // ***************
        // Get WRMS norm of correction; add correction to acor and y
        del = vector_wrms_norm(b, ewt);
        vector_linear_sum(1.0, acor, 1.0, b, acor);
        vector_linear_sum(1.0, zn[0], 1.0, acor, y);

        /* Test for convergence.  If m > 0, an estimate of the convergence
           rate constant is stored in crate, and used in the test.        */
        if (m > 0) {
            crate = MAX(CRDOWN*crate, del / delp);
        }
        dcon = del*MIN(1.0, crate) / tq[4];

        if (dcon <= 1.0) {
            acnrm = (m == 0) ? del : vector_wrms_norm(acor, ewt);
            jcur = false;
            return SOLVED;  // Nonlinear system was solved successfully
        }

        mnewt = ++m;

        /* Stop at maxcor iterations or if iter. seems to be diverging.
           If still not converged and Jacobian data is not current,
           signal to try the solution again                            */
        if ((m == maxcor) || ((m >= 2) && (del > RDIV*delp))) {
            if ((!jcur) && (setupNonNull)) {
                return TRY_AGAIN;
            }
            return CONV_FAIL;
        }

        // Save norm of correction, evaluate f, and loop again
        delp = del;
        f(N, tn, y, ftemp, f_data);
        nfe++;
    }
}

/********************** cv_handle_n_flag *******************************

 This routine takes action on the return value nflag = *nflagPtr
 returned by CVnls, as follows:

 If CVnls succeeded in solving the nonlinear system, then
 cv_handle_n_flag returns the constant DO_ERROR_TEST, which tells CVStep
 to perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented and Nordsieck array zn is restored.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value SETUP_FAILED.

 If it failed due to an unrecoverable failure in solve, then we return
 the value SOLVE_FAILED.

 Otherwise, a recoverable failure occurred when solving the
 nonlinear system (CVnls returned nflag == CONV_FAIL).
   In this case, we return the value REP_CONV_FAIL if ncf is now
   equal to MXNCF or |h| = hmin.
   If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
   PREDICT_AGAIN, telling CVStep to reattempt the step.

*********************************************************************/

static int32
cv_handle_n_flag(CVodeMem cv_mem, int32 *nflagPtr, double saved_t,
                 int32 *ncfPtr) {
    int32 nflag;

    nflag = *nflagPtr;

    if (nflag == SOLVED) {
        return DO_ERROR_TEST;
    }

    // The nonlinear soln. failed; increment ncfn and restore zn
    ncfn++;
    cv_restore(cv_mem, saved_t);

    // Return if lsetup or lsolve failed unrecoverably
    if (nflag == SETUP_FAIL_UNREC) {
        return SETUP_FAILED;
    }

    if (nflag == SOLVE_FAIL_UNREC) {
        return SOLVE_FAILED;
    }

    // At this point, nflag == CONV_FAIL; increment ncf

    (*ncfPtr)++;
    etamax = 1.0;
    // If we had MXNCF failures or |h| = hmin, return REP_CONV_FAIL
    if ((ABS(h) <= hmin*ONEPSM) || (*ncfPtr == MXNCF)) {
        return REP_CONV_FAIL;
    }

    // Reduce step size; return to reattempt the step
    eta = MAX(ETACF, hmin / ABS(h));
    *nflagPtr = PREV_CONV_FAIL;
    cv_rescale(cv_mem);
    return PREDICT_AGAIN;
}

/********************** CVRestore ************************************

 This routine restores the value of tn to saved_t and undoes the
 prediction.  After execution of CVRestore, the Nordsieck array zn has
 the same values as before the call to CVPredict.

********************************************************************/

static void
cv_restore(CVodeMem cv_mem, double saved_t) {
    int32 j;

    tn = saved_t;
    for (int32 k = 1; k <= q; k++) {
        for (j = q; j >= k; j--) {
            vector_linear_sum(1.0, zn[j - 1], -1.0, zn[j], zn[j - 1]);
        }
    }
    return;
}

/******************* cv_do_error_test ********************************

 This routine performs the local error test.
 The weighted local error norm dsm is loaded into *dsmPtr, and
 the test dsm ?<= 1 is made.

 If the test passes, cv_do_error_test returns true.

 If the test fails, we undo the step just taken (call CVRestore),
 set *nflagPtr to PREV_ERR_FAIL, and return false.

 If MXNEF error test failures have occurred or if ABS(h) = hmin,
 we set *kflagPtr = REP_ERR_FAIL. (Otherwise *kflagPtr has the
 value last returned by CVHandleNflag.)

 If more than MXNEF1 error test failures have occurred, an order
 reduction is forced.

******************************************************************/

static bool
cv_do_error_test(CVodeMem cv_mem, int32 *nflagPtr, int32 *kflagPtr,
                 double saved_t, int32 *nefPtr, double *dsmPtr) {
    double dsm;

    dsm = acnrm / tq[2];

    // If est. local error norm dsm passes test, return true
    *dsmPtr = dsm;
    if (dsm <= 1.0) {
        return true;
    }

    // Test failed; increment counters, set nflag, and restore zn array
    (*nefPtr)++;
    netf++;
    *nflagPtr = PREV_ERR_FAIL;
    cv_restore(cv_mem, saved_t);

    // At MXNEF failures or |h| = hmin, return with kflag = REP_ERR_FAIL
    if ((ABS(h) <= hmin*ONEPSM) || (*nefPtr == MXNEF)) {
        *kflagPtr = REP_ERR_FAIL;
        return false;
    }

    // Set etamax = 1 to prevent step size increase at end of this step
    etamax = 1.0;

    // Set h ratio eta from dsm, rescale, and return for retry of step
    if (*nefPtr <= MXNEF1) {
        eta = 1.0 / (llnlmath_rpower_r(BIAS2*dsm, 1.0 / L) + ADDON);
        eta = MAX(ETAMIN, MAX(eta, hmin / ABS(h)));
        if (*nefPtr >= SMALL_NEF) {
            eta = MIN(eta, ETAMXF);
        }
        cv_rescale(cv_mem);
        return false;
    }

    // After MXNEF1 failures, force an order reduction and retry step
    if (q > 1) {
        eta = MAX(ETAMIN, hmin / ABS(h));
        cv_adjust_order(cv_mem, -1);
        L = q;
        q--;
        qwait = L;
        cv_rescale(cv_mem);
        return false;
    }

    // If already at order 1, restart: reload zn from scratch
    eta = MAX(ETAMIN, hmin / ABS(h));
    h *= eta;
    hscale = h;
    qwait = LONG_WAIT;
    f(N, tn, zn[0], tempv, f_data);
    nfe++;
    vector_scale(h, tempv, zn[1]);
    return false;
}

/*************** CVCompleteStep **********************************

 This routine performs various update operations when the solution
 to the nonlinear system has passed the local error test.
 We increment the step counter nst, record the values hu and qu,
 update the tau array, and apply the corrections to the zn array.
 The tau[i] are the last q values of h, with tau[1] the most recent.
 The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 we save acor and tq[5] for a possible order increase.

******************************************************************/

static void
cv_complete_step(CVodeMem cv_mem) {
    int32 i;
    int32 j;

    nst++;
    hu = h;
    qu = q;

    for (i = q; i >= 2; i--) {
        tau[i] = tau[i - 1];
    }
    if ((q == 1) && (nst > 1)) {
        tau[2] = tau[1];
    }
    tau[1] = h;

    for (j = 0; j <= q; j++) {
        vector_linear_sum(l[j], acor, 1.0, zn[j], zn[j]);
    }
    qwait--;
    if ((qwait == 1) && (q != qmax)) {
        vector_scale(1.0, acor, zn[qmax]);
        saved_tq5 = tq[5];
    }
    return;
}

/************* CVPrepareNextStep **********************************

 This routine handles the setting of stepsize and order for the
 next step -- hprime and qprime.  A  with hprime, it sets the
 ratio eta = hprime/h.  It also updates other state variables
 related to a change of step size or order.  Finally, we rescale
 the acor array to be the estimated local error vector.

******************************************************************/

static void
cv_prepare_next_step(CVodeMem cv_mem, double dsm) {
    double etaqm1;
    double etaq;
    double etaqp1;

    // If etamax = 1, defer step size or order changes
    if (etamax == 1.0) {
        qwait = MAX(qwait, 2);
        qprime = q;
        hprime = h;
        eta = 1.0;
        etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;
        vector_scale(1.0 / tq[2], acor, acor);
        return;
    }

    // etaq is the ratio of new to old h at the current order
    etaq = 1.0 / (llnlmath_rpower_r(BIAS2*dsm, 1.0 / L) + ADDON);

    // If no order change, adjust eta and acor in CVSetEta and return
    if (qwait != 0) {
        eta = etaq;
        qprime = q;
        cv_set_eta(cv_mem);
        return;
    }

    /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are
       the ratios of new to old h at orders q-1 and q+1, respectively.
       cv_choose_eta selects the largest; CVSetEta adjusts eta and acor */
    qwait = 2;
    etaqm1 = cv_compute_etaqm1(cv_mem);
    etaqp1 = cv_compute_etaqp1(cv_mem);
    cv_choose_eta(cv_mem, etaqm1, etaq, etaqp1);
    cv_set_eta(cv_mem);
    return;
}

/***************** CVSetEta ***************************************

 This routine adjusts the value of eta according to the various
 heuristic limits and the optional input hmax.  It also resets
 etamax and rescales acor to be the estimated local error vector.

*******************************************************************/

static void
cv_set_eta(CVodeMem cv_mem) {

    // If eta below the threshhold THRESH, reject a change of step size
    if (eta < THRESH) {
        eta = 1.0;
        hprime = h;
    } else {
        // Limit eta by etamax and hmax, then set hprime
        eta = MIN(eta, etamax);
        eta /= MAX(1.0, ABS(h)*hmax_inv*eta);
        hprime = h*eta;
    }

    // Reset etamx for the next step size change, and scale acor
    etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;
    vector_scale(1.0 / tq[2], acor, acor);
    return;
}

/*************** CVComputeEtaqm1 **********************************

 This routine computes and returns the value of etaqm1 for a
 possible decrease in order by 1.

******************************************************************/

static double
cv_compute_etaqm1(CVodeMem cv_mem) {
    double etaqm1;
    double ddn;

    etaqm1 = 0.0;
    if (q > 1) {
        ddn = vector_wrms_norm(zn[q], ewt) / tq[1];
        etaqm1 = 1.0 / (llnlmath_rpower_r(BIAS1*ddn, 1.0 / q) + ADDON);
    }
    return etaqm1;
}

/*************** CVComputeEtaqp1 **********************************

 This routine computes and returns the value of etaqp1 for a
 possible increase in order by 1.

******************************************************************/

static double
cv_compute_etaqp1(CVodeMem cv_mem) {
    double etaqp1;
    double dup;
    double cquot;

    etaqp1 = 0.0;
    if (q != qmax) {
        cquot = (tq[5] / saved_tq5)*llnlmath_rpower_i(h / tau[2], L);
        vector_linear_sum(-cquot, zn[qmax], 1.0, acor, tempv);
        dup = vector_wrms_norm(tempv, ewt) / tq[3];
        etaqp1 = 1.0 / (llnlmath_rpower_r(BIAS3*dup, 1.0 / (L + 1)) + ADDON);
    }
    return etaqp1;
}

/******************* cv_choose_eta **********************************

 Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 q - 1, q, or q + 1, respectively), this routine chooses the
 maximum eta value, sets eta to that value, and sets qprime to the
 corresponding value of q.  If there is a tie, the preference
 order is to (1) keep the same order, then (2) decrease the order,
 and finally (3) increase the order.  If the maximum eta value
 is below the threshhold THRESH, the order is kept unchanged and
 eta is set to 1.

******************************************************************/

static void
cv_choose_eta(CVodeMem cv_mem, double etaqm1, double etaq, double etaqp1) {
    double etam;

    etam = MAX(etaqm1, MAX(etaq, etaqp1));

    if (etam < THRESH) {
        eta = 1.0;
        qprime = q;
        return;
    }

    if (etam == etaq) {
        eta = etaq;
        qprime = q;
    } else if (etam == etaqm1) {
        eta = etaqm1;
        qprime = q - 1;
    } else {
        eta = etaqp1;
        qprime = q + 1;
        vector_scale(1.0, acor, zn[qmax]);
    }
    return;
}

/****************** CVHandleFailure ******************************

 This routine prints error messages for all cases of failure by
 CVStep. It returns to CVode the value that CVode is to return to
 the user.

*****************************************************************/

static int32
cv_handle_failure(CVodeMem cv_mem, int32 kflag) {

    // Set imxer to the index of maximum weighted local error
    vector_prod(acor, ewt, tempv);
    vector_abs(tempv, tempv);

    // Depending on kflag, print error message and return error flag
    switch (kflag) {
    case REP_ERR_FAIL:
        fprintf(errfp, MSG_ERR_FAILS, tn, h);
        return ERR_FAILURE;
    case REP_CONV_FAIL:
        fprintf(errfp, MSG_CONV_FAILS, tn, h);
        return CONV_FAILURE;
    case SETUP_FAILED:
        fprintf(errfp, MSG_SETUP_FAILED, tn);
        return SETUP_FAILURE;
    case SOLVE_FAILED:
        fprintf(errfp, MSG_SOLVE_FAILED, tn);
        return SOLVE_FAILURE;
    default:
        break;
    }
    return ERR_FAILURE;
}

/*******************************************************************/
/********* END Private Helper Functions Implementation *************/
/*******************************************************************/

/***************************************************************/
/************** END CVODE Implementation ***********************/
/***************************************************************/
