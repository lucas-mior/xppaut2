#include "functions.h"

#include "cv2.h"
#include "parserslow.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "integers.h"

#include "phsplan.h"

/*  this is also X free ! */

#define GEAR 5
#define RKQS 8
#define STIFF 9
#define CVODE 10
#define DP5 11
#define DP83 12
#define RB23 13

static struct FitInfo {
    char file[25];
    char varlist[25];
    char collist[25];
    char parlist1[25];
    char parlist2[25];
    int32 dim;
    int32 npars;
    int32 nvars;
    int32 npts;
    int32 maxiter;
    int32 icols[50];
    int32 ipar[50];
    int32 ivar[50];
    double tol;
    double eps;
} fit_info;

static void do_fit_parse_parlist(char *parlist, int32 *ipars, int32 *n);
static void do_fit_parse_varlist(char *varlist, int32 *ivars, int32 *n);
static void do_fit_parse_collist(char *collist, int32 *icols, int32 *n);

void
do_fit_init_info(void) {
    fit_info.tol = .001;
    fit_info.eps = 1e-5;
    fit_info.dim = 0;
    fit_info.npars = 0;
    fit_info.nvars = 0;
    fit_info.varlist[0] = 0;
    fit_info.collist[0] = 0;
    fit_info.parlist1[0] = 0;
    fit_info.parlist2[0] = 0;
    fit_info.npts = 0;
    fit_info.maxiter = 20;
    fit_info.file[0] = 0;
    return;
}

/*
  y     initial condition
  a     initial guesses for the parameters
  t0    vector of output times
  flag  1 for success  0 for failure
  eps   derivative step
  yfit  has y[i1](t0),...,y[im](t0), ..., y[i1](tn),...,y[im](tn)
        which are the values of the test functions at the npts
        times.  yfit is (npts)*nvars int32
  yderv[npar][nvars*(npts)] is the derivative of yfit with rrspect
        to the parameter

  npts   is the number of times to be fitted
  npars    the number of parameters
  nvars    the number of variables
  ipar     the vector of parameters  negative are constants
           positive are initial data
  ivar     the vector of variables

 */
void
do_fit_get_info(double *y, double *a, double *t0, int32 *flag, double eps,
                double *yfit, double **yderv, int32 npts, int32 npars,
                int32 nvars, int32 *ivar, int32 *ipar) {
    int32 iv;
    int32 ip;
    int32 istart = 1;
    int32 l;
    int32 k0;
    int32 ok;
    double yold[MAX_ODE];
    double dp;
    double par;
    *flag = 0;
    // set up all initial data and parameter guesses
    for (l = 0; l < npars; l++) {
        ip = ipar[l];
        if (ip < 0) {
            constants[-ip] = a[l];
        } else {
            y[ip] = a[l];
        }
    }
    for (int32 i = 0; i < NODE; i++) {
        yold[i] = y[i];
    }
    if (DelayFlag) {
        // restart initial data
        if (delay_handle_do_init_delay(DELAY) == 0) {
            return;
        }
    }
    derived_evaluate();
    //   This gets the values at the desired points
    for (int32 i = 0; i < nvars; i++) {
        iv = ivar[i];
        yfit[i] = y[iv];
    }
    for (int32 k = 1; k < npts; k++) {
        k0 = k*nvars;
        ok = do_fit_one_step_int(y, t0[k - 1], t0[k], &istart);
        if (ok == 0) {
            for (int32 i = 0; i < NODE; i++) {
                y[i] = yold[i];
            }

            return;
        }

        for (int32 i = 0; i < nvars; i++) {
            iv = ivar[i];
            yfit[i + k0] = y[iv];
        }
    }
#ifdef CVODE_YES
    if (METHOD == CVODE) {
        cv_end();
    }
#endif
    //  Now we take the derivatives !!
    for (l = 0; l < npars; l++) {
        istart = 1;
        // set up all the initial conditions
        for (int32 j = 0; j < nvars; j++) {
            yderv[l][j] = 0.0;  // no dependence on initial data ...
        }
        for (int32 i = 0; i < NODE; i++) {
            y[i] = yold[i];
        }
        ip = ipar[l];
        if (ip < 0) {
            par = constants[-ip];
            dp = eps*MAX(eps, fabs(par));
            constants[-ip] = par + dp;
        } else {
            par = yold[ip];
            dp = eps*MAX(eps, fabs(par));
            y[ip] = par + dp;
            for (int32 j = 0; j < nvars; j++) {
                if (ip == ivar[j]) {
                    yderv[l][j] =
                        1.0;  // ... except for those ICs that can vary
                }
            }
        }
        if (DelayFlag) {
            // restart initial data
            if (delay_handle_do_init_delay(DELAY) == 0) {
                return;
            }
        }
        derived_evaluate();
        // now loop through all the points
        for (int32 k = 1; k < npts; k++) {
            k0 = k*nvars;
            ok = do_fit_one_step_int(y, t0[k - 1], t0[k], &istart);
            if (ok == 0) {
                for (int32 i = 0; i < NODE; i++) {
                    y[i] = yold[i];
                }

                return;
            }
            for (int32 i = 0; i < nvars; i++) {
                iv = ivar[i];
                yderv[l][i + k0] = (y[iv] - yfit[i + k0]) / dp;
            }
        }
        // Now return the parameter to its old value
        if (ip < 0) {
            constants[-ip] = par;
        }
        derived_evaluate();
#ifdef CVODE_YES
        if (METHOD == CVODE) {
            cv_end();
        }
#endif
    }
    *flag = 1;
    for (int32 i = 0; i < NODE; i++) {
        y[i] = yold[i];
    }

    return;
}

void
do_fit_printem(double **yderv, double *yfit, double *t0, int32 npars,
               int32 nvars, int32 npts) {
    int32 ioff;
    for (int32 i = 0; i < npts; i++) {
        ggets_plintf(" %8.5g ", t0[i]);
        ioff = nvars*i;
        for (int32 j = 0; j < nvars; j++) {
            ggets_plintf(" %g ", yfit[ioff + j]);
            for (int32 k = 0; k < npars; k++) {
                printf(" %g ", yderv[k][ioff + j]);
            }
        }
        ggets_plintf(" \n");
    }
    return;
}

int32
do_fit_one_step_int(double *y, double t0, double t1, int32 *istart) {
    int32 nit;
    int32 kflag;
    double dt = DELTA_T;
    double z;
    double error[MAX_ODE];
    double t = t0;
#ifdef CVODE_YES
    if (METHOD == CVODE) {
        cvode(istart, y, &t, NODE, t1, &kflag, &TOLER, &atoler);
        if (kflag < 0) {
            cvode_err_msg(kflag);
            return 0;
        }
        delay_handle_stor_delay(y);
        return 1;
    }
#endif
    if (METHOD == DP5 || METHOD == DP83) {
        dp(istart, y, &t, NODE, t1, &TOLER, &atoler, METHOD - DP5, &kflag);
        if (kflag != 1) {
            dormpri_dp_err(kflag);
            return 0;
        }
        delay_handle_stor_delay(y);
        return 1;
    }
    if (METHOD == RB23) {
        odesol_rb23(y, &t, t1, istart, NODE, WORK, &kflag);
        if (kflag < 0) {
            ggets_err_msg("Step size too small");
            return 0;
        }
        delay_handle_stor_delay(y);
        return 1;
    }
    if (METHOD == RKQS || METHOD == STIFF) {
        stiff_adaptive(y, NODE, &t, t1, TOLER, &dt, HMIN, WORK, &kflag,
                       NEWT_ERR, METHOD, istart);
        if (kflag) {
            ggets_ping();
            switch (kflag) {
            case 2:
                ggets_err_msg("Step size too small");
                break;
            case 3:
                ggets_err_msg("Too many steps");
                break;
            case -1:
                ggets_err_msg("singular jacobian encountered");
                break;
            case 1:
                ggets_err_msg("stepsize is close to 0");
                break;
            case 4:
                ggets_err_msg("exceeded MAXTRY in stiff");
                break;
            default:
                fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
                exit(EXIT_FAILURE);
            }
            return 0;
        }
        delay_handle_stor_delay(y);
        return 1;
    }
    /* cvode(command,y,t,n,tout,kflag,atol,rtol)
   command =0 continue, 1 is start 2 finish   */
    if (METHOD == GEAR) {
        gear(NODE, &t, t1, y, HMIN, HMAX, TOLER, 2, error, &kflag, istart, WORK,
             IWORK);
        if (kflag < 0) {
            ggets_ping();
            switch (kflag) {
            case -1:
                ggets_err_msg("kflag=-1: minimum step too big");
                break;
            case -2:
                ggets_err_msg("kflag=-2: required order too big");
                break;
            case -3:
                ggets_err_msg("kflag=-3: minimum step too big");
                break;
            case -4:
                ggets_err_msg("kflag=-4: tolerance too small");
                break;
            default:
                fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
                exit(EXIT_FAILURE);
            }

            return 0;
        }
        delay_handle_stor_delay(y);
        return 1;
    }
    if (METHOD == 0) {
        nit = (int32)(fabs(t0 - t1));
        dt = dt / fabs(dt);
        kflag = solver(y, &t, dt, nit, NODE, istart, WORK);

        return 1;
    }
    z = (t1 - t0) / dt;
    nit = (int32)z;
    kflag = solver(y, &t, dt, nit, NODE, istart, WORK);

    if (kflag < 0) {
        return 0;
    }

    if ((dt < 0 && t > t1) || (dt > 0 && t < t1)) {
        dt = t1 - t;
        kflag = solver(y, &t, dt, 1, NODE, istart, WORK);
        if (kflag < 0) {
            return 0;
        }
    }

    return 1;
}

void
do_fit_test(void) {
    double *yfit, a[1000], y0[1000];
    int32 nvars;
    int32 npars;
    int32 ok;
    char collist[30];
    char parlist1[30];
    char parlist2[30];
    char varlist[30];

    static char *n[] = {"File",  "Fitvar", "Params", "Tolerance", "Npts",
                        "NCols", "To Col", "Params", "Epsilon",   "Max iter"};
    int32 status;
    char values[LENGTH(n)][MAX_LEN_SBOX];

    fit_info.nvars = 0;
    fit_info.npars = 0;

    // do fit get params
    snprintf(values[0], sizeof(values[0]), "%s", fit_info.file);
    snprintf(values[1], sizeof(values[1]), "%s", fit_info.varlist);
    snprintf(values[2], sizeof(values[2]), "%s", fit_info.parlist1);
    snprintf(values[3], sizeof(values[3]), "%g", fit_info.tol);
    snprintf(values[4], sizeof(values[4]), "%d", fit_info.npts);
    snprintf(values[5], sizeof(values[5]), "%d", fit_info.dim);
    snprintf(values[6], sizeof(values[6]), "%s", fit_info.collist);
    snprintf(values[7], sizeof(values[7]), "%s", fit_info.parlist2);
    snprintf(values[8], sizeof(values[8]), "%g", fit_info.eps);
    snprintf(values[9], sizeof(values[9]), "%d", fit_info.maxiter);

    status = pop_list_do_string_box(10, 5, 2, "Fit", n, values, 45);
    if (status == 0) {
        return;
    }

    fit_info.tol = atof(values[3]);
    fit_info.npts = atoi(values[4]);
    fit_info.dim = atoi(values[5]);
    fit_info.eps = atof(values[8]);
    fit_info.maxiter = atoi(values[9]);
    strncpy(fit_info.file, values[0], sizeof(fit_info.file));
    strncpy(fit_info.varlist, values[1], sizeof(fit_info.varlist));
    strncpy(fit_info.parlist1, values[2], sizeof(fit_info.varlist));
    strncpy(fit_info.collist, values[6], sizeof(fit_info.varlist));
    strncpy(fit_info.parlist2, values[7], sizeof(fit_info.varlist));

    strncpy(collist, fit_info.collist, sizeof(collist));
    strncpy(varlist, fit_info.varlist, sizeof(varlist));
    strncpy(parlist1, fit_info.parlist1, sizeof(parlist1));
    strncpy(parlist2, fit_info.parlist2, sizeof(parlist2));

    do_fit_parse_collist(collist, fit_info.icols, &nvars);

    if (nvars <= 0) {
        ggets_err_msg("No columns...");
        return;
    }
    fit_info.nvars = nvars;
    nvars = 0;
    do_fit_parse_varlist(varlist, fit_info.ivar, &nvars);

    if (fit_info.nvars != nvars) {
        ggets_err_msg(" # columns != # fitted variables");
        return;
    }
    npars = 0;
    do_fit_parse_parlist(parlist1, fit_info.ipar, &npars);

    do_fit_parse_parlist(parlist2, fit_info.ipar, &npars);

    if (npars <= 0) {
        ggets_err_msg(" No parameters!");
        return;
    }
    fit_info.npars = npars;
    for (int32 i = 0; i < npars; i++) {
        if (fit_info.ipar[i] >= 0) {
            if (fit_info.ipar[i] >= NODE) {
                ggets_err_msg(" Cant vary auxiliary/markov variables! ");
                return;
            }
        }
    }
    for (int32 i = 0; i < nvars; i++) {
        if (fit_info.icols[i] < 2) {
            ggets_err_msg(" Illegal column must be >= 2");
            return;
        }
        if (fit_info.ivar[i] < 0 || fit_info.ivar[i] >= NODE) {
            ggets_err_msg(" Fit only to variables! ");
            return;
        }
    }
    yfit = xmalloc((usize)(fit_info.npts*fit_info.nvars)*sizeof(*yfit));
    for (int32 i = 0; i < NODE; i++) {
        y0[i] = last_ic[i];
    }
    for (int32 i = 0; i < fit_info.npars; i++) {
        if (fit_info.ipar[i] < 0) {
            a[i] = constants[-fit_info.ipar[i]];
        } else {
            a[i] = last_ic[fit_info.ipar[i]];
        }
    }

    // do fit print info
    ggets_plintf("dim=%d maxiter=%d npts=%d file=%s tol=%g eps=%g\n",
                 fit_info.dim, fit_info.maxiter, fit_info.npts, fit_info.file,
                 fit_info.tol, fit_info.eps);

    for (int32 i = 0; i < fit_info.nvars; i++) {
        ggets_plintf(" variable %d to col %d \n", fit_info.ivar[i],
                     fit_info.icols[i]);
    }
    for (int32 i = 0; i < fit_info.npars; i++) {
        ggets_plintf(" P[%d]=%d \n", i, fit_info.ipar[i]);
    }

    ggets_plintf(" Running the fit...\n");
    ok =
        do_fit_run(fit_info.file, fit_info.npts, fit_info.npars, fit_info.nvars,
                   fit_info.maxiter, fit_info.dim, fit_info.eps, fit_info.tol,
                   fit_info.ipar, fit_info.ivar, fit_info.icols, y0, a, yfit);

    free(yfit);
    if (ok == 0) {
        return;
    }

    // get the latest par values ...

    for (int32 i = 0; i < npars; i++) {
        if (fit_info.ipar[i] < 0) {
            constants[-fit_info.ipar[i]] = a[i];
        } else {
            last_ic[fit_info.ipar[i]] = a[i];
        }
    }
    return;
}

/*
   filename is where the data file is -- it is of the form:
   t1 y11 y12 .... y1m
   t2 ....
   ...
   tn yn1 ....     ynm
   icols gives the dependent variable columns -- we assume first col
   is the times
   ndim is the number of y-pts in the a row

*/
int32
do_fit_run(  // double arrays
    char *filename, int32 npts, int32 npars, int32 nvars, int32 maxiter,
    int32 ndim, double eps, double tol, int32 *ipar, int32 *ivar, int32 *icols,
    double *y0, double *a, double *yfit) {
    double *t0, *y, sig[MAX_ODE], *covar, *alpha, chisq, ochisq, alambda,
        **yderv, *work;
    int32 ioff;
    int32 ictrl = 0;
    int32 ok = 0;
    FILE *fp;
    int32 niter = 0;
    int32 good_flag = 0;
    double tol10 = 10*tol;
    double t;
    double ytemp[MAX_ODE];
    /*printf(" %s %d %d %d %d %d \n",
              filename,
            npts,npars,nvars,maxiter,ndim); */

    if ((fp = fopen(filename, "r")) == NULL) {
        ggets_err_msg("No such file...");
        return 0;
    }
    t0 = xmalloc((usize)(npts + 1)*sizeof(*(t0)));
    y = xmalloc((usize)((npts + 1)*nvars)*sizeof(*y));
    // load up the data to fit

    for (int32 i = 0; i < npts; i++) {
        fscanf(fp, "%lg ", &t);

        for (int32 j = 0; j < ndim - 1; j++) {
            fscanf(fp, "%lg ", &ytemp[j]);
        }
        t0[i] = t;

        ioff = nvars*i;
        for (int32 k = 0; k < nvars; k++) {
            y[ioff + k] = ytemp[icols[k] - 2];
        }
    }
    ggets_plintf(" Data loaded ... %f %f ...  %f %f \n", y[0], y[1],
                 y[npts*nvars - 2], y[npts*nvars - 1]);

    work = xmalloc(sizeof(*work)*(usize)(4*npars + npars*npars));
    yderv = xmalloc((usize)npars*sizeof(double *));
    for (int32 i = 0; i < npars; i++) {
        yderv[i] = xmalloc((usize)((npts + 1)*nvars)*sizeof(*(yderv[i])));
    }
    for (int32 i = 0; i < nvars; i++) {
        sig[i] = 1.0;
    }

    covar = xmalloc((usize)(npars*npars)*sizeof(*covar));
    alpha = xmalloc((usize)(npars*npars)*sizeof(*alpha));

    while (good_flag < 3) {  // take 3 good steps after convergence

        ok = do_fit_marlev_step(t0, y0, y, sig, a, npts, nvars, npars, ivar,
                                ipar, covar, alpha, &chisq, &alambda, work,
                                yderv, yfit, &ochisq, ictrl, eps);
        niter++;
        ggets_plintf(" step %d is %d  -- lambda= %g  chisq= %g oldchi= %g\n",
                     niter, ok, alambda, chisq, ochisq);
        ggets_plintf(" params: ");
        for (int32 i = 0; i < npars; i++) {
            ggets_plintf(" %g ", a[i]);
        }
        ggets_plintf("\n");
        if ((ok == 0) || (niter >= maxiter)) {
            break;
        }
        if (ochisq > chisq) {
            if (((ochisq - chisq) < tol10) ||
                (((ochisq - chisq) / MAX(1.0, chisq)) < tol)) {
                good_flag++;
                niter--;  // compensate for good stuff ...
            }
            ochisq = chisq;
        } else {
            chisq = ochisq;
        }

        ictrl = 1;
    }

    if (ok == 0) {
        ggets_err_msg("Error in step...");

        free(work);
        for (int32 i = 0; i < npars; i++) {
            free(yderv[i]);
        }
        free(yderv);
        free(alpha);
        free(covar);
        free(t0);
        free(y);

        return 0;
    }
    if (niter >= maxiter) {
        ggets_err_msg("Max iterations exceeded...");

        free(work);
        for (int32 i = 0; i < npars; i++) {
            free(yderv[i]);
        }
        free(yderv);
        free(alpha);
        free(covar);
        free(t0);
        free(y);

        return 1;
    }
    ictrl = 2;
    do_fit_marlev_step(t0, y0, y, sig, a, npts, nvars, npars, ivar, ipar, covar,
                       alpha, &chisq, &alambda, work, yderv, yfit, &ochisq,
                       ictrl, eps);
    ggets_err_msg(" Success! ");
    // have the covariance matrix -- so what?
    ggets_plintf(" covariance: \n");
    for (int32 i = 0; i < npars; i++) {
        for (int32 j = 0; j < npars; j++) {
            ggets_plintf(" %g ", covar[i + npars*j]);
        }
        ggets_plintf("\n");
    }

    free(work);
    for (int32 i = 0; i < npars; i++) {
        free(yderv[i]);
    }
    free(yderv);
    free(alpha);
    free(covar);
    free(t0);
    free(y);

    return 1;
}

/*   One step of Levenberg-Marquardt

nvars  the number of variables to fit
ivar   their indices
npars  the number of parameters to alter
ipar   their indices
npts   the number of times
ictrl  0 to start  1 to continue  2 to finish up

t0  the npts times
y0  the NODE initial data
y   the (npts)*nvars data points to fit
sig  the nvars  weights
a  the npars initial guesses of the things to be fit
chisq  the chisquare
alpha is work array and also the curvature matrix
covar is the covariance matrix (npars x npars)
alambda is control of step size; start negative 0 to get final value
work is an array of size npar*4+npar*npar
yderv is  npar x (nptts+1)*nvars
yfit  is  (npts)*nvars  on each completed step it has the fitted soln
eps   control numerical derivative
sigma  weights on nvars
*/
int32
do_fit_marlev_step(double *t0, double *y0, double *y, double *sig, double *a,
                   int32 npts, int32 nvars, int32 npars, int32 *ivar,
                   int32 *ipar, double *covar, double *alpha, double *chisq,
                   double *alambda, double *work, double **yderv, double *yfit,
                   double *ochisq, int32 ictrl, double eps) {
    int32 ierr;
    int32 ipivot[1000];

    double *da, *atry, *beta, *oneda;
    da = work;
    atry = work + npars;
    beta = work + 2*npars;
    oneda = work + 3*npars;

    if (ictrl == 0) {
        *alambda = .001;
        if (do_fit_mrqcof(t0, y0, y, sig, a, npts, nvars, npars, ivar, ipar,
                          alpha, chisq, beta, yderv, yfit, eps) == 0) {
            return 0;
        }
        for (int32 i = 0; i < npars; i++) {
            atry[i] = a[i];
        }
        *ochisq = (*chisq);
    }
    for (int32 j = 0; j < npars; j++) {
        for (int32 k = 0; k < npars; k++) {
            covar[j + k*npars] = alpha[j + k*npars];
        }
        covar[j + j*npars] = alpha[j + j*npars]*(1 + (*alambda));
        oneda[j] = beta[j];
    }
    gear_sgefa(covar, npars, npars, ipivot, &ierr);
    if (ierr != -1) {
        ggets_err_msg(" Singular matrix encountered...");
        return 0;
    }

    gear_sgesl(covar, npars, npars, ipivot, oneda, 0);
    for (int32 j = 0; j < npars; j++) {
        da[j] = oneda[j];
    }
    if (ictrl == 2) {  // all done invert alpha to get the covariance
        for (int32 j = 0; j < (npars*npars); j++) {
            alpha[j] = covar[j];
        }
        for (int32 j = 0; j < npars; j++) {
            for (int32 k = 0; k < npars; k++) {
                oneda[k] = 0.0;
            }
            oneda[j] = 1.0;
            gear_sgesl(alpha, npars, npars, ipivot, oneda, 0);
            for (int32 k = 0; k < npars; k++) {
                covar[j + k*npars] = oneda[k];
            }
        }
        return 1;
    }
    for (int32 j = 0; j < npars; j++) {
        atry[j] = a[j] + da[j];
    }
    if (do_fit_mrqcof(t0, y0, y, sig, atry, npts, nvars, npars, ivar, ipar,
                      covar, chisq, da, yderv, yfit, eps) == 0) {
        return 0;
    }

    if (*chisq < *ochisq) {
        // *ochisq=*chisq;
        *alambda *= 0.1;
        for (int32 j = 0; j < npars; j++) {
            for (int32 k = 0; k < npars; k++) {
                alpha[j + k*npars] = covar[j + k*npars];
            }
            beta[j] = da[j];
            a[j] = atry[j];
        }
    } else {
        *alambda *= 10.0;
        // *chisq=*ochisq;
    }
    return 1;
}

int32
do_fit_mrqcof(double *t0, double *y0, double *y, double *sig, double *a,
              int32 npts, int32 nvars, int32 npars, int32 *ivar, int32 *ipar,
              double *alpha, double *chisq, double *beta, double **yderv,
              double *yfit, double eps) {
    int32 flag;
    int32 l;
    int32 k0;
    double sig2i;
    double dy;
    double wt;

    do_fit_get_info(y0, a, t0, &flag, eps, yfit, yderv, npts, npars, nvars,
                    ivar, ipar);
    if (flag == 0) {
        ggets_err_msg(" Integration error ...\n");
        return 0;
    }
    for (int32 i = 0; i < npars; i++) {
        beta[i] = 0.0;
        for (int32 j = 0; j < npars; j++) {
            alpha[i + j*npars] = 0.0;
        }
    }
    *chisq = 0.0;
    for (int32 i = 0; i < nvars; i++) {
        sig2i = 1.0 / (sig[i]*sig[i]);
        for (int32 k = 0; k < npts; k++) {
            k0 = k*nvars + i;
            dy = y[k0] - yfit[k0];
            for (int32 j = 0; j < npars; j++) {
                wt = yderv[j][k0]*sig2i;
                for (l = 0; l < npars; l++) {
                    alpha[j + l*npars] += wt*yderv[l][k0];
                }
                beta[j] += dy*wt;
            }
            (*chisq) += dy*dy*sig2i;

            /* the last loop could be halved because of symmetry, but I am lazy
               and this is an insignificiant amount of the CPU time since
              the evaluation step is really where all the time is used
            */
        }
    }
    return 1;
}

/* gets a list of the data columns to use ... */

void
do_fit_parse_collist(char *collist, int32 *icols, int32 *n) {
    char *item;
    int32 v;
    int32 i = 0;

    item = form_ode_get_first(collist, " ,");

    if (item[0] == 0) {
        return;
    }
    v = atoi(item);
    icols[i] = v;
    i++;
    while ((item = form_ode_do_fit_get_next(" ,")) != NULL) {
        v = atoi(item);
        icols[i] = v;
        i++;
    }
    *n = i;
    return;
}

void
do_fit_parse_varlist(char *varlist, int32 *ivars, int32 *n) {
    char *item;
    int32 v;
    int32 i = 0;

    item = form_ode_get_first(varlist, " ,");
    if (item[0] == 0) {
        return;
    }
    browser_find_variable(item, &v);
    if (v <= 0) {
        return;
    }
    ivars[i] = v - 1;
    i++;
    while ((item = form_ode_do_fit_get_next(" ,")) != NULL) {
        browser_find_variable(item, &v);
        if (v <= 0) {
            return;
        }
        ivars[i] = v - 1;
        i++;
    }
    *n = i;
    return;
}

void
do_fit_parse_parlist(char *parlist, int32 *ipars, int32 *n) {
    char *item;
    int32 v;
    int32 i = 0;
    usize j;
    for (j = 0; j < strlen(parlist); j++) {
        if (parlist[j] != ' ') {
            break;
        }
    }
    if (j == strlen(parlist)) {
        return;
    }
    if (strlen(parlist) == 0) {
        return;
    }
    item = form_ode_get_first(parlist, " ,");
    if (item[0] == 0L) {
        return;
    }

    browser_find_variable(item, &v);
    if (v > 0) {
        ipars[i + *n] = v - 1;
        i++;
    } else {
        v = get_param_index(item);
        if (v <= 0) {
            return;
        }
        ipars[i + *n] = -v;
        i++;
    }
    while ((item = form_ode_do_fit_get_next(" ,")) != NULL) {
        browser_find_variable(item, &v);
        if (v > 0) {
            ipars[i + *n] = v - 1;
            i++;
        } else {
            v = get_param_index(item);
            if (v <= 0) {
                return;
            }
            ipars[i + *n] = -v;
            i++;
        }
    }
    *n = *n + i;
}
