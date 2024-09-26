#include "functions.h"
#include "parserslow.h"

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "integers.h"
#define MAXDAE 400

/*    will have more stuff someday */

static struct DaeWork {
    double *work;
    int32 *iwork;
    int32 status;
} dae_work;

static struct Svar {
    char name[12], *rhs;
    int32 *form;
    int32 index;
    double value;
    double last;
} svar[MAXDAE];

static struct Aeqn {
    char *rhs;
    int32 *form;
} aeqn[MAXDAE];

static int32 nsvar = 0, naeqn = 0;

/* this adds an algebraically defined variable  and a formula
   for the first guess */

static void get_dae_fun(double *y, double *f);

int32
dae_fun_add_svar(char *name, char *rhs) {
    if (nsvar >= MAXDAE) {
        ggets_plintf(" Too many variables\n");
        return 1;
    }

    strcpy(svar[nsvar].name, name);
    svar[nsvar].rhs = xmalloc(80);
    strcpy(svar[nsvar].rhs, rhs);
    ggets_plintf(" Added sol-var[%d] %s = %s \n", nsvar, svar[nsvar].name,
                 svar[nsvar].rhs);
    nsvar++;
    return 0;
}

/* adds algebraically define name to name list */

int32
dae_fun_add_svar_names(void) {
    for (int32 i = 0; i < nsvar; i++) {
        svar[i].index = NVAR;
        if (add_var(svar[i].name, 0.0) == 1)
            return 1;
    }
    return 0;
}

/* adds a right-hand side to slove for zero */

int32
dae_fun_add_aeqn(char *rhs) {
    if (naeqn >= MAXDAE) {
        ggets_plintf(" Too many equations\n");
        return 1;
    }
    aeqn[naeqn].rhs = xmalloc(strlen(rhs) + 5);
    strcpy(aeqn[naeqn].rhs, rhs);
    naeqn++;
    return 0;
}

/* this compiles formulas to set to zero */
int32
dae_fun_compile_svars(void) {
    int32 f[256], n, k;
    if (nsvar != naeqn) {
        ggets_plintf(" #SOL_VAR(%d) must equal #ALG_EQN(%d) ! \n", nsvar,
                     naeqn);
        return 1;
    }

    for (int32 i = 0; i < naeqn; i++) {
        if (add_expr(aeqn[i].rhs, f, &n) == 1) {
            ggets_plintf(" Bad right-hand side for alg-eqn \n");
            return 1;
        }
        aeqn[i].form = xmalloc(sizeof(*(aeqn[i].form))*(usize)(n + 2));
        for (k = 0; k < n; k++)
            aeqn[i].form[k] = f[k];
    }

    for (int32 i = 0; i < nsvar; i++) {
        if (add_expr(svar[i].rhs, f, &n) == 1) {
            ggets_plintf(" Bad initial guess for sol-var \n");
            return 1;
        }
        svar[i].form = xmalloc(100*sizeof(*(svar[i].form)));
        for (k = 0; k < n; k++)
            svar[i].form[k] = f[k];
    }

    /* dae fun init work */
    dae_work.work =
        xmalloc(sizeof(*(dae_work.work))*(usize)(nsvar*nsvar + 10*nsvar));
    dae_work.iwork = xmalloc(sizeof(*(dae_work.iwork))*(usize)nsvar);
    dae_work.status = 1;
    return 0;
}

void
dae_fun_reset_dae(void) {
    dae_work.status = 1;
    return;
}

void
dae_fun_set_init_guess(void) {
    double z;
    dae_work.status = 1;
    if (nsvar == 0)
        return;
    for (int32 i = 0; i < nsvar; i++) {
        z = evaluate(svar[i].form);
        SETVAR(svar[i].index, z);
        svar[i].value = z;
        svar[i].last = z;
    }
    return;
}

void
dae_fun_err_dae(void) {
    switch (dae_work.status) {
    case 2:
        ggets_err_msg(" Warning - no change in Iterates");
        break;
    case -1:
        ggets_err_msg(" Singular jacobian for dae\n");

        break;
    case -2:
        ggets_err_msg(" Maximum iterates exceeded for dae\n");

        break;
    case -3:
        ggets_err_msg(" Newton update out of bounds\n");
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    dae_work.status = 1;
    return;
}

void
get_dae_fun(double *y, double *f) {
    /* better do this in case fixed variables depend on sol_var */
    for (int32 i = 0; i < nsvar; i++)
        SETVAR(svar[i].index, y[i]);
    for (int32 i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    for (int32 i = 0; i < naeqn; i++)
        f[i] = evaluate(aeqn[i].form);
    return;
}

void
dae_fun_do_daes(void) {
    int32 ans;

    /* dae fun solve */
    /* Newton solver for algebraic stuff */
    int32 j, n;
    int32 info;
    double err, del, z, yold;
    double tol = EVEC_ERR, eps = NEWT_ERR;
    int32 maxit = EVEC_ITER, iter = 0;
    double *y, *ynew, *f, *fnew, *jac, *errvec;
    n = nsvar;
    if (nsvar == 0)
        ans = 1;
    if (dae_work.status < 0)
        ans = dae_work.status; /* accepts no change error */
    y = dae_work.work;
    f = y + nsvar;
    fnew = f + nsvar;
    ynew = fnew + nsvar;
    errvec = ynew + nsvar;
    jac = errvec + nsvar;
    for (int32 i = 0; i < n; i++) { /* copy current value as initial guess */
        y[i] = svar[i].last;
        ynew[i] = y[i]; /* keep old guess */
    }
    while (true) {
        get_dae_fun(y, f);
        err = 0.0;
        for (int32 i = 0; i < n; i++) {
            err += fabs(f[i]);
            errvec[i] = f[i];
        }
        if (err < tol) { /* success */
            for (int32 i = 0; i < n; i++) {
                SETVAR(svar[i].index, y[i]);
                svar[i].last = y[i];
            }
            ans = 1;
            break;
        }
        /* compute jacobian */
        for (int32 i = 0; i < n; i++) {
            z = fabs(y[i]);
            if (z < eps)
                z = eps;
            del = eps*z;
            yold = y[i];
            y[i] = y[i] + del;
            get_dae_fun(y, fnew);
            for (j = 0; j < n; j++)
                jac[j*n + i] = (fnew[j] - f[j]) / del;
            y[i] = yold;
        }
        gear_sgefa(jac, n, n, dae_work.iwork, &info);
        if (info != -1) {
            for (int32 i = 0; i < n; i++)
                SETVAR(svar[i].index, ynew[i]);
            ans = -1; /* singular jacobian */
            break;
        }
        gear_sgesl(jac, n, n, dae_work.iwork, errvec, 0); /* get x=J^(-1) f */
        err = 0.0;
        for (int32 i = 0; i < n; i++) {
            y[i] -= errvec[i];
            err += fabs(errvec[i]);
        }
        if (err > (n*BOUND)) {
            for (int32 i = 0; i < n; i++)
                SETVAR(svar[i].index, svar[i].last);
            ans = -3; /* getting too big */
            break;
        }
        if (err < tol) /* not much change */
        {
            for (int32 i = 0; i < n; i++) {
                SETVAR(svar[i].index, y[i]);
                svar[i].last = y[i];
            }
            ans = 2;
            break;
        }
        iter++;
        if (iter > maxit) {
            for (int32 i = 0; i < n; i++)
                SETVAR(svar[i].index, svar[i].last);
            ans = -2; /* too many iterates */
            break;
        }
    }
    dae_work.status = ans;
    if (ans == 1 || ans == 2)
        return; /* accepts a no change error! */
    DelayErr = 1;
    return;
}

/* interface shit -- different for Win95 */

void
dae_fun_get_new_guesses(void) {
    int32 n;
    double z;
    if (nsvar < 1)
        return;
    for (int32 i = 0; i < nsvar; i++) {
        char name[sizeof(svar[i].name) + 23];

        z = svar[i].last;
        snprintf(name, sizeof(name), "Initial %s(%g):", svar[i].name, z);
        ggets_new_string(name, svar[i].rhs);
        if (add_expr(svar[i].rhs, svar[i].form, &n)) {
            ggets_err_msg("Illegal formula");
            return;
        }
        z = evaluate(svar[i].form);
        SETVAR(svar[i].index, z);
        svar[i].value = z;
        svar[i].last = z;
    }
}
