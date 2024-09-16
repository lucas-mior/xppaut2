#include "functions.h"
#include "parserslow.h"

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "getvar.h"
#include "integers.h"
#define MAXDAE 400

extern double variables[];
extern int32 NVAR;

extern int32 DelayErr;

extern double EVEC_ERR, NEWT_ERR, BOUND;
extern int32 EVEC_ITER;

extern int32 NODE, FIX_VAR;
extern int32 *my_ode[];

/*    will have more stuff someday */

typedef struct {
    double *work;
    int32 *iwork;
    int32 status;
} DaeWork;
DaeWork dae_work;

typedef struct {
    char name[12], *rhs;
    int32 *form;
    int32 index;
    double value, last;
} SolVar;

typedef struct {
    char *rhs;
    int32 *form;
} DaeEqn;

SolVar svar[MAXDAE];
DaeEqn aeqn[MAXDAE];

int32 nsvar = 0, naeqn = 0;

/* this adds an algebraically defined variable  and a formula
   for the first guess */

int32
add_svar(char *name, char *rhs) {
    if (nsvar >= MAXDAE) {
        plintf(" Too many variables\n");
        return 1;
    }

    strcpy(svar[nsvar].name, name);
    svar[nsvar].rhs = malloc(80);
    strcpy(svar[nsvar].rhs, rhs);
    plintf(" Added sol-var[%d] %s = %s \n", nsvar, svar[nsvar].name,
           svar[nsvar].rhs);
    nsvar++;
    return 0;
}

/* adds algebraically define name to name list */

int32
add_svar_names(void) {
    int32 i;
    for (i = 0; i < nsvar; i++) {
        svar[i].index = NVAR;
        if (add_var(svar[i].name, 0.0) == 1)
            return 1;
    }
    return 0;
}

/* adds a right-hand side to slove for zero */

int32
add_aeqn(char *rhs) {
    if (naeqn >= MAXDAE) {
        plintf(" Too many equations\n");
        return 1;
    }
    aeqn[naeqn].rhs = malloc(strlen(rhs) + 5);
    strcpy(aeqn[naeqn].rhs, rhs);
    naeqn++;
    return 0;
}

/* this compiles formulas to set to zero */
int32
compile_svars(void) {
    int32 i, f[256], n, k;
    if (nsvar != naeqn) {
        plintf(" #SOL_VAR(%d) must equal #ALG_EQN(%d) ! \n", nsvar, naeqn);
        return 1;
    }

    for (i = 0; i < naeqn; i++) {
        if (add_expr(aeqn[i].rhs, f, &n) == 1) {
            plintf(" Bad right-hand side for alg-eqn \n");
            return 1;
        }
        aeqn[i].form = malloc(sizeof(*(aeqn[i].form))*(n + 2));
        for (k = 0; k < n; k++)
            aeqn[i].form[k] = f[k];
    }

    for (i = 0; i < nsvar; i++) {
        if (add_expr(svar[i].rhs, f, &n) == 1) {
            plintf(" Bad initial guess for sol-var \n");
            return 1;
        }
        svar[i].form = malloc(100*sizeof(*(svar[i].form)));
        for (k = 0; k < n; k++)
            svar[i].form[k] = f[k];
    }
    init_dae_work();
    return 0;
}

void
reset_dae(void) {
    dae_work.status = 1;
    return;
}

void
set_init_guess(void) {
    int32 i;
    double z;
    dae_work.status = 1;
    if (nsvar == 0)
        return;
    for (i = 0; i < nsvar; i++) {
        z = evaluate(svar[i].form);
        SETVAR(svar[i].index, z);
        svar[i].value = z;
        svar[i].last = z;
    }
    return;
}

void
err_dae(void) {
    switch (dae_work.status) {
    case 2:
        err_msg(" Warning - no change in Iterates");
        break;
    case -1:
        err_msg(" Singular jacobian for dae\n");

        break;
    case -2:
        err_msg(" Maximum iterates exceeded for dae\n");

        break;
    case -3:
        err_msg(" Newton update out of bounds\n");
        break;
    }
    dae_work.status = 1;
    return;
}

void
init_dae_work(void) {
    dae_work.work =
        malloc(sizeof(*(dae_work.work))*(nsvar*nsvar + 10*nsvar));
    dae_work.iwork = malloc(sizeof(*(dae_work.iwork))*nsvar);
    dae_work.status = 1;
    return;
}

void
get_dae_fun(double *y, double *f) {
    int32 i;
    /* better do this in case fixed variables depend on sol_var */
    for (i = 0; i < nsvar; i++)
        SETVAR(svar[i].index, y[i]);
    for (i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    for (i = 0; i < naeqn; i++)
        f[i] = evaluate(aeqn[i].form);
    return;
}

void
do_daes(void) {
    int32 ans;
    ans = solve_dae();
    dae_work.status = ans;
    if (ans == 1 || ans == 2)
        return; /* accepts a no change error! */
    DelayErr = 1;
    return;
}

/* Newton solver for algebraic stuff */
int32
solve_dae(void) {
    int32 i, j, n;
    int32 info;
    double err, del, z, yold;
    double tol = EVEC_ERR, eps = NEWT_ERR;
    int32 maxit = EVEC_ITER, iter = 0;
    double *y, *ynew, *f, *fnew, *jac, *errvec;
    n = nsvar;
    if (nsvar == 0)
        return 1;
    if (dae_work.status < 0)
        return dae_work.status; /* accepts no change error */
    y = dae_work.work;
    f = y + nsvar;
    fnew = f + nsvar;
    ynew = fnew + nsvar;
    errvec = ynew + nsvar;
    jac = errvec + nsvar;
    for (i = 0; i < n; i++) { /* copy current value as initial guess */
        y[i] = svar[i].last;
        ynew[i] = y[i]; /* keep old guess */
    }
    while (true) {
        get_dae_fun(y, f);
        err = 0.0;
        for (i = 0; i < n; i++) {
            err += fabs(f[i]);
            errvec[i] = f[i];
        }
        if (err < tol) { /* success */
            for (i = 0; i < n; i++) {
                SETVAR(svar[i].index, y[i]);
                svar[i].last = y[i];
            }
            return 1;
        }
        /* compute jacobian */
        for (i = 0; i < n; i++) {
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
        sgefa(jac, n, n, dae_work.iwork, &info);
        if (info != -1) {
            for (i = 0; i < n; i++)
                SETVAR(svar[i].index, ynew[i]);
            return -1; /* singular jacobian */
        }
        sgesl(jac, n, n, dae_work.iwork, errvec, 0); /* get x=J^(-1) f */
        err = 0.0;
        for (i = 0; i < n; i++) {
            y[i] -= errvec[i];
            err += fabs(errvec[i]);
        }
        if (err > (n*BOUND)) {
            for (i = 0; i < n; i++)
                SETVAR(svar[i].index, svar[i].last);
            return -3; /* getting too big */
        }
        if (err < tol) /* not much change */
        {
            for (i = 0; i < n; i++) {
                SETVAR(svar[i].index, y[i]);
                svar[i].last = y[i];
            }
            return 2;
        }
        iter++;
        if (iter > maxit) {
            for (i = 0; i < n; i++)
                SETVAR(svar[i].index, svar[i].last);
            return -2; /* too many iterates */
        }
    }
    exit(EXIT_FAILURE);
}

/* interface shit -- different for Win95 */

void
get_new_guesses(void) {
    int32 i, n;
    char name[30];
    double z;
    if (nsvar < 1)
        return;
    for (i = 0; i < nsvar; i++) {
        z = svar[i].last;
        sprintf(name, "Initial %s(%g):", svar[i].name, z);
        new_string(name, svar[i].rhs);
        if (add_expr(svar[i].rhs, svar[i].form, &n)) {
            err_msg("Illegal formula");
            return;
        }
        z = evaluate(svar[i].form);
        SETVAR(svar[i].index, z);
        svar[i].value = z;
        svar[i].last = z;
    }
}
