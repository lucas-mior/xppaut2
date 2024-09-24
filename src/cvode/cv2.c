#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "cv2.h"
#include "cvband.h"
#include "cvdense.h"
#include "cvode.h"
#include "functions.h"
#include "integers.h"
#include "llnltyps.h"
#include "vector.h"

static double cv_ropt[OPT_SIZE];
static int32 cv_iopt[OPT_SIZE];

static void cv_func(int64 n, double t, N_Vector y, N_Vector ydot, void *fdata);
static void *cvode_mem;
static N_Vector ycv;

void
start_cv(double *y, double t, int32 n, double *atol, double *rtol) {
    ycv = N_VNew(n);
    for (int32 i = 0; i < n; i++)
        ycv->data[i] = y[i];
    cvode_mem = cvode_malloc(n, cv_func, t, ycv, BDF, NEWTON, SS, rtol, atol,
                             NULL, NULL, FALSE, cv_iopt, cv_ropt);
    if (cv_bandflag == 1)
        CVBand(cvode_mem, cv_bandupper, cv_bandlower, NULL, NULL);
    else
        cv_dense(cvode_mem, NULL, NULL);
    return;
}

void
cv_end(void) {
    N_VFree(ycv);
    cvode_free(cvode_mem);
    return;
}

void
cv_func(int64 n, double t, N_Vector y, N_Vector ydot, void *fdata) {
    (void)fdata;
    my_rhs(t, y->data, ydot->data, (int32)n);
    return;
}

void
cvode_err_msg(int32 kflag) {
    char s[256];
    strcpy(s, "");
    switch (kflag) {
    case 0:
        strcpy(s, "");
        break;
    case -1:
        strcpy(s, "No memory allocated");
        break;
    case -2:
        strcpy(s, "Bad input to CVode");
        break;
    case -3:
        strcpy(s, "Too much work -- try smaller DT");
        break;
    case -4:
        sprintf(s, "Tolerance too low-- try TOL=%g ATOL=%g",
                TOLER*cv_ropt[TOLSF], ATOLER*cv_ropt[TOLSF]);
        break;
    case -5:
        strcpy(s, "Error test failure too frequent ??");
        break;
    case -6:
        strcpy(s, "Converg. failure -- oh well!");
        break;
    case -7:
        strcpy(s, "Setup failed for linsolver in CVODE ???");
        break;
    case -8:
        strcpy(s, "Singular matrix encountered. Hmmm?");
        break;
    case -9:
        strcpy(s, "Flags error...");
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    if (strlen(s) > 0)
        ggets_err_msg(s);
    return;
}

int32
cvode(
    /* command =0 continue, 1 is start 2 finish */
    int32 *command, double *y, double *t, int32 n, double tout, int32 *kflag,
    double *atol, double *rtol) {
    int32 err = 0;
    if (NFlags == 0)
        return ccvode(command, y, t, n, tout, kflag, atol, rtol);
    err = one_flag_step_cvode(command, y, t, n, tout, kflag, atol, rtol);
    if (err == 1)
        *kflag = -9;
    return 1;
}
/* rtol is like our TOLER and atol is something else ?? */
int32
ccvode(
    /* command =0 continue, 1 is start 2 finish */
    int32 *command, double *y, double *t, int32 n, double tout, int32 *kflag,
    double *atol, double *rtol) {
    int32 flag;
    *kflag = 0;
    if (*command == 2) {
        cv_end();
        return 1;
    }
    if (*command == 1) {
        start_cv(y, *t, n, atol, rtol);
        flag = CVode(cvode_mem, tout, ycv, t, NORMAL);
        if (flag != SUCCESS) {
            *kflag = flag;
            cv_end();
            *command = 1;
            return -1;
        }
        *command = 0;
        for (int32 i = 0; i < n; i++)
            y[i] = ycv->data[i];
        return 0;
    }
    flag = CVode(cvode_mem, tout, ycv, t, NORMAL);
    if (flag != SUCCESS) {
        *kflag = flag;
        cv_end();
        *command = 1;

        return -1;
    }
    for (int32 i = 0; i < n; i++)
        y[i] = ycv->data[i];
    return 0;
}
