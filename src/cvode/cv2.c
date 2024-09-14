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

double cv_ropt[OPT_SIZE];
int32 cv_iopt[OPT_SIZE];
extern int32 cv_bandflag, cv_bandupper, cv_bandlower;

static void cvf(int64 n, double t, N_Vector y, N_Vector ydot, void *fdata);
void *cvode_mem;
N_Vector ycv;
extern int32 NFlags;
extern double TOLER, ATOLER;

void
start_cv(double *y, double t, int32 n, double tout, double *atol,
         double *rtol) {
    int32 i;

    ycv = N_VNew(n);
    for (i = 0; i < n; i++)
        ycv->data[i] = y[i];
    cvode_mem = CVodeMalloc(n, cvf, t, ycv, BDF, NEWTON, SS, rtol, atol, NULL,
                            NULL, FALSE, cv_iopt, cv_ropt);
    if (cv_bandflag == 1)
        CVBand(cvode_mem, cv_bandupper, cv_bandlower, NULL, NULL);
    else
        CVDense(cvode_mem, NULL, NULL);
    return;
}

void
end_cv(void) {
    N_VFree(ycv);
    CVodeFree(cvode_mem);
    return;
}

void
cvf(int64 n, double t, N_Vector y, N_Vector ydot, void *fdata) {
    my_rhs(t, y->data, ydot->data, n);
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
                TOLER * cv_ropt[TOLSF], ATOLER * cv_ropt[TOLSF]);
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
    }
    if (strlen(s) > 0)
        err_msg(s);
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
    int32 i, flag;
    *kflag = 0;
    if (*command == 2) {
        end_cv();
        return 1;
    }
    if (*command == 1) {
        start_cv(y, *t, n, tout, atol, rtol);
        flag = CVode(cvode_mem, tout, ycv, t, NORMAL);
        if (flag != SUCCESS) {

            *kflag = flag;
            end_cv();
            *command = 1;
            return -1;
        }
        *command = 0;
        for (i = 0; i < n; i++)
            y[i] = ycv->data[i];
        return 0;
    }
    flag = CVode(cvode_mem, tout, ycv, t, NORMAL);
    if (flag != SUCCESS) {
        *kflag = flag;
        end_cv();
        *command = 1;

        return -1;
    }
    for (i = 0; i < n; i++)
        y[i] = ycv->data[i];
    return 0;
}
