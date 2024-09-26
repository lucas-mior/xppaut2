#include "functions.h"
#include "integers.h"
#include <math.h>
#include "xpplim.h"
#include <stdbool.h>

int32 (*rhs_function)(double t, double *y, double *ydot, int32 neq);

static double coefp[] = {6.875 / 3.00, -7.375 / 3.00, 4.625 / 3.00, -.375};
static double coefc[] = {.375, 2.375 / 3.00, -.625 / 3.00, 0.125 / 3.00};
static double *y_s[4];
static double *y_p[4];
static double *ypred;

static double symp_b[] = {7 / 24., .75, -1. / 24};
static double symp_B[] = {2 / 3., -2. / 3., 1.0};

/* my first symplectic integrator */

static void odesol_bandsol(double *a, double *b, int32 ml, int32 mr, int32 n);
static int32 odesol_bandfac(double *a, int32 ml, int32 mr, int32 n);
static int32 odesol_abmpc(double *y, double *t, double dt, int32 neq);

int32
odesol_symplect3(double *y, double *tim, double dt, int32 nt, int32 neq,
                 int32 *istart, double *work) {
    if (NFlags == 0) {
        for (int32 i = 0; i < nt; i++) {
            odesol2_one_step_symp(y, dt, work, neq, tim);
        }
        delay_handle_stor_delay(y);
        return 0;
    }
    for (int32 i = 0; i < nt; i++) {
        one_flag_step_symp(y, dt, work, neq, tim, istart);
        delay_handle_stor_delay(y);
    }
    return 0;
}

/*   DISCRETE    */

int32
odesol_discrete(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work) {
    if (NFlags == 0) {
        for (int32 i = 0; i < nt; i++) {
            odesol2_one_step_discrete(y, dt, work, neq, tim);
            delay_handle_stor_delay(y);
        }
        return 0;
    }
    for (int32 i = 0; i < nt; i++) {
        one_flag_step_discrete(y, dt, work, neq, tim, istart);
        delay_handle_stor_delay(y);
    }
    return 0;
}

/* Backward Euler  */

int32
odesol_bak_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                 int32 *istart, double *work) {
    int32 j;
    double *jac, *yg, *yp, *yp2, *ytemp, *errvec;
    yp = work;
    yg = yp + neq;
    ytemp = yg + neq;
    errvec = ytemp + neq;
    yp2 = errvec + neq;
    jac = yp2 + neq;
    if (NFlags == 0) {
        for (int32 i = 0; i < nt; i++) {
            if ((j = odesol_one_bak_step(y, tim, dt, neq, yg, yp, yp2, ytemp,
                                         errvec, jac)) != 0)
                return j;
            delay_handle_stor_delay(y);
        }
        return 0;
    }
    for (int32 i = 0; i < nt; i++) {
        if ((j = one_flag_step_backeul(y, tim, dt, neq, yg, yp, yp2, ytemp,
                                       errvec, jac, istart)) != 0)
            return j;
        delay_handle_stor_delay(y);
    }
    return 0;
}

int32
odesol_one_bak_step(double *y, double *t, double dt, int32 neq, double *yg,
                    double *yp, double *yp2, double *ytemp, double *errvec,
                    double *jac) {
    double err = 0.0, err1 = 0.0;

    int32 iter = 0, info, ipivot[MAX_ODE1];
    int32 ml = cv_bandlower, mr = cv_bandupper, mt = ml + mr + 1;
    markov_set_wieners(dt, y, *t);
    *t = *t + dt;
    rhs_function(*t, y, yp2, neq);
    for (int32 i = 0; i < neq; i++)
        yg[i] = y[i];
    while (true) {
        err1 = 0.0;
        err = 0.0;
        rhs_function(*t, yg, yp, neq);
        for (int32 i = 0; i < neq; i++) {
            errvec[i] = yg[i] - .5*dt*(yp[i] + yp2[i]) - y[i];
            err1 += fabs(errvec[i]);
            ytemp[i] = yg[i];
        }
        odesol_get_the_jac(*t, yg, yp, ytemp, jac, neq, NEWT_ERR, -.5*dt);
        if (cv_bandflag) {
            for (int32 i = 0; i < neq; i++)
                jac[i*mt + ml] += 1;
            odesol_bandfac(jac, ml, mr, neq);
            odesol_bandsol(jac, errvec, ml, mr, neq);
        } else {
            for (int32 i = 0; i < neq; i++)
                jac[i*neq + i] += 1.0;
            gear_sgefa(jac, neq, neq, ipivot, &info);
            if (info != -1) {
                return -1;
            }
            gear_sgesl(jac, neq, neq, ipivot, errvec, 0);
        }
        for (int32 i = 0; i < neq; i++) {
            err += fabs(errvec[i]);
            yg[i] -= errvec[i];
        }
        if (err < euler_tol || err1 < euler_tol) {
            for (int32 i = 0; i < neq; i++)
                y[i] = yg[i];
            return 0;
        }
        iter++;
        if (iter > euler_max_iter)
            return -2;
    }
}

void
odesol2_one_step_discrete(double *y, double dt, double *yp, int32 neq,
                          double *t) {
    markov_set_wieners(dt, y, *t);
    rhs_function(*t, y, yp, neq);
    *t = *t + dt;
    for (int32 j = 0; j < neq; j++) {
        y[j] = yp[j];
    }
    return;
}

void
odesol2_one_step_symp(double *y, double h, double *f, int32 n, double *t) {
    int32 s;
    for (s = 0; s < 3; s++) {
        for (int32 j = 0; j < n; j += 2)
            y[j] += (h*symp_b[s]*y[j + 1]);
        rhs_function(*t, y, f, n);
        for (int32 j = 0; j < n; j += 2)
            y[j + 1] += (h*symp_B[s]*f[j + 1]);
    }
    *t += h;
    return;
}

void
odesol2_one_step_euler(double *y, double dt, double *yp, int32 neq, double *t) {
    markov_set_wieners(dt, y, *t);
    rhs_function(*t, y, yp, neq);
    *t += dt;
    for (int32 j = 0; j < neq; j++)
        y[j] = y[j] + dt*yp[j];
    return;
}

void
odesol_one_step_rk4(double *y, double dt, double *yval[3], int32 neq,
                    double *tim) {
    double t = *tim, t1, t2;
    markov_set_wieners(dt, y, t);
    rhs_function(t, y, yval[1], neq);
    for (int32 i = 0; i < neq; i++) {
        yval[0][i] = y[i] + dt*yval[1][i] / 6.00;
        yval[2][i] = y[i] + dt*yval[1][i]*0.5;
    }
    t1 = t + .5*dt;
    rhs_function(t1, yval[2], yval[1], neq);
    for (int32 i = 0; i < neq; i++) {
        yval[0][i] = yval[0][i] + dt*yval[1][i] / 3.00;
        yval[2][i] = y[i] + .5*dt*yval[1][i];
    }
    rhs_function(t1, yval[2], yval[1], neq);
    for (int32 i = 0; i < neq; i++) {
        yval[0][i] = yval[0][i] + dt*yval[1][i] / 3.000;
        yval[2][i] = y[i] + dt*yval[1][i];
    }
    t2 = t + dt;
    rhs_function(t2, yval[2], yval[1], neq);
    for (int32 i = 0; i < neq; i++)
        y[i] = yval[0][i] + dt*yval[1][i] / 6.00;
    *tim = t2;
    return;
}

void
odesol_one_step_heun(double *y, double dt, double *yval[2], int32 neq,
                     double *tim) {
    double t = *tim, t1;
    markov_set_wieners(dt, y, *tim);
    rhs_function(t, y, yval[0], neq);
    for (int32 i = 0; i < neq; i++)
        yval[0][i] = dt*yval[0][i] + y[i];
    t1 = t + dt;
    rhs_function(t1, yval[0], yval[1], neq);
    for (int32 i = 0; i < neq; i++)
        y[i] = .5*(y[i] + yval[0][i] + dt*yval[1][i]);
    *tim = t1;
    return;
}

/*  Euler  */

int32
odesol_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
             int32 *istart, double *work) {
    if (NFlags == 0) {
        for (int32 i = 0; i < nt; i++) {
            odesol2_one_step_euler(y, dt, work, neq, tim);
            delay_handle_stor_delay(y);
        }
        return 0;
    }
    for (int32 i = 0; i < nt; i++) {
        one_flag_step_euler(y, dt, work, neq, tim, istart);
        delay_handle_stor_delay(y);
    }
    return 0;
}

/* Modified Euler  */

int32
odesol_mod_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                 int32 *istart, double *work) {
    double *yval[2];

    yval[0] = work;
    yval[1] = work + neq;
    if (NFlags == 0) {
        for (int32 j = 0; j < nt; j++) {
            odesol_one_step_heun(y, dt, yval, neq, tim);
            delay_handle_stor_delay(y);
        }
        return 0;
    }
    for (int32 j = 0; j < nt; j++) {
        one_flag_step_heun(y, dt, yval, neq, tim, istart);
        delay_handle_stor_delay(y);
    }
    return 0;
}

/*  Runge Kutta    */

int32
odesol_rung_kut(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work) {
    double *yval[3];

    yval[0] = work;
    yval[1] = work + neq;
    yval[2] = work + neq + neq;

    if (NFlags == 0) {
        for (int32 j = 0; j < nt; j++) {
            odesol_one_step_rk4(y, dt, yval, neq, tim);
            delay_handle_stor_delay(y);
        }
        return 0;
    }

    for (int32 j = 0; j < nt; j++) {
        one_flag_step_rk4(y, dt, yval, neq, tim, istart);
        delay_handle_stor_delay(y);
    }
    return 0;
}

/*   ABM   */

int32
odesol_adams(double *y, double *tim, double dt, int32 nstep, int32 neq,
             int32 *ist, double *work) {
    int32 istart = *ist, istpst, ik, n;
    int32 irk;
    double *work1;
    double x0 = *tim, xst = *tim;
    work1 = work;
    if (istart == 1) {
        for (int32 i = 0; i < 4; i++) {
            y_p[i] = work + (4 + i)*neq;
            y_s[i] = work + (8 + i)*neq;
        }
        ypred = work + 3*neq;
        goto n20;
    }
    if (istart > 1)
        goto n350;
    istpst = 0;
    goto n400;

n20:

    x0 = xst;
    rhs_function(x0, y, y_p[3], neq);
    for (int32 k = 1; k < 4; k++) {
        odesol_rung_kut(y, &x0, dt, 1, neq, &irk, work1);
        delay_handle_stor_delay(y);
        for (int32 i = 0; i < neq; i++)
            y_s[3 - k][i] = y[i];
        rhs_function(x0, y, y_p[3 - k], neq);
    }
    istpst = 3;
    if (istpst <= nstep)
        goto n400;
    ik = 4 - nstep;
    for (int32 i = 0; i < neq; i++)
        y[i] = y_s[ik - 1][i];
    xst = xst + nstep*dt;
    istart = ik;
    goto n1000;

n350:

    ik = istart - nstep;
    if (ik <= 1)
        goto n370;
    for (int32 i = 0; i < neq; i++)
        y[i] = y_s[ik - 1][i];
    xst = xst + nstep*dt;
    istart = ik;
    goto n1000;

n370:
    for (int32 i = 0; i < neq; i++)
        y[i] = y_s[0][i];
    if (ik == 1) {
        x0 = xst + dt*nstep;
        goto n450;
    }

    istpst = istart - 1;

n400:

    if (istpst == nstep)
        goto n450;
    for (n = istpst + 1; n < nstep + 1; n++) {
        markov_set_wieners(dt, y, x0);
        odesol_abmpc(y, &x0, dt, neq);
        delay_handle_stor_delay(y);
    }

n450:
    istart = 0;
    xst = x0;

n1000:

    *tim = *tim + nstep*dt;
    *ist = istart;
    return 0;
}

int32
odesol_abmpc(double *y, double *t, double dt, int32 neq) {
    double x1, x0 = *t;
    for (int32 i = 0; i < neq; i++) {
        ypred[i] = 0;
        for (int32 k = 0; k < 4; k++)
            ypred[i] = ypred[i] + coefp[k]*y_p[k][i];
        ypred[i] = y[i] + dt*ypred[i];
    }

    for (int32 i = 0; i < neq; i++)
        for (int32 k = 3; k > 0; k--)
            y_p[k][i] = y_p[k - 1][i];
    x1 = x0 + dt;
    rhs_function(x1, ypred, y_p[0], neq);

    for (int32 i = 0; i < neq; i++) {
        ypred[i] = 0;
        for (int32 k = 0; k < 4; k++)
            ypred[i] = ypred[i] + coefc[k]*y_p[k][i];
        y[i] = y[i] + dt*ypred[i];
    }
    *t = x1;
    rhs_function(x1, y, y_p[0], neq);

    return 1;
}

/* this is rosen  - rosenbock step
    This uses banded routines as well */
int32
odesol_rb23(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
            double *work, int32 *ierr) {
    int32 out = -1;
    if (NFlags == 0) {
        out = odesol_rosen(y, tstart, tfinal, istart, n, work, ierr);
    } else {
        out = one_flag_step_rosen(y, tstart, tfinal, istart, n, work, ierr);
    }
    return out;
}

int32
odesol_rosen(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
             double *work, int32 *ierr) {
    static double htry;
    double epsjac = NEWT_ERR;
    double eps = 1e-15, hmin, hmax;
    double tdir = 1, t0 = *tstart, t = t0;
    double atol = ATOLER, rtol = TOLER;
    double sqrteps = sqrt(eps);
    double thresh = atol / rtol, absh, h = 0;
    double d = 1 / (2. + sqrt(2.)), e32 = 6. + sqrt(2.), tnew;
    /*double ninf;  Is this needed?*/
    int32 n2 = n*n, done = 0, info, ml = cv_bandlower, mr = cv_bandupper,
          mt = ml + mr + 1;
    int32 ipivot[MAX_ODE1], nofailed;
    double temp, err, tdel;
    double *ypnew, *k1, *k2, *k3, *f0, *f1, *f2, *dfdt, *ynew, *dfdy;
    *ierr = 1;
    ypnew = work;
    k1 = ypnew + n;
    k2 = k1 + n;
    k3 = k2 + n;
    f0 = k3 + n;
    f1 = f0 + n;
    f2 = f1 + n;
    dfdt = f2 + n;
    ynew = dfdt + n;
    dfdy = ynew + n;

    if (t0 > tfinal)
        tdir = -1;
    hmax = fabs(tfinal - t);
    if (*istart == 1)
        htry = hmax;
    rhs_function(t0, y, f0, n);
    hmin = 16*eps*fabs(t);
    absh = MIN(hmax, MAX(hmin, htry));
    while (!done) {
        nofailed = 1;
        hmin = 16*eps*fabs(t);
        absh = MIN(hmax, MAX(hmin, absh));
        h = tdir*absh;
        if (1.1*absh >= fabs(tfinal - t)) {
            h = tfinal - t;
            absh = fabs(h);
            done = 1;
        }
        odesol_get_the_jac(t, y, f0, ypnew, dfdy, n, epsjac, 1.0);
        tdel = (t + tdir*MIN(sqrteps*MAX(fabs(t), fabs(t + h)), absh)) - t;
        rhs_function(t + tdel, y, f1, n);
        for (int32 i = 0; i < n; i++)
            dfdt[i] = (f1[i] - f0[i]) / tdel;
        while (true) { /* advance a step  */
            for (int32 i = 0; i < n2; i++)
                dfdy[i] = -h*d*dfdy[i];
            for (int32 i = 0; i < n; i++)
                k1[i] = f0[i] + (h*d)*dfdt[i];
            if (cv_bandflag) {
                for (int32 i = 0; i < n; i++)
                    dfdy[i*mt + ml] += 1;

                odesol_bandfac(dfdy, ml, mr, n);
                odesol_bandsol(dfdy, k1, ml, mr, n);
            } else {
                for (int32 i = 0; i < n; i++)
                    dfdy[i*n + i] += 1;

                gear_sgefa(dfdy, n, n, ipivot, &info);
                gear_sgesl(dfdy, n, n, ipivot, k1, 0);
            }
            for (int32 i = 0; i < n; i++)
                ynew[i] = y[i] + .5*h*k1[i];
            rhs_function(t + .5*h, ynew, f1, n);
            for (int32 i = 0; i < n; i++)
                k2[i] = f1[i] - k1[i];
            if (cv_bandflag)
                odesol_bandsol(dfdy, k2, ml, mr, n);
            else
                gear_sgesl(dfdy, n, n, ipivot, k2, 0);
            for (int32 i = 0; i < n; i++) {
                k2[i] = k2[i] + k1[i];
                ynew[i] = y[i] + h*k2[i];
            }
            tnew = t + h;
            rhs_function(tnew, ynew, f2, n);
            for (int32 i = 0; i < n; i++)
                k3[i] = f2[i] - e32*(k2[i] - f1[i]) - 2*(k1[i] - f0[i]) +
                        (h*d)*dfdt[i];
            if (cv_bandflag)
                odesol_bandsol(dfdy, k3, ml, mr, n);
            else
                gear_sgesl(dfdy, n, n, ipivot, k3, 0);
            /*ninf=0;  This is not used anywhere?
             */
            err = 0.0;
            for (int32 i = 0; i < n; i++) {
                temp = MAX(MAX(fabs(y[i]), fabs(ynew[i])), thresh);
                temp = fabs(k1[i] - 2*k2[i] + k3[i]) / temp;
                if (err < temp)
                    err = temp;
            }
            err = err*(absh / 6);
            if (err > rtol) {
                if (absh < hmin) {
                    *ierr = -1;
                    return -1;
                }
                absh = MAX(hmin,
                           absh*MAX(0.1, pow(0.8*(rtol / err), 1. / 3.)));
                h = tdir*absh;
                nofailed = 0;
                done = 0;
            } else {
                break;
            }
        }
        if (nofailed == 1) {
            temp = 1.25*pow(err / rtol, 1. / 3.);
            if (temp > 0.2)
                absh = absh / temp;
            else
                absh = 5*absh;
        }
        t = tnew;
        for (int32 i = 0; i < n; i++) {
            y[i] = ynew[i];
            f0[i] = f2[i];
        }
    }
    *tstart = t;
    htry = h;
    *istart = 0;
    return 0;
}

/* wait_for_key() {
  char bob[256];
  ggets_plintf(" Pause:");
  gets(bob);
}
*/

/* this assumes that yp is already computed */
void
odesol_get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy,
                   int32 neq, double eps, double scal) {
    double yold, del, dsy;
    if (cv_bandflag)
        odesol_get_band_jac(dfdy, y, t, ypnew, yp, neq, eps, scal);
    else {
        for (int32 i = 0; i < neq; i++) {
            del = eps*MAX(eps, fabs(y[i]));
            dsy = scal / del;
            yold = y[i];
            y[i] = y[i] + del;
            rhs_function(t, y, ypnew, neq);
            for (int32 j = 0; j < neq; j++)
                dfdy[j*neq + i] = dsy*(ypnew[j] - yp[j]);
            y[i] = yold;
        }
    }
    return;
}

void
odesol_get_band_jac(double *a, double *y, double t, double *ypnew,
                    double *ypold, int32 n, double eps, double scal) {
    int32 ml = cv_bandlower, mr = cv_bandupper;
    int32 k, n1 = n - 1, mt = ml + mr + 1;
    double yhat;
    double dy;
    double dsy;
    for (int32 i = 0; i < (n*mt); i++)
        a[i] = 0.0;
    for (int32 i = 0; i < n; i++) {
        yhat = y[i];
        dy = eps*(eps + fabs(yhat));
        dsy = scal / dy;
        y[i] += dy;
        rhs_function(t, y, ypnew, n);
        for (int32 j = -ml; j <= mr; j++) {
            k = i - j;
            if (k < 0 || k > n1)
                continue;
            a[k*mt + j + ml] = dsy*(ypnew[k] - ypold[k]);
        }
        y[i] = yhat;
    }
    return;
}

int32
odesol_bandfac(/*   factors the matrix    */
               double *a, int32 ml, int32 mr, int32 n) {
    int32 n1 = n - 1, mt = ml + mr + 1, row, rowi, m, r0, ri0;
    double al;
    for (row = 0; row < n; row++) {
        r0 = row*mt + ml;
        if ((al = a[r0]) == 0.0)
            return -1 - row;
        al = 1.0 / al;
        m = MIN(mr, n1 - row);
        for (int32 j = 1; j <= m; j++)
            a[r0 + j] = a[r0 + j]*al;
        a[r0] = al;
        for (int32 i = 1; i <= ml; i++) {
            rowi = row + i;
            if (rowi > n1)
                break;
            ri0 = rowi*mt + ml;
            al = a[ri0 - i];
            if (al == 0.0)
                continue;
            for (int32 k = 1; k <= m; k++)
                a[ri0 - i + k] = a[ri0 - i + k] - (al*a[r0 + k]);
            a[ri0 - i] = -al;
        }
    }
    return 0;
}

void
odesol_bandsol(/* requires that the matrix be factored   */
               double *a, double *b, int32 ml, int32 mr, int32 n) {
    int32 r0;
    int32 mt = ml + mr + 1;
    int32 m, n1 = n - 1, row;
    for (int32 i = 0; i < n; i++) {
        r0 = i*mt + ml;
        m = MAX(-ml, -i);
        for (int32 j = m; j < 0; j++)
            b[i] += a[r0 + j]*b[i + j];
        b[i] *= a[r0];
    }
    for (row = n1 - 1; row >= 0; row--) {
        m = MIN(mr, n1 - row);
        r0 = row*mt + ml;
        for (int32 k = 1; k <= m; k++)
            b[row] = b[row] - a[r0 + k]*b[row + k];
    }
    return;
}
