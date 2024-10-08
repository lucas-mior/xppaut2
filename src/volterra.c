#include "functions.h"
#include "integers.h"
#include "xmalloc.h"
#include <stdbool.h>

#include <stdlib.h>
#include "xpplim.h"
#include <math.h>
#include <stdio.h>
#include "parserslow.h"

/* This is an implicit solver for volterra integral and integro-differential
 * equations.  It is based on code found in Peter Linz's book
 * ont Volterra equations.
 * One tries to evaluate:

 *    int_0^t ( (t-t')^-mu K(t,t',u) dt')
 * where  0 <= mu < 1 and K(t,t',u) is cts and Lipschitz.
 * The product method is used combined with the trapezoidal rule for
 * integration. The method is A-stable since it is an implicit scheme.

 * The kernel structure contains the constant mu and the expression for
 * evaluating K(t,t',u) */

#define CONV 2
static double *Memory[MAX_KER];
static int32 current_point;
static int32 kn_flag;

int32 auto_evaluate = 0;

static void volterra_get_kn(double *y, double t);
static double volterra_betnn(double mu, double dt);
static double volterra_alpbetjn(double mu, double dt, int32 l);
static double volterra_alpha1n(double mu, double dt, double t, double t0);

double
volterra_ker_val(int32 in) {
    if (kn_flag) {
        return kernel[in].k_n;
    }
    return kernel[in].k_n1;
}

void
volterra_alloc_memory(void) {
    int32 len;
    int32 formula[256];

    // First parse the kernels   since these were deferred
    for (int32 i = 0; i < nkernel; i++) {
        kernel[i].k_n = 0.0;
        if (parserslow_add_expr(kernel[i].expr, formula, &len)) {
            ggets_plintf("Illegal kernel %s=%s\n", kernel[i].name, kernel[i].expr);
            exit(0);  // fatal error ...
        }
        kernel[i].formula = xmalloc((usize)(len + 2)*sizeof(*(kernel[i].formula)));
        for (int32 j = 0; j < len; j++) {
            kernel[i].formula[j] = formula[j];
        }
        if (kernel[i].flag == CONV) {
            if (parserslow_add_expr(kernel[i].kerexpr, formula, &len)) {
                ggets_plintf("Illegal convolution %s=%s\n", kernel[i].name, kernel[i].kerexpr);
                exit(0);  // fatal error ...
            }
            kernel[i].kerform = xmalloc((usize)(len + 2)*sizeof(*(kernel[i].kerform)));
            for (int32 j = 0; j < len; j++) {
                kernel[i].kerform[j] = formula[j];
            }
        }
    }
    volterra_allocate(max_points, 0);
    return;
}

void
volterra_allocate(int32 npts, int32 flag) {
    int32 i;
    int32 oldmem = max_points;
    int32 ntot = NODE + fix_var + nmarkov;
    npts = ABS(npts);
    max_points = npts;
    // now allocate the memory
    if (nkernel == 0) {
        return;
    }
    if (flag == 1) {
        for (i = 0; i < ntot; i++) {
            free(Memory[i]);
        }
    }
    for (i = 0; i < ntot; i++) {
        Memory[i] = xmalloc(sizeof(*Memory)*(usize)max_points);
        if (Memory[i] == NULL) {
            break;
        }
    }

    if (i < ntot && flag == 0) {
        ggets_plintf("Not enough memory... make Maxpts smaller \n");
        exit(0);
    }
    if (i < ntot) {
        max_points = oldmem;
        for (int32 j = 0; j < i; j++) {
            free(Memory[j]);
        }
        for (i = 0; i < ntot; i++) {
            Memory[i] = xmalloc(sizeof(*Memory)*(usize)max_points);
        }
        ggets_err_msg("Not enough memory...resetting");
    }
    current_point = 0;
    kn_flag = 1;
    volterra_alloc_kernels(flag);
    return;
}

void
volterra_re_evaluate_kernels(void) {
    if (auto_evaluate == 0) {
        return;
    }
    for (int32 i = 0; i < nkernel; i++) {
        if (kernel[i].flag == CONV) {
            for (int32 j = 0; j <= max_points; j++) {
                SETVAR(0, T0 + delta_t*j);
                kernel[i].cnv[j] = evaluate(kernel[i].kerform);
            }
        }
    }
    return;
}

void
volterra_alloc_kernels(int32 flag) {
    int32 n = max_points;
    double mu;

    for (int32 i = 0; i < nkernel; i++) {
        if (kernel[i].flag == CONV) {
            if (flag == 1) {
                free(kernel[i].cnv);
            }
            kernel[i].cnv = xmalloc((usize)(n + 1)*sizeof(*(kernel[i].cnv)));
            for (int32 j = 0; j <= n; j++) {
                SETVAR(0, T0 + delta_t*j);
                kernel[i].cnv[j] = evaluate(kernel[i].kerform);
            }
        }
        // Do the alpha functions here later
        if (kernel[i].mu > 0.0) {
            mu = kernel[i].mu;
            if (flag == 1) {
                free(kernel[i].al);
            }
            kernel[i].al = xmalloc((usize)(n + 1)*sizeof(*(kernel[i].al)));
            for (int32 j = 0; j <= n; j++) {
                kernel[i].al[j] = volterra_alpbetjn(mu, delta_t, j);
            }
        }
    }
    return;
}

/* the following is the main driver for evaluating the sums in the
 * kernel the results here are used in the implicit solver.  The integral
 * up to t_n-1 is evaluated and placed in sum.  Kn-> Kn-1

 * the weights al and bet are computed in general, but specifically
 * for mu=0,.5 since these involve no transcendental functions

   */

/***   FIX THIS TO DO MORE GENERAL STUFF
       K(t,t',u,u') someday...
***/

void
volterra_init_sums(double t0, int32 n, double dt, int32 i0, int32 iend, int32 ishift) {
    double t = t0 + n*dt, tp = t0 + i0*dt;
    double sum[MAX_ODE];
    double al;
    double alpbet;
    double mu;
    int32 nvar = fix_var + NODE + nmarkov;
    int32 l;
    int32 ioff;
    int32 ker;
    SETVAR(0, t);
    SETVAR(prime_start, tp);
    for (l = 0; l < nvar; l++) {
        SETVAR(l + 1, Memory[l][ishift]);
    }
    for (ker = 0; ker < nkernel; ker++) {
        kernel[ker].k_n1 = kernel[ker].k_n;
        mu = kernel[ker].mu;
        if (mu == 0.0) {
            al = .5*dt;
        } else {
            al = volterra_alpha1n(mu, dt, t, tp);
        }
        sum[ker] = al*evaluate(kernel[ker].formula);
        if (kernel[ker].flag == CONV) {
            sum[ker] = sum[ker]*kernel[ker].cnv[n - i0];
        }
    }
    for (int32 i = 1; i <= iend; i++) {
        ioff = (ishift + i) % max_points;
        tp += dt;
        SETVAR(prime_start, tp);
        for (l = 0; l < nvar; l++) {
            SETVAR(l + 1, Memory[l][ioff]);
        }
        for (ker = 0; ker < nkernel; ker++) {
            mu = kernel[ker].mu;
            if (mu == 0.0) {
                alpbet = dt;
            } else {
                alpbet = kernel[ker].al[n - i0 - i];
            }
            if (kernel[ker].flag == CONV) {
                sum[ker] += (alpbet*evaluate(kernel[ker].formula)*kernel[ker].cnv[n - i0 - i]);
            } else {
                sum[ker] += (alpbet*evaluate(kernel[ker].formula));
            }
        }
    }
    for (ker = 0; ker < nkernel; ker++) {
        kernel[ker].sum = sum[ker];
    }
    return;
}

/* the following functions compute integrals for the piecewise
 * -- constant -- product integration rule.  Thus they agree with
 * the trapezoid rule for mu=0 and there is a special case for mu=.5
 * since that involves no transcendentals.  Later I will put in the
 * piecewise --linear-- method
 */

double
volterra_alpha1n(double mu, double dt, double t, double t0) {
    double m1;

    if (mu == .5) {
        return sqrt(fabs(t - t0)) - sqrt(fabs(t - t0 - dt));
    }
    m1 = 1 - mu;
    return .5*(pow(fabs(t - t0), m1) - pow(fabs(t - t0 - dt), m1)) / m1;
}

double
volterra_alpbetjn(double mu, double dt, int32 l) {
    double m1;
    double dif = l*dt;
    if (mu == .5) {
        return sqrt(dif + dt) - sqrt(fabs(dif - dt));
    }
    m1 = 1 - mu;
    return .5*(pow(dif + dt, m1) - pow(fabs(dif - dt), m1)) / m1;
}

double
volterra_betnn(double mu, double dt) {
    double m1;
    if (mu == .5) {
        return sqrt(dt);
    }
    m1 = 1 - mu;
    return .5*pow(dt, m1) / m1;
}

void
volterra_get_kn(double *y, double t) {
    // uses the guessed value y to update Kn

    SETVAR(0, t);
    SETVAR(prime_start, t);
    for (int32 i = 0; i < NODE; i++) {
        SETVAR(i + 1, y[i]);
    }
    for (int32 i = NODE; i < NODE + fix_var; i++) {
        SETVAR(i + 1, evaluate(my_ode[i]));
    }
    for (int32 i = 0; i < nkernel; i++) {
        if (kernel[i].flag == CONV) {
            kernel[i].k_n =
                kernel[i].sum + kernel[i].betnn*evaluate(kernel[i].formula)*kernel[i].cnv[0];
        } else {
            kernel[i].k_n = kernel[i].sum + kernel[i].betnn*evaluate(kernel[i].formula);
        }
    }
    return;
}

int32
volterra(double *y, double *t, double dt, int32 nt, int32 neq, int32 *istart, double *work) {
    double *jac, *yg, *yp, *yp2, *ytemp, *errvec;
    double z;
    double mu;
    double bet;
    int32 j;
    yp = work;
    yg = yp + neq;
    ytemp = yg + neq;
    errvec = ytemp + neq;
    yp2 = errvec + neq;
    jac = yp2 + neq;

    //  Initialization of everything
    if (*istart == 1) {
        current_point = 0;
        kn_flag = 1;
        for (int32 i = 0; i < nkernel; i++) {  // zero the integrals
            kernel[i].k_n = 0.0;
            kernel[i].k_n1 = 0.0;
            mu = kernel[i].mu;  //  compute bet_nn
            if (mu == 0.0) {
                bet = .5*dt;
            } else {
                bet = volterra_betnn(mu, dt);
            }
            kernel[i].betnn = bet;
        }
        SETVAR(0, *t);
        SETVAR(prime_start, *t);
        for (int32 i = 0; i < NODE; i++) {
            if (!eq_type[i]) {
                SETVAR(i + 1, y[i]);  // assign initial data
            }
        }
        for (int32 i = NODE; i < NODE + fix_var; i++) {
            SETVAR(i + 1,
                   evaluate(my_ode[i]));  // set fixed variables  for pass 1
        }
        for (int32 i = 0; i < NODE; i++) {
            if (eq_type[i]) {
                z = evaluate(my_ode[i]);  // reset IC for integral eqns
                SETVAR(i + 1, z);
                y[i] = z;
            }
        }
        for (int32 i = NODE; i < NODE + fix_var; i++) {  // pass 2 for fixed variables
            SETVAR(i + 1, evaluate(my_ode[i]));
        }
        for (int32 i = 0; i < NODE + fix_var + nmarkov; i++) {
            Memory[i][0] = get_ivar(i + 1);  // save everything
        }
        current_point = 1;
        *istart = 0;
    }

    for (int32 i = 0; i < nt; i++)  // the real computation
    {
        *t = *t + dt;
        markov_set_wieners(dt, y, *t);
        if ((j = volterra_step(y, *t, dt, neq, yg, yp, yp2, errvec, jac)) != 0) {
            return j;
        }
        delay_handle_stor_delay(y);
    }
    return 0;
}

int32
volterra_step(double *y, double t, double dt, int32 neq, double *yg, double *yp, double *yp2,
              double *errvec, double *jac) {
    int32 i0;
    int32 iend;
    int32 ishift;
    int32 iter = 0;
    int32 info;
    int32 ipivot[MAX_ODE1];
    int32 ind;
    int32 n1 = NODE + 1;
    double dt2 = .5*dt;
    double err;
    double del;
    double yold;
    double fac;
    double delinv;
    i0 = MAX(0, current_point - max_points);
    iend = MIN(current_point - 1, max_points - 1);
    ishift = i0 % max_points;
    volterra_init_sums(T0, current_point, dt, i0, iend,
                       ishift);  //  initialize all the sums
    kn_flag = 0;
    for (int32 i = 0; i < neq; i++) {
        SETVAR(i + 1, y[i]);
        yg[i] = y[i];
    }
    for (int32 i = NODE; i < NODE + nmarkov; i++) {
        SETVAR(i + 1 + fix_var, y[i]);
    }
    SETVAR(0, t - dt);
    for (int32 i = NODE; i < NODE + fix_var; i++) {
        SETVAR(i + 1, evaluate(my_ode[i]));
    }
    for (int32 i = 0; i < NODE; i++) {
        if (!eq_type[i]) {
            yp2[i] = y[i] + dt2*evaluate(my_ode[i]);
        } else {
            yp2[i] = 0.0;
        }
    }
    kn_flag = 1;
    while (true) {
        volterra_get_kn(yg, t);
        for (int32 i = NODE; i < NODE + fix_var; i++) {
            SETVAR(i + 1, evaluate(my_ode[i]));
        }
        for (int32 i = 0; i < NODE; i++) {
            yp[i] = evaluate(my_ode[i]);
            if (eq_type[i]) {
                errvec[i] = -yg[i] + yp[i];
            } else {
                errvec[i] = -yg[i] + dt2*yp[i] + yp2[i];
            }
        }
        //   Compute Jacobian
        for (int32 i = 0; i < NODE; i++) {
            del = newt_err*MAX(newt_err, fabs(yg[i]));
            yold = yg[i];
            yg[i] += del;
            delinv = 1. / del;
            volterra_get_kn(yg, t);
            for (int32 j = NODE; j < NODE + fix_var; j++) {
                SETVAR(j + 1, evaluate(my_ode[j]));
            }
            for (int32 j = 0; j < NODE; j++) {
                fac = delinv;
                if (!eq_type[j]) {
                    fac *= dt2;
                }
                jac[j*NODE + i] = (evaluate(my_ode[j]) - yp[j])*fac;
            }
            yg[i] = yold;
        }

        for (int32 i = 0; i < NODE; i++) {
            jac[n1*i] -= 1.0;
        }
        gear_sgefa(jac, NODE, NODE, ipivot, &info);
        if (info != -1) {
            return -1;  // Jacobian is singular
        }
        err = 0.0;
        gear_sgesl(jac, NODE, NODE, ipivot, errvec, 0);
        for (int32 i = 0; i < NODE; i++) {
            err = MAX(fabs(errvec[i]), err);
            yg[i] -= errvec[i];
        }
        if (err < euler_tol) {
            break;
        }
        iter++;
        if (iter > euler_max_iter) {
            return -2;  // too many iterates
        }
    }
    // We have a good graphics_point; lets save it
    volterra_get_kn(yg, t);
    for (int32 i = 0; i < NODE; i++) {
        y[i] = yg[i];
    }
    ind = current_point % max_points;
    for (int32 i = 0; i < NODE + fix_var + nmarkov; i++) {
        Memory[i][ind] = GETVAR(i + 1);
    }
    current_point++;

    return 0;
}
