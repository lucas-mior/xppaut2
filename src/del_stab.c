#include <stdlib.h>

#include "functions.h"

#include <math.h>
#include "xpplim.h"
#include "integers.h"
#include "xmalloc.h"

/* this code takes the determinant of a floatcomplex valued matrix
 */

static int32 del_stab_test_sign(double old, double new);
static double del_stab_get_arg(double *delay2, double *coef, int32 m, int32 n, COMPLEX lambda);
static void del_stab_process_root(double real, double im);
static COMPLEX del_stab_z_determ(COMPLEX *z, int32 n);
static double del_stab_z_abs(COMPLEX z);
static COMPLEX del_stab_rtoc(double x, double y);
static void del_stab_switch_rows(COMPLEX *z, int32 i1, int32 i2, int32 n);
static COMPLEX del_stab_c_exp2(COMPLEX z);
static COMPLEX del_stab_c_div(COMPLEX z, COMPLEX w);
static COMPLEX del_stab_z_mult(COMPLEX z, COMPLEX w);
static COMPLEX del_stab_z_dif(COMPLEX z, COMPLEX w);
static int32 del_stab_plot_args(double *coef, double *delay2, int32 n, int32 m, int32 npts,
                                double almax, double wmax);
static void del_stab_z_make(COMPLEX *z, double *delay2, int32 n, int32 m, double *coef,
                            COMPLEX lambda);

#define Z(a, b, NN) z[(a) + NN*(b)]

/*typedef struct{
  double r,i;
}COMPLEX;
*/

/* The
 code here replaces the do_sing code if the equation is
   a delay2 differential equation.
*/

void
del_stab_do_delay_sing(double *x, double eps, double err, double big, int32 maxit, int32 n,
                       int32 *ierr, double *stabinfo) {
    double rr[2];

    double colnorm = 0;
    double colmax;
    double colsum;
    double *work, old_x[MAX_ODE], sign;
    double *coef, yp[MAX_ODE], y[MAX_ODE], xp[MAX_ODE], dx;
    int32 kmem = n*(2*n + 5) + 50, i, okroot;

    double *ev;
    ev = xmalloc((usize)(2*n)*sizeof(*ev));
    for (i = 0; i < (2*n); i++) {
        ev[i] = 0.0;
    }
    // first we establish how many delays there are
    del_stab_flag = 0;
    for (i = 0; i < n; i++) {
        old_x[i] = x[i];
    }
    work = xmalloc((usize)kmem*sizeof(*work));
    gear_rooter(x, err, eps, big, work, ierr, maxit, n);
    if (*ierr != 0) {
        del_stab_flag = 1;
        free(work);
        ggets_err_msg("Could not converge to root");
        for (i = 0; i < n; i++) {
            x[i] = old_x[i];
        }
        return;
    }
    // OKAY -- we have the root
    NDelay = 0;
    rhs_function(0.0, x, y, n);  // one more evaluation to get delays
    for (i = 0; i < n; i++) {
        variable_shift[0][i] = x[i];  // unshifted
        variable_shift[1][i] = x[i];
    }
    free(work);
    coef = xmalloc((usize)(n*n*(NDelay + 1))*sizeof(*coef));

    // now we must compute a bunch of jacobians
    // first the normal one
    del_stab_flag = -1;
    which_delay = -1;
    colmax = 0.0;

    for (i = 0; i < n; i++) {
        colsum = 0.0;
        for (int32 j = 0; j < n; j++) {
            xp[j] = x[j];
        }
        dx = eps*gear_amax(eps, fabs(x[i]));
        xp[i] = xp[i] + dx;
        rhs_function(0.0, xp, yp, n);
        for (int32 j = 0; j < n; j++) {
            coef[j*n + i] = (yp[j] - y[j]) / dx;
            colsum += fabs(coef[j*n + i]);
        }
        if (colsum > colmax) {
            colmax = colsum;
        }
    }
    colnorm = colmax;
    for (int32 j = 0; j < n; j++) {
        xp[j] = x[j];
    }
    // now the jacobians for the delays
    for (int32 k = 0; k < NDelay; k++) {
        which_delay = k;
        colmax = 0.0;
        for (i = 0; i < n; i++) {
            colsum = 0.0;
            for (int32 j = 0; j < n; j++) {
                variable_shift[1][j] = variable_shift[0][j];
            }
            dx = eps*gear_amax(eps, fabs(x[i]));
            variable_shift[1][i] = x[i] + dx;
            rhs_function(0.0, x, yp, n);
            variable_shift[1][i] = x[i];
            for (int32 j = 0; j < n; j++) {
                coef[j*n + i + n*n*(k + 1)] = (yp[j] - y[j]) / dx;
                colsum += fabs(coef[j*n + i + n*n*(k + 1)]);
            }
            if (colsum > colmax) {
                colmax = colsum;
            }
        }
        colnorm += colmax;
    }
    sign = del_stab_plot_args(coef, delay_list, n, NDelay, delay_grid, colnorm, colnorm);

    okroot = del_stab_find_positive_root(coef, delay_list, n, NDelay, err, eps, big, maxit, rr);
    if (okroot > 0) {
        ev[0] = rr[0];
        ev[1] = rr[1];
    }
    free(coef);
    *stabinfo = (double)fabs(sign);
    // if(*stabinfo>0)
    i = (int32)sign;
    if (i == 0 && okroot == 1 && alpha_max > 0) {
        i = 2;
    }

    eig_list_create_eq_box(ABS(i), 2, 0, 0, 0, x, n);
    del_stab_flag = 1;
    free(ev);
    if (okroot == 1) {
        *stabinfo = alpha_max;
    }
    return;
}

COMPLEX
del_stab_z_dif(COMPLEX z, COMPLEX w) {
    COMPLEX sum;
    sum.r = z.r - w.r;
    sum.i = z.i - w.i;
    return sum;
}

COMPLEX
del_stab_z_mult(COMPLEX z, COMPLEX w) {
    COMPLEX sum;
    sum.r = z.r*w.r - z.i*w.i;
    sum.i = z.r*w.i + z.i*w.r;
    return sum;
}

COMPLEX
del_stab_c_div(COMPLEX z, COMPLEX w) {
    COMPLEX sum;
    double amp = w.r*w.r + w.i*w.i;
    sum.r = (z.r*w.r + z.i*w.i) / amp;
    sum.i = (z.i*w.r - z.r*w.i) / amp;
    return sum;
}

COMPLEX
del_stab_c_exp2(COMPLEX z) {
    COMPLEX sum;
    double ex = exp(z.r);
    sum.r = ex*cos(z.i);
    sum.i = ex*sin(z.i);
    return sum;
}

void
del_stab_switch_rows(COMPLEX *z, int32 i1, int32 i2, int32 n) {
    COMPLEX zt;
    for (int32 j = 0; j < n; j++) {
        zt = Z(i1, j, n);
        Z(i1, j, n) = Z(i2, j, n);
        Z(i2, j, n) = zt;
    }
    return;
}

COMPLEX
del_stab_rtoc(double x, double y) {
    COMPLEX sum;
    sum.i = y;
    sum.r = x;
    return sum;
}

double
del_stab_z_abs(COMPLEX z) {
    return sqrt(z.i*z.i + z.r*z.r);
}

COMPLEX
del_stab_z_determ(COMPLEX *z, int32 n) {
    int32 imax = 0;
    double q;
    double qmax;
    COMPLEX sign = del_stab_rtoc(1.0, 0.0), mult, sum, zd;
    for (int32 j = 0; j < n; j++) {
        qmax = 0.0;
        for (int32 i = j; i < n; i++) {
            q = del_stab_z_abs(Z(i, j, n));
            if (q > qmax) {
                qmax = q;
                imax = i;
            }
        }
        if (qmax == 0.0) {
            return del_stab_rtoc(0.0, 0.0);
        }
        del_stab_switch_rows(z, imax, j, n);
        if (imax > j) {
            sign = del_stab_z_mult(del_stab_rtoc(-1.0, 0.0), sign);
        }
        zd = Z(j, j, n);
        for (int32 i = j + 1; i < n; i++) {
            mult = del_stab_c_div(Z(i, j, n), zd);
            for (int32 k = j + 1; k < n; k++) {
                Z(i, k, n) = del_stab_z_dif(Z(i, k, n), del_stab_z_mult(mult, Z(j, k, n)));
            }
        }
    }
    sum = sign;
    for (int32 j = 0; j < n; j++) {
        sum = del_stab_z_mult(sum, Z(j, j, n));
    }
    return sum;
}

void
del_stab_z_make(COMPLEX *z, double *delay2, int32 n, int32 m, double *coef, COMPLEX lambda) {
    int32 km;
    COMPLEX temp;
    COMPLEX eld;

    for (int32 j = 0; j < n; j++) {
        for (int32 i = 0; i < n; i++) {
            if (i == j) {
                temp = lambda;
            } else {
                temp = del_stab_rtoc(0.0, 0.0);
            }
            z[i + j*n] =
                del_stab_z_dif(temp, del_stab_rtoc(coef[i + j*n], 0.0));  // initialize the array
        }
    }
    for (int32 k = 0; k < m; k++) {
        km = (k + 1)*n*n;
        temp = del_stab_rtoc(-delay2[k],
                             0.0);  // convert delay2 to floatcomplex number
        eld = del_stab_c_exp2(del_stab_z_mult(temp, lambda));  // compute exp(-lambda*tau)
        for (int32 j = 0; j < n; j++) {
            for (int32 i = 0; i < n; i++) {
                z[i + j*n] = del_stab_z_dif(
                    z[i + j*n], del_stab_z_mult(eld, del_stab_rtoc(coef[km + i + n*j], 0.0)));
            }
        }
    }
    return;
}

int32
del_stab_find_positive_root(double *coef, double *delay2, int32 n, int32 m, double err, double eps,
                            double big, int32 maxit, double *rr) {
    COMPLEX lambda;
    COMPLEX lambdap;
    COMPLEX det, *z, detp;
    double jac[4];
    double xl;
    double yl;
    double r;
    double xlp;
    double ylp;

    lambda.r = alpha_max;
    lambda.i = omega_max;

    z = xmalloc(sizeof(*z)*(usize)(n*n));

    // now Newtons Method for maxit times
    for (int32 k = 0; k < maxit; k++) {
        del_stab_z_make(z, delay2, n, m, coef, lambda);
        det = del_stab_z_determ(z, n);

        r = del_stab_z_abs(det);
        if (r < err) {  // within the tolerance
            del_stab_process_root(lambda.r, lambda.i);
            alpha_max = lambda.r;
            omega_max = lambda.i;
            return 1;
        }
        xl = lambda.r;
        yl = lambda.i;

        // compute the Jacobian
        if (fabs(xl) > eps) {
            r = eps*fabs(xl);
        } else {
            r = eps*eps;
        }
        xlp = xl + r;
        lambdap = del_stab_rtoc(xlp, yl);
        del_stab_z_make(z, delay2, n, m, coef, lambdap);
        detp = del_stab_z_determ(z, n);
        jac[0] = (detp.r - det.r) / r;
        jac[2] = (detp.i - det.i) / r;
        if (fabs(yl) > eps) {
            r = eps*fabs(yl);
        } else {
            r = eps*eps;
        }
        ylp = yl + r;
        lambdap = del_stab_rtoc(xl, ylp);
        del_stab_z_make(z, delay2, n, m, coef, lambdap);
        detp = del_stab_z_determ(z, n);
        jac[1] = (detp.r - det.r) / r;
        jac[3] = (detp.i - det.i) / r;
        r = jac[0]*jac[3] - jac[1]*jac[2];
        if (r == 0) {
            ggets_plintf(" singular jacobian \n");
            return -1;
        }
        xlp = (jac[3]*det.r - jac[1]*det.i) / r;
        ylp = (-jac[2]*det.r + jac[0]*det.i) / r;
        xl = xl - xlp;
        yl = yl - ylp;
        r = fabs(xlp) + fabs(ylp);
        lambda.r = xl;
        lambda.i = yl;
        if (r < err) {  // within the tolerance
            del_stab_process_root(lambda.r, lambda.i);
            alpha_max = lambda.r;
            omega_max = lambda.i;
            rr[0] = alpha_max;
            rr[1] = omega_max;
            return 1;
        }
        if (r > big) {
            ggets_plintf("Failed to converge \n");
            return -1;
        }
    }

    ggets_plintf("Max iterates exceeded \n");
    return -1;
}

void
del_stab_process_root(double real, double im) {
    ggets_plintf("Root: %g + I %g \n", real, im);
    return;
}

double
del_stab_get_arg(double *delay2, double *coef, int32 m, int32 n, COMPLEX lambda) {
    int32 km;
    COMPLEX *z;
    COMPLEX temp;
    COMPLEX eld;
    double arg;

    if (m == 0) {
        return 0;  // no delays so don't use this!
    }
    z = xmalloc(sizeof(*z)*(usize)(n*n));
    for (int32 j = 0; j < n; j++) {
        for (int32 i = 0; i < n; i++) {
            if (i == j) {
                temp = lambda;
            } else {
                temp = del_stab_rtoc(0.0, 0.0);
            }
            z[i + j*n] =
                del_stab_z_dif(temp, del_stab_rtoc(coef[i + j*n], 0.0));  // initialize the array
        }
    }
    for (int32 k = 0; k < m; k++) {
        km = (k + 1)*n*n;
        temp = del_stab_rtoc(-delay2[k],
                             0.0);  // convert delay2 to floatcomplex number
        eld = del_stab_c_exp2(del_stab_z_mult(temp, lambda));  // compute exp(-lambda*tau)
        for (int32 j = 0; j < n; j++) {
            for (int32 i = 0; i < n; i++) {
                z[i + j*n] = del_stab_z_dif(
                    z[i + j*n], del_stab_z_mult(eld, del_stab_rtoc(coef[km + i + n*j], 0.0)));
            }
        }
    }
    //  the array is done
    temp = del_stab_z_determ(z, n);
    free(z);
    arg = atan2(temp.i, temp.r);
    return arg;
}

int32
del_stab_test_sign(double old, double new) {
    if (old > 0.0 && new < 0.0) {
        if (old > 2.9 && new < -2.9) {
            return 1;
        }
        return 0;  // doesnt pass threshold
    }
    if (old < 0.0 && new > 0.0) {
        if (old < -2.9 && new > 2.9) {
            return -1;
        }
        return 0;
    }
    return 0;
}

/* code for establishing delay2 stability
   sign=del_stab_plot_args(coef,delay2,n,m,npts,amax,wmax)
    coef is a real array of length  (m+1)*n^2
    each n^2 block is the jacobian with respect to the mth delay2
    m total delays
    n is size of system
    npts is number of pts on each part of contour
    contour is
      i wmax -----<---------    amax+i wmax
       |                            |
       V                            ^
       |                            |
     -i wmax ----->-----------  amax-i wmax

     sign is the number of roots in the contour using the argument
     principle
*/

int32
del_stab_plot_args(double *coef, double *delay2, int32 n, int32 m, int32 npts, double almax,
                   double wmax) {
    int32 sign = 0;
    COMPLEX lambda;
    double x;
    double y;
    double arg;
    double oldarg = 0.0;
    double ds;  // steplength
                // first the contour from i wmax -- -i wmax
    ds = 2*wmax / npts;
    x = 0.0;
    for (int32 i = 0; i < npts; i++) {
        y = wmax - i*ds;
        lambda = del_stab_rtoc(x, y);
        arg = del_stab_get_arg(delay2, coef, m, n, lambda);
        sign = sign + del_stab_test_sign(oldarg, arg);
        oldarg = arg;
    }
    // lower contour
    y = -wmax;
    ds = almax / npts;
    for (int32 i = 0; i < npts; i++) {
        x = i*ds;
        lambda = del_stab_rtoc(x, y);
        arg = del_stab_get_arg(delay2, coef, m, n, lambda);
        sign = sign + del_stab_test_sign(oldarg, arg);
        oldarg = arg;
    }
    // right contour
    x = almax;
    ds = 2*wmax / npts;
    for (int32 i = 0; i < npts; i++) {
        y = -wmax + i*ds;
        lambda = del_stab_rtoc(x, y);
        arg = del_stab_get_arg(delay2, coef, m, n, lambda);
        sign = sign + del_stab_test_sign(oldarg, arg);
        oldarg = arg;
    }

    // top contour
    y = wmax;
    ds = almax / npts;
    for (int32 i = 0; i < npts; i++) {
        x = almax - i*ds;
        lambda = del_stab_rtoc(x, y);
        arg = del_stab_get_arg(delay2, coef, m, n, lambda);
        sign = sign + del_stab_test_sign(oldarg, arg);
        oldarg = arg;
    }
    return sign;
}
