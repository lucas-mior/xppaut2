#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "functions.h"
#include "xmalloc.h"
#include "parserslow.h"
#include "integers.h"
#include "xpplim.h"

/* this has a bunch of numerical routines
 * averaging
 * adjoints
 * transpose
 * maximal liapunov exponent */

#define READEM 1

static double **my_adj;
static int32 adj_len;
static double **my_h;
static double *my_liap[2];

static struct MyTrans {
    int32 here;
    int32 col0;
    int32 ncol;
    int32 colskip;
    int32 row0;
    int32 nrow;
    int32 rowskip;
    double **data;
    char firstcol[MAX_ODE_NAME_LENGTH];
} my_trans;

static int32 LIAP_FLAG = 0;
static int32 LIAP_N;
static int32 LIAP_I;
static int32 LIAP_N;
static int32 LIAP_I;
static double adj_eps = 1.e-8;
static double adj_err = 1.e-3;
static int32 adj_maxit = 20;
static int32 adj_here = 0;
static int32 H_HERE = 0;
static int32 h_len;
static int32 hodd_ev = 0;
bool adj_range = false;
static int32 *coup_fun[MAX_ODE];
static char *coup_string[MAX_ODE];

static void adjoints_h_back(void);
static void adjoints_back(void);
static int32 adjoints_make_h(double **orb, double **adj, int32 nt, int32 node, int32 silent2);
static int32 adjoints_step_eul(double **jac, int32 k, int32 k2, double *yold, double *work,
                               int32 node, double dt);
static int32 adjoints_hrw_liapunov(double *liap, int32 batch, double eps);
static int32 adjoints_adjoint(double **orbit, double **adjnt, int32 nt, double dt, double eps,
                              double minerr, int32 maxit, int32 node);

void
adjoints_init_trans(void) {
    my_trans.here = 0;
    strncpy(my_trans.firstcol, uvar_names[0], sizeof(my_trans.firstcol));
    my_trans.ncol = 2;
    my_trans.nrow = 1;
    my_trans.rowskip = 1;
    my_trans.colskip = 1;
    my_trans.row0 = 1;
    my_trans.col0 = 2;
    return;
}

void
adjoints_dump_transpose_info(FILE *fp, int32 f) {
    char bob[256];

    if (f == READEM) {
        fgets(bob, 255, fp);
    } else {
        fprintf(fp, "# Transpose variables etc\n");
    }

    lunch_io_string(my_trans.firstcol, 11, fp, f);
    lunch_io_int(&my_trans.ncol, fp, f, "n columns");
    lunch_io_int(&my_trans.nrow, fp, f, "n rows");
    lunch_io_int(&my_trans.rowskip, fp, f, "row skip");
    lunch_io_int(&my_trans.colskip, fp, f, "col skip");
    lunch_io_int(&my_trans.row0, fp, f, "row 0");
    return;
}

int32
adjoints_do_transpose(void) {
    int32 ii;
    int32 status;
    static char *strings[] = {"*0Column 1", "NCols", "ColSkip", "Row 1", "NRows", "RowSkip"};
    char values[LENGTH(strings)][MAX_LEN_SBOX];
    int32 inrow;
    int32 incol;

    snprintf(values[0], sizeof(values[0]), "%s", my_trans.firstcol);
    snprintf(values[1], sizeof(values[1]), "%d", my_trans.ncol);
    snprintf(values[2], sizeof(values[2]), "%d", my_trans.colskip);
    snprintf(values[3], sizeof(values[3]), "%d", my_trans.row0);
    snprintf(values[4], sizeof(values[4]), "%d", my_trans.nrow);
    snprintf(values[5], sizeof(values[5]), "%d", my_trans.rowskip);

    if (my_trans.here) {
        for (int32 i = 0; i <= my_trans.nrow; i++) {
            free(my_trans.data[i]);
        }
        free(my_trans.data);
        my_trans.here = 0;
        adjoints_data_back();
    }

    status = pop_list_do_string_box(LENGTH(strings), LENGTH(strings), 1, "Transpose Data", strings,
                                    values, 33);
    if (status == 0) {
        return 0;
    }

    browser_find_variable(values[0], &ii);
    if (ii > -1) {
        my_trans.col0 = ii + 1;
    } else {
        ggets_err_msg("No such columns");
        return 0;
    }
    strncpy(my_trans.firstcol, values[0], sizeof(my_trans.firstcol));
    ii = atoi(values[4]);
    if (ii >= NEQ) {
        ii = NEQ - 1;
    }
    my_trans.nrow = ii;
    my_trans.ncol = atoi(values[1]);
    my_trans.colskip = atoi(values[2]);
    my_trans.row0 = atoi(values[3]);
    my_trans.rowskip = atoi(values[5]);

    my_trans.data = xmalloc(sizeof(*(my_trans.data))*(usize)(NEQ + 1));
    for (int32 i = 0; i <= my_trans.nrow; i++) {
        my_trans.data[i] = xmalloc(sizeof(my_trans.data[i])*(usize)my_trans.ncol);
    }
    for (int32 i = my_trans.nrow + 1; i <= NEQ; i++) {
        my_trans.data[i] = storage[i];
    }
    for (int32 j = 0; j < my_trans.ncol; j++) {
        my_trans.data[0][j] = (double)(j + 1);
    }

    for (int32 i = 0; i < my_trans.ncol; i++) {
        incol = my_trans.col0 - 1 + i*my_trans.colskip;
        if (incol > NEQ) {
            incol = NEQ;
        }
        for (int32 j = 0; j < my_trans.nrow; j++) {
            inrow = my_trans.row0 + j*my_trans.rowskip;
            if (inrow > storind) {
                inrow = storind;
            }
            my_trans.data[j + 1][i] = storage[incol][inrow];
        }
    }

    browser_set_data(my_trans.data, 1);
    browser_refresh(my_trans.ncol);
    my_trans.here = 1;
    return 1;
}

void
adjoints_alloc_h_stuff(void) {
    for (int32 i = 0; i < NODE; i++) {
        coup_fun[i] = xmalloc(100*sizeof(*coup_fun));
        coup_string[i] = xmalloc(80);
        memcpy(coup_string[i], "0", strlen("0"));
    }
    return;
}

void
adjoints_data_back(void) {
    four_here = 0;
    browser_set_data(storage, 1);
    browser_refresh(storind);
    return;
}

void
adjoints_back(void) {
    if (adj_here) {
        browser_set_data(my_adj, 1);
        browser_refresh(adj_len);
    }
    return;
}

void
adjoints_h_back(void) {
    if (H_HERE) {
        browser_set_data(my_h, 1);
        browser_refresh(h_len);
    }
    return;
}

/* Here is how to do the range over adjoints and h functions
 * unfortunately, h functions are always computed even if you dont want them
 * they will just be zeros

 * Step 1. Compute a singel orbit, adjoint, and H function (to load the
 * program with the correct right-hand sides for H function. Or just load in
 * set file where it was done
 * Step 2. Set transient to some reasonable number to assure convergence
 * onto the limit cycle as you change parameters and total to be at least
 * 2 periods beyond the transient
 * Step 3. Se up Poincare map - period - stop on section. This lets you
 * get the period
 * Step 4. In numerics averaging - click on adjrange
 * Step 5. Initconds range over the parameter. It should find the periodic
 * orbit, adjoint, and H function and save. Files are of the form
 * orbit.parname_parvalue.dat etc
 */
void
adjoints_make_adj_com(int32 com) {
    static char key[] = "nmaohpr";
    switch (key[com]) {
    case 'n':
        adjoints_new_adjoint();
        break;
    case 'm':
        adjoints_new_h_fun(0);
        break;
    case 'a':
        adjoints_back();
        break;
    case 'o':
        adjoints_data_back();
        break;
    case 'h':
        adjoints_h_back();
        break;
    case 'p':
        // adj2 adjoint parameters
        ggets_new_int("Maximum iterates :", &adj_maxit);
        ggets_new_float("Adjoint error tolerance :", &adj_err);
        break;
    case 'r':
        adj_range = true;
        break;
    default:
        break;
    }
    return;
}

void
adjoints_new_h_fun(int32 silent2) {
    int32 n = 2;

    if (!adj_here) {
        ggets_err_msg("Must compute adjoint first!");
        return;
    }
    if (storind != adj_len) {
        ggets_err_msg("incompatible data and adjoint");
        return;
    }
    if (H_HERE) {
        free(my_h[0]);
        free(my_h[1]);
        if (hodd_ev) {
            free(my_h[2]);
            free(my_h[3]);
        }
        free(my_h);
        H_HERE = 0;
        hodd_ev = 0;
    }
    if (NEQ > 2) {
        hodd_ev = 1;
        n = 4;
    }
    h_len = storind;
    adjoints_data_back();
    my_h = xmalloc(sizeof(*my_h)*(usize)(NEQ + 1));
    for (int32 i = 0; i < n; i++) {
        my_h[i] = xmalloc(sizeof(*my_h)*(usize)h_len);
    }
    for (int32 i = n; i <= NEQ; i++) {
        my_h[i] = storage[i];
    }
    if (adjoints_make_h(storage, my_adj, h_len, NODE, silent2)) {
        H_HERE = 1;
        adjoints_h_back();
    }
    ggets_ping();
    return;
}

void
adjoints_dump_h_stuff(FILE *fp, int32 f) {
    char bob[256];
    if (f == READEM) {
        fgets(bob, 255, fp);
    } else {
        fprintf(fp, "# Coupling stuff for H funs\n");
    }
    for (int32 i = 0; i < NODE; i++) {
        lunch_io_string(coup_string[i], 79, fp, f);
    }
    return;
}

int32
adjoints_make_h(double **orb, double **adj, int32 nt, int32 node, int32 silent2) {
    int32 n;
    int32 rval = 0;
    double sum;
    int32 n0 = node + 1 + fix_var, k2;
    if (silent2 == 0) {
        for (int32 i = 0; i < NODE; i++) {
            char name[sizeof(uvar_names[i]) + 18];
            snprintf(name, sizeof(name), "Coupling for %s eqn:", uvar_names[i]);
            ggets_new_string(name, coup_string[i]);
            if (parserslow_add_expr(coup_string[i], coup_fun[i], &n)) {
                ggets_err_msg("Illegal formula");
                goto bye;
            }
        }
    }
    // formulae are fine .. lets do it ...
    for (int32 j = 0; j < nt; j++) {
        // j is phi variable
        sum = 0.0;

        for (int32 k = 0; k < nt; k++) {
            k2 = k + j;
            if (k2 >= nt) {
                k2 = k2 - nt + 1;
            }
            for (int32 i = 0; i < node; i++) {
                set_ivar(i + 1, (double)orb[i + 1][k]);
                set_ivar(i + n0 + 1, (double)orb[i + 1][k2]);
            }
            main_rhs_update_based_on_current();

            for (int32 i = 0; i < node; i++) {
                double z = evaluate(coup_fun[i]);

                sum = sum + z*adj[i + 1][k];
            }
        }
        my_h[0][j] = orb[0][j];
        my_h[1][j] = sum / (double)nt;
    }
    if (hodd_ev) {
        for (int32 k = 0; k < nt; k++) {
            k2 = nt - k - 1;
            my_h[2][k] = .5*(my_h[1][k] - my_h[1][k2]);
            my_h[3][k] = .5*(my_h[1][k] + my_h[1][k2]);
        }
    }
    rval = 1;

bye:
    NSYM = NSYM_START;
    NCON = NCON_START;
    return rval;
}

void
adjoints_new_adjoint(void) {
    int32 n = NODE + 1;
    if (adj_here) {
        adjoints_data_back();
        for (int32 i = 0; i < n; i++) {
            free(my_adj[i]);
        }
        free(my_adj);
        adj_here = 0;
    }
    adj_len = storind;
    my_adj = xmalloc((usize)(NEQ + 1)*sizeof(*my_adj));
    for (int32 i = 0; i < n; i++) {
        my_adj[i] = xmalloc(sizeof(*my_adj)*(usize)adj_len);
    }
    for (int32 i = n; i <= NEQ; i++) {
        my_adj[i] = storage[i];
    }
    if (adjoints_adjoint(storage, my_adj, adj_len, delta_t*NJMP, adj_eps, adj_err, adj_maxit,
                         NODE)) {
        adj_here = 1;
        adjoints_back();
    }
    ggets_ping();
    return;
}

/* ADJOINT ROUTINE
 *
 * This assumes that you have already computed the periodic orbit
 * and have stored in in an array **orbit
 * including time in the first column

 * The righthand sides of the equations are
 *   rhs_function(t,y,yp,n)
 * and the coupling function for ``H'' functions is
 *   couple(y,yhat,f,n)

 * where yhat is presynaptic and y is postynaptic
 * variable. f returns the coupling vector.

 * adjoint is the same size as orbit and when returned has
 * t in the first column. */

int32
adjoints_adjoint(double **orbit, double **adjnt, int32 nt, double dt, double eps, double minerr,
                 int32 maxit, int32 node) {
    double **jac, *yold, ytemp, *fold, *fdev;
    double *yprime;
    double *work;
    double t;
    double prod;
    double del;
    int32 l;
    int32 k2;
    int32 rval = 0;
    int32 n2 = node*node;
    double error;

    work = xmalloc((usize)(n2 + 4*node)*sizeof(*work));
    yprime = xmalloc((usize)node*sizeof(*yprime));
    yold = xmalloc((usize)node*sizeof(*yold));
    fold = xmalloc((usize)node*sizeof(*fold));
    fdev = xmalloc((usize)node*sizeof(*fdev));
    jac = xmalloc((usize)n2*sizeof(*jac));

    for (int32 i = 0; i < n2; i++) {
        jac[i] = xmalloc((usize)nt*sizeof(*jac));
        if (jac[i] == NULL) {
            ggets_err_msg("Insufficient storage");
            return 0;
        }
    }

    /* Now we compute the
     * transpose time reversed jacobian -- this is floatcomplex !! */
    for (int32 k = 0; k < nt; k++) {
        l = nt - 1 - k;  // reverse the limit cycle
        for (int32 i = 0; i < node; i++) {
            yold[i] = (double)orbit[i + 1][l];
        }
        rhs_function(0.0, yold, fold, node);
        for (int32 j = 0; j < node; j++) {
            ytemp = yold[j];
            del = eps*fabs(ytemp);
            if (del < eps) {
                del = eps;
            }

            yold[j] += del;
            rhs_function(0.0, yold, fdev, node);
            yold[j] = ytemp;
            for (int32 i = 0; i < node; i++) {
                jac[i + node*j][k] = (fdev[i] - fold[i]) / del;
            }
        }
    }

    // now we iterate to get a good adjoint using implicit Euler's method
    ytemp = 0.0;
    for (int32 i = 0; i < node; i++) {
        yold[i] = 1. + .01*(markov_ndrand48() - .5);  // random initial data

        ytemp += fabs(yold[i]);
    }
    for (int32 i = 0; i < node; i++) {
        yold[i] = yold[i] / ytemp;
        fdev[i] = yold[i];
    }

    ggets_plintf("%f %f \n", yold[0], yold[1]);

    for (l = 0; l < maxit; l++) {
        for (int32 k = 0; k < nt - 1; k++) {
            k2 = k + 1;
            if (k2 >= nt) {
                k2 = k2 - nt;
            }
            if (adjoints_step_eul(jac, k, k2, yold, work, node, dt) == 0) {
                rval = 0;
                goto bye;
            }
        }
        ytemp = 0.0;
        error = 0.0;

        for (int32 i = 0; i < node; i++) {
            if (fabs(yold[i]) > bound) {
                rval = 0;
                ggets_err_msg("Out of bounds");
                goto bye;
            }
            error += fabs(yold[i] - fdev[i]);
            ytemp += fabs(yold[i]);
        }

        for (int32 i = 0; i < node; i++) {
            yold[i] = yold[i] / ytemp;
            fdev[i] = yold[i];
        }
        printf("%f %f \n", yold[0], yold[1]);
        ggets_plintf("err=%f \n", error);
        if (error < minerr) {
            break;
        }
    }
    // onelast time to compute the adjoint
    prod = 0.0;  // for normalization
    t = 0.0;
    for (int32 k = 0; k < nt; k++) {
        l = nt - k - 1;
        t += dt;
        for (int32 i = 0; i < node; i++) {
            fdev[i] = (double)orbit[i + 1][l];
        }
        rhs_function(0.0, fdev, yprime, node);
        for (int32 j = 0; j < node; j++) {
            adjnt[j + 1][l] = (double)yold[j];
            prod += yold[j]*yprime[j]*dt;
        }
        k2 = k + 1;
        if (k2 >= nt) {
            k2 -= nt;
        }
        if (adjoints_step_eul(jac, k, k2, yold, work, node, dt) == 0) {
            rval = 0;
            goto bye;
        }
    }

    prod = prod / t;
    ggets_plintf(" Multiplying the adjoint by 1/%g to normalize\n", prod);
    for (int32 k = 0; k < nt; k++) {
        for (int32 j = 0; j < node; j++) {
            adjnt[j + 1][k] = adjnt[j + 1][k] / (double)prod;
        }
        adjnt[0][k] = orbit[0][k];
    }
    rval = 1;

bye:
    free(work);
    free(yprime);
    free(yold);
    free(fold);
    free(fdev);
    for (int32 i = 0; i < n2; i++) {
        free(jac[i]);
    }
    free(jac);
    return rval;
}

int32
adjoints_step_eul(double **jac, int32 k, int32 k2, double *yold, double *work, int32 node,
                  double dt) {
    int32 n2 = node*node;
    int32 info;
    int32 ipvt[MAX_ODE];
    double *mat;
    double *fold;
    fold = work;
    mat = work + node;

    for (int32 j = 0; j < node; j++) {
        fold[j] = 0.0;
        for (int32 i = 0; i < node; i++) {
            fold[j] = fold[j] + jac[i + j*node][k]*yold[i];
        }
    }
    for (int32 j = 0; j < node; j++) {
        yold[j] = yold[j] + .5*dt*fold[j];
    }
    for (int32 i = 0; i < n2; i++) {
        mat[i] = -jac[i][k2]*dt*.5;
    }
    for (int32 i = 0; i < node; i++) {
        mat[i + i*node] = 1. + mat[i + i*node];
    }
    gear_sgefa(mat, node, node, ipvt, &info);
    if (info != -1) {
        ggets_err_msg("Univertible Jacobian");
        return 0;
    }
    gear_sgesl(mat, node, node, ipvt, yold, 0);
    return 1;
}

/* this is some code for the maximal liapunov exponent
 * I assume you have computed an orbit and it is in storage

 * at each time point, I use y+dy as an initial condition
 * I then integrate for one time step
 * I subtract this from y(t+dt) and divide by the norm of dy.
 * I take the log of this and sum up the logs dividing by Ndt
 * to get an approximation */

void
adjoints_do_liapunov(void) {
    double z;
    double *x;
    ggets_new_int("Range over parameters?(0/1)", &LIAP_FLAG);
    if (LIAP_FLAG != 1) {
        adjoints_hrw_liapunov(&z, 0, NEWT_ERR);
        return;
    }
    x = &my_data[0];
    integrate_do_range(x, 0);
    // done the range
    for (int32 i = 0; i < LIAP_I; i++) {
        storage[0][i] = my_liap[0][i];
        storage[1][i] = my_liap[1][i];
    }
    storind = LIAP_I;
    browser_refresh(storind);
    LIAP_FLAG = 0;
    free(my_liap[0]);
    free(my_liap[1]);
    return;
}

void
adjoints_alloc_liap(int32 n) {
    if (LIAP_FLAG == 0) {
        return;
    }
    my_liap[0] = xmalloc(sizeof(*my_liap)*(usize)(n + 1));
    my_liap[1] = xmalloc(sizeof(*my_liap)*(usize)(n + 1));
    LIAP_N = (n + 1);
    LIAP_I = 0;
    return;
}

void
adjoints_do_this_liaprun(int32 i, double p) {
    double liap;
    if (LIAP_FLAG == 0) {
        return;
    }
    my_liap[0][i] = p;
    if (adjoints_hrw_liapunov(&liap, 1, NEWT_ERR)) {
        my_liap[1][i] = liap;
        LIAP_I++;
    }
    return;
}

int32
adjoints_hrw_liapunov(double *liap, int32 batch, double eps) {
    double y[MAX_ODE];
    double yp[MAX_ODE];
    double dy[MAX_ODE];
    double nrm;
    double t0;
    double t1;
    double sum = 0.0;
    char bob[256];
    int32 istart = 1;

    if (storind < 2) {
        if (batch == 0) {
            ggets_err_msg("You need to compute an orbit first");
        }
        return 0;
    }

    // lets make an initial random perturbation
    for (int32 i = 0; i < NODE; i++) {
        dy[i] = 0;
    }
    dy[0] = eps;

    for (int32 j = 0; j < (storind - 1); j++) {
        t0 = storage[0][j];
        t1 = storage[0][j + 1];
        istart = 1;
        for (int32 i = 0; i < NODE; i++) {
            y[i] = storage[i + 1][j] + dy[i];
        }
        do_fit_one_step_int(y, t0, t1, &istart);
        for (int32 i = 0; i < NODE; i++) {
            yp[i] = (y[i] - storage[i + 1][j + 1]);
        }

        {
            // adj2 norm vec
            double *v = yp;
            double *mu = &nrm;
            int32 n = NODE;

            double sum2 = 0.0;
            for (int32 i2 = 0; i2 < n; i2++) {
                sum += (v[i2]*v[i2]);
            }
            sum = sqrt(sum);
            if (sum > 0) {
                for (int32 i2 = 0; i2 < n; i2++) {
                    v[i2] = v[i2] / sum;
                }
            }
            *mu = sum2;
        }

        nrm = nrm / eps;
        if (nrm == 0.0) {
            if (batch == 0) {
                ggets_err_msg("Liapunov:-infinity exponent!");
            }
            fprintf(stderr, "Something wrong here: %s\n", __func__);
            exit(EXIT_FAILURE);
        }
        sum = sum + log(nrm);
        for (int32 i = 0; i < NODE; i++) {
            dy[i] = eps*yp[i];
        }
    }
    t1 = storage[0][storind - 1] - storage[0][0];
    if (fabs(t1) > 1e-12) {
        sum = sum / t1;
    }
    *liap = sum;

    if (batch == 0) {
        snprintf(bob, sizeof(bob), "Maximal exponent is %g", sum);
        ggets_err_msg(bob);
    }

    return 1;
}
