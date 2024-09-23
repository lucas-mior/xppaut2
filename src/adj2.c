#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include "xpplim.h"

/* this has a bunch of numerical routines
 * averaging
 * adjoints
 * transpose
 * maximal liapunov exponent */

#define READEM 1

extern double MyData[MAX_ODE];
extern int32 (*rhs)(double t, double *y, double *ydot, int32 neq);
extern double **storage;
extern int32 storind;
extern int32 FOUR_HERE;
extern int32 NODE;
extern int32 INFLAG;
extern int32 NEQ;
extern int32 NJMP;
extern int32 FIX_VAR;
extern int32 NMarkov;
extern int32 nvec;
extern double TEND;
static double **my_adj;
static int32 adj_len;
static double **my_h;
static double *my_liap[2];

extern char *info_message;
static struct {
    int32 here, col0, ncol, colskip;
    int32 row0, nrow, rowskip;
    double **data;
    char firstcol[11];
} my_trans;

static int32 LIAP_FLAG = 0;
static int32 LIAP_N;
static int32 LIAP_I;
static int32 LIAP_N;
static int32 LIAP_I;
extern double NEWT_ERR;
static double ADJ_EPS = 1.e-8, ADJ_ERR = 1.e-3;
static int32 ADJ_MAXIT = 20, ADJ_HERE = 0, H_HERE = 0, h_len, HODD_EV = 0;
int32 AdjRange = 0;
extern double DELTA_T;
extern double BOUND;
static int32 *coup_fun[MAX_ODE];
static char *coup_string[MAX_ODE];

extern int32 *my_ode[];
extern int32 NSYM;
extern int32 NSYM_START;
extern int32 NCON;
extern int32 NCON_START;
extern int32 DCURY;

static void adj2_h_back(void);
static void adj_back(void);
static void adj2_adjoint_parameters(void);
static int32 make_h(double **orb, double **adj, int32 nt, int32 node,
                    int32 silent);
static int32 adj2_step_eul(double **jac, int32 k, int32 k2, double *yold,
                           double *work, int32 node, double dt);
static void adj2_norm_vec(double *v, double *mu, int32 n);
static int32 adj2_hrw_liapunov(double *liap, int32 batch, double eps);

void
adj2_init_trans(void) {
    my_trans.here = 0;
    strcpy(my_trans.firstcol, uvar_names[0]);
    my_trans.ncol = 2;
    my_trans.nrow = 1;
    my_trans.rowskip = 1;
    my_trans.colskip = 1;
    my_trans.row0 = 1;
    my_trans.col0 = 2;
    return;
}

void
adj2_dump_transpose_info(FILE *fp, int32 f) {
    char bob[256];

    if (f == READEM)
        fgets(bob, 255, fp);
    else
        fprintf(fp, "# Transpose variables etc\n");

    io_string(my_trans.firstcol, 11, fp, f);
    io_int(&my_trans.ncol, fp, f, "n columns");
    io_int(&my_trans.nrow, fp, f, "n rows");
    io_int(&my_trans.rowskip, fp, f, "row skip");
    io_int(&my_trans.colskip, fp, f, "col skip");
    io_int(&my_trans.row0, fp, f, "row 0");
    return;
}

int32
adj2_do_transpose(void) {
    int32 ii;
    int32 status;
    static char *n[] = {"*0Column 1", "NCols", "ColSkip",
                        "Row 1",      "NRows", "RowSkip"};
    char values[LENGTH(n)][MAX_LEN_SBOX];

    snprintf(values[0], sizeof(values[0]), "%s", my_trans.firstcol);
    snprintf(values[1], sizeof(values[1]), "%d", my_trans.ncol);
    snprintf(values[2], sizeof(values[2]), "%d", my_trans.colskip);
    snprintf(values[3], sizeof(values[3]), "%d", my_trans.row0);
    snprintf(values[4], sizeof(values[4]), "%d", my_trans.nrow);
    snprintf(values[5], sizeof(values[5]), "%d", my_trans.rowskip);

    if (my_trans.here) {
        for (int32 i = 0; i <= my_trans.nrow; i++)
            free(my_trans.data[i]);
        free(my_trans.data);
        my_trans.here = 0;
        adj_data_back();
    }

    status = do_string_box(6, 6, 1, "Transpose Data", n, values, 33);
    if (status != 0) {
        int32 inrow, incol;

        find_variable(values[0], &ii);
        if (ii > -1)
            my_trans.col0 = ii + 1;
        else {
            err_msg("No such columns");
            return 0;
        }
        strcpy(my_trans.firstcol, values[0]);
        ii = atoi(values[4]);
        if (ii >= NEQ)
            ii = NEQ - 1;
        my_trans.nrow = ii;
        my_trans.ncol = atoi(values[1]);
        my_trans.colskip = atoi(values[2]);
        my_trans.row0 = atoi(values[3]);
        my_trans.rowskip = atoi(values[5]);

        my_trans.data = xmalloc(sizeof(*(my_trans.data))*(usize)(NEQ + 1));
        for (int32 i = 0; i <= my_trans.nrow; i++)
            my_trans.data[i] =
                xmalloc(sizeof(my_trans.data[i])*(usize)my_trans.ncol);
        for (int32 i = my_trans.nrow + 1; i <= NEQ; i++)
            my_trans.data[i] = storage[i];
        for (int32 j = 0; j < my_trans.ncol; j++)
            my_trans.data[0][j] = (double)(j + 1);

        for (int32 i = 0; i < my_trans.ncol; i++) {
            incol = my_trans.col0 - 1 + i*my_trans.colskip;
            if (incol > NEQ)
                incol = NEQ;
            for (int32 j = 0; j < my_trans.nrow; j++) {
                inrow = my_trans.row0 + j*my_trans.rowskip;
                if (inrow > storind)
                    inrow = storind;
                my_trans.data[j + 1][i] = storage[incol][inrow];
            }
        }

        set_browser_data(my_trans.data, 1);
        refresh_browser(my_trans.ncol);
        my_trans.here = 1;
        return 1;
    }
    return 0;
}

void
adj2_alloc_h_stuff(void) {
    for (int32 i = 0; i < NODE; i++) {
        coup_fun[i] = xmalloc(100*sizeof(*coup_fun));
        coup_string[i] = xmalloc(80);
        strcpy(coup_string[i], "0");
    }
    return;
}

void
adj_data_back(void) {
    FOUR_HERE = 0;
    set_browser_data(storage, 1);
    refresh_browser(storind);
    return;
}

void
adj_back(void) {
    if (ADJ_HERE) {
        set_browser_data(my_adj, 1);
        refresh_browser(adj_len);
    }
    return;
}

void
adj2_h_back(void) {
    if (H_HERE) {
        set_browser_data(my_h, 1);
        refresh_browser(h_len);
    }
    return;
}

/*  Here is how to do the range over adjoints and h functions
    unfortunately, h functions are always computed even if you dont want them
    they will just be zeros

    Step 1. Compute a singel orbit, adjoint, and H function (to load the
    program with the correct right-hand sides for H function. Or just load in
    set file where it was done
    Step 2.  Set transient to some reasonable number to assure convergence
    onto the limit cycle as you change parameters and total to be at least
    2 periods beyond the transient
    Step 3. Se up Poincare map - period - stop on section. This lets you
    get the period
    Step 4. In numerics averaging - click on adjrange
    Step 5. Initconds range over the parameter. It should find the periodic
    orbit, adjoint, and H function and save. Files are of the form
    orbit.parname_parvalue.dat etc
*/
void
adj2_make_adj_com(int32 com) {
    static char key[] = "nmaohpr";
    switch (key[com]) {
    case 'n':
        adj2_new_adjoint();
        break;
    case 'm':
        adj2_new_h_fun(0);
        break;
    case 'a':
        adj_back();
        break;
    case 'o':
        adj_data_back();
        break;
    case 'h':
        adj2_h_back();
        break;
    case 'p':
        adj2_adjoint_parameters();
        break;
    case 'r':
        AdjRange = 1;
        break;
    default:
        break;
    }
    return;
}

void
adj2_adjoint_parameters(void) {
    new_int("Maximum iterates :", &ADJ_MAXIT);
    new_float("Adjoint error tolerance :", &ADJ_ERR);
    return;
}

void
adj2_new_h_fun(int32 silent) {
    int32 n = 2;
    if (!ADJ_HERE) {
        err_msg("Must compute adjoint first!");
        return;
    }
    if (storind != adj_len) {
        err_msg("incompatible data and adjoint");
        return;
    }
    if (H_HERE) {
        free(my_h[0]);
        free(my_h[1]);
        if (HODD_EV) {
            free(my_h[2]);
            free(my_h[3]);
        }
        free(my_h);
        H_HERE = 0;
        HODD_EV = 0;
    }
    if (NEQ > 2) {
        HODD_EV = 1;
        n = 4;
    }
    h_len = storind;
    adj_data_back();
    my_h = xmalloc(sizeof(*my_h)*(usize)(NEQ + 1));
    for (int32 i = 0; i < n; i++)
        my_h[i] = xmalloc(sizeof(*my_h)*(usize)h_len);
    for (int32 i = n; i <= NEQ; i++)
        my_h[i] = storage[i];
    if (make_h(storage, my_adj, h_len, NODE, silent)) {
        H_HERE = 1;
        adj2_h_back();
    }
    ping();
    return;
}

void
adj2_dump_h_stuff(FILE *fp, int32 f) {
    char bob[256];
    if (f == READEM)
        fgets(bob, 255, fp);
    else
        fprintf(fp, "# Coupling stuff for H funs\n");
    for (int32 i = 0; i < NODE; i++)
        io_string(coup_string[i], 79, fp, f);
    return;
}

int32
make_h(double **orb, double **adj, int32 nt, int32 node, int32 silent) {
    int32 j, rval = 0;
    double sum;
    double z;
    int32 n0 = node + 1 + FIX_VAR, k2, k;
    if (silent == 0) {
        for (int32 i = 0; i < NODE; i++) {
            char name[sizeof(uvar_names[i]) + 18];
            snprintf(name, sizeof(name), "Coupling for %s eqn:", uvar_names[i]);
            new_string(name, coup_string[i]);
            if (add_expr(coup_string[i], coup_fun[i], &j)) {
                err_msg("Illegal formula");
                goto bye;
            }
        }
    }
    /*  formulae are fine .. lets do it ... */
    for (j = 0; j < nt; j++) {
        /* j is phi variable  */
        sum = 0.0;

        for (k = 0; k < nt; k++) {
            k2 = k + j;
            if (k2 >= nt)
                k2 = k2 - nt + 1;
            for (int32 i = 0; i < node; i++) {
                set_ivar(i + 1, (double)orb[i + 1][k]);
                set_ivar(i + n0 + 1, (double)orb[i + 1][k2]);
            }
            z = 0.0;
            update_based_on_current();

            for (int32 i = 0; i < node; i++) {
                z = evaluate(coup_fun[i]);

                sum = sum + (double)z*adj[i + 1][k];
            }
        }
        my_h[0][j] = orb[0][j];
        my_h[1][j] = sum / (double)nt;
    }
    if (HODD_EV) {
        for (k = 0; k < nt; k++) {
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
adj2_new_adjoint(void) {
    int32 n = NODE + 1;
    if (ADJ_HERE) {
        adj_data_back();
        for (int32 i = 0; i < n; i++)
            free(my_adj[i]);
        free(my_adj);
        ADJ_HERE = 0;
    }
    adj_len = storind;
    my_adj = xmalloc((usize)(NEQ + 1)*sizeof(*my_adj));
    for (int32 i = 0; i < n; i++)
        my_adj[i] = xmalloc(sizeof(*my_adj)*(usize)adj_len);
    for (int32 i = n; i <= NEQ; i++)
        my_adj[i] = storage[i];
    if (adj2_adjoint(storage, my_adj, adj_len, DELTA_T*NJMP, ADJ_EPS, ADJ_ERR,
                     ADJ_MAXIT, NODE)) {
        ADJ_HERE = 1;
        adj_back();
    }
    ping();
    return;
}

/*    ADJOINT ROUTINE
 *
 *    This assumes that you have already computed the periodic orbit
 *      and have stored in in an array **orbit
 *    including time in the first column

 *    The righthand sides of the equations are
 *      rhs(t,y,yp,n)
 *    and the coupling function for ``H'' functions is
 *      couple(y,yhat,f,n)

 *      where yhat is presynaptic and y is postynaptic
 *   variable.  f returns the coupling vector.

 *  adjoint is the same size as orbit and when returned has
 *  t in the first column.  */

int32
adj2_adjoint(double **orbit, double **adjnt, int32 nt, double dt, double eps,
             double minerr, int32 maxit, int32 node) {
    double **jac, *yold, ytemp, *fold, *fdev;
    double *yprime;
    double *work;
    double t, prod, del;
    int32 j, l, k2, rval = 0;
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
            err_msg("Insufficient storage");
            return 0;
        }
    }

    /*  Now we compute the
          transpose time reversed jacobian  --  this is floatcomplex !! */
    for (int32 k = 0; k < nt; k++) {
        l = nt - 1 - k; /* reverse the limit cycle  */
        for (int32 i = 0; i < node; i++)
            yold[i] = (double)orbit[i + 1][l];
        rhs(0.0, yold, fold, node);
        for (j = 0; j < node; j++) {
            ytemp = yold[j];
            del = eps*fabs(ytemp);
            if (del < eps)
                del = eps;

            yold[j] += del;
            rhs(0.0, yold, fdev, node);
            yold[j] = ytemp;
            for (int32 i = 0; i < node; i++)
                jac[i + node*j][k] = (fdev[i] - fold[i]) / del;
        }
    }

    /* now we iterate to get a good adjoint using implicit Euler's method */
    ytemp = 0.0;
    for (int32 i = 0; i < node; i++) {
        yold[i] = 1. + .01*(ndrand48() - .5); /* random initial data */

        ytemp += fabs(yold[i]);
    }
    for (int32 i = 0; i < node; i++) {
        yold[i] = yold[i] / ytemp;
        fdev[i] = yold[i];
    }

    plintf("%f %f \n", yold[0], yold[1]);

    for (l = 0; l < maxit; l++) {
        for (int32 k = 0; k < nt - 1; k++) {
            k2 = k + 1;
            if (k2 >= nt)
                k2 = k2 - nt;
            if (adj2_step_eul(jac, k, k2, yold, work, node, dt) == 0) {
                rval = 0;
                goto bye;
            }
        }
        ytemp = 0.0;
        error = 0.0;

        for (int32 i = 0; i < node; i++) {
            if (fabs(yold[i]) > BOUND) {
                rval = 0;
                err_msg("Out of bounds");
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
        plintf("err=%f \n", error);
        if (error < minerr)
            break;
    }
    /*  onelast time to compute the adjoint  */
    prod = 0.0; /* for normalization   */
    t = 0.0;
    for (int32 k = 0; k < nt; k++) {
        l = nt - k - 1;
        t += dt;
        for (int32 i = 0; i < node; i++)
            fdev[i] = (double)orbit[i + 1][l];
        rhs(0.0, fdev, yprime, node);
        for (j = 0; j < node; j++) {
            adjnt[j + 1][l] = (double)yold[j];
            prod += yold[j]*yprime[j]*dt;
        }
        k2 = k + 1;
        if (k2 >= nt)
            k2 -= nt;
        if (adj2_step_eul(jac, k, k2, yold, work, node, dt) == 0) {
            rval = 0;
            goto bye;
        }
    }

    prod = prod / t;
    plintf(" Multiplying the adjoint by 1/%g to normalize\n", prod);
    for (int32 k = 0; k < nt; k++) {
        for (j = 0; j < node; j++)
            adjnt[j + 1][k] = adjnt[j + 1][k] / (double)prod;
        adjnt[0][k] = orbit[0][k];
    }
    rval = 1;

bye:
    free(work);
    free(yprime);
    free(yold);
    free(fold);
    free(fdev);
    for (int32 i = 0; i < n2; i++)
        free(jac[i]);
    free(jac);
    return rval;
}

int32
adj2_step_eul(double **jac, int32 k, int32 k2, double *yold, double *work,
              int32 node, double dt) {
    int32 j, i, n2 = node*node, info;
    int32 ipvt[MAX_ODE];
    double *mat;
    double *fold;
    fold = work;
    mat = work + node;

    for (j = 0; j < node; j++) {
        fold[j] = 0.0;
        for (i = 0; i < node; i++)
            fold[j] = fold[j] + jac[i + j*node][k]*yold[i];
    }
    for (j = 0; j < node; j++)
        yold[j] = yold[j] + .5*dt*fold[j];
    for (i = 0; i < n2; i++)
        mat[i] = -jac[i][k2]*dt*.5;
    for (i = 0; i < node; i++)
        mat[i + i*node] = 1. + mat[i + i*node];
    sgefa(mat, node, node, ipvt, &info);
    if (info != -1) {
        err_msg("Univertible Jacobian");
        return 0;
    }
    sgesl(mat, node, node, ipvt, yold, 0);
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
adj2_do_liapunov(void) {
    double z;
    int32 i;
    double *x;
    new_int("Range over parameters?(0/1)", &LIAP_FLAG);
    if (LIAP_FLAG != 1) {
        adj2_hrw_liapunov(&z, 0, NEWT_ERR);
        return;
    }
    x = &MyData[0];
    do_range(x, 0);
    /* done the range */
    for (i = 0; i < LIAP_I; i++) {
        storage[0][i] = my_liap[0][i];
        storage[1][i] = my_liap[1][i];
    }
    storind = LIAP_I;
    refresh_browser(storind);
    LIAP_FLAG = 0;
    free(my_liap[0]);
    free(my_liap[1]);
    return;
}

void
adj2_alloc_liap(int32 n) {
    if (LIAP_FLAG == 0)
        return;
    my_liap[0] = xmalloc(sizeof(*my_liap)*(usize)(n + 1));
    my_liap[1] = xmalloc(sizeof(*my_liap)*(usize)(n + 1));
    LIAP_N = (n + 1);
    LIAP_I = 0;
    return;
}

void
adj2_do_this_liaprun(int32 i, double p) {
    double liap;
    if (LIAP_FLAG == 0)
        return;
    my_liap[0][i] = p;
    adj2_hrw_liapunov(&liap, 1, NEWT_ERR);
    my_liap[1][i] = liap;
    LIAP_I++;
    return;
}

void
adj2_norm_vec(double *v, double *mu, int32 n) {
    int32 i;
    double sum = 0.0;
    for (i = 0; i < n; i++)
        sum += (v[i]*v[i]);
    sum = sqrt(sum);
    if (sum > 0)
        for (i = 0; i < n; i++)
            v[i] = v[i] / sum;
    *mu = sum;
    return;
}

int32
adj2_hrw_liapunov(double *liap, int32 batch, double eps) {
    double y[MAX_ODE];
    double yp[MAX_ODE], nrm, dy[MAX_ODE];
    double t0;
    double t1;
    double sum = 0.0;
    char bob[256];
    int32 istart = 1;
    int32 i;
    int32 j;
    if (storind < 2) {
        if (batch == 0)
            err_msg("You need to compute an orbit first");
        return 0;
    }

    /* lets make an initial random perturbation */
    for (i = 0; i < NODE; i++)
        dy[i] = 0;
    dy[0] = eps;

    for (j = 0; j < (storind - 1); j++) {
        t0 = storage[0][j];
        t1 = storage[0][j + 1];
        istart = 1;
        for (i = 0; i < NODE; i++)
            y[i] = storage[i + 1][j] + dy[i];
        one_step_int(y, t0, t1, &istart);
        for (i = 0; i < NODE; i++)
            yp[i] = (y[i] - storage[i + 1][j + 1]);
        adj2_norm_vec(yp, &nrm, NODE);
        nrm = nrm / eps;
        if (nrm == 0.0) {
            if (batch == 0)
                err_msg("Liapunov:-infinity exponent!");
            return 0; /* something wrong here */
        }
        sum = sum + log(nrm);
        for (i = 0; i < NODE; i++)
            dy[i] = eps*yp[i];
    }
    t1 = storage[0][storind - 1] - storage[0][0];
    if (fabs(t1) > 1e-12)
        sum = sum / t1;
    *liap = sum;
    if (batch == 0) {
        snprintf(bob, sizeof(bob), "Maximal exponent is %g", sum);
        err_msg(bob);
    }

    return 1;
}
