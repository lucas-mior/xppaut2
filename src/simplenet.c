#include "functions.h"
#include "xmalloc.h"
#include "integers.h"
#include <stdbool.h>

#include "parserslow.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * n is the number of values to return
 * ncon is the number of connections each guy gets
 * index[n][ncon] is the list of indices to which it connects
 * weight[n][ncon] is the list of weights
 * root is the address of lowest entry in the variable array
 *
 * To wit:
 * for sparse type
 * name(i) produces values[i]
 * value[i] = simplenet_sum(k=0..ncon-1)of(w[i][k]*root[index[i][k]])
 *
 * ODE FILE CALL:
 *
 * special f=sparse( n, ncon, w,index,rootname)
 * special f=sconv(type,n,ncon,w,rootname)
 * index,w are tables that are loaded by the "tabular"
 * command.
 * rootname here is the name of the root quantity.
 * type is either ep0
 * for conv.
 *
 * conv0 conve convp type
 *
 * value[i]=    simplenet_sum(j=-ncon;j<=ncon;i++){
 *               k=i+j;
 *         0 type      if(k>=0 && k<n)
 *               value[i]+=wgt[j+ncon]*rootname[k]
 *         e type k=ABS(i+j); if(k<2n)
 *                                if(k>=n)k=n-k;
 *                                ...
 *         p type
 *            k=mod(2*n+i+j,n)
 *
 * for example  discretized diffusion
 * tabular dd % 3 -1 1 3*ABS(t)-2
 * special diff=conv(even, 51, 2, dd,v0)
 * v[0..50]'=f(v[j],w[j])+d*diff([j])
 *
 * another example
 * nnetwork
 *
 * tabular wgt % 51 -25 25 .032*cos(.5*pi*t/25)
 * special stot=conv(0, 51, 25, wgt,v0)
 * v[0..50]'=f(v[j],w[j])-gsyn*stot([j])*(v([j])-vsyn)
 *
 * last example -- random sparse network 51 cells 5 connections each
 *
 * tabular w % 251 0 250 .2*rand(1)
 * tabular con % 251 0 250 flr(rand(1)*5)
 * special stot=sparse(51,5,w,con,v0)
 * v[0..50]'=f(v[j],w[j])-gsyn*stot([j])*(v([j])-vsyn)
 *
 * more stuff:
 * special k= delmmult(n,m,w,tau,root)
 *           w has n*m values and tau has n*m delays
 *           n is the length of root (cols) and m is number of rows
 *           k(i) = simplenet_sum(j=0,n-1) w[i*n+j] delay(root(j),tau[i*n+j])
 *           Note delays can only work on variables and not on
 *           fixed  values
 * special k=delsparse(m,nc,w,index,tau,root)
 *           m is the number of return values
 *           nc is number of connections per node. Must be the sam
 *           w is  n*nc is weights
 *           index is  m*nc that has the indices of connections
 *
 *           tau is m*nc is list of delays
 *
 *        k[i]=simplenet_sum( 0<=j < nc)
 * w[j+nc*i]*delay(root[c[j+nc*i]],tau[j+nc*i])
 *
 * special f=mmult(n,m,w,root)  -- note that here m is the number of return
 *                           values and n is the size of input root
 * f[j] = simplenet_sum(i=0,n-1)of(w(i+n*j)*root(i)) j=0,..,m-1
 * special f=fmmult(n,m,w,root1,root2,fname)
 * f[j] = simplenet_sum(i=0,n-1)of(w(i+n*j)*fname(root1[i],root2[j])
 * special f=fsparse( n, ncon,w,index,root1,root2,fname)
 * special f=fconv(type,n,ncon,w,root1,root2,fname)
 * simplenet_sum(j=-ncon,ncon)w(j)*fname(root1[i+j],root2[i])
 * similarly for fsparse
 *
 * special k=fftcon(type,n,wgt,v0)
 *
 * uses the fft to convolve v0 with the wgt which must be of
 * length n if type=periodic
 *        2n if type=0
 *
 * special f=interp(type,n,root)
 *         this produces an interpolated value of x for
 *         x \in [0,n)
 *         type =0 for linear (all that works now)
 *         root0,....rootn-1  are values at the integers
 *         Like a table, but the int64 values are variables
 *
 * special f=findext(type,n,skip,root)
 * if type=1  mx(0)=maximum mx(1)=index
 * if type=-1 mx(2)=minimum mx(3)=index
 * if type=0 mx(0)=maximum mx(1)=index,mx(2)=minimum mx(3)=index
 *
 * special ydot=import(soname,sofun,nret,root,w1,w2,...wm)
 *
 *  run right hand side in C
 *  soname is shared object library say mlnet.so
 *  sofun  is shared object function
 *  nret is the number of return values
 *  root is the name of the first variable
 *
 * sofun(int32 nret, int32 root, double *con, double *var, double *z[50],double
 * *ydot)
 *
 * *z[50] contains a list of pointers  z[0] -> w1, .... 50 is hard coded
 *
 * NOTE that the user-defined parameters start at #6 and are in order
 * including derived parameters but XPP takes care of this so start at 0 */

#include <stdio.h>

#define EVEN 0
#define PERIODIC 2
#define MAXW 50

/* simple network stuff */

#define MAXVEC 100

typedef struct Vectorizer {
    char name[20];
    int32 root;
    int32 length;
    int32 il;
    int32 ir;
} Vectorizer;

static Vectorizer my_vec[MAXVEC];
static int32 n_vector = 0;

typedef struct Network {
    int32 type;
    int32 ncon;
    int32 n;
    char name[20];
    char soname[256];
    char sofun[256];

    int32 root;
    int32 root2;
    int32 f[20];
    int32 iwgt;
    int32 *gcom;  // for group commands

    double *values, *weight, *index, *taud;  // for delays
    double *fftr, *ffti, *dr, *di;
    double *wgtlist[MAXW];
} Network;

#define CONVE 0
#define CONV0 1
#define CONVP 2
#define SPARSE 3
#define FCONVE 10
#define FCONV0 11
#define FCONVP 12
#define FSPARSE 13
#define FFTCON0 4
#define FFTCONP 5
#define MMULT 6
#define FMMULT 7
#define GILLTYPE 25
#define INTERP 30
#define FINDEXT 35   // find extrema in list of variables
#define DEL_MUL 40   // for delayed coupled networks  - global coupling
#define DEL_SPAR 41  // sparse with unequal in degree and delays
#define IMPORT 50    // not really a network type

static int32 simplenet_g_namelist(char *s, char *root, int32 *flag, int32 *i1, int32 *i2);
static int32 simplenet_gil_parse(char *s, int32 *ind, int32 *nn);
static void simplenet_update_fft(int32 ind);
static int32 simplenet_is_network(char *s);
static void simplenet_init(double *v, int32 n);
static double simplenet_interp(double x, int32 i);
static int32 simplenet_parse_import(char *s, char *soname, char *sofun, int32 *n, char *vname,
                                    int32 *m, char *tname[MAXW]);
static int32 simplenet_get_vector_info(char *str, char *name, int32 *root, int32 *length, int32 *il,
                                       int32 *ir);
static int32 simplenet_get_imp_str(char *in, int32 *i, char *out);
static int32 simplenet_import_error(void);

static Network my_net[MAX_NET];
static int32 n_network = 0;
double
simplenet_interp(double x, int32 i) {
    int32 jlo = (int32)x;
    double *y;
    int32 n = my_net[i].n;
    double dx = x - (double)jlo;
    y = &variables[my_net[i].root];
    if (jlo < 0 || jlo > (n - 1)) {
        return 0.0;  // out of range
    }
    return (1 - dx)*y[jlo] + dx*y[jlo + 1];
}

int32
simplenet_add_vectorizer(char *name, char *rhs) {
    int32 i;
    int32 ivar;
    int32 il;
    int32 ir;
    int32 ind;
    int32 len;
    int32 flag;

    for (i = 0; i < n_vector; i++) {
        if (strcmp(name, my_vec[i].name) == 0) {
            break;
        }
    }

    ind = i;
    flag = simplenet_get_vector_info(rhs, name, &ivar, &len, &il, &ir);

    if (flag == 0) {
        return 0;
    }

    my_vec[ind].root = ivar;
    my_vec[ind].length = len;
    my_vec[ind].il = il;
    my_vec[ind].ir = ir;
    ggets_plintf("adding vector %s based on variable %d of length %d ends %d %d\n", name, ivar, len,
                 il, ir);

    return 1;
}

void
simplenet_add_vectorizer_name(char *name, char *rhs) {
    (void)rhs;
    if (n_vector >= MAXVEC) {
        ggets_plintf("Too many vectors \n");
        exit(0);
    }
    strcpy(my_vec[n_vector].name, name);
    if (parserslow_add_vector_name(n_vector, name)) {
        exit(0);
    }
    n_vector++;
    return;
}

double
simplenet_vector_value(double x, int32 i) {
    int32 il = my_vec[i].il, ir = my_vec[i].ir, n = my_vec[i].length, k = (int32)x;
    int32 root = my_vec[i].root;

    if ((k >= 0) && (k < n)) {
        return variables[root + k];
    }
    if (il == PERIODIC) {
        return variables[root + ((k + n) % n)];
    }
    if (k < 0) {
        if (il == 0.0) {
            return 0.0;
        }
        return variables[root - k - 1];
    }
    if (k >= n) {
        if (ir == 0.0) {
            return 0.0;
        }
        return variables[2*n - k - 1 + root];
    }
    return 0;
}

double
simplenet_network_value(double x, int32 i) {
    int32 j = (int32)x;
    if (my_net[i].type == INTERP) {
        return simplenet_interp(x, i);
    }
    if (j >= 0 && j < my_net[i].n) {
        return my_net[i].values[j];
    }
    return 0.0;
}

void
simplenet_init(double *v, int32 n) {
    for (int32 i = 0; i < n; i++) {
        v[i] = 0.0;
    }
    return;
}

int32
simplenet_add_spec_fun(char *name, char *rhs) {
    int32 i;
    int32 ind;
    int32 elen;
    int32 type;
    int32 iwgt;
    int32 itau;
    int32 iind;
    int32 ivar;
    int32 ivar2;
    int32 ntype;
    int32 ntot;
    int32 ncon;
    int32 ntab;
    char *str;
    char junk[256];
    char rootname[20];
    char wgtname[20];
    char tauname[20];
    char indname[20];
    char root2name[20];
    char fname[20];
    char sofun[256], soname[256], *tname[MAXW];
    type = simplenet_is_network(rhs);
    if (type == 0) {
        return 0;
    }
    ggets_plintf("type=%d \n", type);
    for (i = 0; i < n_network; i++) {
        if (strcmp(name, my_net[i].name) == 0) {
            break;
        }
    }
    ind = i;
    if (ind >= n_network) {
        ggets_plintf(" No such name %s ?? \n", name);
        return 0;
    }
    switch (type) {
    case 1:  // convolution
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = -1;
        if (str[0] == 'E') {
            ntype = CONVE;
        }
        if (str[0] == '0' || str[0] == 'Z') {
            ntype = CONV0;
        }
        if (str[0] == 'P') {
            ntype = CONVP;
        }
        if (ntype == -1) {
            ggets_plintf(" No such convolution type %s \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ntot = atoi(str);
        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);
        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);
        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        ggets_plintf(" Added net %s type %d len=%d x %d using %s var[%d] \n", name, ntype, ntot,
                     ncon, wgtname, ivar);

        return 1;
    case 2:  // sparse
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = SPARSE;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(indname, str);
        iind = find_lookup(indname);

        if (iind < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, indname);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].index = my_table[iind].y;

        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        ggets_plintf(" Added sparse %s len=%d x %d using %s var[%d]  and %s\n", name, ntot, ncon,
                     wgtname, ivar, indname);
        return 1;
    case 3:  // convolution
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = -1;
        if (str[0] == 'E') {
            ntype = FCONVE;
        }
        if (str[0] == '0' || str[0] == 'Z') {
            ntype = FCONV0;
        }
        if (str[0] == 'P') {
            ntype = FCONVP;
        }
        if (ntype == -1) {
            ggets_plintf(" No such convolution type %s \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ntot = atoi(str);
        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);
        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);
        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(root2name, str);
        ivar2 = get_var_index(root2name);
        if (ivar2 < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, root2name);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(fname, str);
        snprintf(junk, sizeof(junk), "%s(%s,%s)", fname, rootname, root2name);
        if (parserslow_add_expr(junk, my_net[ind].f, &elen)) {
            ggets_plintf(" bad function %s \n", fname);
            return 0;
        }
        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].type = ntype;
        my_net[ind].root = my_net[ind].f[0];  // this is strange - I am adding the compiled names
        my_net[ind].root2 = my_net[ind].f[1];
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        ggets_plintf(" Added net %s type %d len=%d x %d using %s %s(var[%d],var[%d]) \n", name,
                     ntype, ntot, ncon, wgtname, fname, ivar, ivar2);
        return 1;
    case 4:  // sparse
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = FSPARSE;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(indname, str);
        iind = find_lookup(indname);

        if (iind < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, indname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(root2name, str);
        ivar2 = get_var_index(root2name);
        if (ivar2 < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, root2name);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(fname, str);
        snprintf(junk, sizeof(junk), "%s(%s,%s)", fname, rootname, root2name);
        if (parserslow_add_expr(junk, my_net[ind].f, &elen)) {
            ggets_plintf(" bad function %s \n", fname);
            return 0;
        }

        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].index = my_table[iind].y;

        my_net[ind].type = ntype;
        my_net[ind].root = my_net[ind].f[0];  // this is strange - I am adding the compiled names
        my_net[ind].root2 = my_net[ind].f[1];
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        ggets_plintf(" Sparse %s len=%d x %d using %s %s(var[%d],var[%d]) and %s\n", name, ntot,
                     ncon, wgtname, fname, ivar, ivar2, indname);
        return 1;
    case 5:  // fft convolution
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = -1;
        if (str[0] == '0' || str[0] == 'Z') {
            ntype = FFTCON0;
        }
        if (str[0] == 'P') {
            ntype = FFTCONP;
        }
        if (ntype == -1) {
            ggets_plintf(" No such fft convolution type %s \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ntot = atoi(str);
        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);
        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }
        ntab = tabular_get_lookup_len(iwgt);
        if (type == FFTCONP && ntab < ntot) {
            ggets_plintf(" In %s, weight is length %d < %d \n", name, ntab, ntot);
            return 0;
        }
        if (type == FFTCON0 && ntab < (2*ntot)) {
            ggets_plintf(" In %s, weight is length %d < %d \n", name, ntab, 2*ntot);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        if (ntype == FFTCON0) {
            ncon = 2*ntot;
        } else {
            ncon = ntot;
        }
        my_net[ind].fftr = xmalloc((usize)(ncon + 2)*sizeof(*(my_net[ind].fftr)));
        my_net[ind].ffti = xmalloc((usize)(ncon + 2)*sizeof(*(my_net[ind].ffti)));
        my_net[ind].dr = xmalloc((usize)(ncon + 2)*sizeof(*(my_net[ind].dr)));
        my_net[ind].di = xmalloc((usize)(ncon + 2)*sizeof(*(my_net[ind].di)));
        my_net[ind].iwgt = iwgt;
        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        simplenet_update_fft(ind);

        ggets_plintf(" Added net %s type %d len=%d x %d using %s var[%d] \n", name, ntype, ntot,
                     ncon, wgtname, ivar);
        return 1;
    case 6:  // MMULT    ntot=n,ncon=m
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = MMULT;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        my_net[ind].values = xmalloc((usize)(ncon + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ncon);
        my_net[ind].weight = my_table[iwgt].y;

        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ncon;
        my_net[ind].ncon = ntot;
        ggets_plintf(" Added mmult %s len=%d x %d using %s var[%d]\n", name, ntot, ncon, wgtname,
                     ivar, indname);
        return 1;
    case 7:  // FMMULT
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = FMMULT;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(root2name, str);
        ivar2 = get_var_index(root2name);
        if (ivar2 < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, root2name);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(fname, str);
        snprintf(junk, sizeof(junk), "%s(%s,%s)", fname, rootname, root2name);
        if (parserslow_add_expr(junk, my_net[ind].f, &elen)) {
            ggets_plintf(" bad function %s \n", fname);
            return 0;
        }
        my_net[ind].values = xmalloc((usize)(ncon + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ncon);
        my_net[ind].weight = my_table[iwgt].y;

        my_net[ind].type = ntype;
        my_net[ind].root = my_net[ind].f[0];  // this is strange - I am adding the compiled names
        my_net[ind].root2 = my_net[ind].f[1];
        my_net[ind].n = ncon;
        my_net[ind].ncon = ntot;
        ggets_plintf(" Added fmmult %s len=%d x %d using %s %s(var[%d],var[%d])\n", name, ntot,
                     ncon, wgtname, fname, ivar, ivar2);
        return 1;

    case FINDEXT:
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = atoi(str);
        if (ntype > 1 || ntype < (-1)) {
            ggets_plintf("In %s,  type =-1,0,1 not %s \n", name, ntype);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ntot = atoi(str);
        if (ntot <= 0) {
            ggets_plintf("In %s,  n>0 not %s \n", name, ntot);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);
        if (ncon <= 0) {
            ggets_plintf("In %s,  skip>=1 not %s \n", name, ncon);
            return 0;
        }
        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        my_net[ind].values = xmalloc(6*sizeof(*(my_net[ind].values)));
        my_net[ind].type = FINDEXT;
        my_net[ind].root = ivar;
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        my_net[ind].iwgt = ntype;
        ggets_plintf(" Added findextr %s: type=%d len=%d  skip= %d using var[%d] \n", name, ntype,
                     ntot, ncon, ivar);
        return 1;

    case 30:
        /* interpolation array
           z=INTERP(meth,n,root)
        */
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ivar = atoi(str);
        my_net[ind].type = INTERP;
        my_net[ind].iwgt = ivar;
        str = form_ode_do_fit_get_next(",");
        ivar = atoi(str);
        if (ivar < 1) {
            ggets_plintf("Need more than 1 entry for interpolate\n");
            return 0;
        }
        my_net[ind].n = ivar;  // # entries in array
        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        my_net[ind].root = ivar;
        ggets_plintf("Added interpolator %s length %d on %s \n", name, my_net[ind].n, rootname);
        return 1;

    case IMPORT:
        ntype = IMPORT;
        for (i = 0; i < MAXW; i++) {
            tname[i] = xmalloc(25);
        }
        simplenet_parse_import(rhs, soname, sofun, &ncon, rootname, &ntab, tname);
        my_net[ind].values = xmalloc((usize)(ncon + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ncon);
        my_net[ind].n = ncon;
        ivar = get_var_index(rootname);
        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }
        strcpy(my_net[ind].soname, soname);
        strcpy(my_net[ind].sofun, sofun);
        my_net[ind].root = ivar;
        my_net[ind].type = ntype;
        my_net[ind].ncon = 0;
        for (i = 0; i < ntab; i++) {
            iwgt = find_lookup(tname[i]);
            ggets_plintf("Found %s\n", tname[i]);
            if (iwgt < 0) {
                ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
                return 0;
            }
            my_net[ind].wgtlist[i] = my_table[iwgt].y;
        }
        for (i = 0; i < MAXW; i++) {
            free(tname[i]);
        }
        ggets_plintf(" Added import %s len=%d  with %s %s var[%d] %d weights\n", name,
                     my_net[ind].n, soname, sofun, ivar, ntab);

        return 1;
    case DEL_MUL:

        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = DEL_MUL;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(tauname, str);
        itau = find_lookup(tauname);

        if (itau < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, tauname);
            return 0;
        }

        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        my_net[ind].values = xmalloc((usize)(ncon + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ncon);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].taud = my_table[itau].y;

        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ncon;
        my_net[ind].ncon = ntot;
        ggets_plintf(" Added del_mul %s len=%d x %d using %s var[%d] with delay %s\n", name, ntot,
                     ncon, wgtname, ivar, indname, tauname);
        NDELAYS = 1;
        return 1;
    case DEL_SPAR:
        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ntype = DEL_SPAR;
        ntot = atoi(str);

        if (ntot <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        ncon = atoi(str);

        if (ncon <= 0) {
            ggets_plintf(" %s must be positive int32 \n", str);
            return 0;
        }
        str = form_ode_do_fit_get_next(",");
        strcpy(wgtname, str);
        iwgt = find_lookup(wgtname);

        if (iwgt < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, wgtname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(indname, str);
        iind = find_lookup(indname);

        if (iind < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, indname);
            return 0;
        }

        str = form_ode_do_fit_get_next(",");
        strcpy(tauname, str);
        itau = find_lookup(tauname);

        if (itau < 0) {
            ggets_plintf("in network %s,  %s is not a table \n", name, tauname);
            return 0;
        }

        str = form_ode_do_fit_get_next(")");
        strcpy(rootname, str);
        ivar = get_var_index(rootname);

        if (ivar < 0) {
            ggets_plintf(" In %s , %s is not valid variable\n", name, rootname);
            return 0;
        }

        my_net[ind].values = xmalloc((usize)(ntot + 1)*sizeof(*(my_net[ind].values)));
        simplenet_init(my_net[ind].values, ntot);
        my_net[ind].weight = my_table[iwgt].y;
        my_net[ind].index = my_table[iind].y;
        my_net[ind].taud = my_table[itau].y;
        my_net[ind].type = ntype;
        my_net[ind].root = ivar;
        my_net[ind].n = ntot;
        my_net[ind].ncon = ncon;
        ggets_plintf(" Added sparse %s len=%d x %d using %s var[%d]  and %s with "
                     "dely %s\n",
                     name, ntot, ncon, wgtname, ivar, indname, tauname);
        NDELAYS = 1;
        return 1;
    case 10:
        /*
           z=GILL(meth,rxn list)
           e.g
           z=GILL(meth,r{1-15})
           GILL is different -
           iwgt=evaluation method - 0 is standard
                                    1 - tau-leap
           root=number of reactions
           values[0]=time of next reaction
           values[1..root]=number of times this rxn took place
           gcom contains list of all the fixed holding the reactions
        */

        form_ode_get_first(rhs, "(");
        str = form_ode_do_fit_get_next(",");
        ivar = atoi(str);
        str = form_ode_do_fit_get_next(")");
        my_net[ind].type = GILLTYPE;
        if (ivar > 0) {
            ggets_plintf(" Tau leaping not implemented yet. Changing to 0\n");
            ivar = 0;
        }
        my_net[ind].iwgt = ivar;
        my_net[ind].gcom = xmalloc(1000*sizeof(*(my_net[ind].gcom)));
        if (simplenet_gil_parse(str, my_net[ind].gcom, &ivar2) == 0) {
            return 0;
        }
        my_net[ind].root = ivar2;
        my_net[ind].n = ivar2 + 1;
        my_net[ind].ncon = -1;
        my_net[ind].values = xmalloc((usize)(ivar2 + 2)*sizeof(*(my_net[ind].values)));
        ggets_plintf("Added gillespie chain with %d reactions \n", ivar2);
        return 1;

        /*  case 8:
        form_ode_get_first(rhs,"(");
        str=form_ode_do_fit_get_next(",");
        ntot=atoi(str);
        str=form_ode_do_fit_get_next("{");
        i=0;
        elen=strlen(str);

        while(1){
          cc=str[i];
          if(cc=='}'){junk[i]=0;
                       break;
          }
          junk[i]=cc;
          i++;
          if(i==elen){
            printf("Illegal syntax for GROUP %s \n",str);
            return 0;
          }

        }
        ggets_plintf("total=%d str=%s\n",ntot,junk);

        return 0; */
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
}

void
simplenet_add_special_name(char *name, char *rhs) {
    if (simplenet_is_network(rhs)) {
        ggets_plintf(" netrhs = |%s| \n", rhs);
        if (n_network >= MAX_NET) {
            return;
        }
        strcpy(my_net[n_network].name, name);
        parserslow_add_net_name(n_network, name);
        n_network++;
    } else {
        ggets_plintf(" No such special type ...\n");
    }
    return;
}

int32
simplenet_is_network(char *s) {
    ani_de_space(s);
    strupr(s);
    if (s[0] == 'C' && s[1] == 'O' && s[2] == 'N' && s[3] == 'V') {
        return 1;
    }
    if (s[0] == 'S' && s[1] == 'P' && s[2] == 'A' && s[3] == 'R') {
        return 2;
    }
    if (s[0] == 'F' && s[1] == 'C' && s[2] == 'O' && s[3] == 'N' && s[4] == 'V') {
        return 3;
    }
    if (s[0] == 'F' && s[1] == 'S' && s[2] == 'P' && s[3] == 'A' && s[4] == 'R') {
        return 4;
    }
    if (s[0] == 'F' && s[1] == 'F' && s[2] == 'T' && s[3] == 'C') {
        return 5;
    }
    if (s[0] == 'M' && s[1] == 'M' && s[2] == 'U' && s[3] == 'L') {
        return 6;
    }
    if (s[0] == 'F' && s[1] == 'M' && s[2] == 'M' && s[3] == 'U' && s[4] == 'L') {
        return 7;
    }
    if (s[0] == 'G' && s[1] == 'I' && s[2] == 'L' && s[3] == 'L') {
        return 10;
    }
    if (s[0] == 'I' && s[1] == 'N' && s[2] == 'T' && s[3] == 'E' && s[4] == 'R') {
        return INTERP;
    }
    if (s[0] == 'F' && s[1] == 'I' && s[2] == 'N' && s[3] == 'D' && s[4] == 'E') {
        return FINDEXT;
    }
    if (s[0] == 'D' && s[1] == 'E' && s[2] == 'L' && s[3] == 'M') {
        return DEL_MUL;
    }
    if (s[0] == 'D' && s[1] == 'E' && s[2] == 'L' && s[3] == 'S') {
        return DEL_SPAR;
    }
    if (s[0] == 'I' && s[1] == 'M' && s[2] == 'P' && s[3] == 'O') {
        return IMPORT;
    }
    return 0;
}

void
simplenet_eval_all_nets(void) {
    for (int32 ind = 0; ind < n_network; ind++) {
        // evaluate network
        int32 i;
        int32 k;
        int32 ij;
        int32 imin;
        int32 imax;
        double ymin;
        double ymax;
        int32 skip;
        int32 mmt;
        int32 in0;
        double sum;
        double z;
        int32 n = my_net[ind].n, *f;
        int32 ncon = my_net[ind].ncon;
        double *w, *y, *cc, *values, *tau;
        int32 twon = 2*n;
        int32 root = my_net[ind].root;
        int32 root2 = my_net[ind].root2;
        cc = my_net[ind].index;
        w = my_net[ind].weight;
        values = my_net[ind].values;
        switch (my_net[ind].type) {
        case FINDEXT:
            y = &variables[root];
            mmt = my_net[ind].iwgt;
            skip = ncon;
            imax = 0;
            ymax = y[0];
            imin = 0;
            ymin = y[0];
            i = 0;

            if (mmt == 1 || mmt == 0) {  // get max
                while (i < n) {
                    if (y[i] > ymax) {
                        imax = i;
                        ymax = y[i];
                    }
                    i += skip;
                }
            }
            if (mmt == (-1) || mmt == 0) {  // get min
                i = 0;
                while (i < n) {
                    if (y[i] < ymin) {
                        imin = i;
                        ymin = y[i];
                    }
                    i += skip;
                }
            }
            values[1] = (double)imax;
            values[0] = ymax;
            values[3] = (double)imin;
            values[2] = ymin;
            break;
        case INTERP:  // do nothing!
            break;
        case GILLTYPE:
            if (my_net[ind].ncon == -1 && my_net[ind].iwgt > 0) {
                my_net[ind].weight = xmalloc((usize)(my_net[ind].root*NODE)*sizeof(double));
                markov_make_gill_nu(my_net[ind].weight, NODE, my_net[ind].root, my_net[ind].values);
                my_net[ind].ncon = 0;
            }
            markov_one_gill_step(my_net[ind].iwgt, my_net[ind].root, my_net[ind].gcom,
                                 my_net[ind].values);
            break;
        case CONVE:
            y = &variables[root];
            for (i = 0; i < n; i++) {
                sum = 0.0;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = ABS(i + j);
                    if (k < twon) {
                        if (k >= n) {
                            k = ABS(twon - 2 - k);
                        }
                        sum += (w[j + ncon]*y[k]);
                    }
                }
                values[i] = sum;
            }
            break;
        case CONV0:
            y = &variables[root];
            for (i = 0; i < n; i++) {
                sum = 0.0;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = i + j;
                    if (k < n && k >= 0) {
                        sum += (w[j + ncon]*y[k]);
                    }
                }
                values[i] = sum;
            }
            break;
        case CONVP:
            y = &variables[root];
            for (i = 0; i < n; i++) {
                sum = 0.0;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = ((twon + i + j) % n);
                    sum += (w[j + ncon]*y[k]);
                }
                values[i] = sum;
            }
            break;
        case FFTCONP:
            y = &variables[root];
            simplenet_fft_conv(0, n, values, y, my_net[ind].fftr, my_net[ind].ffti, my_net[ind].dr,
                               my_net[ind].di);
            break;

        case FFTCON0:
            y = &variables[root];
            simplenet_fft_conv(1, n, values, y, my_net[ind].fftr, my_net[ind].ffti, my_net[ind].dr,
                               my_net[ind].di);
            break;

        case IMPORT:
            extra_get_import_values(n, values, my_net[ind].soname, my_net[ind].sofun,
                                    my_net[ind].root, my_net[ind].wgtlist, variables,
                                    &constants[6]);
            break;
        case DEL_MUL:
            tau = my_net[ind].taud;
            in0 = my_net[ind].root;
            for (int32 j = 0; j < n; j++) {
                sum = 0.0;
                for (i = 0; i < ncon; i++) {
                    ij = j*ncon + i;
                    sum += (w[ij]*delay_handle_get_delay(i + in0, tau[ij]));
                }
                values[j] = sum;
            }
            break;
        case MMULT:
            y = &variables[root];
            for (int32 j = 0; j < n; j++) {
                sum = 0.0;
                for (i = 0; i < ncon; i++) {
                    ij = j*ncon + i;
                    sum += (w[ij]*y[i]);
                }
                values[j] = sum;
            }
            break;
        case DEL_SPAR:
            tau = my_net[ind].taud;
            in0 = my_net[ind].root;
            for (i = 0; i < n; i++) {
                sum = 0.0;
                for (int32 j = 0; j < ncon; j++) {
                    ij = i*ncon + j;
                    k = (int32)cc[ij];
                    if (k >= 0) {
                        sum += (w[ij]*delay_handle_get_delay(k + in0, tau[ij]));
                    }
                }
                values[i] = sum;
            }
            break;
        case SPARSE:
            y = &variables[root];
            for (i = 0; i < n; i++) {
                sum = 0.0;
                for (int32 j = 0; j < ncon; j++) {
                    ij = i*ncon + j;
                    k = (int32)cc[ij];
                    if (k >= 0) {
                        sum += (w[ij]*y[k]);
                    }
                }
                values[i] = sum;
            }
            break;

            //     f stuff
        case FCONVE:
            f = my_net[ind].f;

            for (i = 0; i < n; i++) {
                sum = 0.0;
                f[1] = root + i;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = ABS(i + j);
                    if (k < twon) {
                        if (k >= n) {
                            k = ABS(twon - 2 - k);
                        }
                        f[0] = root2 + k;

                        z = evaluate(f);
                        sum += (w[j + ncon]*z);
                    }
                }
                values[i] = sum;
            }
            break;
        case FCONV0:
            f = my_net[ind].f;

            for (i = 0; i < n; i++) {
                sum = 0.0;
                f[1] = root + i;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = i + j;
                    if (k < n && k >= 0) {
                        f[0] = root2 + k;
                        z = evaluate(f);
                        sum += (w[j + ncon]*z);
                    }
                }
                values[i] = sum;
            }
            break;
        case FCONVP:
            f = my_net[ind].f;

            for (i = 0; i < n; i++) {
                f[1] = root + i;
                sum = 0.0;
                for (int32 j = -ncon; j <= ncon; j++) {
                    k = ((twon + i + j) % n);
                    f[0] = root2 + k;
                    z = evaluate(f);
                    sum += (w[j + ncon]*z);
                }
                values[i] = sum;
            }
            break;
        case FSPARSE:
            f = my_net[ind].f;

            for (i = 0; i < n; i++) {
                f[1] = root + i;
                sum = 0.0;
                for (int32 j = 0; j < ncon; j++) {
                    ij = i*ncon + j;
                    k = (int32)cc[ij];
                    if (k >= 0) {
                        f[0] = root2 + k;
                        z = evaluate(f);
                        sum += (w[ij]*z);
                    }
                }
                values[i] = sum;
            }
            break;
        case FMMULT:
            f = my_net[ind].f;

            for (int32 j = 0; j < n; j++) {
                f[1] = root + j;

                sum = 0.0;
                for (i = 0; i < ncon; i++) {
                    ij = j*ncon + i;

                    f[0] = root2 + i;

                    z = evaluate(f);

                    sum += (w[ij]*z);
                }

                values[j] = sum;
            }
            break;
        default:
            fprintf(stderr, "Unexpected case in %s.\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
    return;
}

void
simplenet_update_all_ffts(void) {
    for (int32 i = 0; i < n_network; i++) {
        if (my_net[i].type == FFTCON0 || my_net[i].type == FFTCONP) {
            simplenet_update_fft(i);
        }
    }
    return;
}
/*
 tabular weights are of size 2k+1
 and go from -k ... k
 for FFT's
 they are reordered as follows
 fftr[i]=wgt[i+k] i = 0.. k
 fftr[i+k]=wgt[i] i=1 .. k-1
*/
void
simplenet_update_fft(int32 ind) {
    int32 dims[2];
    double *w = my_net[ind].weight;
    double *fftr = my_net[ind].fftr;
    double *ffti = my_net[ind].ffti;
    int32 n;
    int32 n2;
    int32 type = my_net[ind].type;
    if (type == FFTCONP) {
        n = my_net[ind].n;
        n2 = n / 2;
        for (int32 i = 0; i < n; i++) {
            ffti[i] = 0.0;
        }
        for (int32 i = 0; i <= n2; i++) {
            fftr[i] = w[i + n2];
        }
        for (int32 i = 0; i < n2; i++) {
            fftr[n2 + i + 1] = w[i];
        }
        dims[0] = n;
        fftn(1, dims, fftr, ffti, 1, 1.);
    }
    if (type == FFTCON0) {
        n = 2*my_net[ind].n;
        n2 = n / 2;
        for (int32 i = 0; i < n; i++) {
            ffti[i] = 0.0;
        }
        for (int32 i = 0; i <= n2; i++) {
            fftr[i] = w[i + n2];
        }
        for (int32 i = 1; i < n2; i++) {
            fftr[n2 + i] = w[i];
        }
        dims[0] = n;
        fftn(1, dims, fftr, ffti, 1, 1.);
    }
    return;
}

void
simplenet_fft_conv(int32 it, int32 n, double *values, double *yy, double *fftr, double *ffti,
                   double *dr, double *di) {
    int32 dims[2];
    double x;
    double y;
    int32 n2 = 2*n;
    switch (it) {
    case 0:
        dims[0] = n;
        for (int32 i = 0; i < n; i++) {
            di[i] = 0.0;
            dr[i] = yy[i];
        }

        fftn(1, dims, dr, di, 1, -2.0);

        for (int32 i = 0; i < n; i++) {
            x = dr[i]*fftr[i] - di[i]*ffti[i];
            y = dr[i]*ffti[i] + di[i]*fftr[i];
            dr[i] = x;
            di[i] = y;
        }

        fftn(1, dims, dr, di, -1, -2.0);
        for (int32 i = 0; i < n; i++) {
            values[i] = dr[i];
        }

        return;
    case 1:
        dims[0] = n2;
        for (int32 i = 0; i < n2; i++) {
            di[i] = 0.0;
            if (i < n) {
                dr[i] = yy[i];
            } else {
                dr[i] = 0.0;
            }
        }
        fftn(1, dims, dr, di, 1, -2.0);
        for (int32 i = 0; i < n2; i++) {
            x = dr[i]*fftr[i] - di[i]*ffti[i];
            y = dr[i]*ffti[i] + di[i]*fftr[i];
            dr[i] = x;
            di[i] = y;
        }
        fftn(1, dims, dr, di, -1, -2.0);
        for (int32 i = 0; i < n; i++) {
            values[i] = dr[i];
        }
        return;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
}

/* parsing stuff to get gillespie code quickly */

int32
simplenet_gil_parse(char *s, int32 *ind, int32 *nn) {
    int32 i = 0;
    int32 n = (int32)strlen(s);
    char piece[50];
    char b[20];
    char bn[25];
    char c;
    int32 i1;
    int32 i2;
    int32 jp = 0;
    int32 f;
    int32 k = 0;
    int32 iv;
    int32 id;
    int32 m;
    ggets_plintf("s=|%s|", s);
    while (true) {
        c = s[i];
        if (c == ',' || i > (n - 1)) {
            piece[jp] = 0;
            if (simplenet_g_namelist(piece, b, &f, &i1, &i2) == 0) {
                printf("Bad gillespie list %s\n", s);
                return 0;
            }
            if (f == 0) {
                ggets_plintf("added %s\n", b);
                iv = get_var_index(b);
                if (iv < 0) {
                    ggets_plintf("No such name %s\n", b);
                    return 0;
                }
                ind[k] = iv;
                k++;
            } else {
                ggets_plintf("added %s{%d-%d}\n", b, i1, i2);
                m = i2 - i1 + 1;
                for (id = 0; id < m; id++) {
                    snprintf(bn, sizeof(bn), "%s%d", b, id + i1);
                    iv = get_var_index(bn);
                    if (iv < 0) {
                        ggets_plintf("No such name %s\n", bn);
                        return 0;
                    }
                    ind[k] = iv;
                    k++;
                }
            }
            if (i > (n - 1)) {
                *nn = k;
                return 1;
            }
            jp = 0;
        } else {
            piece[jp] = c;
            jp++;
        }
        i++;
    }
}

/* plucks info out of  xxx{aa-bb}  or returns string */
int32
simplenet_g_namelist(char *s, char *root, int32 *flag, int32 *i1, int32 *i2) {
    int32 i;
    int32 n = (int32)strlen(s);
    int32 ir = -1;
    int32 j = 0;
    char c;
    char num[20];
    *flag = 0;
    for (i = 0; i < n; i++) {
        if (s[i] == '{') {
            ir = i;
        }
    }
    if (ir < 0) {
        strcpy(root, s);
        return 1;
    }
    for (i = 0; i < ir; i++) {
        root[i] = s[i];
    }
    root[ir] = 0;
    *flag = 1;
    j = 0;
    for (i = ir + 1; i < n; i++) {
        c = s[i];
        if (c == '-') {
            break;
        }
        num[j] = c;
        j++;
    }
    if (i == n) {
        ggets_plintf("Illegal syntax %s\n", s);
        return 0;
    }
    num[j] = 0;
    *i1 = atoi(num);
    ir = i + 1;
    j = 0;
    for (i = ir; i < n; i++) {
        c = s[i];
        if (c == '}') {
            break;
        }
        num[j] = c;
        j++;
    }
    num[j] = 0;
    *i2 = atoi(num);
    return 1;
}

int32
simplenet_get_imp_str(char *in, int32 *i, char *out) {
    int32 j = 0;
    int32 done = 1;
    char c;
    int32 k = 0;
    while (done > 0) {
        c = in[*i];

        if ((c == ',') || (c == ')')) {
            out[j] = 0;
            *i = *i + 1;
            done = 0;
            if (c == ')') {
                k = 1;
            } else {
                k = 0;
            }

        } else {
            out[j] = c;
            j++;
            *i = *i + 1;
        }
    }

    return k;
}

int32
simplenet_import_error(void) {
    printf("k=import(soname,sofun,nret,var0,w1,...,wm)");
    return 0;
}

int32
simplenet_parse_import(char *s, char *soname, char *sofun, int32 *n, char *vname, int32 *m,
                       char *tname[MAXW]) {
    char temp[256];
    int32 j;
    char c;

    int32 i = 0;
    int32 done = 1;

    while (done > 0) {
        c = s[i];
        i++;
        if (c == '(') {
            done = 0;
        }
    }

    j = simplenet_get_imp_str(s, &i, temp);
    strcpy(soname, temp);
    if (j == 1) {
        return simplenet_import_error();
    }

    j = simplenet_get_imp_str(s, &i, temp);
    strcpy(sofun, temp);
    if (j == 1) {
        return simplenet_import_error();
    }

    j = simplenet_get_imp_str(s, &i, temp);
    *n = atoi(temp);
    if (j == 1 || *n <= 0) {
        return simplenet_import_error();
    }

    j = simplenet_get_imp_str(s, &i, temp);
    strcpy(vname, temp);
    *m = 0;
    if (j == 1) {
        printf("No weights....\n");
        return 1;
    }

    done = 1;

    while (done > 0) {
        j = simplenet_get_imp_str(s, &i, temp);
        strcpy(tname[*m], temp);
        *m = *m + 1;
        if (j == 1) {
            done = 0;
        }
    }
    return 1;
}

int32
simplenet_get_vector_info(char *str, char *name, int32 *root, int32 *length, int32 *il, int32 *ir) {
    int32 i = 0;
    int32 ivar;
    int32 n = (int32)strlen(str);
    char c;
    int32 j;
    char temp[100];
    ani_de_space(str);
    for (i = 0; i < n; i++) {
        if (str[i] == '(') {
            break;
        }
    }
    i++;
    j = 0;
    while (true) {
        c = str[i];
        if (c == ',') {
            i++;
            break;
        }
        temp[j] = c;
        i++;
        j++;
    }
    temp[j] = 0;
    ivar = get_var_index(temp);

    if (ivar < 0) {
        ggets_plintf(" In vector %s , %s is not valid variable\n", name, temp);
        return 0;
    }
    *root = ivar;
    j = 0;
    while (true) {
        c = str[i];
        if (c == ',') {
            i++;
            break;
        }
        temp[j] = c;
        i++;
        j++;
    }
    temp[j] = 0;
    n = atoi(temp);

    *length = n;
    *il = PERIODIC;
    if (str[i] == 'e' || str[i] == 'E') {
        *il = EVEN;
    }
    if (str[i] == 'z' || str[i] == 'Z') {
        *il = 0.0;
    }
    i++;
    i++;
    *ir = PERIODIC;

    if (str[i] == 'e' || str[i] == 'E') {
        *ir = EVEN;
    }
    if (str[i] == 'z' || str[i] == 'Z') {
        *ir = 0.0;
    }

    return 1;
}
