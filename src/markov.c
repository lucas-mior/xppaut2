#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>

#include <strings.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "xpplim.h"
#include "parserslow.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-12
#define RNMX (1.0 - EPS)
#define PI 3.1415926

static int64 myrandomseed = -1;

typedef struct Markov {
    int32 **command;
    char **trans;
    double *fixed;
    int32 nstates;
    double *states;
    int32 type; /* 0 is default and state dependent.  1 is fixed for all time */
    char name[12];
} Markov;

static Markov markov[MAX_MARK];

static double *my_mean[MAX_ODE];
static double *my_variance[MAX_ODE];
static int32 stoch_len;

int32 STOCH_FLAG;
static int32 STOCH_HERE;
static int32 N_TRIALS;
static int32 Wiener[MAX_PAR];
int32 NWiener;

static double markov_ran1(long *idum);
static double markov_gammln(double xx);
static void markov_init_stoch(int32 len);
static double markov_new_state(double old, int32 index, double dt);
static void update_markov(double *x, double t, double dt);
static int32 compile_markov(int32 index, int32 j, int32 k);
static void add_markov_entry(int32 index, int32 j, int32 k, char *expr);
static void create_markov(int32 nstates, double *st, int32 type, char *name);
static void markov_extract_expr(char *source, char *dest, int32 *i0);

void
markov_add_wiener(int32 index) {
    Wiener[NWiener] = index;
    NWiener++;
    return;
}

void
markov_set_wieners(double dt, double *x, double t) {
    int32 i;
    update_markov(x, t, fabs(dt));
    for (i = 0; i < NWiener; i++)
        constants[Wiener[i]] = normal(0.00, 1.00) / sqrt(fabs(dt));
    return;
}

void
markov_add(int32 nstate, char *name) {
    double st[50];
    int32 i;
    for (i = 0; i < 50; i++)
        st[i] = (double)i;
    create_markov(nstate, st, 0, name);
    return;
}

int32
build_markov(
    /*   FILE *fptr; */
    char **ma, char *name) {
    /*int32 nn;
     */
    int32 len = 0, ll;
    char line[256], expr[256];
    int32 istart;

    int32 i, j, nstates, index;
    index = -1;
    /* find it -- if not defined, then abort  */
    for (i = 0; i < NMarkov; i++) {
        ll = (int32)strlen(markov[i].name);
        if (strncasecmp(name, markov[i].name, (usize)ll) == 0) {
            if (len < ll) {
                index = i;
                len = ll;
            }
        }
    }
    if (index == -1) {
        ggets_plintf(" Markov variable |%s| not found \n", name);
        exit(0);
    }
    /* get number of states  */
    nstates = markov[index].nstates;
    if (ConvertStyle)
        fprintf(convertf, "markov %s %d\n", name, nstates);
    ggets_plintf(" Building %s %d states...\n", name, nstates);
    for (i = 0; i < nstates; i++) {
        /* fgets(line,256,fptr); */
        snprintf(line, sizeof(line), "%s", ma[i]);
        if (ConvertStyle)
            fprintf(convertf, "%s", line);
        /*nn=strlen(line)+1;*/
        /* if((save_eqn[NLINES]=xmalloc(nn))==NULL){
          ggets_plintf("saveeqn-prob\n");exit(0);}
          strncpy(save_eqn[NLINES++],line,nn); */
        istart = 0;
        for (j = 0; j < nstates; j++) {
            markov_extract_expr(line, expr, &istart);
            ggets_plintf("%s ", expr);
            add_markov_entry(index, i, j, expr);
        }
        ggets_plintf("\n");
    }
    return index;
}

int32
markov_old_build(FILE *fptr, char *name) {
    /*int32 nn;*/
    int32 len = 0, ll;
    char line[256], expr[256];
    int32 istart;

    int32 i, j, nstates, index;
    index = -1;
    /* find it -- if not defined, then abort  */
    for (i = 0; i < NMarkov; i++) {
        ll = (int32)strlen(markov[i].name);
        if (strncasecmp(name, markov[i].name, (usize)ll) == 0) {
            if (len < ll) {
                index = i;
                len = ll;
            }
        }
    }
    if (index == -1) {
        ggets_plintf(" Markov variable |%s| not found \n", name);
        exit(0);
    }
    /* get number of states  */
    nstates = markov[index].nstates;
    if (ConvertStyle)
        fprintf(convertf, "markov %s %d\n", name, nstates);
    ggets_plintf(" Building %s ...\n", name);
    for (i = 0; i < nstates; i++) {
        fgets(line, 256, fptr);

        if (ConvertStyle)
            fprintf(convertf, "%s", line);
        /*nn=strlen(line)+1;*/
        /* if((save_eqn[NLINES]=xmalloc(nn))==NULL)exit(0);
           strncpy(save_eqn[NLINES++],line,nn); */
        istart = 0;
        for (j = 0; j < nstates; j++) {
            markov_extract_expr(line, expr, &istart);
            ggets_plintf("%s ", expr);
            add_markov_entry(index, i, j, expr);
        }
        ggets_plintf("\n");
    }
    return index;
}

void
markov_extract_expr(char *source, char *dest, int32 *i0) {
    char ch;
    int32 len = 0;
    int32 flag = 0;
    while (true) {
        ch = source[*i0];
        *i0 = *i0 + 1;
        if (ch == '}')
            break;
        if (ch == '{')
            flag = 1;
        else {
            if (flag) {
                dest[len] = ch;
                len++;
            }
        }
    }
    dest[len] = 0;
    return;
}

void
create_markov(int32 nstates, double *st, int32 type, char *name) {
    int32 i;
    int32 n2 = nstates*nstates;
    int32 j = NMarkov;
    if (j >= MAX_MARK) {
        ggets_plintf("Too many Markov chains...\n");
        exit(0);
    }

    markov[j].nstates = nstates;
    markov[j].states = xmalloc((usize)nstates*sizeof(*(markov[j].states)));
    if (type == 0) {
        markov[j].trans = xmalloc((usize)n2*sizeof(char *));
        markov[j].command = xmalloc((usize)n2*sizeof(int32 *));
    } else {
        markov[j].fixed = xmalloc((usize)n2*sizeof(*(markov[j].fixed)));
    }

    for (i = 0; i < nstates; i++)
        markov[j].states[i] = st[i];
    strcpy(markov[j].name, name);
    NMarkov++;
    return;
}

void
add_markov_entry(int32 index, int32 j, int32 k, char *expr) {
    int32 l0 = markov[index].nstates*j + k;
    int32 type = markov[index].type;
    if (type == 0) {
        markov[index].trans[l0] =
            xmalloc(sizeof(*(markov[index].trans[l0]))*(strlen(expr) + 1));
        strcpy(markov[index].trans[l0], expr);
        /*  compilation step -- can be delayed */
        /*
         if(add_expr(expr,com,&leng)){
           ggets_plintf("Illegal expression %s\n",expr);
           exit(0);
         }
         markov[index].command[l0]=xmalloc(sizeof(int32)*(leng+2));
         for(i=0;i<leng;i++){
           markov[index].command[l0][i]=com[i];

         }
        */
        /*  end of compilation   */

    } else {
        markov[index].fixed[l0] = atof(expr);
    }
    return;
}

void
markov_compile_all(void) {
    int32 index, j, k, ns, l0;
    if (NMarkov == 0)
        return;
    for (index = 0; index < NMarkov; index++) {
        ns = markov[index].nstates;
        for (j = 0; j < ns; j++) {
            for (k = 0; k < ns; k++) {
                l0 = ns*j + k;
                if (compile_markov(index, j, k) == -1) {
                    ggets_plintf("Bad expression %s[%d][%d] = %s \n",
                                 markov[index].name, j, k,
                                 markov[index].trans[l0]);
                    exit(0);
                }
            }
        }
    }
    return;
}

int32
compile_markov(int32 index, int32 j, int32 k) {
    char *expr;
    int32 l0 = markov[index].nstates*j + k;
    int32 leng2;
    int32 i;
    int32 com[256];
    expr = markov[index].trans[l0];

    if (add_expr(expr, com, &leng2))
        return -1;
    markov[index].command[l0] =
        xmalloc(sizeof(*(markov[index].command[l0]))*(usize)(leng2 + 2));
    for (i = 0; i < leng2; i++) {
        markov[index].command[l0][i] = com[i];
    }

    return 1;
}

void
update_markov(double *x, double t, double dt) {
    int32 i;
    double yp[MAX_ODE];
    if (NMarkov == 0)
        return;
    set_ivar(0, t);
    for (i = 0; i < NODE; i++)
        set_ivar(i + 1, x[i]);
    for (i = NODE + FIX_VAR; i < NODE + FIX_VAR + NMarkov; i++)
        set_ivar(i + 1, x[i - FIX_VAR]);
    for (i = NODE; i < NODE + FIX_VAR; i++)
        set_ivar(i + 1, evaluate(my_ode[i]));
    for (i = 0; i < NMarkov; i++)
        yp[i] = markov_new_state(x[NODE + i], i, dt);
    for (i = 0; i < NMarkov; i++) {
        x[NODE + i] = yp[i];
        set_ivar(i + NODE + FIX_VAR + 1, yp[i]);
    }
    return;
}

double
markov_new_state(double old, int32 index, double dt) {
    double prob;
    double sum;
    double coin = markov_ndrand48();
    int32 row = -1, rns;
    double *st;
    int32 i, ns = markov[index].nstates;
    int32 type = markov[index].type;
    st = markov[index].states;
    for (i = 0; i < ns; i++)
        if (fabs(st[i] - old) < .0001) {
            row = i;
            break;
        }
    if (row == -1)
        return old;
    rns = row*ns;
    sum = 0.0;
    if (type == 0) {
        for (i = 0; i < ns; i++) {
            if (i != row) {
                prob = evaluate(markov[index].command[rns + i])*dt;
                sum = sum + prob;
                if (coin <= sum) {
                    return st[i];
                }
            }
        }
    } else {
        for (i = 0; i < ns; i++) {
            if (i != row) {
                prob = markov[index].fixed[rns + i]*dt;
                sum = sum + prob;
                if (coin <= sum) {
                    return st[i];
                }
            }
        }
    }

    return old;
}

void
markov_make_gill_nu(double *nu, int32 n, int32 m, double *v) {
    /* nu[j+m*i] = nu_{i,j} i=1,n-1 -- assume first eqn is tr'=tr+z(0)
       i species j reaction
      need this for improved tau stepper
     */
    double *y, *yp, *yold;
    int32 ir;
    int32 iy;

    y = xmalloc((usize)n*sizeof(*y));
    yold = xmalloc((usize)n*sizeof(*yold));
    yp = xmalloc((usize)n*sizeof(*yp));
    for (ir = 0; ir < m; ir++)
        v[ir + 1] = 0;
    rhs_only(yold);
    for (ir = 0; ir < m; ir++) {
        v[ir + 1] = 1;
        rhs_only(yp);
        for (iy = 0; iy < n; iy++) {
            nu[ir + m*iy] = yp[iy];
            ggets_plintf("ir=%d iy=%d nu=%g\n", ir + 1, iy, yp[iy] - yold[iy]);
        }
        v[ir + 1] = 0;
    }

    free(y);
    free(yp);
    free(yold);
    return;
}

void
markov_one_gill_step(int32 meth, int32 nrxn, int32 *rxn, double *v) {
    double rate = 0, test;
    double r[1000];
    /*double rold[1000]; Not used*/

    int32 i;
    switch (meth) {
    case 0: /* std gillespie method */
        for (i = 0; i < nrxn; i++) {
            v[i + 1] = 0.0;
            r[i] = get_ivar(rxn[i]);
            rate += r[i];
        }
        if (rate <= 0.0)
            return;
        v[0] = -log(markov_ndrand48()) / rate; /* next step */
        test = rate*markov_ndrand48();
        rate = r[0];
        for (i = 0; i < nrxn; i++) {
            if (test < rate) {
                v[i + 1] = 1.0;
                break;
            }
            rate += r[i + 1];
        }
        break;
    case 1: /* tau stepping method  */
        perror("Tau stepping method not implemented yet.");
        /*for(i=0;i<nrxn;i++)
          rold[i]=get_ivar(rxn[i]);
            */
        break;
    default:
        break;
    }
    return;
}

void
markov_do_stochast_com(int32 i) {
    static char key[] = "ncdmvhofpislaxe2";
    char ch = key[i];

    if (ch == 27)
        return;
    switch (ch) {
    case 'n':
        ggets_new_int("Seed:", &RandSeed);
        markov_nsrand48(RandSeed);
        break;
    case 'd':
        adj2_data_back();
        break;
    case 'm':
        markov_mean_back();
        break;
    case 'v':
        markov_variance_back();
        break;
    case 'c':
        /* markov compute em */
        /* markov free stoch */
        if (STOCH_HERE) {
            adj2_data_back();
            for (int32 i2 = 0; i2 < (NEQ + 1); i2++) {
                free(my_mean[i2]);
                free(my_variance[i2]);
            }
            STOCH_HERE = 0;
        }
        STOCH_FLAG = 1;
        integrate_do_range(&MyData[0], 0);
        init_conds_redraw_ics();

        STOCH_FLAG = 0;
        break;
    case 'h':
        histogram_compute();
        break;
    case 'o':
        histogram_back();
        break;
    case 'f':
        histogram_compute_fourier();
        break;
    case 'p':
        histogram_compute_power();
        break;
    case 'i':
        do_fit_test();
        init_conds_redraw_params();
        init_conds_redraw_ics();
        break;
    case 's':
        histogram_column_mean();
        break;
    case 'l':
        adj2_do_liapunov();
        break;
    case 'a':
        histogram_compute_stacor();
        break;
    case 'x':
        histogram_compute_correl();
        break;
    case 'e':
        histogram_compute_sd();
        break;
    case '2':
        histogram_new_2d();
        break;
    default:
        break;
    }
    return;
}

void
markov_mean_back(void) {
    if (STOCH_HERE) {
        set_browser_data(my_mean, 1);
        /*    my_browser.data=my_mean;
              my_browser.col0=1; */
        refresh_browser(stoch_len);
        storind = stoch_len;
    }
    return;
}

void
markov_variance_back(void) {
    if (STOCH_HERE) {
        set_browser_data(my_variance, 1);
        /*    my_browser.data=my_variance;
              my_browser.col0=1; */
        refresh_browser(stoch_len);
        storind = stoch_len;
    }
    return;
}



void
markov_init_stoch(int32 len) {
    int32 i;
    int32 j;
    N_TRIALS = 0;
    stoch_len = len;
    for (i = 0; i < (NEQ + 1); i++) {
        my_mean[i] = xmalloc(sizeof(*(my_mean[i]))*(usize)stoch_len);
        my_variance[i] = xmalloc(sizeof(*(my_variance[i]))*(usize)stoch_len);
        for (j = 0; j < stoch_len; j++) {
            my_mean[i][j] = 0.0;
            my_variance[i][j] = 0.0;
        }
    }
    for (j = 0; j < stoch_len; j++) {
        my_mean[0][j] = storage[0][j];
        my_variance[0][j] = storage[0][j];
    }
    STOCH_HERE = 1;
    return;
}

void
markov_append_stoch(int32 first, int32 length) {
    int32 i;
    int32 j;
    double z;
    if (first == 0)
        markov_init_stoch(length);
    if (length != stoch_len || !STOCH_HERE)
        return;
    for (i = 0; i < stoch_len; i++) {
        for (j = 1; j <= NEQ; j++) {
            z = storage[j][i];
            my_mean[j][i] = my_mean[j][i] + z;
            my_variance[j][i] = my_variance[j][i] + z*z;
        }
    }
    N_TRIALS++;
    return;
}

void
markov_do_stats(int32 ierr) {
    int32 i;
    int32 j;
    double ninv, mean;
    /*  STOCH_FLAG=0; */
    if (ierr != -1 && N_TRIALS > 0) {
        ninv = 1. / (double)(N_TRIALS);
        for (i = 0; i < stoch_len; i++) {
            for (j = 1; j <= NEQ; j++) {
                mean = my_mean[j][i]*ninv;
                my_mean[j][i] = mean;
                my_variance[j][i] = (my_variance[j][i]*ninv - mean*mean);
            }
        }
    }
    return;
}

double
markov_gammln(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146,     -86.50532032941677,
                            24.01409824083091,     -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5};
    int32 j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5)*log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
        ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005*ser / x);
}

double
markov_poidev(double xm) {
    static double sq, alxm, g, oldm = (-1.0);

    double em, t, y;

    if (xm < 12.0) {
        if (xm != oldm) {
            oldm = xm;
            g = exp(-xm);
        }
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= markov_ndrand48();
        } while (t > g);
    } else {
        if (xm != oldm) {
            oldm = xm;
            sq = sqrt(2.0*xm);
            alxm = log(xm);
            g = xm*alxm - markov_gammln(xm + 1.0);
        }
        do {
            do {
                y = tan(PI*markov_ndrand48());
                em = sq*y + xm;
            } while (em < 0.0);
            em = floor(em);
            t = 0.9*(1.0 + y*y) *
                exp(em*alxm - markov_gammln(em + 1.0) - g);
        } while (markov_ndrand48() > t);
    }
    return em;
}

double
markov_ndrand48(void) {
    return markov_ran1(&myrandomseed);
}

void
markov_nsrand48(int32 seed) {
    myrandomseed = -seed;
    return;
}

double
markov_ran1(long *idum) {
    int32 j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA*(*idum - k*IQ) - IR*k;
            if (*idum < 0)
                *idum += IM;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA*(*idum - k*IQ) - IR*k;
    if (*idum < 0)
        *idum += IM;
    j = (int32)(iy / NDIV);
    iy = iv[j];
    iv[j] = *idum;
    if ((temp = AM*(double)iy) > RNMX)
        return RNMX;
    else
        return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
