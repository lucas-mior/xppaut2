#include "functions.h"
#include <strings.h>

#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>

extern Window main_win, info_pop;
extern Display *display;
extern int32 DCURY, NDELAYS;
extern int32 RandSeed;
#include "struct.h"
extern GRAPH *MyGraph;
#define VOLTERRA 6
#define BACKEUL 7
#define RKQS 8
#define STIFF 9
#define CVODE 10
#define GEAR 5
#define DP5 11
#define DP83 12
#define RB23 13
#define SYMPLECT 14

extern int32 NKernel, MyStart, MaxPoints;
extern int32 NFlags;
extern double STOL;
extern double MyTime;
extern char *info_message, *meth_hint[];
extern int32 DelayGrid;
extern double OmegaMax, AlphaMax;

/*   This is numerics.c
 *   The input is primitive and eventually, I want to make it so
        that it uses nice windows for input.
        For now, I just will let it remain command driven
*/

typedef struct {
    double tmod;
    int32 maxvar, sos, type, sign;
    char section[256];
    int32 formula[256];
} POINCARE_MAP;

POINCARE_MAP my_pmap;

extern int32 (*solver)(double *y, double *tim, double dt, int32 nt, int32 neq,
                       int32 *istart, double *work);
extern double DELTA_T, TEND, T0, TRANS, NULL_ERR, EVEC_ERR, NEWT_ERR;
extern double BOUND, DELAY, TOLER, ATOLER, HMIN, HMAX;
double *fft_data, *hist_data;
extern double POIPLN;

extern double BVP_TOL, BVP_EPS;

extern int32 NMESH, NJMP, METHOD, NC_ITER;
extern int32 EVEC_ITER;
extern int32 BVP_MAXIT, BVP_NL, BVP_NR;

extern int32 POIMAP, POIVAR, POISGN, SOS;

extern int32 HIST, HVAR, hist_ind, FOREVER, INFLAG;
extern int32 MaxEulIter;
extern double EulTol;

extern int32 AutoEvaluate;

int32 cv_bandflag = 0, cv_bandupper = 1, cv_bandlower = 1;
extern Window command_pop;

/*   This is the input for the various functions */

/*   I will need access to storage  */

extern double **storage;
extern int32 storind;

extern int32 NODE, NEQ; /* as well as the number of odes etc  */

void
chk_volterra(void) {
    if (NKernel > 0)
        METHOD = VOLTERRA;
    return;
}

void
check_pos(int32 *j) {
    if (*j <= 0)
        *j = 1;
    return;
}

void
quick_num(int32 com) {
    char key[] = "tsrdnviobec";
    if (com >= 0 && com < 11)
        get_num_par(key[com]);
    return;
}

void
set_total(double total) {
    int32 n;
    n = (int32)(total / fabs(DELTA_T)) + 1;
    TEND = n*fabs(DELTA_T);
    return;
}

void
get_num_par(int32 ch)

{
    double temp;
    int32 tmp;
    switch (ch) {
    case 'a':
        make_adj();
        break;
    case 't':
        /* total */
        new_float("total :", &TEND);
        FOREVER = 0;
        if (TEND < 0) {
            FOREVER = 1;
            TEND = -TEND;
        }
        break;
    case 's':
        /* start */
        new_float("start time :", &T0);
        break;
    case 'r':
        /* transient */
        new_float("transient :", &TRANS);
        break;
    case 'd':
        /* DT */
        temp = DELTA_T;
        new_float("Delta t :", &DELTA_T);
        if (DELTA_T == 0.0)
            DELTA_T = temp;
        if (DELAY > 0.0) {
            free_delay();
            if (alloc_delay(DELAY)) {
                INFLAG = 0; /*  Make sure no last ics allowed */
            }
        } else
            free_delay();
        if (NKernel > 0) {
            INFLAG = 0;
            MyStart = 1;
            alloc_kernels(1);
        }
        break;
    case 'n':
        /* ncline */
        new_int("ncline mesh :", &NMESH);
        check_pos(&NMESH);
        break;
    case 'v':
        new_int("Maximum iterates :", &BVP_MAXIT);
        check_pos(&BVP_MAXIT);
        new_float("Tolerance :", &BVP_TOL);
        new_float("Epsilon :", &BVP_EPS);
        reset_bvp();
        break;
    case 'i':
        /* sing pt */
        new_int("Maximum iterates :", &EVEC_ITER);
        check_pos(&EVEC_ITER);
        new_float("Newton tolerance :", &EVEC_ERR);
        new_float("Jacobian epsilon :", &NEWT_ERR);
        if (NFlags > 0)
            new_float("SMIN :", &STOL);

        break;
    case 'o':
        /* noutput */
        new_int("n_out :", &NJMP);
        check_pos(&NJMP);

        break;
    case 'b':
        /* bounds */
        new_float("Bounds :", &BOUND);
        BOUND = fabs(BOUND);
        break;
    case 'm':
        /* method */
        get_method();
        if (METHOD == VOLTERRA && NKernel == 0) {
            err_msg("Volterra only for integral eqns");
            METHOD = 4;
        }
        if (NKernel > 0)
            METHOD = VOLTERRA;
        if (METHOD == GEAR || METHOD == RKQS || METHOD == STIFF) {
            new_float("Tolerance :", &TOLER);
            new_float("minimum step :", &HMIN);
            new_float("maximum step :", &HMAX);
        }
        if (METHOD == CVODE || METHOD == DP5 || METHOD == DP83 ||
            METHOD == RB23) {
            new_float("Relative tol:", &TOLER);
            new_float("Abs. Toler:", &ATOLER);
        }
        if (METHOD == BACKEUL || METHOD == VOLTERRA) {
            new_float("Tolerance :", &EulTol);
            new_int("MaxIter :", &MaxEulIter);
        }
        if (METHOD == VOLTERRA) {
            tmp = MaxPoints;
            new_int("MaxPoints:", &tmp);
            new_int("AutoEval(1=yes) :", &AutoEvaluate);
            allocate_volterra(tmp, 1);
        }

        if (METHOD == CVODE || METHOD == RB23) {
            new_int("Banded system(0/1)?", &cv_bandflag);
            if (cv_bandflag == 1) {
                new_int("Lower band:", &cv_bandlower);
                new_int("Upper band:", &cv_bandupper);
            }
        }
        if (METHOD == SYMPLECT) {
            if ((NODE % 2) != 0) {
                err_msg("Symplectic is only for even dimensions");
                METHOD = 4;
            }
        }
        break;
    case 'e':
        /* delay */
        if (NDELAYS == 0)
            break;
        new_float("Maximal delay :", &DELAY);
        new_float("real guess :", &AlphaMax);
        new_float("imag guess :", &OmegaMax);
        new_int("DelayGrid :", &DelayGrid);
        if (DELAY > 0.0) {
            free_delay();
            if (alloc_delay(DELAY)) {
                INFLAG = 0; /*  Make sure no last ics allowed */
            }
        } else
            free_delay();
        break;
    case 'c':
        /* color */
        if (COLOR == 0)
            break;
        set_col_par();
        break;
    case 'h':
        do_stochast();
        break;
    case 'f':
        /* FFT */
        break;
    case 'p':
        /*Poincare map */
        get_pmap_pars();
        break;
    case 'u':
        /* ruelle */
        ruelle();
        break;
    case 'k':
        /*lookup table */
        new_lookup();
        break;
    case 27:
        do_meth();
        TEND = fabs(TEND);
        alloc_meth();
        help();
        break;
    default:
        break;
    }
    return;
}

void
chk_delay(void) {
    if (DELAY > 0.0) {
        free_delay();
        if (alloc_delay(DELAY)) {
            INFLAG = 0; /*  Make sure no last ics allowed */
        }
    } else
        free_delay();
    return;
}

void
set_delay(void) {
    if (NDELAYS == 0)
        return;
    if (DELAY > 0.0) {
        free_delay();
        if (alloc_delay(DELAY)) {
            INFLAG = 0;
        }
    }
    return;
}

void
ruelle(void) {
    new_int("x-axis shift ", &(MyGraph->xshft));
    new_int("y-axis shift ", &(MyGraph->yshft));
    new_int("z-axis shift", &(MyGraph->zshft));
    if (MyGraph->xshft < 0)
        MyGraph->xshft = 0;
    if (MyGraph->yshft < 0)
        MyGraph->yshft = 0;
    if (MyGraph->zshft < 0)
        MyGraph->zshft = 0;
    return;
}

void
init_numerics(void)
/*    these are the default values of the numerical parameters   */
{
    DELTA_T = .05;
    TEND = 20.0;
    T0 = 0.0;
    TRANS = 0.0;
    NULL_ERR = .001;
    EVEC_ERR = .001;
    NEWT_ERR = .001;
    BOUND = 100.0;
    DELAY = 0.0;
    TOLER = .00001;
    HMIN = .001;
    HMAX = 1.0;

    POIPLN = 0.0;
    NMESH = 50;
    NJMP = 1;
    METHOD = 4;
    NC_ITER = 100;
    EVEC_ITER = 100;

    /* new improved poincare map */

    my_pmap.maxvar = 1;
    my_pmap.type = 0;
    my_pmap.sos = 0;
    my_pmap.sign = 1;
    my_pmap.tmod = 8.*atan(1.0);
    snprintf(my_pmap.section, sizeof(my_pmap.section), " ");

    POIMAP = 0;
    POIVAR = 1;
    POISGN = 1;
    SOS = 0;
    return;
}

void
meth_dialog(void) {
    /*static char *n[]={"*6Method","Abs tol","Rel Tol","DtMin","DtMax",
                      "Banded(y/n)","UpperBand","LowerBand"};*/
    char values[8][MAX_LEN_SBOX];
    snprintf(values[0], sizeof(values[0]), "%d", METHOD);
    snprintf(values[1], sizeof(values[1]), "%g", ATOLER);
    snprintf(values[2], sizeof(values[2]), "%g", TOLER);
    return;
}

void
compute_one_period(double period, double *x, char *name) {
    int32 opm = POIMAP;
    char filename[256];
    double ot = TRANS, ote = TEND;
    FILE *fp;
    TRANS = 0;
    T0 = 0;
    MyTime = 0;
    TEND = period;
    POIMAP = 0; /* turn off poincare map */
    reset_browser();

    usual_integrate_stuff(x);
    snprintf(filename, sizeof(filename), "orbit.%s.dat", name);
    fp = fopen(filename, "w");
    if (fp != NULL) {
        write_mybrowser_data(fp);
        fclose(fp);
    } else {
        TRANS = ot;
        POIMAP = opm;
        TEND = ote;

        return;
    }
    adj2_new_adjoint();
    snprintf(filename, sizeof(filename), "adjoint.%s.dat", name);
    fp = fopen(filename, "w");
    if (fp != NULL) {
        write_mybrowser_data(fp);
        fclose(fp);
        adj_data_back();
    }
    adj2_new_h_fun(1);
    snprintf(filename, sizeof(filename), "hfun.%s.dat", name);
    fp = fopen(filename, "w");
    if (fp != NULL) {
        write_mybrowser_data(fp);
        fclose(fp);
        adj_data_back();
    }

    reset_browser();

    TRANS = ot;
    POIMAP = opm;
    TEND = ote;
    return;
}

void
get_pmap_pars_com(int32 l) {
    static char mkey[] = "nsmp";
    char ch;
    static char *n[] = {"*0Variable", "Section", "Direction (+1,-1,0)",
                        "Stop on sect(y/n)"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    static char *yn[] = {"N", "Y"};
    int32 status;
    char n1[15];
    int32 i1 = POIVAR;

    ch = mkey[l];

    POIMAP = 0;
    if (ch == 's')
        POIMAP = 1;
    if (ch == 'm')
        POIMAP = 2;
    if (ch == 'p')
        POIMAP = 3;

    if (POIMAP == 0)
        return;

    ind_to_sym(i1, n1);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%.16g", POIPLN);
    snprintf(values[2], sizeof(values[2]), "%d", POISGN);
    snprintf(values[3], sizeof(values[3]), "%s", yn[SOS]);
    status = do_string_box(4, 4, 1, "Poincare map", n, values, 45);
    if (status != 0) {
        find_variable(values[0], &i1);
        if (i1 < 0) {
            POIMAP = 0;
            err_msg("No such section");
            return;
        }
        POIVAR = i1;
        POISGN = atoi(values[2]);
        if (values[3][0] == 'Y' || values[3][0] == 'y')
            SOS = 1;
        else
            SOS = 0;
        POIPLN = atof(values[1]);
    }
    return;
}

void
get_method(void) {
    char ch;
    int32 i;
    int32 nmeth;

    Window temp = main_win;
    static char *n[] = {"(D)iscrete",    "(E)uler",   "(M)od. Euler",
                        "(R)unge-Kutta", "(A)dams",   "(G)ear",
                        "(V)olterra",    "(B)ackEul", "(Q)ualst.RK4",
                        "(S)tiff",       "(C)Vode",   "DoPri(5)",
                        "DoPri(8)3",     "Rosen(2)3", "sYmplectic"};
    static char key[] = "demragvbqsc582y";

#ifdef CVODE_YES
    nmeth = 15;
#else
    nmeth = 15;
#endif
    ch = (char)pop_up_list(&temp, "Method", n, key, nmeth, 15, METHOD, 10,
                           DCURY + 8, meth_hint, info_pop, info_message);
    for (i = 0; i < nmeth; i++)
        if (ch == key[i])
            METHOD = i;
    if (i > (nmeth - 1))
        i = nmeth - 1;
    return;
}

void
user_set_color_par(int32 flag, char *via, double lo, double hi) {
    int32 ivar;
    MyGraph->min_scale = lo;
    if (hi > lo)
        MyGraph->color_scale = (hi - lo);
    else
        MyGraph->color_scale = 1;

    if (strncasecmp("speed", via, 5) == 0) {
        MyGraph->ColorFlag = 1;
    } else {
        find_variable(via, &ivar);
        if (ivar >= 0) {
            MyGraph->ColorValue = ivar;
            MyGraph->ColorFlag = 2;
        } else {
            MyGraph->ColorFlag = 0; /* no valid colorizing */
        }
    }
    if (flag == 0) { /* force overwrite  */
        MyGraph->ColorFlag = 0;
    }
    return;
}

void
set_col_par_com(int32 i) {
    int32 j, ivar;
    double temp[2];
    double maxder = 0.0, minder = 0.0, sum = 0.0;
    char ch, name[20];
    MyGraph->ColorFlag = i;
    if (MyGraph->ColorFlag == 0) {
        /* set color to black/white */
        return;
    }
    if (MyGraph->ColorFlag == 2) {
        ind_to_sym(MyGraph->ColorValue, name);
        new_string("Color via:", name);
        find_variable(name, &ivar);

        if (ivar >= 0)
            MyGraph->ColorValue = ivar;
        else {

            err_msg("No such quantity!");
            MyGraph->ColorFlag = 0;
            return;
        }
    }

    /*   This will be uncommented    ..... */
    ch = (char)TwoChoice("(O)ptimize", "(C)hoose", "Color", "oc");

    if (ch == 'c') {
        temp[0] = MyGraph->min_scale;
        temp[1] = MyGraph->min_scale + MyGraph->color_scale;
        new_float("Min :", &temp[0]);
        new_float("Max :", &temp[1]);
        if (temp[1] > temp[0] &&
            ((MyGraph->ColorFlag == 2) ||
             (MyGraph->ColorFlag == 1 && temp[0] >= 0.0))) {
            MyGraph->min_scale = temp[0];
            MyGraph->color_scale = (temp[1] - temp[0]);
        } else {
            err_msg("Min>=Max or Min<0 error");
        }
        return;
    }
    if (MyGraph->ColorFlag == 1) {
        if (storind < 2)
            return;
        maxder = 0.0;
        minder = 1.e20;
        for (i = 1; i < my_browser.maxrow; i++) {
            sum = 0.0;
            for (j = 0; j < NODE; j++)
                sum += (double)fabs((double)(my_browser.data[1 + j][i] -
                                             my_browser.data[1 + j][i - 1]));
            if (sum < minder)
                minder = sum;
            if (sum > maxder)
                maxder = sum;
        }
        if (minder >= 0.0 && maxder > minder) {
            MyGraph->color_scale = (maxder - minder) / (fabs(DELTA_T*NJMP));
            MyGraph->min_scale = minder / (fabs(DELTA_T*NJMP));
        }
    } else {
        get_max(MyGraph->ColorValue, &temp[0], &temp[1]);
        MyGraph->min_scale = temp[0];
        MyGraph->color_scale = (temp[1] - temp[0]);
        if (MyGraph->color_scale == 0.0)
            MyGraph->color_scale = 1.0;
    }
    return;
}

void
do_meth(void) {
    if (NKernel > 0)
        METHOD = VOLTERRA;
    switch (METHOD) {
    case 0:
        solver = discrete;
        DELTA_T = 1;
        break;
    case 1:
        solver = euler;
        break;
    case 2:
        solver = mod_euler;
        break;
    case 3:
        solver = rung_kut;
        break;
    case 4:
        solver = adams;
        break;
    case 5:
        NJMP = 1;
        break;
    case 6:
        solver = volterra;
        break;
    case SYMPLECT:
        solver = symplect3;
        break;
    case BACKEUL:
        solver = bak_euler;
        break;
    case RKQS:
    case STIFF:
    case CVODE:
    case DP5:
    case DP83:
    case RB23:
        NJMP = 1;
        break;
    default:
        solver = rung_kut;
    }
    return;
}
