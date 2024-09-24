#include "functions.h"
#include "parserslow.h"
#include <strings.h>

#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>

#include "struct.h"
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

/*   This is numerics.c
 *   The input is primitive and eventually, I want to make it so
        that it uses nice windows for input.
        For now, I just will let it remain command driven
*/

extern int32 (*solver)(double *y, double *tim, double dt, int32 nt, int32 neq,
                       int32 *istart, double *work);
double *fft_data;
double *hist_data;

int32 cv_bandflag = 0;
int32 cv_bandupper = 1;
int32 cv_bandlower = 1;

/*   This is the input for the various functions */

/*   I will need access to storage  */

extern double **storage;

extern int32 NODE, NEQ; /* as well as the number of odes etc  */

static void get_method(void);
static void ruelle(void);
static void check_pos(int32 *j);

void
numerics_chk_volterra(void) {
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
numerics_quick_num(int32 com) {
    char key[] = "tsrdnviobec";
    if (com >= 0 && com < 11)
        get_num_par(key[com]);
    return;
}

void
numerics_set_total(double total) {
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
        menudrive_make_adj();
        break;
    case 't':
        /* total */
        ggets_new_float("total :", &TEND);
        FOREVER = 0;
        if (TEND < 0) {
            FOREVER = 1;
            TEND = -TEND;
        }
        break;
    case 's':
        /* start */
        ggets_new_float("start time :", &T0);
        break;
    case 'r':
        /* transient */
        ggets_new_float("transient :", &TRANS);
        break;
    case 'd':
        /* DT */
        temp = DELTA_T;
        ggets_new_float("Delta t :", &DELTA_T);
        if (DELTA_T == 0.0)
            DELTA_T = temp;
        if (DELAY > 0.0) {
            delay_handle_free_delay();
            if (delay_handle_alloc_delay(DELAY)) {
                INFLAG = 0; /*  Make sure no last ics allowed */
            }
        } else
            delay_handle_free_delay();
        if (NKernel > 0) {
            INFLAG = 0;
            MyStart = 1;
            volterra_alloc_kernels(1);
        }
        break;
    case 'n':
        /* ncline */
        ggets_new_int("ncline mesh :", &NMESH);
        check_pos(&NMESH);
        break;
    case 'v':
        ggets_new_int("Maximum iterates :", &BVP_MAXIT);
        check_pos(&BVP_MAXIT);
        ggets_new_float("Tolerance :", &BVP_TOL);
        ggets_new_float("Epsilon :", &BVP_EPS);
        pp_shoot_reset_bvp();
        break;
    case 'i':
        /* sing pt */
        ggets_new_int("Maximum iterates :", &EVEC_ITER);
        check_pos(&EVEC_ITER);
        ggets_new_float("Newton tolerance :", &EVEC_ERR);
        ggets_new_float("Jacobian epsilon :", &NEWT_ERR);
        if (NFlags > 0)
            ggets_new_float("SMIN :", &STOL);

        break;
    case 'o':
        /* noutput */
        ggets_new_int("n_out :", &NJMP);
        check_pos(&NJMP);

        break;
    case 'b':
        /* bounds */
        ggets_new_float("Bounds :", &BOUND);
        BOUND = fabs(BOUND);
        break;
    case 'm':
        /* method */
        get_method();
        if (METHOD == VOLTERRA && NKernel == 0) {
            ggets_err_msg("Volterra only for integral eqns");
            METHOD = 4;
        }
        if (NKernel > 0)
            METHOD = VOLTERRA;
        if (METHOD == GEAR || METHOD == RKQS || METHOD == STIFF) {
            ggets_new_float("Tolerance :", &TOLER);
            ggets_new_float("minimum step :", &HMIN);
            ggets_new_float("maximum step :", &HMAX);
        }
        if (METHOD == CVODE || METHOD == DP5 || METHOD == DP83 ||
            METHOD == RB23) {
            ggets_new_float("Relative tol:", &TOLER);
            ggets_new_float("Abs. Toler:", &ATOLER);
        }
        if (METHOD == BACKEUL || METHOD == VOLTERRA) {
            ggets_new_float("Tolerance :", &euler_tol);
            ggets_new_int("MaxIter :", &euler_max_iter);
        }
        if (METHOD == VOLTERRA) {
            tmp = MaxPoints;
            ggets_new_int("MaxPoints:", &tmp);
            ggets_new_int("AutoEval(1=yes) :", &AutoEvaluate);
            volterra_allocate(tmp, 1);
        }

        if (METHOD == CVODE || METHOD == RB23) {
            ggets_new_int("Banded system(0/1)?", &cv_bandflag);
            if (cv_bandflag == 1) {
                ggets_new_int("Lower band:", &cv_bandlower);
                ggets_new_int("Upper band:", &cv_bandupper);
            }
        }
        if (METHOD == SYMPLECT) {
            if ((NODE % 2) != 0) {
                ggets_err_msg("Symplectic is only for even dimensions");
                METHOD = 4;
            }
        }
        break;
    case 'e':
        /* delay */
        if (NDELAYS == 0)
            break;
        ggets_new_float("Maximal delay :", &DELAY);
        ggets_new_float("real guess :", &AlphaMax);
        ggets_new_float("imag guess :", &OmegaMax);
        ggets_new_int("DelayGrid :", &DelayGrid);
        if (DELAY > 0.0) {
            delay_handle_free_delay();
            if (delay_handle_alloc_delay(DELAY)) {
                INFLAG = 0; /*  Make sure no last ics allowed */
            }
        } else
            delay_handle_free_delay();
        break;
    case 'c':
        /* color */
        if (COLOR == 0)
            break;
        menudrive_set_col_par();
        break;
    case 'h':
        menudrive_do_stochast();
        break;
    case 'f':
        /* FFT */
        break;
    case 'p':
        /*Poincare map */
        menudrive_get_pmap_pars();
        break;
    case 'u':
        /* ruelle */
        ruelle();
        break;
    case 'k':
        /*lookup table */
        menudrive_new_lookup();
        break;
    case 27:
        numerics_do_meth();
        TEND = fabs(TEND);
        storage_alloc_meth();
        menu_help();
        break;
    default:
        break;
    }
    return;
}

void
numerics_chk_delay(void) {
    if (DELAY > 0.0) {
        delay_handle_free_delay();
        if (delay_handle_alloc_delay(DELAY)) {
            INFLAG = 0; /*  Make sure no last ics allowed */
        }
    } else
        delay_handle_free_delay();
    return;
}

void
numerics_set_delay(void) {
    if (NDELAYS == 0)
        return;
    if (DELAY > 0.0) {
        delay_handle_free_delay();
        if (delay_handle_alloc_delay(DELAY)) {
            INFLAG = 0;
        }
    }
    return;
}

void
ruelle(void) {
    ggets_new_int("x-axis shift ", &(MyGraph->xshft));
    ggets_new_int("y-axis shift ", &(MyGraph->yshft));
    ggets_new_int("z-axis shift", &(MyGraph->zshft));
    if (MyGraph->xshft < 0)
        MyGraph->xshft = 0;
    if (MyGraph->yshft < 0)
        MyGraph->yshft = 0;
    if (MyGraph->zshft < 0)
        MyGraph->zshft = 0;
    return;
}

void
numerics_compute_one_period(double period, double *x, char *name) {
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
        write_my_browser_data(fp);
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
        write_my_browser_data(fp);
        fclose(fp);
        adj2_data_back();
    }
    adj2_new_h_fun(1);
    snprintf(filename, sizeof(filename), "hfun.%s.dat", name);
    fp = fopen(filename, "w");
    if (fp != NULL) {
        write_my_browser_data(fp);
        fclose(fp);
        adj2_data_back();
    }

    reset_browser();

    TRANS = ot;
    POIMAP = opm;
    TEND = ote;
    return;
}

void
numerics_get_pmap_pars_com(int32 l) {
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

    graf_par_ind_to_sym(i1, n1);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%.16g", POIPLN);
    snprintf(values[2], sizeof(values[2]), "%d", POISGN);
    snprintf(values[3], sizeof(values[3]), "%s", yn[SOS]);
    status = do_string_box(4, 4, 1, "Poincare map", n, values, 45);
    if (status != 0) {
        browse_find_variable(values[0], &i1);
        if (i1 < 0) {
            POIMAP = 0;
            ggets_err_msg("No such section");
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
numerics_user_set_color_par(int32 flag, char *via, double lo, double hi) {
    int32 ivar;
    MyGraph->min_scale = lo;
    if (hi > lo)
        MyGraph->color_scale = (hi - lo);
    else
        MyGraph->color_scale = 1;

    if (strncasecmp("speed", via, 5) == 0) {
        MyGraph->ColorFlag = 1;
    } else {
        browse_find_variable(via, &ivar);
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
numerics_set_col_par_com(int32 i) {
    int32 j;
    int32 ivar;
    double temp[2];
    double maxder = 0.0, minder = 0.0, sum = 0.0;
    char ch, name[20];
    MyGraph->ColorFlag = i;
    if (MyGraph->ColorFlag == 0) {
        /* set color to black/white */
        return;
    }
    if (MyGraph->ColorFlag == 2) {
        graf_par_ind_to_sym(MyGraph->ColorValue, name);
        ggets_new_string("Color via:", name);
        browse_find_variable(name, &ivar);

        if (ivar >= 0)
            MyGraph->ColorValue = ivar;
        else {

            ggets_err_msg("No such quantity!");
            MyGraph->ColorFlag = 0;
            return;
        }
    }

    /*   This will be uncommented    ..... */
    ch = (char)menudrive_two_choice("(O)ptimize", "(C)hoose", "Color", "oc");

    if (ch == 'c') {
        temp[0] = MyGraph->min_scale;
        temp[1] = MyGraph->min_scale + MyGraph->color_scale;
        ggets_new_float("Min :", &temp[0]);
        ggets_new_float("Max :", &temp[1]);
        if (temp[1] > temp[0] &&
            ((MyGraph->ColorFlag == 2) ||
             (MyGraph->ColorFlag == 1 && temp[0] >= 0.0))) {
            MyGraph->min_scale = temp[0];
            MyGraph->color_scale = (temp[1] - temp[0]);
        } else {
            ggets_err_msg("Min>=Max or Min<0 error");
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
        graf_par_get_max(MyGraph->ColorValue, &temp[0], &temp[1]);
        MyGraph->min_scale = temp[0];
        MyGraph->color_scale = (temp[1] - temp[0]);
        if (MyGraph->color_scale == 0.0)
            MyGraph->color_scale = 1.0;
    }
    return;
}

void
numerics_do_meth(void) {
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
