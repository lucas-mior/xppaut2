#include <string.h>
#include <stdbool.h>
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "functions.h"
#include "auto_nox.h"
#include "autlim.h"
#include "parserslow.h"
#include "x_auto.h"

#include "autevd.h"
#include "integers.h"
#include "auto_c.h"
#include "autlib.h"

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define MAXLINELENGTH 100000
#define PARAM_BOX 1

#define RUBBOX 0

#define ESC 27

#define UPT 6
#define SPT 7

#define OPEN_3 1
#define NO_OPEN_3 0
#define OVERWRITE 0
#define APPEND 1

/* calculation types */

int32 TypeOfCalc = 0;
#define LPE2 1
#define LPP2 2
#define HB2 3
#define TR2 4
#define BR2 5
#define PD2 6
#define FP2 7
#define BV1 8
#define EQ1 9
#define PE1 10
#define DI1 11
#define HO2 12

#define HI_P 0 /* uhi vs par */
#define NR_P 1 /* norm vs par */
#define HL_P 2 /* Hi and Lo vs par  periodic only */
#define PE_P 3 /* period vs par   */
#define P_P 4  /* param vs param  */

#define FR_P 10 /* freq vs par   */
#define AV_P 11 /* ubar vs par */
#define SPER 3
#define UPER 4
#define CSEQ 1
#define CUEQ 2

#define DISCRETE 0

int32 SEc = 20;
int32 UEc = 0;
int32 SPc = 26;
int32 UPc = 28;

/* two parameter colors  need to do this
 * LP is 20 (red)
 * HB is 28 (blue)
 * TR is 26 (green)
 * PD is 24 (orange)
 * BR is 27 (turquoise)
 * FP is 25 (olive) */

static int32 LPP_color = 0;
static int32 LPE_color = 20;
static int32 HB_color = 28;
static int32 TR_color = 26;
static int32 PD_color = 23;
static int32 BR_color = 27;
static int32 FP_color = 25;

int32 RestartLabel = 0;

int32 auto_ntst = 15;
int32 auto_nmx = 200;
int32 auto_npr = 50;
int32 auto_ncol = 4;
double auto_ds = .02;
double auto_dsmax = .5;
double auto_dsmin = .001;

double auto_rl0 = 0.0;
double auto_rl1 = 2;
double auto_a0 = 0.0;
double auto_a1 = 1000.;
double auto_xmax = 2.5;
double auto_xmin = -.5;
double auto_ymax = 3.0;
double auto_ymin = -3.0;
double auto_epsl = 1e-4;
double auto_epsu = 1e-4;
double auto_epss = 1e-4;
int32 auto_var = 0;

int32 load_all_labeled_orbits = 0;

static int32 SuppressBP = 0;
Rotchk blrtn;

GrabPoint grabpt;

int32 AutoTwoParam = 0;
int32 NAutoPar = 8;
int32 Auto_index_to_array[8];
int32 AutoPar[8];

double outperiod[20];
int64 UzrPar[20];
int32 NAutoUzr;

static char this_auto_file[200];
char fort3[200];
char fort7[200];
char fort8[200];
char fort9[200];
static char TMPSWAP[200];

static double XfromAuto;
static double YfromAuto;
static double XfromAuto;
static double YfromAuto;
static int32 FromAutoFlag = 0;

int32 HomoFlag = 0;
double homo_l[100];
double homo_r[100];
static double HOMO_SHIFT = 0.0;

Bifurcation Auto;
AdvAuto aauto;

int32 NewPeriodFlag;

static AutoAX Old1p;
static AutoAX Old2p;

static void auto_nox_start_at_homoclinic(void);
static void auto_nox_2p_hopf(void);
static void auto_nox_load_orbitx(int32 ibr, int32 flag, int32 lab, double per);
static void auto_nox_pscolset2(int32 flag2);

/* color plot stuff */
void
auto_nox_colset(int32 type) {
    switch (type) {
    case CSEQ:
        auto_x11_col(SEc);
        break;
    case CUEQ:
        auto_x11_col(UEc);
        break;
    case SPER:
        auto_x11_col(SPc);
        break;
    case UPER:
        auto_x11_col(UPc);
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
auto_nox_pscolset2(int32 flag2) {
    switch (flag2) {
    case LPE2:
        graphics_set_linestyle(LPE_color - 19);
        break;
    case LPP2:
        graphics_set_linestyle(LPP_color);
        break;
    case HB2:
        graphics_set_linestyle(HB_color - 19);
        break;
    case TR2:
        graphics_set_linestyle(TR_color - 19);
        break;
    case BR2:
        graphics_set_linestyle(BR_color - 19);
        break;
    case PD2:
        graphics_set_linestyle(PD_color - 19);
        break;
    case FP2:
        graphics_set_linestyle(FP_color - 19);
        break;
    default:
        graphics_set_linestyle(0);
    }
    return;
}

void
auto_nox_colset2(int32 flag2) {
    auto_x11_line_width(2);
    switch (flag2) {
    case LPE2:
        auto_x11_col(LPE_color);
        break;
    case LPP2:
        auto_x11_col(LPP_color);
        break;
    case HB2:
        auto_x11_col(HB_color);
        break;
    case TR2:
        auto_x11_col(TR_color);
        break;
    case BR2:
        auto_x11_col(BR_color);
        break;
    case PD2:
        auto_x11_col(PD_color);
        break;
    case FP2:
        auto_x11_col(FP_color);
        break;
    default:
        auto_x11_col(0);
    }
    return;
}

void
auto_nox_store_point(double x, double y) {
    if (Auto.plot == P_P) {
        XfromAuto = x;
        YfromAuto = y;
        FromAutoFlag = 1;
    }
    return;
}

void
auto_nox_get_str(char *xlabel, char *ylabel) {
    sprintf(xlabel, "%s", upar_names[AutoPar[Auto.icp1]]);
    switch (Auto.plot) {
    case HI_P:
    case HL_P:
        sprintf(ylabel, "%s", uvar_names[Auto.var]);
        break;
    case NR_P:
        sprintf(ylabel, "Norm");
        break;
    case PE_P:
        sprintf(ylabel, "Period");
        break;
    case FR_P:
        sprintf(ylabel, "Frequency");
        break;
    case P_P:
        sprintf(ylabel, "%s", upar_names[AutoPar[Auto.icp2]]);
        break;
    case AV_P:
        sprintf(ylabel, "%s_bar", uvar_names[Auto.var]);
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
auto_nox_draw_ps_axes(void) {
    char sx[20], sy[20];
    graphics_set_scale(Auto.xmin, Auto.ymin, Auto.xmax, Auto.ymax);
    auto_nox_get_str(sx, sy);
    axes2_box(Auto.xmin, Auto.xmax, Auto.ymin, Auto.ymax, sx, sy, 0);
    return;
}

void
auto_nox_draw_svg_axes(void) {
    char sx[20], sy[20];
    graphics_set_scale(Auto.xmin, Auto.ymin, Auto.xmax, Auto.ymax);
    auto_nox_get_str(sx, sy);
    axes2_box(Auto.xmin, Auto.xmax, Auto.ymin, Auto.ymax, sx, sy, 0);
    return;
}

void
auto_nox_draw_bix_axes(void) {
    int32 x0 = Auto.x0;
    int32 y0 = Auto.y0;
    int32 ii;
    int32 i0;
    int32 x1 = x0 + Auto.wid, y1 = y0 + Auto.hgt;
    char junk[20], xlabel[20], ylabel[20];

    auto_x11_clear_plot();
    auto_x11_line(x0, y0, x1, y0);
    auto_x11_line(x1, y0, x1, y1);
    auto_x11_line(x1, y1, x0, y1);
    auto_x11_line(x0, y1, x0, y0);
    sprintf(junk, "%g", Auto.xmin);
    auto_x11_text(x0, y1 + DCURYs + 2, junk);
    sprintf(junk, "%g", Auto.xmax);
    ii = (int32)strlen(junk)*DCURXs;
    auto_x11_text(x1 - ii, y1 + DCURYs + 2, junk);
    sprintf(junk, "%g", Auto.ymin);
    ii = (int32)strlen(junk);
    i0 = 9 - ii;
    if (i0 < 0)
        i0 = 0;
    auto_x11_text(i0*DCURXs, y1, junk);
    sprintf(junk, "%g", Auto.ymax);
    ii = (int32)strlen(junk);
    i0 = 9 - ii;
    if (i0 < 0)
        i0 = 0;
    auto_x11_text(i0*DCURXs, y0 + DCURYs, junk);
    auto_nox_get_str(xlabel, ylabel);
    auto_x11_text((x0 + x1) / 2, y1 + DCURYs + 2, xlabel);
    auto_x11_text(10*DCURXs, DCURYs, ylabel);
    auto_x11_refresh_display();
    return;
}

int32
auto_nox_ix_val(double x) {
    double temp = (double)Auto.wid*(x - Auto.xmin) / (Auto.xmax - Auto.xmin);
    return (int32)temp + Auto.x0;
}

int32
auto_nox_iy_val(double y) {
    double temp = (double)Auto.hgt*(y - Auto.ymin) / (Auto.ymax - Auto.ymin);
    return Auto.hgt - (int32)temp + Auto.y0;
}

int32
auto_nox_check_bnds(int32 ix, int32 iy) {
    int32 x1 = Auto.x0, x2 = Auto.x0 + Auto.wid;
    int32 y1 = Auto.y0, y2 = Auto.y0 + Auto.hgt;
    if ((ix >= x1) && (ix < x2) && (iy >= y1) && (iy < y2))
        return 1;
    return 0;
}

void
auto_nox_renamef(char *old, char *new) {
    rename(old, new);
}

void
auto_nox_copyf(char *old, char *new) {
    FILE *fo;
    FILE *fn;
    int32 c;
    fo = fopen(old, "r");
    fn = fopen(new, "w");

    while ((c = getc(fo)) != EOF) {
        putc(c, fn);
    }
    fclose(fo);
    fclose(fn);
    return;
}

void
auto_nox_appendf(char *old, char *new) {
    FILE *fo;
    FILE *fn;
    FILE *ft;
    int32 c;

    fo = fopen(old, "r");
    fn = fopen(new, "r");
    if (fn == NULL) {
        fclose(fo);

        auto_nox_copyf(old, new);
        return;
    }
    ft = fopen(TMPSWAP, "w");
    if (ft == NULL) {
        printf("Can't open %s \n", TMPSWAP);
        return;
    }
    while ((c = getc(fo)) != EOF)
        putc(c, ft);
    fclose(fo);
    while ((c = getc(fn)) != EOF)
        putc(c, ft);
    fclose(fn);
    fclose(ft);
    auto_nox_copyf(TMPSWAP, new);
    auto_nox_deletef(TMPSWAP);
    return;
}

void
auto_nox_deletef(char *old) {
    remove(old);
}

void
auto_nox_close(int32 flg) {
    /* labels compatible with A2K  */
    char string[1000];
    if (flg == 0) { /*Overwrite*/
        sprintf(string, "%s.b", this_auto_file);
        auto_nox_renamef(fort7, string);
        sprintf(string, "%s.d", this_auto_file);
        auto_nox_renamef(fort9, string);

        sprintf(string, "%s.s", this_auto_file);
        auto_nox_renamef(fort8, string);
    } else { /*APPEND*/
        sprintf(string, "%s.b", this_auto_file);
        auto_nox_appendf(fort7, string);
        sprintf(string, "%s.d", this_auto_file);
        auto_nox_appendf(fort9, string);
        sprintf(string, "%s.s", this_auto_file);
        auto_nox_appendf(fort8, string);
    }

    auto_nox_deletef(fort8);

    fp8_is_open = 0;
    auto_nox_deletef(fort7);
    auto_nox_deletef(fort9);
    auto_nox_deletef(fort3);
    return;
}

void
auto_nox_open(int32 flg) {
    /* compatible with new auto */
    char *basec, *bname, *dirc, *dname;
    char *HOME;

    basec = strdup(this_file);
    dirc = strdup(this_file);
    bname = (char *)basename(basec);
    dname = (char *)dirname(dirc);

    HOME = getenv("HOME");
    if (HOME == NULL)
        HOME = dname;

    sprintf(this_auto_file, "%s/%s", HOME, bname);
    sprintf(fort3, "%s/%s", HOME, "fort.3");
    sprintf(fort7, "%s/%s", HOME, "fort.7");
    sprintf(fort8, "%s/%s", HOME, "fort.8");
    sprintf(fort9, "%s/%s", HOME, "fort.9");
    sprintf(TMPSWAP, "%s/%s", HOME, "__tmp__");

    if (flg == 1) {
        char string[sizeof(this_auto_file) + 2];
        snprintf(string, sizeof(string), "%s.s", this_auto_file);
        auto_nox_copyf(string, fort3);
    }
    return;
}

/* MAIN Running routine  Assumes that Auto structure is set up */

void
auto_nox_do(int32 iold, int32 isave) {
    auto_x11_redraw_menus();

    /* auto_nox_set_auto */
    /* this sets up all the continuation initialization
     * it is equivalent to reading in auto parameters
     * and running init in auto */

    /* Caution - need to include NICP here */
    NAutoUzr = Auto.nper;
    autevd_init_auto(NODE, Auto.nfpar, Auto.ips, Auto.irs, Auto.ilp, Auto.ntst,
                     Auto.isp, Auto.isw, Auto.nmx, Auto.npr, Auto.ds,
                     Auto.dsmin, Auto.dsmax, Auto.rl0, Auto.rl1, Auto.a0,
                     Auto.a1, Auto.icp1, Auto.icp2, Auto.icp3, Auto.icp4,
                     Auto.icp5, Auto.epsl, Auto.epsu, Auto.epss, Auto.ncol);

    auto_nox_open(iold); /* this copies the relevant files .s  to fort.3 */
    go_go_auto();        /* this complets the initialization and calls the
                             main routines
                         */
    /*     run_aut(Auto.nfpar,itp); THIS WILL CHANGE TO gogoauto stuff */
    auto_nox_close(isave); /* this copies fort.8 to the .s file and other
                          irrelevant stuff
                       */

    if (RestartLabel != 0) {
        printf("RestartLabel=%d itp=%d ips=%d nfpar=%d ilp=%d isw=%d isp=%d "
               "A2p=%d \n",
               RestartLabel, Auto.itp, Auto.ips, Auto.nfpar, Auto.ilp, Auto.isw,
               Auto.isp, AutoTwoParam);
        Auto.irs = RestartLabel;
        RestartLabel = 0;
        auto_nox_do(iold, isave);
    }
    ggets_ping();
    init_conds_redraw_params();
}

int32
auto_nox_name_to_index(char *s) {
    int32 i;
    int32 in;
    browse_find_variable(s, &in);
    if (in == 0)
        return 10;
    in = init_conds_find_user_name(PARAM_BOX, s);
    for (i = 0; i < NAutoPar; i++)
        if (AutoPar[i] == in)
            return i;
    return -1;
}

int32
auto_nox_par_to_name(int64 index, char *s) {
    if (index == 10) {
        sprintf(s, "T");
        return 1;
    }
    if (index < 0 || index > 8)
        return 0;
    sprintf(s, "%s", upar_names[AutoPar[index]]);
    return 1;
}

void
auto_nox_per_par(void) {
    static char *m[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    static char key[] = "0123456789";
    char values[LENGTH(m)][MAX_LEN_SBOX];
    char bob[100], *ptr;
    static char *n[] = {"Uzr1", "Uzr2", "Uzr3", "Uzr4", "Uzr5",
                        "Uzr6", "Uzr7", "Uzr8", "Uzr9"};
    int32 status, i, in;
    char ch;
    ch = (char)auto_x11_pop_up_list("Number", m, key, 10, 12, Auto.nper, 10, 10,
                                    no_hint, Auto.hinttxt);
    for (i = 0; i < 10; i++)
        if (ch == key[i])
            Auto.nper = i;
    NAutoUzr = Auto.nper;
    if (Auto.nper > 0) {
        for (i = 0; i < 9; i++) {
            auto_nox_par_to_name(Auto.uzrpar[i], bob);

            sprintf(values[i], "%s=%g", bob, Auto.period[i]);
        }
        status = pop_list_do_string_box(9, 5, 2, "AutoPer", n, values, 45);
        if (status != 0)
            for (i = 0; i < 9; i++) {
                ptr = form_ode_get_first(values[i], "=");
                in = auto_nox_name_to_index(ptr);
                if (in >= 0) {
                    Auto.uzrpar[i] = in;
                    ptr = form_ode_do_fit_get_next("@");
                    Auto.period[i] = atof(ptr);
                }
            }
    }
    for (i = 0; i < 9; i++) {
        outperiod[i] = Auto.period[i];
        UzrPar[i] = Auto.uzrpar[i];
    }
    return;
}

/* auto parameters are 1-8 (0-7) and since there are only 8, need to associate
 * them with real xpp parameters for which there may be many */
void
auto_nox_params(void) {
    static char *n[] = {"*2Par1", "*2Par2", "*2Par3", "*2Par4",
                        "*2Par5", "*2Par6", "*2Par7", "*2Par8"};
    int32 status;
    int32 in;
    char values[LENGTH(n)][MAX_LEN_SBOX];
    for (int32 i = 0; i < 8; i++) {
        if (i < NAutoPar)
            sprintf(values[i], "%s", upar_names[AutoPar[i]]);
        else
            values[i][0] = '\0'; /*sprintf(values[i],"");*/
    }
    status = pop_list_do_string_box(8, 8, 1, "Parameters", n, values, 38);
    if (status != 0) {
        for (int32 i = 0; i < 8; i++) {
            if (i < NAutoPar) {
                in = init_conds_find_user_name(PARAM_BOX, values[i]);
                if (in >= 0) {
                    AutoPar[i] = in;
                    in = get_param_index(values[i]);
                    Auto_index_to_array[i] = in;
                }
            }
        }
    }
    return;
}

void
auto_nox_num_par(void) {
    static char *n[] = {"Ntst",     "Nmax",     "NPr",   "Ds",      "Dsmin",
                        "Ncol",     "EPSL",     "Dsmax", "Par Min", "Par Max",
                        "Norm Min", "Norm Max", "EPSU",  "EPSS",    "IAD",
                        "MXBF",     "IID",      "ITMX",  "ITNW",    "NWTN",
                        "IADS",     "SuppBP"};
    int32 status;
    char values[LENGTH(n)][MAX_LEN_SBOX];
    sprintf(values[0], "%d", Auto.ntst);
    sprintf(values[1], "%d", Auto.nmx);
    sprintf(values[2], "%d", Auto.npr);
    sprintf(values[3], "%g", Auto.ds);
    sprintf(values[4], "%g", Auto.dsmin);
    sprintf(values[7], "%g", Auto.dsmax);
    sprintf(values[8], "%g", Auto.rl0);
    sprintf(values[9], "%g", Auto.rl1);
    sprintf(values[10], "%g", Auto.a0);
    sprintf(values[11], "%g", Auto.a1);
    sprintf(values[5], "%d", Auto.ncol);
    sprintf(values[6], "%g", Auto.epsl);
    sprintf(values[12], "%g", Auto.epsu);
    sprintf(values[13], "%g", Auto.epss);
    sprintf(values[14], "%d", aauto.iad);
    sprintf(values[15], "%d", aauto.mxbf);
    sprintf(values[16], "%d", aauto.iid);
    sprintf(values[17], "%d", aauto.itmx);
    sprintf(values[18], "%d", aauto.itnw);
    sprintf(values[19], "%d", aauto.nwtn);
    sprintf(values[20], "%d", aauto.iads);
    sprintf(values[21], "%d", SuppressBP);

    status = pop_list_do_string_box(22, 7, 4, "AutoNum", n, values, 25);
    if (status != 0) {
        Auto.ntst = atoi(values[0]);
        Auto.nmx = atoi(values[1]);
        Auto.npr = atoi(values[2]);
        Auto.ds = atof(values[3]);
        Auto.dsmin = atof(values[4]);
        Auto.dsmax = atof(values[7]);
        Auto.rl0 = atof(values[8]);
        Auto.rl1 = atof(values[9]);
        Auto.a0 = atof(values[10]);
        Auto.a1 = atof(values[11]);
        Auto.ncol = atoi(values[5]);
        Auto.epsl = atof(values[6]);
        Auto.epsu = atof(values[12]);
        Auto.epss = atof(values[13]);
        aauto.iad = atoi(values[14]);
        aauto.mxbf = atoi(values[15]);
        aauto.iid = atoi(values[16]);
        aauto.itmx = atoi(values[17]);
        aauto.itnw = atoi(values[18]);
        aauto.nwtn = atoi(values[19]);
        aauto.iads = atoi(values[20]);
        SuppressBP = atoi(values[21]);
    }
    return;
}

void
auto_nox_plot_par(void) {
    static char *m[] = {"Hi",         "Norm",      "hI-lo",      "Period",
                        "Two par",    "(Z)oom in", "Zoom (O)ut", "last 1 par",
                        "last 2 par", "Fit",       "fRequency",  "Average",
                        "Default",    "Scroll"};
    static char key[] = "hniptzo12frads";
    char ch;

    static char *n[] = {"*1Y-axis", "*2Main Parm", "*2Secnd Parm", "Xmin",
                        "Ymin",     "Xmax",        "Ymax"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    int32 i;
    int32 ii1, ii2, ji1, ji2;
    int32 i1 = Auto.var + 1;
    char n1[15];
    ch = (char)auto_x11_pop_up_list("Plot Type", m, key, 14, 10, Auto.plot, 10,
                                    50, aaxes_hint, Auto.hinttxt);
    if (ch == ESC)
        return;
    for (i = 0; i < 5; i++)
        if (ch == key[i])
            Auto.plot = i;
    if (ch == key[10])
        Auto.plot = 10;
    if (ch == key[11])
        Auto.plot = 11;
    if (ch == key[5]) {
        if (auto_x11_rubber(&ii1, &ji1, &ii2, &ji2, RUBBOX) != 0) {
            auto_nox_zoom_in(ii1, ji1, ii2, ji2);
            diagram_redraw();
        }
        return;
    }

    if (ch == key[6]) {
        if (auto_x11_rubber(&ii1, &ji1, &ii2, &ji2, RUBBOX) != 0) {
            auto_nox_zoom_out(ii1, ji1, ii2, ji2);

            diagram_redraw();
        }
        return;
    }

    if (ch == key[7]) {
        auto_nox_load_last_plot(1);
        auto_nox_draw_bix_axes();
        return;
    }

    if (ch == key[8]) {
        auto_nox_load_last_plot(2);
        auto_nox_draw_bix_axes();
        return;
    }

    if (ch == key[9]) {
        /* auto fit */
        double xlo = Auto.xmin;
        double xhi = Auto.xmax;
        double ylo = Auto.ymin;
        double yhi = Auto.ymax;
        diagram_bound(&xlo, &xhi, &ylo, &yhi);
        Auto.xmin = xlo;
        Auto.xmax = xhi;
        Auto.ymin = ylo;
        Auto.ymax = yhi;

        diagram_redraw();
        return;
    }
    if (ch == key[12]) {
        /* auto default */
        Auto.xmin = auto_xmin;
        Auto.xmax = auto_xmax;
        Auto.ymin = auto_ymin;
        Auto.ymax = auto_ymax;

        diagram_redraw();
        return;
    }
    if (ch == key[13]) {
        auto_x11_scroll();
        diagram_redraw();
        /* printf("I am done scrolling!!"); */
        return;
    }
    graf_par_ind_to_sym(i1, n1);
    sprintf(values[0], "%s", n1);
    sprintf(values[1], "%s", upar_names[AutoPar[Auto.icp1]]);
    sprintf(values[2], "%s", upar_names[AutoPar[Auto.icp2]]);
    sprintf(values[3], "%g", Auto.xmin);
    sprintf(values[4], "%g", Auto.ymin);
    sprintf(values[5], "%g", Auto.xmax);
    sprintf(values[6], "%g", Auto.ymax);
    status = pop_list_do_string_box(7, 7, 1, "AutoPlot", n, values, 31);
    if (status != 0) {
        /*  get variable names  */
        browse_find_variable(values[0], &i);
        if (i > 0)
            Auto.var = i - 1;
        /*  Now check the parameters  */
        i1 = init_conds_find_user_name(PARAM_BOX, values[1]);
        if (i1 >= 0) {
            for (i = 0; i < NAutoPar; i++) {
                if (i1 == AutoPar[i]) {
                    Auto.icp1 = i;
                }
            }
        }
        i1 = init_conds_find_user_name(PARAM_BOX, values[2]);
        if (i1 >= 0) {
            for (i = 0; i < NAutoPar; i++) {
                if (i1 == AutoPar[i]) {
                    Auto.icp2 = i;
                }
            }
        }

        Auto.xmin = atof(values[3]);
        Auto.ymin = atof(values[4]);
        Auto.xmax = atof(values[5]);
        Auto.ymax = atof(values[6]);
        auto_nox_draw_bix_axes();
        if (Auto.plot < 4)
            auto_nox_keep_last_plot(1);
        if (Auto.plot == 4)
            auto_nox_keep_last_plot(2);
    }
    return;
}

void
auto_nox_zoom_in(int32 i1, int32 j1, int32 i2, int32 j2) {
    double x1, y1, x2, y2;
    int32 temp;
    double dx = (Auto.xmax - Auto.xmin);
    double dy = (Auto.ymax - Auto.ymin);

    if (i1 > i2) {
        temp = i1;
        i1 = i2;
        i2 = temp;
    }
    if (j2 > j1) {
        temp = j1;
        j1 = j2;
        j2 = temp;
    }

    x1 = Auto.xmin + (double)(i1 - Auto.x0)*(dx) / (double)Auto.wid;
    x2 = Auto.xmin + (double)(i2 - Auto.x0)*(dx) / (double)Auto.wid;
    y1 =
        Auto.ymin + (double)(Auto.hgt + Auto.y0 - j1)*(dy) / (double)Auto.hgt;
    y2 =
        Auto.ymin + (double)(Auto.hgt + Auto.y0 - j2)*(dy) / (double)Auto.hgt;

    if ((i1 == i2) || (j1 == j2)) {
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        dx = dx / 2;
        dy = dy / 2;
        /*Shrink by thirds and center (track) about the point clicked*/
        Auto.xmin = x1 - dx / 2;
        Auto.xmax = x1 + dx / 2;
        Auto.ymin = y1 - dy / 2;
        Auto.ymax = y1 + dy / 2;
    } else {
        Auto.xmin = x1;
        Auto.ymin = y1;
        Auto.xmax = x2;
        Auto.ymax = y2;
    }
    return;
}

void
auto_nox_zoom_out(int32 i1, int32 j1, int32 i2, int32 j2) {
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    int32 temp;
    double dx = (Auto.xmax - Auto.xmin);
    double dy = (Auto.ymax - Auto.ymin);
    double a1, a2, b1, b2;

    if (i1 > i2) {
        temp = i1;
        i1 = i2;
        i2 = temp;
    }
    if (j2 > j1) {
        temp = j1;
        j1 = j2;
        j2 = temp;
    }
    a1 = (double)(i1 - Auto.x0) / (double)Auto.wid;
    a2 = (double)(i2 - Auto.x0) / (double)Auto.wid;
    b1 = (double)(Auto.hgt + Auto.y0 - j1) / (double)Auto.hgt;
    b2 = (double)(Auto.hgt + Auto.y0 - j2) / (double)Auto.hgt;

    if ((i1 == i2) || (j1 == j2)) {
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        dx = dx*2;
        dy = dy*2;
        /*Shrink by thirds and center (track) about the point clicked*/
        Auto.xmin = x1 - dx / 2;
        Auto.xmax = x1 + dx / 2;
        Auto.ymin = y1 - dy / 2;
        Auto.ymax = y1 + dy / 2;
    } else {
        x1 = (a1*Auto.xmax - a2*Auto.xmin) / (a1 - a2);
        x2 = (Auto.xmin - Auto.xmax + a1*Auto.xmax - a2*Auto.xmin) /
             (a1 - a2);
        y1 = (b1*Auto.ymax - b2*Auto.ymin) / (b1 - b2);
        y2 = (Auto.ymin - Auto.ymax + b1*Auto.ymax - b2*Auto.ymin) /
             (b1 - b2);
        Auto.xmin = x1;
        Auto.ymin = y1;
        Auto.xmax = x2;
        Auto.ymax = y2;
    }
    return;
}

void
auto_nox_xy_plot(double *x, double *y1, double *y2, double par1, double par2,
             double per, double *uhigh, double *ulow, double *ubar, double a) {
    switch (Auto.plot) {
    case HI_P:
        *x = par1;
        *y1 = uhigh[Auto.var];
        *y2 = *y1;
        break;
    case NR_P:
        *x = par1;
        *y1 = a;
        *y2 = *y1;
        break;
    case HL_P:
        *x = par1;
        *y1 = uhigh[Auto.var];
        *y2 = ulow[Auto.var];
        break;
    case AV_P:
        *x = par1;
        *y1 = ubar[Auto.var];
        *y2 = *y1;
        break;
    case PE_P:
        *x = par1;
        *y1 = per;
        *y2 = *y1;
        break;
    case FR_P:
        *x = par1;
        if (per > 0)
            *y1 = 1. / per;
        else
            *y1 = 0.0;
        *y2 = *y1;
        break;
    case P_P:
        *x = par1;
        *y1 = par2;
        *y2 = *y1;
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

int32
auto_nox_plot_point(int32 flag2, int32 icp1, int32 icp2) {
    int32 j = 1;
    if (icp1 != Auto.icp1)
        j = 0;
    if (flag2 > 0 && icp2 != Auto.icp2)
        j = 0;
    return j;
}

void
auto_nox_add_ps_point(double *par, double per, double *uhigh, double *ulow,
                  double *ubar, double a, int32 type, int32 flag, int32 icp1,
                  int32 icp2, int32 flag2) {
    double x, y1, y2, par1, par2 = 0;
    int32 type1 = type;
    par1 = par[icp1];
    if (icp2 < NAutoPar)
        par2 = par[icp2];
    auto_nox_xy_plot(&x, &y1, &y2, par1, par2, per, uhigh, ulow, ubar, a);
    if (flag == 0) {
        Auto.lastx = x;
        Auto.lasty = y1;
    }

    if (flag2 == 0 && Auto.plot == P_P)
        return;
    if (flag2 > 0 && Auto.plot != P_P)
        return;

    if ((flag2 > 0) && (Auto.plot == P_P))
        type1 = CSEQ;

    switch (type1) {
    case CSEQ:
        if (Auto.plot == PE_P || Auto.plot == FR_P)
            break;
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;

        if (PS_Color) {
            graphics_set_linestyle(1);
            if (flag2 > 0)
                auto_nox_pscolset2(flag2);
        } else
            graphics_set_linestyle(8);
        graphics_line_abs((double)x, (double)y1, (double)Auto.lastx,
                          (double)Auto.lasty);
        break;
    case CUEQ:
        if (Auto.plot == PE_P || Auto.plot == FR_P)
            break;
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        if (Auto.plot != P_P) {
            if (PS_Color)
                graphics_set_linestyle(0);
            else
                graphics_set_linestyle(4);
        } else {
            auto_nox_pscolset2(flag2);
        }
        graphics_line_abs((double)x, (double)y1, (double)Auto.lastx,
                          (double)Auto.lasty);
        break;
    case UPER:
        if (PS_Color)
            graphics_set_linestyle(9);
        else
            graphics_set_linestyle(0);
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        PointType = UPT;
        graphics_point_abs((double)x, (double)y1);
        graphics_point_abs((double)x, (double)y2);
        break;
    case SPER:
        if (PS_Color)
            graphics_set_linestyle(7);
        else
            graphics_set_linestyle(0);
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        PointType = SPT;
        graphics_point_abs((double)x, (double)y1);
        graphics_point_abs((double)x, (double)y2);
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }

    Auto.lastx = x;
    Auto.lasty = y1;
    return;
}

void
auto_nox_line(double x1i, double y1i, double x2i, double y2i) {
    double xmin, ymin, xmax, ymax;
    double x1 = x1i, x2 = x2i, y1 = y1i, y2 = y2i;
    double x1d, x2d, y1d, y2d;
    double x1_out, y1_out, x2_out, y2_out;

    graphics_get_scale(&xmin, &ymin, &xmax, &ymax);
    graphics_set_scale(Auto.xmin, Auto.ymin, Auto.xmax, Auto.ymax);
    if (graphics_clip(x1, x2, y1, y2, &x1_out, &y1_out, &x2_out, &y2_out)) {
        x1d = x1_out;
        x2d = x2_out;
        y1d = y1_out;
        y2d = y2_out;
        auto_x11_line_trans(x1d, y1d, x2d, y2d);
    }

    graphics_set_scale(xmin, ymin, xmax, ymax);
}
/* this bit of code is for writing points - it only saves what is
   in the current view

*/
int32
auto_nox_check_plot_type(int32 flag2, int32 icp1, int32 icp2) {
    if (flag2 == 0 && Auto.plot == P_P)
        return 0;
    if (flag2 > 0 && Auto.plot != P_P)
        return 0;
    if (icp1 != Auto.icp1)
        return 0;
    if (flag2 > 0 && icp1 != icp2)
        return 0;
    return 1;
}
/* main plotting code  */
void
auto_nox_add_point(double *par, double per, double *uhigh, double *ulow,
               double *ubar, double a, int32 type, int32 flg, int32 lab,
               int32 icp1, int32 icp2, int32 flag2, double *evr, double *evi) {
    double x, y1, y2, par1, par2 = 0;
    int32 ix, iy1, iy2, type1 = type;
    char bob[5];
    sprintf(bob, "%d", lab);
    par1 = par[icp1];
    if (icp2 < NAutoPar)
        par2 = par[icp2];
    auto_nox_xy_plot(&x, &y1, &y2, par1, par2, per, uhigh, ulow, ubar,
                 a); /* figure out who sits on axes */
    if (flg == 0) {
        Auto.lastx = x;
        Auto.lasty = y1;
    }
    ix = auto_nox_ix_val(x);
    iy1 = auto_nox_iy_val(y1);
    iy2 = auto_nox_iy_val(y2);
    auto_x11_bw();
    if (flag2 == 0 && Auto.plot == P_P) /* if the point was a 1 param run and we
                                           are in 2 param plot, skip */
    {
        auto_nox_plot_stab(evr, evi, NODE);
        auto_x11_refresh_display();
        return;
    }
    if (flag2 > 0 && Auto.plot != P_P) { /* two parameter and not in two
                                            parameter plot, just skip it */
        auto_nox_plot_stab(evr, evi, NODE);
        auto_x11_refresh_display();
        return;
    }

    if ((flag2 > 0) && (Auto.plot == P_P))
        type1 = CSEQ;
    switch (type1) {
    case CSEQ:
        if (Auto.plot == PE_P || Auto.plot == FR_P)
            break;
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        auto_x11_line_width(2);
        auto_nox_colset(type);
        if (flag2 > 0)
            auto_nox_colset2(flag2);
        auto_nox_line(x, y1, Auto.lastx, Auto.lasty);
        auto_x11_bw();
        break;
    case CUEQ:
        if (Auto.plot == PE_P || Auto.plot == FR_P)
            break;
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        auto_x11_line_width(1);
        auto_nox_colset(type);
        if (flag2 > 0)
            auto_nox_colset2(flag2);
        auto_nox_line(x, y1, Auto.lastx, Auto.lasty);
        auto_x11_bw();
        break;
    case UPER:
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        auto_x11_line_width(1);
        auto_nox_colset(type);
        if (flag2 > 0)
            auto_nox_colset2(flag2);
        if (auto_nox_check_bnds(ix, iy1))
            auto_x11_circle(ix, iy1, 3);
        if (auto_nox_check_bnds(ix, iy2))
            auto_x11_circle(ix, iy2, 3);
        auto_x11_bw();
        break;
    case SPER:
        if (icp1 != Auto.icp1)
            break;
        if (flag2 > 0 && Auto.icp2 != icp2)
            break;
        auto_x11_line_width(1);
        auto_nox_colset(type);
        if (flag2 > 0)
            auto_nox_colset2(flag2);
        if (auto_nox_check_bnds(ix, iy1))
            auto_x11_fill_circle(ix, iy1, 3);
        if (auto_nox_check_bnds(ix, iy2))
            auto_x11_fill_circle(ix, iy2, 3);
        auto_x11_bw();
        break;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (lab != 0) {
        if (icp1 == Auto.icp1) {
            if (flag2 == 0 || (flag2 > 0 && Auto.icp2 == icp2)) {
                auto_x11_line_width(1);
                if (auto_nox_check_bnds(ix, iy1)) {
                    auto_x11_line(ix - 4, iy1, ix + 4, iy1);
                    auto_x11_line(ix, iy1 - 4, ix, iy1 + 4);
                }
                if (auto_nox_check_bnds(ix, iy2)) {
                    auto_x11_line(ix - 4, iy2, ix + 4, iy2);
                    auto_x11_line(ix, iy2 - 4, ix, iy2 + 4);
                }
                if (auto_nox_check_bnds(ix, iy1))
                    auto_x11_text(ix + 8, iy1 + 8, bob);
            }
        }
    }

    Auto.lastx = x;
    Auto.lasty = y1;
    auto_nox_plot_stab(evr, evi, NODE);
    auto_x11_refresh_display();
}

void
auto_nox_get_bif_sym(char *at, int32 itp) {
    int32 i = itp % 10;
    switch (i) {
    case 1:
    case 6:
        sprintf(at, "BP");
        break;
    case 2:
    case 5:
        sprintf(at, "LP");
        break;
    case 3:
        sprintf(at, "HB");
        break;
    case -4:
        sprintf(at, "UZ");
        break;
    case 7:
        sprintf(at, "PD");
        break;
    case 8:
        sprintf(at, "TR");
        break;
    case 9:
        sprintf(at, "EP");
        break;
    case -9:
        sprintf(at, "MX");
        break;
    default:
        sprintf(at, "  ");
        break;
    }
    return;
}

void
auto_nox_info_header(int32 icp1, int32 icp2) {
    char bob[80];
    char p1name[12], p2name[12];

    strncpy(p1name, upar_names[AutoPar[icp1]], sizeof(p1name));
    if (icp2 < NAutoPar)
        strncpy(p2name, upar_names[AutoPar[icp2]], sizeof(p2name));
    else
        strncpy(p2name, "   ", sizeof(p2name));
    many_pops_small_base();
    sprintf(bob, "  Br  Pt Ty  Lab %10s %10s       norm %10s     period",
            p1name, p2name, uvar_names[Auto.var]);
    auto_x11_draw_info(bob, 10, DCURYs + 1);
    return;
}

void
auto_nox_new_info(int32 ibr, int32 pt, char *ty, int32 lab, double *par,
                  double norm, double u0, double per, int32 icp1, int32 icp2) {
    char bob[80];
    double p1, p2 = 0.0;
    auto_x11_clear_info();
    auto_nox_info_header(icp1, icp2);
    p1 = par[icp1];
    if (icp2 < NAutoPar)
        p2 = par[icp2];
    sprintf(bob, "%4d %4d %2s %4d %10.4g %10.4g %10.4g %10.4g %10.4g", ibr, pt,
            ty, lab, p1, p2, norm, u0, per);
    auto_x11_draw_info(bob, 10, 2*DCURYs + 2);
    auto_x11_refresh_display();
}

void
auto_nox_traverse_out(Diagram *d, int32 *ix, int32 *iy, int32 dodraw) {
    double norm, per, *par, par1, par2 = 0, *evr, *evi;
    int32 pt, itp, ibr, lab, icp1, icp2;
    double x, y1, y2;
    char symb[3];
    if (d == NULL) {
        /*ggets_err_msg("Can not traverse to NULL diagram.");*/
        return;
    }
    norm = d->norm;
    par = d->par;

    per = d->per;
    lab = d->lab;
    itp = d->itp;
    ibr = d->ibr;
    icp1 = d->icp1;
    icp2 = d->icp2;
    pt = d->ntot;

    evr = d->evr;
    evi = d->evi;

    auto_nox_get_bif_sym(symb, itp);
    par1 = par[icp1];
    if (icp2 < NAutoPar)
        par2 = par[icp2];
    auto_nox_xy_plot(&x, &y1, &y2, par1, par2, per, d->uhi, d->ulo, d->ubar, norm);

    *ix = auto_nox_ix_val(x);
    *iy = auto_nox_iy_val(y1);
    if (dodraw == 1) {
        auto_x11_xor_cross(*ix, *iy);
        auto_nox_plot_stab(evr, evi, NODE);
        auto_nox_new_info(ibr, pt, symb, lab, par, norm, d->u0[Auto.var], per,
                          icp1, icp2);
    }
    if (lab > 0 && load_all_labeled_orbits > 0)
        auto_nox_load_orbitx(ibr, 1, lab, per);
    return;
}

void
auto_nox_win(void) {
    char bob[256];
    if (Auto.exist == 0) {
        if (NODE > NAUTO) {
            sprintf(bob, "Auto restricted to less than %d variables", NAUTO);
            ggets_err_msg(bob);
            return;
        }
        auto_x11_make("It's AUTO man!", "AUTO");
        Auto.exist = 1;
    }
    return;
}

void
auto_nox_load_last_plot(int32 flg) {
    if (flg == 1) { /* one parameter */
        Auto.xmin = Old1p.xmin;
        Auto.xmax = Old1p.xmax;
        Auto.ymin = Old1p.ymin;
        Auto.ymax = Old1p.ymax;
        Auto.icp1 = Old1p.icp1;
        Auto.icp2 = Old1p.icp2;
        Auto.plot = Old1p.plot;
        Auto.var = Old1p.var;
    }
    if (flg == 2) { /* two parameter */
        Auto.xmin = Old2p.xmin;
        Auto.xmax = Old2p.xmax;
        Auto.ymin = Old2p.ymin;
        Auto.ymax = Old2p.ymax;
        Auto.icp1 = Old2p.icp1;
        Auto.icp2 = Old2p.icp2;
        Auto.plot = Old2p.plot;
        Auto.var = Old2p.var;
    }
    return;
}

void
auto_nox_keep_last_plot(int32 flg) {
    if (flg == 1) { /* one parameter */
        Old1p.xmin = Auto.xmin;
        Old1p.xmax = Auto.xmax;
        Old1p.ymin = Auto.ymin;
        Old1p.ymax = Auto.ymax;
        Old1p.icp1 = Auto.icp1;
        Old1p.icp2 = Auto.icp2;
        Old1p.plot = Auto.plot;
        Old1p.var = Auto.var;
    }
    if (flg == 2) {
        Old2p.xmin = Auto.xmin;
        Old2p.xmax = Auto.xmax;
        Old2p.ymin = Auto.ymin;
        Old2p.ymax = Auto.ymax;
        Old2p.icp1 = Auto.icp1;
        Old2p.icp2 = Auto.icp2;
        Old2p.plot = P_P;
        Old2p.var = Auto.var;
    }
    return;
}

void
auto_nox_init_win(void) {
    int32 i;
    if (NODE > NAUTO)
        return;
    diagram_start(NODE);
    for (i = 0; i < 10; i++) {
        Auto.period[i] = 11. + 3.*i;
        Auto.uzrpar[i] = 10;
        outperiod[i] = Auto.period[i];
        UzrPar[i] = 10;
    }
    NAutoPar = 8;
    if (NUPAR < 8)
        NAutoPar = NUPAR;
    for (i = 0; i < NAutoPar; i++)
        AutoPar[i] = i;
    for (i = 0; i < NAutoPar; i++) {
        Auto_index_to_array[i] = get_param_index(upar_names[AutoPar[i]]);
    }
    Auto.nper = 0;
    grabpt.flag = 0; /*  no point in buffer  */
    Auto.exist = 0;
    blrtn.irot = 0;
    for (i = 0; i < NODE; i++)
        blrtn.nrot[i] = 0;
    blrtn.torper = TOR_PERIOD;

    {
        /* auto nox create file name */
        char *basec, *bname, *dirc, *dname;
        char *HOME;

        basec = strdup(this_file);
        dirc = strdup(this_file);
        bname = (char *)basename(basec);
        dname = (char *)dirname(dirc);

        HOME = getenv("HOME");
        if (HOME == NULL)
            HOME = dname;

        sprintf(this_auto_file, "%s/%s", HOME, bname);
    }

    /*  Control -- done automatically   */
    Auto.irs = 0;
    Auto.ips = 1;
    Auto.isp = 1;
    if (SuppressBP == 1)
        Auto.isp = 0;
    Auto.ilp = 1;
    Auto.isw = 1;
    Auto.nbc = NODE;
    Auto.nfpar = 1;
    HomoFlag = 0;
    /*  User controls this      */
    Auto.ncol = auto_ncol;
    Auto.ntst = auto_ntst;
    Auto.nmx = auto_nmx;
    Auto.npr = auto_npr;
    Auto.ds = auto_ds;
    Auto.dsmax = auto_dsmax;
    Auto.dsmin = auto_dsmin;
    Auto.rl0 = auto_rl0;
    Auto.rl1 = auto_rl1;
    Auto.a0 = auto_a0;
    Auto.a1 = auto_a1;

    Auto.epsl = auto_epsl;
    Auto.epsu = auto_epsu;
    Auto.epss = auto_epss;

    /* The diagram plotting stuff    */

    Auto.xmax = auto_xmax;
    Auto.xmin = auto_xmin;
    Auto.ymax = auto_ymax;
    Auto.ymin = auto_ymin;
    Auto.plot = HL_P;
    Auto.var = auto_var;

    /* xpp parameters    */

    Auto.icp1 = 0;
    Auto.icp2 = 1;
    Auto.icp3 = 1;
    Auto.icp4 = 1;
    Auto.icp5 = 1;
    auto_nox_keep_last_plot(1);
    auto_nox_keep_last_plot(2);
    aauto.iad = 3;
    aauto.mxbf = 5;
    aauto.iid = 2;
    aauto.itmx = 8;
    aauto.itnw = 7;
    aauto.nwtn = 3;
    aauto.iads = 1;
    x_auto.nunstab = 1;
    x_auto.nstab = NODE - 1;
    return;
}

void
auto_nox_plot_stab(double *evr, double *evi, int32 n) {
    int32 i, ix, iy;
    int32 r = Auto.st_wid;

    double x;
    double y;
    auto_x11_line_width(0);
    auto_x11_clr_stab();
    for (i = 0; i < n; i++) {
        x = evr[i];
        if (x < -1.95)
            x = -1.95;
        if (x > 1.95)
            x = 1.95;
        y = evi[i];
        if (y < -1.95)
            y = -1.95;
        if (y > 1.95)
            y = 1.95;
        x = r*(x + 2.0) / 4.0;
        y = r - r*(y + 2.0) / 4.0;
        ix = (int32)x;
        iy = (int32)y;
        auto_x11_stab_line(ix - 2, iy, ix + 2, iy);
        auto_x11_stab_line(ix, iy - 2, ix, iy + 2);
    }
    return;
}

int32
auto_nox_yes_reset(void) {
    char string[256];
    if (NBifs <= 1)
        return 0;
    kill_diagrams();
    FromAutoFlag = 0;
    NBifs = 1;
    grabpt.flag = 0;
    sprintf(string, "%s.b", this_auto_file);
    auto_nox_deletef(string);
    sprintf(string, "%s.d", this_auto_file);
    auto_nox_deletef(string);
    sprintf(string, "%s.s", this_auto_file);
    auto_nox_deletef(string);
    mark_flag = 0;
    return 1;
}

int32
auto_nox_reset(void) {
    char ch;
    if (NBifs <= 1)
        return 0;
    ch = (char)menudrive_two_choice("YES", "NO", "Destroy AUTO diagram & files",
                                    "yn");
    if (ch != 'y')
        return 0;

    return auto_nox_yes_reset();
}

void
auto_nox_grab(void) {
    auto_x11_traverse_diagram();
    return;
}

void
auto_nox_get_start_period(double *p) {
    *p = storage[0][storind - 1];
    return;
}

void
auto_nox_find_best_homo_shift(int32 n) {
    /* this code looks for the best value
     * of the shift to be close as possible to the saddle
     * point of the homoclinic when starting from a
     * long periodic orbit */
    int32 i;
    int32 j;
    double dmin = 10000.0;
    double d;
    double tshift = 0.0;
    for (i = 0; i < storind; i++) {
        d = 0.0;
        for (j = 0; j < n; j++) {
            d += fabs(storage[j + 1][i] - homo_l[j]);
        }
        if (d < dmin) {
            dmin = d;
            tshift = storage[0][i];
        }
    }
    HOMO_SHIFT = tshift;
    printf("shifting %g\n", HOMO_SHIFT);
    return;
}

void
auto_nox_get_shifted_orbit(double *u, double t, double p, int32 n) {
    double ts;
    int32 i, i1, i2, ip, j;
    double lam;
    if (t > 1.0)
        t -= 1.0;
    if (t < 0.0)
        t += 1.0;
    ts = fmod(t*p + HOMO_SHIFT, p);
    for (i = 0; i < storind; i++) {
        ip = (i + 1) % storind;
        if ((ts >= storage[0][i]) && (ts < storage[0][ip])) {
            i1 = i;
            i2 = ip;
            lam = ts - storage[0][i];
            for (j = 0; j < n; j++)
                u[j] =
                    (1.0 - lam)*storage[j + 1][i1] + lam*storage[j + 1][i2];
            break;
        }
    }
    return;
}

void
auto_nox_get_start_orbit(double *u, double t, int32 n) {
    double tnorm;
    double lam;
    int32 i1, i2, j;
    if (t > 1.0)
        t -= 1.0;
    if (t < 0.0)
        t += 1.0;
    tnorm = t*(storind - 1);
    i1 = (int32)tnorm;
    i2 = i1 + 1;
    if (i2 >= storind)
        i2 -= storind;
    lam = (tnorm - (double)i1);

    for (j = 0; j < n; j++)
        u[j] = (1.0 - lam)*storage[j + 1][i1] + lam*storage[j + 1][i2];
    return;
}

void
auto_nox_run(void) {
    int32 itp1, itp2, itp, ips;
    char ch;
    if (grabpt.flag == 0) { /* the first call to AUTO   */
        /* auto start choice */
        static char *m[] = {"Steady state", "Periodic", "Bdry Value",
                            "Homoclinic", "hEteroclinic"};
        static char key[] = "spbhe";
        char ch2;
        HomoFlag = 0;
        if (METHOD == DISCRETE) {
            /* auto new discrete */
            int32 opn = NO_OPEN_3, cls = OVERWRITE;
            NewPeriodFlag = 0;

            if (NBifs > 1)
                auto_nox_reset();

            TypeOfCalc = DI1;
            Auto.ips = -1;
            Auto.irs = 0;
            Auto.itp = 0;
            Auto.ilp = 1;
            Auto.isw = 1;
            Auto.isp = 1;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.nfpar = 1;
            AutoTwoParam = 0;
            auto_nox_do(opn, cls);
            return;
        }
        ch2 = (char)auto_x11_pop_up_list("Start", m, key, 5, 13, 0, 10, 10,
                                         arun_hint, Auto.hinttxt);
        if (ch2 == 's') {
            auto_nox_new_ss();
            return;
        }
        if (ch2 == 'p') {
            /* auto start at per */
            int32 opn = NO_OPEN_3, cls = OVERWRITE;

            TypeOfCalc = PE1;
            Auto.ips = 2;
            Auto.irs = 0;
            Auto.itp = 0;
            Auto.ilp = 1;
            Auto.isw = 1;

            Auto.isp = 2;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.nfpar = 1;
            AutoTwoParam = 0;
            NewPeriodFlag = 1;
            auto_nox_do(opn, cls);
            return;
        }
        if (ch2 == 'b') {
            int32 opn, cls;
            Auto.nbc = NODE;

            /* auto start at bvp */
            opn = NO_OPEN_3;
            cls = OVERWRITE;
            pp_shoot_compile_bvp();
            if (BVP_FLAG == 0)
                return;
            TypeOfCalc = BV1;
            Auto.ips = 4;
            Auto.irs = 0;
            Auto.itp = 0;
            Auto.ilp = 1;
            Auto.isw = 1;

            Auto.isp = 2;
            if (SuppressBP == 1)
                Auto.isp = 0;

            Auto.nfpar = 1;
            AutoTwoParam = 0;
            NewPeriodFlag = 2;
            auto_nox_do(opn, cls);
            return;
        }
        if (ch2 == 'h') {
            HomoFlag = 1;
            auto_nox_start_at_homoclinic();
            return;
        }

        if (ch2 == 'e') {
            HomoFlag = 2;
            auto_nox_start_at_homoclinic();
            return;
        }

        auto_x11_redraw_menus();
        ggets_ping();
        return;
    }
    if (grabpt.lab == 0) {
        ch = (char)menudrive_two_choice("YES", "NO",
                                        "Not Labeled Pt: New Start?", "y");
        if (ch == 'y') {
            /* auto start diff ss */
            TypeOfCalc = EQ1;
            Auto.ips = 1;
            if (METHOD == DISCRETE)
                Auto.ips = -1;
            Auto.irs = 0;
            Auto.itp = 0;
            Auto.ilp = 1;
            Auto.isw = 1;
            Auto.isp = 1;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.nfpar = 1;
            AutoTwoParam = 0;
            auto_nox_do(NO_OPEN_3, APPEND);
        }

        ggets_ping();
        return;
    }

    itp = grabpt.itp;
    itp1 = itp % 10;
    itp2 = itp / 10;
    ips = Auto.ips;
    /*  printf(" ips=%d itp=%d itp1= %d itp2=%d\n",ips,itp,itp1,itp2); */
    if (itp1 == 3 || itp2 == 3) {
        /* its a HOPF Point  */
        /* auto nox hopf choice */
        static char *m[] = {"Periodic", "Extend", "New Point", "Two Param"};
        static char key[] = "pent";
        char ch2;
        if (METHOD == DISCRETE) {
            auto_nox_2p_hopf();
            return;
        }

        ch2 = (char)auto_x11_pop_up_list("Hopf Pt", m, key, 4, 10, 0, 10, 10,
                                         no_hint, Auto.hinttxt);
        if (ch2 == 'p') {
            auto_nox_new_per();
            return;
        }
        if (ch2 == 'e') {
            auto_nox_extend_ss();
            return;
        }
        if (ch2 == 'n') {
            auto_nox_new_ss();
            return;
        }
        if (ch2 == 't') {
            auto_nox_2p_hopf();
            return;
        }
        auto_x11_redraw_menus();

        ggets_ping();
        return;
    }
    if (itp1 == 7 || itp2 == 7) {
        /* period doubling */

        /* auto nox per doub choice */
        static char *m[] = {"Doubling", "Two Param", "Fixed period", "Extend"};
        static char key[] = "dtfe";
        char ch2;
        ch2 = (char)auto_x11_pop_up_list("Per. Doub.", m, key, 4, 10, 0, 10, 10,
                                         no_hint, Auto.hinttxt);
        if (ch2 == 'd') {
            /* auto period double */
            blrtn.torper = grabpt.torper;
            Auto.ntst = 2*Auto.ntst;
            Auto.irs = grabpt.lab;
            Auto.nfpar = 1; /* grabpt.nfpar; */

            Auto.itp = grabpt.itp;
            Auto.ilp = 1;
            Auto.isw = -1;
            TypeOfCalc = PE1;
            Auto.isp = 2;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.ips = 2;
            AutoTwoParam = 0;
            auto_nox_do(OPEN_3, APPEND);
            return;
        }
        if (ch2 == 'e') {
            auto_nox_new_per();
            return;
        }
        if (ch2 == 'f') {
            auto_nox_2p_fixper();
            return;
        }
        if (ch2 == 't') {
            /* auto twopar double */
            blrtn.torper = grabpt.torper;
            Auto.irs = grabpt.lab;
            Auto.itp = grabpt.itp;
            Auto.nfpar = 2;
            AutoTwoParam = PD2;
            TypeOfCalc = PD2;
            Auto.ips = 2;
            Auto.ilp = 0;
            Auto.isw = 2;
            Auto.isp = 0;
            auto_nox_do(OPEN_3, APPEND);
            return;
        }
        auto_x11_redraw_menus();

        ggets_ping();
        return;
    }
    if (ips == 9) {
        auto_nox_homo_choice(itp);
        ggets_ping();
        return;
    }
    if (itp1 == 2 || itp2 == 2) { /* limit point */
        Auto.ips = 1;
        auto_nox_2p_limit(Auto.ips);
        ggets_ping();
        return;
    }
    if (itp1 == 5 || itp2 == 5) { /* limit pt of periodic or BVP */
        if (Auto.ips != 4)
            Auto.ips = 2; /* this is a bit dangerous - the idea is that
                             if you are doing BVPs, then that is all you are
                             doing
                          */
        auto_nox_2p_limit(Auto.ips);
        ggets_ping();
        return;
    }
    if (itp1 == 6 || itp2 == 6 || itp1 == 1 || itp2 == 1) { /* branch point  */

        auto_nox_branch_choice(grabpt.ibr, ips);
        ggets_ping();
        return;
    }
    if (itp1 == 8 || itp2 == 8) { /* Torus 2 parameter */
        /* auto_nox_torus_choice */
        static char *m[] = {"Two Param", "Fixed period", "Extend"};
        /*static char *m[]={"Fixed period","Extend"}; */
        static char key[] = "tfe";
        char ch2;
        ch2 = (char)auto_x11_pop_up_list("Torus", m, key, 3, 10, 0, 10, 10,
                                         no_hint, Auto.hinttxt);
        if (ch2 == 'e') {
            auto_nox_new_per();
            return;
        }
        if (ch2 == 'f') {
            auto_nox_2p_fixper();
            return;
        }
        if (ch2 == 't') {
            /* auto torus */
            blrtn.torper = grabpt.torper;
            Auto.irs = grabpt.lab;
            Auto.itp = grabpt.itp;
            Auto.nfpar = 2;
            AutoTwoParam = TR2;
            TypeOfCalc = TR2;
            Auto.ips = 2;
            Auto.ilp = 0;
            Auto.isw = 2;
            Auto.isp = 0;
            auto_nox_do(OPEN_3, APPEND);
            return;
        }
        auto_x11_redraw_menus();

        ggets_ping();
        return;
    }
    if (grabpt.ibr < 0) {
        /* its a periodic -- just extend it  */
        /* auto_nox_periodic_choice */
        static char *m[] = {"Extend", "Fixed Period"};
        static char key[] = "ef";
        char ch2;
        ch2 = (char)auto_x11_pop_up_list("Periodic ", m, key, 2, 14, 0, 10, 10,
                                         no_hint, Auto.hinttxt);
        if (ch2 == 'e') {
            auto_nox_new_per();
            return;
        }
        if (ch2 == 'f') {
            auto_nox_2p_fixper();
            return;
        }

        auto_x11_redraw_menus();

        ggets_ping();
        return;
    }
    if (grabpt.ibr > 0 && ips != 4) { /*  old steady state -- just extend it  */
        auto_nox_extend_ss();
        ggets_ping();
        return;
    }
    if (grabpt.ibr > 0 && ips == 4) {
        /* auto extend bvp */
        TypeOfCalc = BV1;
        Auto.irs = grabpt.lab;
        Auto.itp = grabpt.itp;
        Auto.nfpar = grabpt.nfpar;
        Auto.ilp = 1;
        Auto.isw = 1;
        Auto.isp = 2;
        if (SuppressBP == 1)
            Auto.isp = 0;
        Auto.ips = 4;
        AutoTwoParam = 0;
        auto_nox_do(OPEN_3, APPEND);
        ggets_ping();
        return;
    }
    return;
}

void
auto_nox_homo_choice(int32 itp) {
    if (itp != 5) {
        /* auto extend homoclinic */
        Auto.irs = grabpt.lab;
        Auto.itp = grabpt.itp;

        TypeOfCalc = HO2;
        AutoTwoParam = HO2;
        NewPeriodFlag = 1;
        Auto.ips = 9;

        Auto.nfpar = 2;
        Auto.ilp = 1;
        Auto.isw = 1;
        Auto.isp = 0;
        Auto.nbc = 0;

        if (HomoFlag == 1)
            x_auto.iequib = 1;
        if (HomoFlag == 2)
            x_auto.iequib = -2;

        auto_nox_do(OPEN_3, APPEND);
    }
    return;
}

void
auto_nox_branch_choice(int32 ibr, int32 ips) {
    static char *m[] = {"Switch", "Extend", "New Point", "Two Param"};
    static char key[] = "sent";
    char ch;
    int32 ipsuse;
    ch = (char)auto_x11_pop_up_list("Branch Pt", m, key, 4, 10, 0, 10, 10,
                                    no_hint, Auto.hinttxt);

    if (ch == 's') {
        if (ibr < 0 && ips == 2) {
            /* auto switch per */
            TypeOfCalc = PE1;
            blrtn.torper = grabpt.torper;
            Auto.irs = grabpt.lab;
            Auto.itp = grabpt.itp;
            Auto.nfpar = 1; /*grabpt.nfpar;*/
            Auto.ilp = 1;
            Auto.isw = -1;
            Auto.isp = 2;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.ips = 2;
            AutoTwoParam = 0;
            auto_nox_do(OPEN_3, APPEND);
        } else if (ips == 4) {
            /* auto switch bvp */
            TypeOfCalc = BV1;
            Auto.irs = grabpt.lab;
            Auto.itp = grabpt.itp;
            Auto.nfpar = grabpt.nfpar;
            Auto.ilp = 1;
            Auto.isw = -1;
            Auto.isp = 2;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.ips = 4;
            AutoTwoParam = 0;
            auto_nox_do(OPEN_3, APPEND);
        } else {
            /* auto switch ss */
            TypeOfCalc = EQ1;
            Auto.irs = grabpt.lab;
            Auto.itp = grabpt.itp;
            Auto.nfpar = grabpt.nfpar;
            Auto.ilp = 1;
            Auto.isw = -1;
            Auto.isp = 1;
            if (SuppressBP == 1)
                Auto.isp = 0;
            Auto.ips = 1;
            if (METHOD == DISCRETE)
                Auto.ips = -1;
            AutoTwoParam = 0;
            auto_nox_do(OPEN_3, APPEND);
        }
        return;
    }
    if (ch == 'e') {
        auto_nox_extend_ss();
        return;
    }
    if (ch == 'n') {
        auto_nox_new_ss();
        return;
    }
    if (ch == 't') {
        ipsuse = 1;
        if (ips == 4)
            ipsuse = 4;
        if (ibr < 0)
            ipsuse = 2;
        auto_nox_2p_branch(ipsuse);
        /* auto_nox_2p_limit(ips); */
        return;
    }
    auto_x11_redraw_menus();
}

/*  RUN AUTO HERE */
/*  these are for setting the parameters to run for different choices    */

/*  Just a short recall of the AUTO parameters
   NBC = 0 unless it really is a BVP problem (not periodics or heteroclinics)
   NICP = 1 for 1 parameter and 2 fro 2 parameter and the rest will
            be taken care of in AUTO  and I think 2 for hetero?
   ILP = 1 (0) detection (no) of folds usually 1
   ISP = 2  detect all special points! but I think maybe set to 0, 1 for
            BVP I think
            for 2 parameter continuation ?
   ISW = -1 branch switching 1 is for normal 2 for two parameter of folds,
tori,HB, PD!!

   IPS   1 - std for steady states of ODEs
         -1 maps
         2 periodic orbits
         4 BVP  (set NBC=NODE)
         9 Homoclinic

for example   2 P continuation of HB
              IPS=1 ILP=1 NICP=2 ISP=0 ISW=2
BVP problem   IPS=4, NICP=1 NBC=NODE ISP=1 ISW=1 ILP=1

discrete dynamical system with two par of Hopf
first IPS=-1 ISP=ISW=1  then
NICP=2, ISW=2 at Hopf

*/

/* Start a new point for bifurcation diagram   */

void
auto_nox_new_ss(void) {
    int32 opn = NO_OPEN_3, cls = OVERWRITE;
    NewPeriodFlag = 0;

    if (NBifs > 1)
        auto_nox_reset();

    TypeOfCalc = EQ1;
    Auto.ips = 1;
    Auto.irs = 0;
    Auto.itp = 0;
    Auto.ilp = 1;
    Auto.isw = 1;
    Auto.isp = 1;
    if (SuppressBP == 1)
        Auto.isp = 0;

    Auto.nfpar = 1;
    AutoTwoParam = 0;
    auto_nox_do(opn, cls);
    return;
}

void
auto_nox_extend_ss(void) {
    /*Prevent crash on hopf of infinite period. here

    Typical abort message after crash is currently something like:

    fmt: read unexpected character
    apparent state: unit 3 named ~/fort.3
    last format: (4x,1p7e18.10)
    lately reading sequential formatted external IO

    */

    if (isinf(grabpt.per)) {
        pop_list_respond_box("Okay", "Can't continue infinite period Hopf!");
        return;
    }

    TypeOfCalc = EQ1;
    Auto.irs = grabpt.lab;
    Auto.itp = grabpt.itp;
    Auto.nfpar = grabpt.nfpar;
    Auto.ilp = 1;
    Auto.isw = 1;
    Auto.ips = 1;
    if (METHOD == DISCRETE)
        Auto.ips = -1;
    Auto.isp = 1;
    if (SuppressBP == 1)
        Auto.isp = 0;

    AutoTwoParam = 0;
    auto_nox_do(OPEN_3, APPEND);
    return;
}

int32
auto_nox_get_homo_info(int32 *nun, int32 *nst, double *ul, double *ur) {
    char **s;
    char v[100][MAX_LEN_SBOX];
    int32 n = 2 + 2*NODE;
    int32 i;
    int32 flag = 0;
    s = xmalloc((usize)n*sizeof(*s));
    for (i = 0; i < n; i++) {
        s[i] = xmalloc(100);
    }
    sprintf(s[0], "dim unstable");
    sprintf(v[0], "%d", *nun);
    sprintf(s[NODE + 1], "dim stable");
    sprintf(v[NODE + 1], "%d", *nst);
    for (i = 0; i < NODE; i++) {
        sprintf(s[i + 1], "%s_L", uvar_names[i]);
        sprintf(v[i + 1], "%g", ul[i]);
        sprintf(s[i + 2 + NODE], "%s_R", uvar_names[i]);
        sprintf(v[i + 2 + NODE], "%g", ur[i]);
    }

    flag = pop_list_do_string_box(n, n / 2, 2, "Homoclinic info", s, v, 16);
    if (flag != 0) {
        *nun = atoi(v[0]);
        *nst = atoi(v[NODE + 1]);
        for (i = 0; i < NODE; i++) {
            ul[i] = atof(v[i + 1]);
            if (HomoFlag == 2)
                ur[i] = atof(v[i + 2 + NODE]);
        }
    }
    for (i = 0; i < n; i++) {
        free(s[i]);
    }
    free(s);

    return flag;
}

void
auto_nox_start_at_homoclinic(void) {
    int32 opn = NO_OPEN_3;
    int32 flag;
    Auto.irs = 0;
    Auto.itp = 0;
    TypeOfCalc = HO2;

    AutoTwoParam = HO2;
    NewPeriodFlag = 1;
    Auto.ips = 9;

    Auto.nfpar = 2;
    Auto.ilp = 1; /* maybe 1 someday also in extend homo, but for now, no 3
                     param allowed    */
    Auto.isw = 1;
    Auto.isp = 0;
    Auto.nbc = 0;

    if (HomoFlag == 1) {
        x_auto.iequib = 1;
        auto_nox_find_best_homo_shift(NODE);
    }
    if (HomoFlag == 2)
        x_auto.iequib = -2;
    flag =
        auto_nox_get_homo_info(&x_auto.nunstab, &x_auto.nstab, homo_l, homo_r);
    if (flag) {
        /* TODO: for some reason, the second argument was `close`, which maps
         * to the libc function with this name. That does not make any sense
         * so I changed it to 1 */
        auto_nox_do(opn, 1);
    }
    return;
}

void
auto_nox_new_per(void) {
    /* same for extending periodic  */
    blrtn.torper = grabpt.torper;

    /*Prevent crash on hopf of infinite period. here

    Typical abort message after crash is currently something like:

    fmt: read unexpected character
    apparent state: unit 3 named ~/fort.3
    last format: (4x,1p7e18.10)
    lately reading sequential formatted external IO

    */

    if (isinf(grabpt.per)) {
        pop_list_respond_box("Okay", "Can't continue infinite period Hopf.");
        return;
    }
    TypeOfCalc = PE1;
    Auto.irs = grabpt.lab;
    Auto.itp = grabpt.itp;
    /* Auto.nfpar=grabpt.nfpar; */
    Auto.nfpar = 1;
    Auto.ilp = 1;
    Auto.isw = 1; /* -1 */
    Auto.isp = 2;
    if (SuppressBP == 1)
        Auto.isp = 0;
    Auto.ips = 2;
    AutoTwoParam = 0;
    auto_nox_do(OPEN_3, APPEND);
    return;
}

void
auto_nox_2p_limit(int32 ips) {
    int32 ipsuse = 1;
    int32 itp1;
    int32 itp2;
    blrtn.torper = grabpt.torper;
    Auto.irs = grabpt.lab;
    itp1 = (grabpt.itp) % 10;
    itp2 = abs(grabpt.itp) / 10;
    Auto.itp = grabpt.itp;
    Auto.nfpar = 2;
    Auto.ilp = 0; /* was 1 */
    Auto.isw = 2;
    Auto.isp = 0; /* was 2 */
    /* fix ips now */
    if (ips == 4)
        ipsuse = 4;
    else {
        if ((itp1 == 5) || (itp2 == 5))
            ipsuse = 2;
    }

    Auto.ips = ipsuse;
    AutoTwoParam = LPP2;
    if (ipsuse == 1) {
        TypeOfCalc = LPE2;
        AutoTwoParam = LPE2;
    } else {
        TypeOfCalc = LPP2;
        AutoTwoParam = LPP2;
    }
    /* printf("ips=%d  itp=%d \n",Auto.ips,Auto.itp); */
    auto_nox_do(OPEN_3, APPEND);
    return;
}

void
auto_nox_2p_branch(int32 ips) {
    int32 ipsuse = 1;
    int32 itp1;
    int32 itp2;
    blrtn.torper = grabpt.torper;
    Auto.irs = grabpt.lab;
    itp1 = (grabpt.itp) % 10;
    itp2 = abs(grabpt.itp) / 10;
    Auto.itp = grabpt.itp;
    Auto.nfpar = 2;
    Auto.ilp = 0; /* was 1 */
    Auto.isw = 2;
    Auto.isp = 0; /* was 2 */
    if (ips == 4)
        ipsuse = 4;
    else {
        if ((itp1 == 6) || (itp2 == 6))
            ipsuse = 2;
    }

    Auto.ips = ipsuse;
    if (METHOD == DISCRETE)
        Auto.ips = -1;
    AutoTwoParam = BR2;
    TypeOfCalc = BR2;
    auto_nox_do(OPEN_3, APPEND);
    return;
}

void
auto_nox_2p_fixper(void) {
    Auto.irs = grabpt.lab;
    Auto.itp = grabpt.itp;
    Auto.nfpar = 2;
    Auto.ilp = 1; /* was1 */
    Auto.isw = 1;
    Auto.isp = 0;
    Auto.ips = 2;
    AutoTwoParam = FP2;
    TypeOfCalc = FP2;
    auto_nox_do(OPEN_3, APPEND);
    return;
}

void
auto_nox_2p_hopf(void) {
    /*Prevent crash on hopf of infinite period. here

    Typical abort message after crash is currently something like:

    fmt: read unexpected character
    apparent state: unit 3 named ~/fort.3
    last format: (4x,1p7e18.10)
    lately reading sequential formatted external IO

    */

    if (isinf(grabpt.per)) {
        pop_list_respond_box("Okay", "Can't continue infinite period Hopf.");
        return;
    }

    Auto.irs = grabpt.lab;
    Auto.itp = grabpt.itp;
    Auto.nfpar = 2;
    Auto.ilp = 0; /* was 1 */
    Auto.isw = 2;
    Auto.isp = 0;
    Auto.ips = 1;
    if (METHOD == DISCRETE)
        Auto.ips = -1;
    AutoTwoParam = HB2;
    TypeOfCalc = HB2;
    auto_nox_do(OPEN_3, APPEND);
    return;
}

/**********   END RUN AUTO *********************/

void
auto_nox_err(char *s) {
    pop_list_respond_box("OK", s);
}

void
auto_nox_load_orbit(void) {
    auto_nox_load_orbitx(grabpt.ibr, grabpt.flag, grabpt.lab, grabpt.per);
    return;
}

void
auto_nox_load_orbitx(int32 ibr, int32 flag, int32 lab, double per) {
    FILE *fp;
    double *x;
    int32 i, j, nstor;
    double u[NAUTO], t;
    double period;
    char string[256];
    int32 nrow, ndim, label, flg;

    if ((ibr > 0 && (Auto.ips != 4) && (Auto.ips != 3) && (Auto.ips != 9)) ||
        flag == 0)
        return;
    // either nothing grabbed or just a fixed point and that is already loaded
    sprintf(string, "%s.s", this_auto_file);
    fp = fopen(string, "r");
    if (fp == NULL) {
        auto_nox_err("No such file");
        return;
    }
    label = lab;
    period = per;
    flg = auto_nox_move_to_label(label, &nrow, &ndim, fp);
    nstor = ndim;
    if (ndim > NODE)
        nstor = NODE;
    if (flg == 0) {
        printf("Could not find label %d in file %s \n", label, string);
        auto_nox_err("Cant find labeled pt");
        fclose(fp);
        return;
    }
    x = &MyData[0];
    for (i = 0; i < nrow; i++) {
        auto_nox_get_a_row(u, &t, ndim, fp);
        if (Auto.ips != 4)
            storage[0][i] = t*period;
        else
            storage[0][i] = t;

        for (j = 0; j < nstor; j++) {
            storage[j + 1][i] = u[j];
            x[j] = u[j];
        }
        main_rhs_extra(x, (double)storage[0][i], nstor, NEQ);
        for (j = nstor; j < NEQ; j++)
            storage[j + 1][i] = (double)x[j];
    }
    storind = nrow;
    refresh_browser(nrow);
    /* insert auxiliary stuff here */
    if (load_all_labeled_orbits == 2)
        menudrive_clr_all_scrns();
    menudrive_drw_all_scrns();
    fclose(fp);
    return;
}

void
auto_nox_save(void) {
    int32 ok;
    FILE *fp;
    /*char filename[256];*/
    char filename[XPP_MAX_NAME];
    int32 status;
    /* XGetInputFocus(display,&w,&rev); */

    sprintf(filename, "%s.auto", basename(this_auto_file));
    /* status=dialog_box_get("Save Auto","Filename",filename,"Ok","Cancel",60);
    XSetInputFocus(display,w,rev,CurrentTime);
    */
    status = init_conds_file_selector("Save Auto", filename, "*.auto");
    if (status == 0)
        return;
    browse_open_write_file(&fp, filename, &ok);
    if (!ok)
        return;
    auto_nox_save_numerics(fp);
    auto_nox_save_graph(fp);
    status = diagram_save(fp, NODE);
    if (status != 1) {
        fclose(fp);
        return;
    }
    auto_nox_save_q_file(fp);
    fclose(fp);
    return;
}

void
auto_nox_save_numerics(FILE *fp) {
    int32 i;
    fprintf(fp, "%d ", NAutoPar);
    for (i = 0; i < NAutoPar; i++)
        fprintf(fp, "%d ", AutoPar[i]);
    fprintf(fp, "%d\n", NAutoUzr);
    for (i = 0; i < 9; i++)
        fprintf(fp, "%g %ld\n", outperiod[i], UzrPar[i]);
    fprintf(fp, "%d %d %d \n", Auto.ntst, Auto.nmx, Auto.npr);
    fprintf(fp, "%g %g %g \n", Auto.ds, Auto.dsmin, Auto.dsmax);
    fprintf(fp, "%g %g %g %g\n", Auto.rl0, Auto.rl1, Auto.a0, Auto.a1);
    fprintf(fp, "%d %d %d %d %d %d %d\n", aauto.iad, aauto.mxbf, aauto.iid,
            aauto.itmx, aauto.itnw, aauto.nwtn, aauto.iads);
    return;
}

void
auto_nox_load_numerics(FILE *fp) {
    int32 in;
    fscanf(fp, "%d ", &NAutoPar);
    for (int64 i = 0; i < NAutoPar; i++) {
        fscanf(fp, "%d ", &AutoPar[i]);
        in = get_param_index(upar_names[AutoPar[i]]);
        Auto_index_to_array[i] = in;
    }
    fscanf(fp, "%d ", &NAutoUzr);
    for (int32 i = 0; i < 9; i++) {
        Auto.nper = NAutoUzr;
        fscanf(fp, "%lg %ld\n", &outperiod[i], &UzrPar[i]);
        Auto.period[i] = outperiod[i];
        Auto.uzrpar[i] = UzrPar[i];
    }

    fscanf(fp, "%d %d %d \n", &Auto.ntst, &Auto.nmx, &Auto.npr);
    fscanf(fp, "%lg %lg %lg \n", &Auto.ds, &Auto.dsmin, &Auto.dsmax);
    fscanf(fp, "%lg %lg %lg %lg\n", &Auto.rl0, &Auto.rl1, &Auto.a0, &Auto.a1);
    fscanf(fp, "%d %d %d %d %d %d %d\n", &aauto.iad, &aauto.mxbf, &aauto.iid,
           &aauto.itmx, &aauto.itnw, &aauto.nwtn, &aauto.iads);
    return;
}

void
auto_nox_save_graph(FILE *fp) {
    fprintf(fp, "%g %g %g %g %d %d \n", Auto.xmin, Auto.ymin, Auto.xmax,
            Auto.ymax, Auto.var, Auto.plot);
    return;
}

void
auto_nox_load_graph(FILE *fp) {
    fscanf(fp, "%lg %lg %lg %lg %d %d \n", &Auto.xmin, &Auto.ymin, &Auto.xmax,
           &Auto.ymax, &Auto.var, &Auto.plot);
    return;
}

void
auto_nox_save_q_file(/* I am keeping the name q_file even though they are
                        s_files */
                     FILE *fp) {
    char string[500];
    FILE *fq;
    sprintf(string, "%s.s", this_auto_file);
    fq = fopen(string, "r");
    if (fq == NULL) {
        auto_nox_err("Couldnt open s-file");
        return;
    }
    while (!feof(fq)) {
        fgets(string, 500, fq);
        fputs(string, fp);
        /* break; */
    }
    fclose(fq);
    return;
}

void
auto_nox_make_q_file(FILE *fp) {
    char string[500];
    FILE *fq;
    sprintf(string, "%s.s", this_auto_file);
    fq = fopen(string, "w");

    if (fq == NULL) {
        auto_nox_err("Couldnt open s-file");
        return;
    }

    while (!feof(fp)) {
        fgets(string, 500, fp);
        if (!auto_nox_no_info_noinfo(string)) {
            fputs(string, fq);
        }
    }
    fclose(fq);
    return;
}

int32
auto_nox_no_info_noinfo(char *s) {
    /* get rid of any blank lines */
    int32 n = (int32)strlen(s);
    int32 i;
    if (n == 0)
        return 1;
    for (i = 0; i < n; i++) {
        if (!isspace(s[i]))
            return 0;
    }
    return 1;
}

void
auto_nox_load(void) {
    int32 ok;
    FILE *fp;
    /*char filename[256];*/
    char filename[XPP_MAX_NAME];
    int32 status;
    if (NBifs > 1) {
        ok = auto_nox_reset();
        if (ok == 0)
            return;
    }

    sprintf(filename, "%s.auto", basename(this_auto_file));

    status = init_conds_file_selector("Load Auto", filename, "*.auto");
    if (status == 0)
        return;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        auto_nox_err("Cannot open file");
        return;
    }

    auto_nox_load_numerics(fp);
    auto_nox_load_graph(fp);
    status = diagram_load(fp, NODE);
    if (status != 1) {
        fclose(fp);
        return;
    }
    auto_nox_make_q_file(fp);
    fclose(fp);
    return;
}

int32
auto_nox_move_to_label(int32 mylab, int32 *nrow, int32 *ndim, FILE *fp) {
    int32 ibr, ntot, itp, lab, nfpar, isw, ntpl, nar, nskip;
    int32 i;
    char line[MAXLINELENGTH];
    while (true) {
        fgets(line, MAXLINELENGTH, fp);
        sscanf(line, "%d%d %d %d %d %d %d %d %d", &ibr, &ntot, &itp, &lab,
               &nfpar, &isw, &ntpl, &nar, &nskip);
        if (mylab == lab) {
            *nrow = ntpl;
            *ndim = nar - 1;
            return 1;
        }
        for (i = 0; i < nskip; i++)
            fgets(line, MAXLINELENGTH, fp);
        if (feof(fp))
            break;
    }
    return 0;
}

void
auto_nox_get_a_row(double *u, double *t, int32 n, FILE *fp) {
    int32 i;
    fscanf(fp, "%lg ", t);
    for (i = 0; i < n; i++)
        fscanf(fp, "%lg ", &u[i]);
    return;
}

void
auto_nox_file(void) {
    static char *m[] = {"Import orbit",   "Save diagram",  "Load diagram",
                        "Postscript",     "SVG",           "Reset diagram",
                        "Clear grab",     "Write pts",     "All info",
                        "init Data",      "Toggle redraw", "auto raNge",
                        "sElect 2par pt", "draw laBled",   "lOad branch"};
    static char key[] = "islpvrcwadtnebo";
    char ch;
    ch = (char)auto_x11_pop_up_list("File", m, key, 15, 15, 0, 10, 10,
                                    afile_hint, Auto.hinttxt);
    if (ch == 'i') {
        auto_nox_load_orbit();
        return;
    }
    if (ch == 's') {
        auto_nox_save();
        return;
    }
    if (ch == 'l') {
        auto_nox_load();
        return;
    }
    if (ch == 'r')
        auto_nox_reset();
    if (ch == 'c')
        grabpt.flag = 0;
    if (ch == 'p') {
        NoBreakLine = 1;
        diagram_post_auto();
        NoBreakLine = 0;
    }
    if (ch == 'v') {
        NoBreakLine = 1;
        diagram_svg_auto();
        NoBreakLine = 0;
    }
    if (ch == 'w')
        diagram_write_pts();
    if (ch == 'a')
        diagram_write_info_out();
    if (ch == 'd')
        diagram_write_init_data_file();
    if (ch == 't') {
        auto_redraw_flag = 1 - auto_redraw_flag;
        if (auto_redraw_flag == 1)
            ggets_err_msg("Redraw is ON");
        else
            ggets_err_msg("Redraw is OFF");
    }
    if (ch == 'o') {
        if (mark_flag < 2)
            ggets_err_msg("Mark a branch first using S and E");
        else
            diagram_load_browser_with_branch(mark_ibrs, mark_ipts, mark_ipte);
    }
    if (ch == 'n') {
        if (mark_flag < 2)
            ggets_err_msg("Mark a branch first using S and E");
        else
            auto_x11_do_range();
    }
    if (ch == 'e') {
        if (Auto.plot != P_P) {
            ggets_err_msg("Must be in 2 parameter plot");
            return;
        }
        if (FromAutoFlag) {
            /* auto nox set point */
            FromAutoFlag = 0;
            set_val(upar_names[AutoPar[Auto.icp1]], XfromAuto);
            set_val(upar_names[AutoPar[Auto.icp2]], YfromAuto);
            derived_evaluate();
            tabular_redo_all_fun_tables();
            init_conds_redraw_params();
        }
    }
    if (ch == 'b') {
        if (load_all_labeled_orbits == 0) {
            load_all_labeled_orbits = 1;
            ggets_err_msg("Draw orbits - no erase");
            return;
        }
        if (load_all_labeled_orbits == 1) {
            load_all_labeled_orbits = 2;
            ggets_err_msg("Draw orbits - erase first");
            return;
        }
        if (load_all_labeled_orbits == 2) {
            load_all_labeled_orbits = 0;
            ggets_err_msg("Draw orbits off");
            return;
        }
    }
    return;
}
