#include "functions.h"
#include "parserslow.h"
#include "integers.h"

#include "read_dir.h"

/*#include "macdirent.h"
 */

#include <dirent.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "xpplim.h"
#include "auto_nox.h"

#define Param 1
#define IC 2

#define DFNORMAL 1
#define MAXOPT 1000
#define READEM 1

static char *interopt[MAXOPT];
static int32 Nopts = 0;
int32 RunImmediately = 0;

int32 IX_PLT[10];
int32 IY_PLT[10];
int32 IZ_PLT[10];
int32 NPltV;
int32 MultiWin = 0;
double X_LO[10];
double Y_LO[10];
double X_HI[10];
double Y_HI[10];
int32 START_LINE_TYPE = 1;
InternSet intern_set[MAX_INTERN_SET];
int32 Nintern_set = 0;

/*   this file has all of the phaseplane parameters defined
     and created.  All other files should use external stuff
    to use them. (Except eqn forming stuff)
 */

double last_ic[MAX_ODE];

char delay_string[MAX_ODE][80];
int32 itor[MAX_ODE];
char this_file[XPP_MAX_NAME];
char this_internset[XPP_MAX_NAME];
int32 mov_ind;
int32 storind;
int32 STORFLAG;
int32 INFLAG;
int32 MAXSTOR;
double x_3d[2];
double y_3d[2];
double z_3d[2];
int32 IXPLT;
int32 IYPLT;
int32 IZPLT;
int32 AXES;
int32 TIMPLOT;
int32 PLOT_3D;
double MY_XLO;
double MY_YLO;
double MY_XHI;
double MY_YHI;
double TOR_PERIOD = 6.2831853071795864770;
int32 TORUS = 0;
int32 NEQ;
char options[100];

/*   Numerical stuff ....   */

double DELTA_T;
double TEND;
double T0;
double TRANS;
double NULL_ERR;
double EVEC_ERR;
double NEWT_ERR;
double BOUND;
double DELAY;
double TOLER;
double ATOLER;
double HMIN;
double HMAX;
double BVP_EPS;
double BVP_TOL;

double POIPLN;

int32 euler_max_iter;
double euler_tol;
int32 NMESH;
int32 NJMP;
int32 METHOD;
int32 color_flag;
int32 NC_ITER;
int32 EVEC_ITER;
int32 BVP_MAXIT;
int32 BVP_FLAG;

int32 POIMAP;
int32 POIVAR;
int32 POISGN;
int32 SOS;
int32 FFT;
int32 NULL_HERE;
int32 POIEXT;
int32 HIST;
int32 HVAR;
int32 hist_ind;
int32 FOREVER;

/*  control of range stuff  */

int32 ENDSING;
int32 SHOOT;
int32 PAR_FOL;

/*  custon color stuff  */

int32 xorfix;
int32 silent;
int32 got_file;

static void load_eqn_split_apart(char *bob, char *name, char *value);
static void load_eqn_do_intern_set(char *name1, char *value);
static void load_eqn_fil_int(FILE *fpt, int32 *val);
static void load_eqn_fil_flt(FILE *fpt, double *val);
static void load_eqn_read_defaults(FILE *fp);

void
load_eqn_dump_torus(FILE *fp, int32 f) {
    char bob[256];
    if (f == READEM) {
        fgets(bob, 255, fp);
    } else {
        fprintf(fp, "# Torus information \n");
    }
    lunch_io_int(&TORUS, fp, f, " Torus flag 1=ON");
    lunch_io_double(&TOR_PERIOD, fp, f, "Torus period");
    if (TORUS) {
        for (int32 i = 0; i < NEQ; i++) {
            lunch_io_int(&itor[i], fp, f, uvar_names[i]);
        }
    }
    return;
}

void
load_eqn(void) {
    int32 no_eqn = 1;
    int32 okay = 0;
    int32 std = 0;
    FILE *fptr;
    struct dirent *dp;

    integrate_init_ar_ic();
    for (int32 i = 0; i < MAX_ODE; i++) {
        itor[i] = 0;
        strcpy(delay_string[i], "0.0");
    }
    if (strcmp(this_file, "/dev/stdin") == 0) {
        std = 1;
    }
    if (got_file == 1 && (std == 0) &&
        (dp = (struct dirent *)opendir(this_file)) != NULL) {
        no_eqn = 1;
        okay = 0;
        read_dir_change_dir(this_file);
        okay = form_ode_make_eqn();
        return;
    } else {
        if (got_file == 1 && (fptr = fopen(this_file, "r")) != NULL) {
            if (std == 1) {
                sprintf(this_file, "console");
            }
            okay = form_ode_get_eqn(fptr);
            if (std == 0) {
                fclose(fptr);
            }

            if (okay == 1) {
                no_eqn = 0;
            }
        }
    }
    if (no_eqn) {
        while (okay == 0) {
            struct dirent *dp2;
            char odeclassrm[256];
            if (getenv("XPPSTART") != NULL) {
                sprintf(odeclassrm, "%s", getenv("XPPSTART"));

                if ((dp2 = (struct dirent *)opendir(odeclassrm)) != NULL) {
                    read_dir_change_dir(odeclassrm);
                }
            }

            okay = form_ode_make_eqn();
        }
    }
    return;
}

void
load_eqn_set_x_vals(void) {
    /*
    Set up the default look here.
    */

    tfBell = 1;
    /*
    No gradients tends to look cleaner but some
    may prefer gradients improved contrast/readability.
    */
    /*fixed is the new X11 default fixed font. 9x15 is dead and gone.
     */
    if (strlen(font_name_big) == 0) {
        strcpy(font_name_big, "fixed");
    }

    if (strlen(font_name_small) == 0) {
        strcpy(font_name_small, "6x13");
    }

    if (strlen(user_black) == 0) {
        sprintf(user_black, "#%s", "000000");
    }

    if (strlen(user_white) == 0) {
        sprintf(user_white, "#%s", "EDE9E3");
    }

    if (strlen(user_main_win_color) == 0) {
        sprintf(user_main_win_color, "#%s", "808080");
    }

    if (strlen(UserDrawWinColor) == 0) {
        sprintf(UserDrawWinColor, "#%s", "FFFFFF");
    }

    if (UserGradients < 0) {
        UserGradients = 1;
    }
    return;
}

void
load_eqn_set_all_vals(void) {
    FILE *fp;

    if (notAlreadySet.TIMEPLOT) {
        TIMPLOT = 1;
        notAlreadySet.TIMEPLOT = 0;
    }
    if (notAlreadySet.FOREVER) {
        FOREVER = 0;
        notAlreadySet.FOREVER = 0;
    }
    if (notAlreadySet.BVP_TOL) {
        BVP_TOL = 1.e-5;
        notAlreadySet.BVP_TOL = 0;
    }
    if (notAlreadySet.BVP_EPS) {
        BVP_EPS = 1.e-5;
        notAlreadySet.BVP_EPS = 0;
    }
    if (notAlreadySet.BVP_MAXIT) {
        BVP_MAXIT = 20;
        notAlreadySet.BVP_MAXIT = 0;
    }
    if (notAlreadySet.BVP_FLAG) {
        BVP_FLAG = 0;
        notAlreadySet.BVP_FLAG = 0;
    }
    if (notAlreadySet.NMESH) {
        NMESH = 40;
        notAlreadySet.NMESH = 0;
    }
    if (notAlreadySet.NOUT) {
        NJMP = 1;
        notAlreadySet.NOUT = 0;
    }
    if (notAlreadySet.SOS) {
        SOS = 0;
        notAlreadySet.SOS = 0;
    }
    if (notAlreadySet.FFT) {
        FFT = 0;
        notAlreadySet.FFT = 0;
    }
    if (notAlreadySet.HIST) {
        HIST = 0;
        notAlreadySet.HIST = 0;
    }
    if (notAlreadySet.PltFmtFlag) {
        PltFmtFlag = 0;
        notAlreadySet.PltFmtFlag = 0;
    }
    if (notAlreadySet.AXES) {
        AXES = 0;
        notAlreadySet.AXES = 0;
    }
    if (notAlreadySet.TOLER) {
        TOLER = 0.001;
        notAlreadySet.TOLER = 0;
    }
    if (notAlreadySet.ATOLER) {
        ATOLER = 0.001;
        notAlreadySet.ATOLER = 0;
    }
    if (notAlreadySet.euler_max_iter) {
        euler_max_iter = 10;
        notAlreadySet.euler_max_iter = 0;
    }
    if (notAlreadySet.euler_tol) {
        euler_tol = 1.e-7;
        notAlreadySet.euler_tol = 0;
    }
    if (notAlreadySet.DELAY) {
        DELAY = 0.0;
        notAlreadySet.DELAY = 0;
    }
    if (notAlreadySet.DTMIN) {
        HMIN = 1e-12;
        notAlreadySet.DTMIN = 0;
    }
    if (notAlreadySet.EVEC_ITER) {
        EVEC_ITER = 100;
        notAlreadySet.EVEC_ITER = 0;
    }
    if (notAlreadySet.EVEC_ERR) {
        EVEC_ERR = .001;
        notAlreadySet.EVEC_ERR = 0;
    }
    if (notAlreadySet.NULL_ERR) {
        NULL_ERR = .001;
        notAlreadySet.NULL_ERR = 0;
    }
    if (notAlreadySet.NEWT_ERR) {
        NEWT_ERR = .001;
        notAlreadySet.NEWT_ERR = 0;
    }
    if (notAlreadySet.NULL_HERE) {
        NULL_HERE = 0;
        notAlreadySet.NULL_HERE = 0;
    }
    del_stab_flag = DFNORMAL;
    if (notAlreadySet.DTMAX) {
        HMAX = 1.000;
        notAlreadySet.DTMAX = 0;
    }
    if (notAlreadySet.POIMAP) {
        POIMAP = 0;
        notAlreadySet.POIMAP = 0;
    }
    if (notAlreadySet.POIVAR) {
        POIVAR = 1;
        notAlreadySet.POIVAR = 0;
    }
    if (notAlreadySet.POIEXT) {
        POIEXT = 0;
        notAlreadySet.POIEXT = 0;
    }
    if (notAlreadySet.POISGN) {
        POISGN = 1;
        notAlreadySet.POISGN = 0;
    }
    if (notAlreadySet.POIPLN) {
        POIPLN = 0.0;
        notAlreadySet.POIPLN = 0;
    }

    storind = 0;
    mov_ind = 0;

    STORFLAG = 0;

    INFLAG = 0;
    solver = odesol_rung_kut;
    PLOT_3D = 0;
    if (notAlreadySet.METHOD) {
        METHOD = 3;
        notAlreadySet.METHOD = 0;
    }
    if (notAlreadySet.XLO) {
        MY_XLO = 0.0;
        x_3d[0] = MY_XLO;
        notAlreadySet.XLO = 0;
        notAlreadySet.XMIN = 0;
    }
    if (notAlreadySet.XHI) {
        MY_XHI = 20.0;
        x_3d[1] = MY_XHI;
        notAlreadySet.XHI = 0;
        notAlreadySet.XMAX = 0;
    }
    if (notAlreadySet.YLO) {
        MY_YLO = -1;
        y_3d[0] = MY_YLO;
        notAlreadySet.YLO = 0;
        notAlreadySet.YMIN = 0;
    }
    if (notAlreadySet.YHI) {
        MY_YHI = 1;
        y_3d[0] = MY_YHI;
        notAlreadySet.YHI = 0;
        notAlreadySet.YMAX = 0;
    }

    if (notAlreadySet.BOUND) {
        BOUND = 100;
        notAlreadySet.BOUND = 0;
    }
    if (notAlreadySet.MAXSTOR) {
        MAXSTOR = 5000;
        notAlreadySet.MAXSTOR = 0;
    }

    if (notAlreadySet.T0) {
        T0 = 0.0;
        notAlreadySet.T0 = 0;
    }
    if (notAlreadySet.TRANS) {
        TRANS = 0.0;
        notAlreadySet.TRANS = 0;
    }
    if (notAlreadySet.DT) {
        DELTA_T = .05;
        notAlreadySet.DT = 0;
    }

    if (notAlreadySet.XMIN) {
        x_3d[0] = -12;
        notAlreadySet.XMIN = 0;
        notAlreadySet.XLO = 0;
    }
    if (notAlreadySet.XMAX) {
        x_3d[1] = 12;
        notAlreadySet.XMAX = 0;
        notAlreadySet.XHI = 0;
    }
    if (notAlreadySet.YMIN) {
        y_3d[0] = -12;
        notAlreadySet.YMIN = 0;
        notAlreadySet.YLO = 0;
    }
    if (notAlreadySet.YMAX) {
        y_3d[1] = 12;
        notAlreadySet.YMAX = 0;
        notAlreadySet.YHI = 0;
    }
    if (notAlreadySet.ZMIN) {
        z_3d[0] = -12;
        notAlreadySet.ZMIN = 0;
    }
    if (notAlreadySet.ZMAX) {
        z_3d[1] = 12;
        notAlreadySet.ZMAX = 0;
    }

    if (notAlreadySet.TEND) {
        TEND = 20.00;
        notAlreadySet.TEND = 0;
    }
    if (notAlreadySet.IXPLT) {
        IXPLT = 0;
        notAlreadySet.IXPLT = 0;
    }
    if (notAlreadySet.IYPLT) {
        IYPLT = 1;
        notAlreadySet.IYPLT = 0;
    }
    if (notAlreadySet.IZPLT) {
        IZPLT = 1;
        notAlreadySet.IZPLT = 0;
    }

    if (notAlreadySet.NPLOT) {
        if (NEQ > 2) {
            if (notAlreadySet.IZPLT) {
                IZPLT = 2;
            }
        }
        NPltV = 1;
        for (int32 i = 0; i < 10; i++) {
            IX_PLT[i] = IXPLT;
            IY_PLT[i] = IYPLT;
            IZ_PLT[i] = IZPLT;
            X_LO[i] = 0;
            Y_LO[i] = -1;
            X_HI[i] = 20;
            Y_HI[i] = 1;
        }
        notAlreadySet.NPLOT = 0;
    }
    // internal options go here
    load_eqn_set_internopts(NULL);

    if ((fp = fopen(options, "r")) != NULL) {
        load_eqn_read_defaults(fp);
        fclose(fp);
    }

    integrate_init_range();
    adjoints_init_trans();
    array_plot_init_my();
    txt_init_view();

    numerics_chk_volterra();

    //

    if (IZPLT > NEQ) {
        IZPLT = NEQ;
    }
    if (IYPLT > NEQ) {
        IYPLT = NEQ;
    }
    if (IXPLT == 0 || IYPLT == 0) {
        TIMPLOT = 1;
    } else {
        TIMPLOT = 0;
    }
    if (x_3d[0] >= x_3d[1]) {
        x_3d[0] = -1;
        x_3d[1] = 1;
    }
    if (y_3d[0] >= y_3d[1]) {
        y_3d[0] = -1;
        y_3d[1] = 1;
    }
    if (z_3d[0] >= z_3d[1]) {
        z_3d[0] = -1;
        z_3d[1] = 1;
    }
    if (MY_XLO >= MY_XHI) {
        MY_XLO = -2.0;
        MY_XHI = 2.0;
    }
    if (MY_YLO >= MY_YHI) {
        MY_YLO = -2.0;
        MY_YHI = 2.0;
    }
    if (AXES < 5) {
        x_3d[0] = MY_XLO;
        y_3d[0] = MY_YLO;
        x_3d[1] = MY_XHI;
        y_3d[1] = MY_YHI;
    }
    storage_init_stor(MAXSTOR, NEQ + 1);
    if (AXES >= 5) {
        PLOT_3D = 1;
    }
    numerics_chk_delay();  // check for delay allocation
    adjoints_alloc_h_stuff();

    volterra_alloc_memory();  // allocate stuff for volterra equations
    storage_alloc_meth();
    integrate_arr_ic_start();  // take care of all predefined array ics
    return;
}

void
load_eqn_read_defaults(FILE *fp) {
    char bob[100];
    char *ptr;
    fgets(bob, 80, fp);
    ptr = form_ode_get_first(bob, " ");
    if (notAlreadySet.BIG_FONT_NAME) {
        strcpy(font_name_big, ptr);
        notAlreadySet.BIG_FONT_NAME = 0;
    }

    fgets(bob, 80, fp);
    ptr = form_ode_get_first(bob, " ");
    if (notAlreadySet.SMALL_FONT_NAME) {
        strcpy(font_name_small, ptr);
        notAlreadySet.SMALL_FONT_NAME = 0;
    }

    if (notAlreadySet.paper_white) {
        load_eqn_fil_int(fp, &paper_white);
        notAlreadySet.paper_white = 0;
    }
    if (notAlreadySet.IXPLT) {
        load_eqn_fil_int(fp, &IXPLT);
        notAlreadySet.IXPLT = 0;
    }
    if (notAlreadySet.IYPLT) {
        load_eqn_fil_int(fp, &IYPLT);
        notAlreadySet.IYPLT = 0;
    }
    if (notAlreadySet.IZPLT) {
        load_eqn_fil_int(fp, &IZPLT);
        notAlreadySet.IZPLT = 0;
    }
    if (notAlreadySet.AXES) {
        load_eqn_fil_int(fp, &AXES);
        notAlreadySet.paper_white = 0;
    }
    if (notAlreadySet.NOUT) {
        load_eqn_fil_int(fp, &NJMP);
        notAlreadySet.NOUT = 0;
    }
    if (notAlreadySet.NMESH) {
        load_eqn_fil_int(fp, &NMESH);
        notAlreadySet.NMESH = 0;
    }
    if (notAlreadySet.METHOD) {
        load_eqn_fil_int(fp, &METHOD);
        notAlreadySet.METHOD = 0;
    }

    if (notAlreadySet.TIMEPLOT) {
        load_eqn_fil_int(fp, &TIMPLOT);
        notAlreadySet.TIMEPLOT = 0;
    }
    if (notAlreadySet.MAXSTOR) {
        load_eqn_fil_int(fp, &MAXSTOR);
        notAlreadySet.MAXSTOR = 0;
    }
    if (notAlreadySet.TEND) {
        load_eqn_fil_flt(fp, &TEND);
        notAlreadySet.TEND = 0;
    }
    if (notAlreadySet.DT) {
        load_eqn_fil_flt(fp, &DELTA_T);
        notAlreadySet.DT = 0;
    }
    if (notAlreadySet.T0) {
        load_eqn_fil_flt(fp, &T0);
        notAlreadySet.T0 = 0;
    }
    if (notAlreadySet.TRANS) {
        load_eqn_fil_flt(fp, &TRANS);
        notAlreadySet.TRANS = 0;
    }
    if (notAlreadySet.BOUND) {
        load_eqn_fil_flt(fp, &BOUND);
        notAlreadySet.BOUND = 0;
    }
    if (notAlreadySet.DTMIN) {
        load_eqn_fil_flt(fp, &HMIN);
        notAlreadySet.DTMIN = 0;
    }
    if (notAlreadySet.DTMAX) {
        load_eqn_fil_flt(fp, &HMAX);
        notAlreadySet.DTMIN = 0;
    }
    if (notAlreadySet.TOLER) {
        load_eqn_fil_flt(fp, &TOLER);
        notAlreadySet.TOLER = 0;
    }
    if (notAlreadySet.DELAY) {
        load_eqn_fil_flt(fp, &DELAY);
        notAlreadySet.DELAY = 0;
    }
    if (notAlreadySet.XLO) {
        load_eqn_fil_flt(fp, &MY_XLO);
        notAlreadySet.XLO = 0;
    }
    if (notAlreadySet.XHI) {
        load_eqn_fil_flt(fp, &MY_XHI);
        notAlreadySet.XHI = 0;
    }
    if (notAlreadySet.YLO) {
        load_eqn_fil_flt(fp, &MY_YLO);
        notAlreadySet.YLO = 0;
    }
    if (notAlreadySet.YHI) {
        load_eqn_fil_flt(fp, &MY_YHI);
        notAlreadySet.YHI = 0;
    }
    return;
}

void
load_eqn_fil_flt(FILE *fpt, double *val) {
    char bob[80];
    fgets(bob, 80, fpt);
    *val = atof(bob);
    return;
}

void
load_eqn_fil_int(FILE *fpt, int32 *val) {
    char bob[80];
    fgets(bob, 80, fpt);
    *val = atoi(bob);
    return;
}

/* here is some new code for internal set files:
 * format of the file is a long string of the form:
 * { x=y, z=w, q=p , .... }
 */

void
load_eqn_add_intern_set(char *name, char *does) {
    char bob[1024];
    char ch;
    int32 n;
    int32 j = Nintern_set;
    int32 k = 0;

    if (Nintern_set >= MAX_INTERN_SET) {
        ggets_plintf(" %s not added -- too many must be less than %d \n", name,
                     MAX_INTERN_SET);
        return;
    }
    intern_set[j].use = 1;
    n = (int32)strlen(name);
    intern_set[j].name = xmalloc((usize)n + 1);
    strcpy(intern_set[j].name, name);
    n = (int32)strlen(does);
    bob[0] = '$';
    bob[1] = ' ';
    k = 2;
    for (int32 i = 0; i < n; i++) {
        ch = does[i];
        if (ch == ',') {
            bob[k] = ' ';
            k++;
        }
        if (ch == '}' || ch == '{') {
            continue;
        }
        if (ch != ',') {
            bob[k] = ch;
            k++;
        }
    }
    bob[k] = 0;
    intern_set[j].does = xmalloc((usize)n + 3);
    strcpy(intern_set[j].does, bob);
    ggets_plintf(" added %s doing %s \n", intern_set[j].name,
                 intern_set[j].does);
    Nintern_set++;
    return;
}

void
load_eqn_extract_action(char *ptr) {
    char name[256];
    char value[256];
    char tmp[2048];
    char *junk;
    char *mystring;
    strcpy(tmp, ptr);
    junk = form_ode_get_first(tmp, " ");
    if (junk == NULL) {
        // No more tokens--should this throw an error?
    }

    while ((mystring = form_ode_do_fit_get_next(" ,;\n")) != NULL) {
        load_eqn_split_apart(mystring, name, value);
        if (strlen(name) > 0 && strlen(value) > 0) {
            load_eqn_do_intern_set(name, value);
        }
    }
    return;
}

void
load_eqn_extract_internset(int32 j) {
    load_eqn_extract_action(intern_set[j].does);
    return;
}

void
load_eqn_do_intern_set(char *name1, char *value) {
    int32 i;
    char name[20];
    convert(name1, name);

    i = init_conds_find_user_name(IC, name);
    if (i > -1) {
        last_ic[i] = atof(value);
    } else {
        i = init_conds_find_user_name(Param, name);
        if (i > -1) {
            set_val(name, atof(value));
        } else {
            load_eqn_set_option(name, value, 1, NULL);
        }
    }
    storage_alloc_meth();
    numerics_do_meth();
}
/*  ODE options stuff  here !!   */

int32
load_eqn_msc(char *s1, char *s2) {
    usize n = strlen(s1);
    if (n > strlen(s2)) {
        return 0;
    }
    for (usize i = 0; i < n; i++) {
        if (s1[i] != s2[i]) {
            return 0;
        }
    }
    return 1;
}

void
load_eqn_set_internopts(OptionsSet *mask) {
    char *ptr, name[20], value[80], *junk, *mystring;
    if (Nopts == 0) {
        return;
    }
    //  parsem here
    for (int32 i = 0; i < Nopts; i++) {
        ptr = interopt[i];
        junk = form_ode_get_first(ptr, " ,");
        if (junk == NULL) {
            // No more tokens.  Should this throw an error?
        }
        while ((mystring = form_ode_do_fit_get_next(" ,\n\r")) != NULL) {
            load_eqn_split_apart(mystring, name, value);
            if (strlen(name) > 0 && strlen(value) > 0) {
                /*
                if (strcmp("mwcolor",name)==0)
                {
                        if (strlen(user_main_win_color)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("dwcolor",name)==0)
                {
                        if (strlen(UserDrawWinColor)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("forecolor",name)==0)
                {
                        if (strlen(user_white)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("backcolor",name)==0)
                {
                        if (strlen(user_black)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("backimage",name)==0)
                {
                        if (strlen(user_bg_bitmap)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("smallfont",name)==0)
                {
                        if (strlen(font_name_small)!=0)
                        {
                                continue;
                        }
                }

                if (strcmp("bigfont",name)==0)
                {
                        if (strlen(font_name_big)!=0)
                        {
                                continue;
                        }
                }
                */
                load_eqn_set_option(name, value, 0, mask);
            }
        }
    }

    for (int32 i = 0; i < Nopts; i++) {
        free(interopt[i]);
    }
    Nopts = 0;
    return;
}

void
load_eqn_set_internopts_xpprc_and_comline(void) {
    char *ptr, name[20], value[80], *junk, *mystring;
    OptionsSet *tempNAS;
    //  parsem here
    // Check for QUIET and LOGFILE options first...
    char intrnoptcpy[255]; /*Must use copy to avoid side effects of strtok used
                              in get_first below*/
    if (Nopts == 0) {
        return;
    }
    for (int32 i = 0; i < Nopts; i++) {
        strcpy(intrnoptcpy, interopt[i]);
        ptr = intrnoptcpy;
        junk = form_ode_get_first(ptr, " ,");
        if (junk == NULL) {
            // No more tokens.  Should this throw an error?
        }
        while ((mystring = form_ode_do_fit_get_next(" ,\n\r")) != NULL) {
            load_eqn_split_apart(mystring, name, value);
            strupr(name);

            if (strlen(name) == 5) {
                strupr(name);
                if (strcmp(name, "QUIET") == 0) {
                    load_eqn_set_option(name, value, 0, NULL);
                }
            } else if (strlen(name) == 7) {
                strupr(name);

                if (strcmp(name, "LOGFILE") == 0) {
                    load_eqn_set_option(name, value, 0, NULL);
                }
            }
        }
    }

    // We make a BOOLEAN MASK using the current OptionsSet
    /*This allows options to be overwritten multiple times within .xpprc
    but prevents overwriting across comline, .xpprc etc.
    */
    tempNAS = xmalloc(sizeof(*tempNAS));
    *tempNAS = notAlreadySet;

    for (int32 i = 0; i < Nopts; i++) {
        ptr = interopt[i];
        junk = form_ode_get_first(ptr, " ,");
        while ((mystring = form_ode_do_fit_get_next(" ,\n\r")) != NULL) {
            load_eqn_split_apart(mystring, name, value);
            if (strlen(name) > 0 && strlen(value) > 0) {
                load_eqn_set_option(name, value, 0, tempNAS);
            }
        }
    }
    free(tempNAS);

    /*
    We leave a fresh start for options specified in the ODE file.
    */
    for (int32 i = 0; i < Nopts; i++) {
        free(interopt[i]);
    }

    Nopts = 0;
    return;
}

void
load_eqn_split_apart(char *bob, char *name, char *value) {
    int32 k;
    int32 l;

    l = (int32)strlen(bob);
    k = (int32)strcspn(bob, "=");
    if (k == l) {
        value[0] = 0;
        strcpy(name, bob);
    } else {
        strncpy(name, bob, (usize)k);
        name[k] = '\0';
        for (int32 i = k + 1; i < l; i++) {
            value[i - k - 1] = bob[i];
        }
        value[l - k - 1] = '\0';
    }
    return;
}

void
load_eqn_check_for_xpprc(void) {
    FILE *fp;
    char rc[256];
    char bob[256];
    sprintf(rc, "%s/.xpprc", getenv("HOME"));
    fp = fopen(rc, "r");
    if (fp == NULL) {
        return;
    }
    while (!feof(fp)) {
        bob[0] = '\0';
        fgets(bob, 255, fp);
        if (bob[0] == '@') {
            load_eqn_stor_internopts(bob);
        }
    }
    fclose(fp);
    return;
}

void
load_eqn_stor_internopts(char *s1) {
    int32 n = (int32)strlen(s1);
    if (Nopts > MAXOPT) {
        ggets_plintf("WARNING -- to many options set %s ignored\n", s1);
        return;
    }
    interopt[Nopts] = xmalloc((usize)n + 1);
    sprintf(interopt[Nopts], "%s", s1);
    Nopts++;
    return;
}

void
load_eqn_set_option(char *s1, char *s2, int32 force, OptionsSet *mask) {
    int32 i;
    int32 f;
    char xx[4];
    char yy[4];
    char zz[4];
    char xxl[6];
    char xxh[6];
    char yyl[6];
    char yyh[6];
    static char mkey[] = "demragvbqsc582y";
    static char Mkey[] = "DEMRAGVBQSC582Y";
    strupr(s1);
    if (load_eqn_msc("QUIET", s1)) {
        if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
            ggets_plintf("QUIET option must be 0 or 1.\n");
            exit(-1);
        }
        if (OVERRIDE_QUIET ==
            0)  // Will be 1 if -quiet was specified on the command line.
        {
            XPPVERBOSE = (atoi(s2) == 0);
        }
        return;
    }
    if (load_eqn_msc("LOGFILE", s1)) {
        if (OVERRIDE_LOGFILE ==
            0)  // Will be 1 if -logfile was specified on the command line.
        {
            if (logfile != NULL) {
                fclose(logfile);
            }
            logfile = fopen(s2, "w");
        }
        return;
    }
    if (load_eqn_msc("BELL", s1)) {
        if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
            ggets_plintf("BELL option must be 0 or 1.\n");
            exit(-1);
        }
        tfBell = atoi(s2);
        return;
    }
    if ((load_eqn_msc("BIGFONT", s1)) || (load_eqn_msc("BIG", s1))) {
        if ((notAlreadySet.BIG_FONT_NAME || force) ||
            ((mask != NULL) && (mask->BIG_FONT_NAME == 1))) {
            strcpy(font_name_big, s2);
            notAlreadySet.BIG_FONT_NAME = 0;
        }
        return;
    }
    if ((load_eqn_msc("SMALLFONT", s1)) || (load_eqn_msc("SMALL", s1))) {
        if ((notAlreadySet.SMALL_FONT_NAME || force) ||
            ((mask != NULL) && (mask->SMALL_FONT_NAME == 1))) {
            strcpy(font_name_small, s2);
            notAlreadySet.SMALL_FONT_NAME = 0;
        }
        return;
    }
    if (load_eqn_msc("FORECOLOR", s1)) {
        if ((notAlreadySet.user_black || force) ||
            ((mask != NULL) && (mask->user_black == 1))) {
            sprintf(user_black, "#%s", s2);
            notAlreadySet.user_black = 0;
        }
        return;
    }
    if (load_eqn_msc("BACKCOLOR", s1)) {
        if ((notAlreadySet.user_white || force) ||
            ((mask != NULL) && (mask->user_white == 1))) {
            sprintf(user_white, "#%s", s2);
            notAlreadySet.user_white = 0;
        }
        return;
    }
    if (load_eqn_msc("MWCOLOR", s1)) {
        if ((notAlreadySet.user_main_win_color || force) ||
            ((mask != NULL) && (mask->user_main_win_color == 1))) {
            sprintf(user_main_win_color, "#%s", s2);
            notAlreadySet.user_main_win_color = 0;
        }
        return;
    }
    if (load_eqn_msc("DWCOLOR", s1)) {
        if ((notAlreadySet.UserDrawWinColor || force) ||
            ((mask != NULL) && (mask->UserDrawWinColor == 1))) {
            sprintf(UserDrawWinColor, "#%s", s2);
            notAlreadySet.UserDrawWinColor = 0;
        }
        return;
    }
    if (load_eqn_msc("GRADS", s1)) {
        if ((notAlreadySet.UserGradients || force) ||
            ((mask != NULL) && (mask->UserGradients == 1))) {
            if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
                ggets_plintf("GRADS option must be 0 or 1.\n");
                exit(-1);
            }
            UserGradients = atoi(s2);
            notAlreadySet.UserGradients = 0;
        }
        return;
    }

    if (load_eqn_msc("PLOTFMT", s1)) {
        if ((notAlreadySet.PLOTFORMAT || force) ||
            ((mask != NULL) && (mask->PLOTFORMAT == 1))) {
            strcpy(plot_format, s2);
            notAlreadySet.PLOTFORMAT = 0;
        }
        return;
    }

    if (load_eqn_msc("BACKIMAGE", s1)) {
        if ((notAlreadySet.user_bg_bitmap || force) ||
            ((mask != NULL) && (mask->user_bg_bitmap == 1))) {
            strcpy(user_bg_bitmap, s2);
            notAlreadySet.user_bg_bitmap = 0;
        }
        return;
    }
    if (load_eqn_msc("WIDTH", s1)) {
        if ((notAlreadySet.UserMinWidth || force) ||
            ((mask != NULL) && (mask->UserMinWidth == 1))) {
            UserMinWidth = atoi(s2);
            notAlreadySet.UserMinWidth = 0;
        }
        return;
    }
    if (load_eqn_msc("HEIGHT", s1)) {
        if ((notAlreadySet.UserMinHeight || force) ||
            ((mask != NULL) && (mask->UserMinHeight == 1))) {
            UserMinHeight = atoi(s2);
            notAlreadySet.UserMinHeight = 0;
        }
        return;
    }
    if (load_eqn_msc("YNC", s1)) {
        if ((notAlreadySet.YNullColor || force) ||
            ((mask != NULL) && (mask->YNullColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                YNullColor = i;
            }
            notAlreadySet.YNullColor = 0;
        }
        return;
    }
    if (load_eqn_msc("XNC", s1)) {
        if ((notAlreadySet.XNullColor || force) ||
            ((mask != NULL) && (mask->XNullColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                XNullColor = i;
                notAlreadySet.XNullColor = 0;
            }
        }
        return;
    }

    if (load_eqn_msc("SMC", s1)) {
        if ((notAlreadySet.StableManifoldColor || force) ||
            ((mask != NULL) && (mask->StableManifoldColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                StableManifoldColor = i;
                notAlreadySet.StableManifoldColor = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("UMC", s1)) {
        if ((notAlreadySet.UnstableManifoldColor || force) ||
            ((mask != NULL) && (mask->UnstableManifoldColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                UnstableManifoldColor = i;
                notAlreadySet.UnstableManifoldColor = 0;
            }
        }
        return;
    }

    if (load_eqn_msc("LT", s1)) {
        if ((notAlreadySet.START_LINE_TYPE || force) ||
            ((mask != NULL) && (mask->START_LINE_TYPE == 1))) {
            i = atoi(s2);
            if (i < 2 && i > -6) {
                START_LINE_TYPE = i;
                graphics_reset_all_line_type();
                notAlreadySet.START_LINE_TYPE = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("SEED", s1)) {
        if ((notAlreadySet.RandSeed || force) ||
            ((mask != NULL) && (mask->RandSeed == 1))) {
            i = atoi(s2);
            if (i >= 0) {
                RandSeed = i;
                markov_nsrand48(RandSeed);
                notAlreadySet.RandSeed = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("BACK", s1)) {
        if ((notAlreadySet.paper_white || force) ||
            ((mask != NULL) && (mask->paper_white == 1))) {
            if (s2[0] == 'w' || s2[0] == 'W') {
                paper_white = 1;
            } else {
                paper_white = 0;
            }
            notAlreadySet.paper_white = 0;
        }
        return;
    }
    if (load_eqn_msc("COLORMAP", s1)) {
        if ((notAlreadySet.COLORMAP || force) ||
            ((mask != NULL) && (mask->COLORMAP == 1))) {
            i = atoi(s2);
            if (i < 7) {
                custom_color = i;
            }
            notAlreadySet.COLORMAP = 0;
        }
        return;
    }
    if (load_eqn_msc("NPLOT", s1)) {
        if ((notAlreadySet.NPLOT || force) ||
            ((mask != NULL) && (mask->NPLOT == 1))) {
            NPltV = atoi(s2);
            notAlreadySet.NPLOT = 0;
        }
        return;
    }

    if (load_eqn_msc("DLL_LIB", s1)) {
        if ((notAlreadySet.DLL_LIB || force) ||
            ((mask != NULL) && (mask->DLL_LIB == 1))) {
            sprintf(dll_lib, "%s", s2);
            dll_flag += 1;
            notAlreadySet.DLL_LIB = 0;
        }
        return;
    }
    if (load_eqn_msc("DLL_FUN", s1)) {
        if ((notAlreadySet.DLL_FUN || force) ||
            ((mask != NULL) && (mask->DLL_FUN == 1))) {
            sprintf(dll_fun, "%s", s2);
            dll_flag += 2;
            notAlreadySet.DLL_FUN = 0;
        }
        return;
    }
    // can now initialize several plots
    if (load_eqn_msc("SIMPLOT", s1)) {
        SimulPlotFlag = 1;
        return;
    }
    if (load_eqn_msc("MULTIWIN", s1)) {
        MultiWin = 1;
        return;
    }
    for (int32 j = 2; j <= 8; j++) {
        sprintf(xx, "XP%d", j);
        sprintf(yy, "YP%d", j);
        sprintf(zz, "ZP%d", j);
        sprintf(xxh, "XHI%d", j);
        sprintf(xxl, "XLO%d", j);
        sprintf(yyh, "YHI%d", j);
        sprintf(yyl, "YLO%d", j);
        if (load_eqn_msc(xx, s1)) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IX_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(yy, s1)) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IY_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(zz, s1)) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IZ_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(xxh, s1)) {
            X_HI[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(xxl, s1)) {
            X_LO[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(yyh, s1)) {
            Y_HI[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(yyl, s1)) {
            Y_LO[j] = atof(s2);
            return;
        }
    }
    if (load_eqn_msc("XP", s1)) {
        if ((notAlreadySet.XP || force) ||
            ((mask != NULL) && (mask->XP == 1))) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IXPLT = i;
            }
            notAlreadySet.XP = 0;
            notAlreadySet.IXPLT = 0;
        }
        return;
    }
    if (load_eqn_msc("YP", s1)) {
        if ((notAlreadySet.YP || force) ||
            ((mask != NULL) && (mask->YP == 1))) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IYPLT = i;
            }
            notAlreadySet.YP = 0;
            notAlreadySet.IYPLT = 0;
        }
        return;
    }
    if (load_eqn_msc("ZP", s1)) {
        if ((notAlreadySet.ZP || force) ||
            ((mask != NULL) && (mask->ZP == 1))) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                IZPLT = i;
            }

            notAlreadySet.ZP = 0;
            notAlreadySet.IZPLT = 0;
        }
        return;
    }
    if (load_eqn_msc("AXES", s1)) {
        if ((notAlreadySet.AXES || force) ||
            ((mask != NULL) && (mask->AXES == 1))) {
            if (s2[0] == '3') {
                AXES = 5;
            } else {
                AXES = 0;
            }

            notAlreadySet.AXES = 0;
        }
        return;
    }

    if (load_eqn_msc("NJMP", s1)) {
        if ((notAlreadySet.NOUT || force) ||
            ((mask != NULL) && (mask->NOUT == 1))) {
            NJMP = atoi(s2);
            notAlreadySet.NOUT = 0;
        }
        return;
    }
    if (load_eqn_msc("NOUT", s1)) {
        if ((notAlreadySet.NOUT || force) ||
            ((mask != NULL) && (mask->NOUT == 1))) {
            NJMP = atoi(s2);
            notAlreadySet.NOUT = 0;
        }
        return;
    }
    if (load_eqn_msc("NMESH", s1)) {
        if ((notAlreadySet.NMESH || force) ||
            ((mask != NULL) && (mask->NMESH == 1))) {
            NMESH = atoi(s2);
            notAlreadySet.NMESH = 0;
        }
        return;
    }
    if (load_eqn_msc("METH", s1)) {
        if ((notAlreadySet.METHOD || force) ||
            ((mask != NULL) && (mask->METHOD == 1))) {
            for (i = 0; i < 15; i++) {
                if (s2[0] == mkey[i] || s2[0] == Mkey[i]) {
                    METHOD = i;
                }
            }

            notAlreadySet.METHOD = 0;
        }
        return;
    }
    if (load_eqn_msc("VMAXPTS", s1)) {
        if ((notAlreadySet.VMAXPTS || force) ||
            ((mask != NULL) && (mask->VMAXPTS == 1))) {
            MaxPoints = atoi(s2);
            notAlreadySet.VMAXPTS = 0;
        }
        return;
    }
    if (load_eqn_msc("MAXSTOR", s1)) {
        if ((notAlreadySet.MAXSTOR || force) ||
            ((mask != NULL) && (mask->MAXSTOR == 1))) {
            MAXSTOR = atoi(s2);
            notAlreadySet.MAXSTOR = 0;
        }
        return;
    }
    if (load_eqn_msc("TOR_PER", s1)) {
        if ((notAlreadySet.TOR_PER || force) ||
            ((mask != NULL) && (mask->TOR_PER == 1))) {
            TOR_PERIOD = atof(s2);
            TORUS = 1;
            notAlreadySet.TOR_PER = 0;
        }
        return;
    }
    if (load_eqn_msc("JAC_EPS", s1)) {
        if ((notAlreadySet.JAC_EPS || force) ||
            ((mask != NULL) && (mask->JAC_EPS == 1))) {
            NEWT_ERR = atof(s2);
            notAlreadySet.JAC_EPS = 0;
        }
        return;
    }
    if (load_eqn_msc("NEWT_TOL", s1)) {
        if ((notAlreadySet.NEWT_TOL || force) ||
            ((mask != NULL) && (mask->NEWT_TOL == 1))) {
            EVEC_ERR = atof(s2);
            notAlreadySet.NEWT_TOL = 0;
        }
        return;
    }
    if (load_eqn_msc("NEWT_ITER", s1)) {
        if ((notAlreadySet.NEWT_ITER || force) ||
            ((mask != NULL) && (mask->NEWT_ITER == 1))) {
            EVEC_ITER = atoi(s2);
            notAlreadySet.NEWT_ITER = 0;
        }
        return;
    }
    if (load_eqn_msc("FOLD", s1)) {
        if ((notAlreadySet.FOLD || force) ||
            ((mask != NULL) && (mask->FOLD == 1))) {
            browse_find_variable(s2, &i);
            if (i > 0) {
                itor[i - 1] = 1;
                TORUS = 1;
            }
        }
        return;
    }
    if (load_eqn_msc("TOTAL", s1)) {
        if ((notAlreadySet.TEND || force) ||
            ((mask != NULL) && (mask->TEND == 1))) {
            TEND = atof(s2);
            notAlreadySet.TEND = 0;
        }
        return;
    }
    if (load_eqn_msc("DTMIN", s1)) {
        if ((notAlreadySet.DTMIN || force) ||
            ((mask != NULL) && (mask->DTMIN == 1))) {
            HMIN = atof(s2);
            notAlreadySet.DTMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("DTMAX", s1)) {
        if ((notAlreadySet.DTMAX || force) ||
            ((mask != NULL) && (mask->DTMAX == 1))) {
            HMAX = atof(s2);
            notAlreadySet.DTMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("DT", s1)) {
        if ((notAlreadySet.DT || force) ||
            ((mask != NULL) && (mask->DT == 1))) {
            DELTA_T = atof(s2);
            notAlreadySet.DT = 0;
        }
        return;
    }
    if (load_eqn_msc("T0", s1)) {
        if ((notAlreadySet.T0 || force) ||
            ((mask != NULL) && (mask->T0 == 1))) {
            T0 = atof(s2);
            notAlreadySet.T0 = 0;
        }
        return;
    }
    if (load_eqn_msc("TRANS", s1)) {
        if ((notAlreadySet.TRANS || force) ||
            ((mask != NULL) && (mask->TRANS == 1))) {
            TRANS = atof(s2);
            notAlreadySet.TRANS = 0;
        }
        return;
    }
    if (load_eqn_msc("BOUND", s1)) {
        if ((notAlreadySet.BOUND || force) ||
            ((mask != NULL) && (mask->BOUND == 1))) {
            BOUND = atof(s2);
            notAlreadySet.BOUND = 0;
        }
        return;
    }
    if (load_eqn_msc("ATOL", s1)) {
        if ((notAlreadySet.ATOLER || force) ||
            ((mask != NULL) && (mask->ATOLER == 1))) {
            ATOLER = atof(s2);
            notAlreadySet.ATOLER = 0;
        }
        return;
    }
    if (load_eqn_msc("TOL", s1)) {
        if ((notAlreadySet.TOLER || force) ||
            ((mask != NULL) && (mask->TOLER == 1))) {
            TOLER = atof(s2);
            notAlreadySet.TOLER = 0;
        }
        return;
    }

    if (load_eqn_msc("DELAY", s1)) {
        if ((notAlreadySet.DELAY || force) ||
            ((mask != NULL) && (mask->DELAY == 1))) {
            DELAY = atof(s2);
            notAlreadySet.DELAY = 0;
        }
        return;
    }
    if (load_eqn_msc("BANDUP", s1)) {
        if ((notAlreadySet.BANDUP || force) ||
            ((mask != NULL) && (mask->BANDUP == 1))) {
            cv_bandflag = 1;
            cv_bandupper = atoi(s2);
            notAlreadySet.BANDUP = 0;
        }
        return;
    }
    if (load_eqn_msc("BANDLO", s1)) {
        if ((notAlreadySet.BANDLO || force) ||
            ((mask != NULL) && (mask->BANDLO == 1))) {
            cv_bandflag = 1;
            cv_bandlower = atoi(s2);
            notAlreadySet.BANDLO = 0;
        }
        return;
    }

    if (load_eqn_msc("PHI", s1)) {
        if ((notAlreadySet.PHI || force) ||
            ((mask != NULL) && (mask->PHI == 1))) {
            PHI0 = atof(s2);
            notAlreadySet.PHI = 0;
        }
        return;
    }
    if (load_eqn_msc("THETA", s1)) {
        if ((notAlreadySet.THETA || force) ||
            ((mask != NULL) && (mask->THETA == 1))) {
            THETA0 = atof(s2);
            notAlreadySet.THETA = 0;
        }
        return;
    }
    if (load_eqn_msc("XLO", s1)) {
        if ((notAlreadySet.XLO || force) ||
            ((mask != NULL) && (mask->XLO == 1))) {
            MY_XLO = atof(s2);
            notAlreadySet.XLO = 0;
        }
        return;
    }
    if (load_eqn_msc("YLO", s1)) {
        if ((notAlreadySet.YLO || force) ||
            ((mask != NULL) && (mask->YLO == 1))) {
            MY_YLO = atof(s2);
            notAlreadySet.YLO = 0;
        }
        return;
    }

    if (load_eqn_msc("XHI", s1)) {
        if ((notAlreadySet.XHI || force) ||
            ((mask != NULL) && (mask->XHI == 1))) {
            MY_XHI = atof(s2);
            notAlreadySet.XHI = 0;
        }
        return;
    }
    if (load_eqn_msc("YHI", s1)) {
        if ((notAlreadySet.YHI || force) ||
            ((mask != NULL) && (mask->YHI == 1))) {
            MY_YHI = atof(s2);
            notAlreadySet.YHI = 0;
        }
        return;
    }
    if (load_eqn_msc("XMAX", s1)) {
        if ((notAlreadySet.XMAX || force) ||
            ((mask != NULL) && (mask->XMAX == 1))) {
            x_3d[1] = atof(s2);
            notAlreadySet.XMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("YMAX", s1)) {
        if ((notAlreadySet.YMAX || force) ||
            ((mask != NULL) && (mask->YMAX == 1))) {
            y_3d[1] = atof(s2);
            notAlreadySet.YMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("ZMAX", s1)) {
        if ((notAlreadySet.ZMAX || force) ||
            ((mask != NULL) && (mask->ZMAX == 1))) {
            z_3d[1] = atof(s2);
            notAlreadySet.ZMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("XMIN", s1)) {
        if ((notAlreadySet.XMIN || force) ||
            ((mask != NULL) && (mask->XMIN == 1))) {
            x_3d[0] = atof(s2);
            notAlreadySet.XMIN = 0;
            if ((notAlreadySet.XLO || force) ||
                ((mask != NULL) && (mask->XLO == 1))) {
                MY_XLO = atof(s2);
                notAlreadySet.XLO = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("YMIN", s1)) {
        if ((notAlreadySet.YMIN || force) ||
            ((mask != NULL) && (mask->YMIN == 1))) {
            y_3d[0] = atof(s2);
            notAlreadySet.YMIN = 0;
            if ((notAlreadySet.YLO || force) ||
                ((mask != NULL) && (mask->YLO == 1))) {
                MY_YLO = atof(s2);
                notAlreadySet.YLO = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("ZMIN", s1)) {
        if ((notAlreadySet.ZMIN || force) ||
            ((mask != NULL) && (mask->ZMIN == 1))) {
            z_3d[0] = atof(s2);
            notAlreadySet.ZMIN = 0;
        }
        return;
    }

    if (load_eqn_msc("POIMAP", s1)) {
        if ((notAlreadySet.POIMAP || force) ||
            ((mask != NULL) && (mask->POIMAP == 1))) {
            if (s2[0] == 'm' || s2[0] == 'M') {
                POIMAP = 2;
            }
            if (s2[0] == 's' || s2[0] == 'S') {
                POIMAP = 1;
            }
            if (s2[0] == 'p' || s2[0] == 'P') {
                POIMAP = 3;
            }
            notAlreadySet.POIMAP = 0;
        }
        return;
    }

    if (load_eqn_msc("POIVAR", s1)) {
        if ((notAlreadySet.POIVAR || force) ||
            ((mask != NULL) && (mask->POIVAR == 1))) {
            browse_find_variable(s2, &i);
            if (i > -1) {
                POIVAR = i;
            }

            notAlreadySet.POIVAR = 0;
        }
        return;
    }
    if (load_eqn_msc("OUTPUT", s1)) {
        if ((notAlreadySet.OUTPUT || force) ||
            ((mask != NULL) && (mask->OUTPUT == 1))) {
            strcpy(batch_out, s2);
            notAlreadySet.OUTPUT = 0;
        }
        return;
    }

    if (load_eqn_msc("POISGN", s1)) {
        if ((notAlreadySet.POISGN || force) ||
            ((mask != NULL) && (mask->POISGN == 1))) {
            POISGN = atoi(s2);
            notAlreadySet.POISGN = 0;
        }
        return;
    }

    if (load_eqn_msc("POISTOP", s1)) {
        if ((notAlreadySet.POISTOP || force) ||
            ((mask != NULL) && (mask->POISTOP == 1))) {
            SOS = atoi(s2);
            notAlreadySet.POISTOP = 0;
        }
        return;
    }
    if (load_eqn_msc("STOCH", s1)) {
        if ((notAlreadySet.STOCH || force) ||
            ((mask != NULL) && (mask->STOCH == 1))) {
            STOCH_FLAG = atoi(s2);
            notAlreadySet.STOCH = 0;
        }
        return;
    }
    if (load_eqn_msc("POIPLN", s1)) {
        if ((notAlreadySet.POIPLN || force) ||
            ((mask != NULL) && (mask->POIPLN == 1))) {
            POIPLN = atof(s2);
            notAlreadySet.POIPLN = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGEOVER", s1)) {
        if ((notAlreadySet.RANGEOVER || force) ||
            ((mask != NULL) && (mask->RANGEOVER == 1))) {
            strcpy(range.item, s2);
            notAlreadySet.RANGEOVER = 0;
        }

        return;
    }
    if (load_eqn_msc("RANGESTEP", s1)) {
        if ((notAlreadySet.RANGESTEP || force) ||
            ((mask != NULL) && (mask->RANGESTEP == 1))) {
            range.steps = atoi(s2);
            notAlreadySet.RANGESTEP = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGELOW", s1)) {
        if ((notAlreadySet.RANGELOW || force) ||
            ((mask != NULL) && (mask->RANGELOW == 1))) {
            range.plow = atof(s2);
            notAlreadySet.RANGELOW = 0;
        }

        return;
    }

    if (load_eqn_msc("RANGEHIGH", s1)) {
        if ((notAlreadySet.RANGEHIGH || force) ||
            ((mask != NULL) && (mask->RANGEHIGH == 1))) {
            range.phigh = atof(s2);
            notAlreadySet.RANGEHIGH = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGERESET", s1)) {
        if ((notAlreadySet.RANGERESET || force) ||
            ((mask != NULL) && (mask->RANGERESET == 1))) {
            if (s2[0] == 'y' || s2[0] == 'Y') {
                range.reset = 1;
            } else {
                range.reset = 0;
            }
            notAlreadySet.RANGERESET = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGEOLDIC", s1)) {
        if ((notAlreadySet.RANGEOLDIC || force) ||
            ((mask != NULL) && (mask->RANGEOLDIC == 1))) {
            if (s2[0] == 'y' || s2[0] == 'Y') {
                range.oldic = 1;
            } else {
                range.oldic = 0;
            }

            notAlreadySet.RANGEOLDIC = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGE", s1)) {
        if ((notAlreadySet.RANGE || force) ||
            ((mask != NULL) && (mask->RANGE == 1))) {
            batch_range = atoi(s2);
            notAlreadySet.RANGE = 0;
        }
        return;
    }

    if (load_eqn_msc("NTST", s1)) {
        if ((notAlreadySet.NTST || force) ||
            ((mask != NULL) && (mask->NTST == 1))) {
            auto_ntst = atoi(s2);
            notAlreadySet.NTST = 0;
        }
        return;
    }
    if (load_eqn_msc("NMAX", s1)) {
        if ((notAlreadySet.NMAX || force) ||
            ((mask != NULL) && (mask->NMAX == 1))) {
            auto_nmx = atoi(s2);
            notAlreadySet.NMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("NPR", s1)) {
        if ((notAlreadySet.NPR || force) ||
            ((mask != NULL) && (mask->NPR == 1))) {
            auto_npr = atoi(s2);
            notAlreadySet.NPR = 0;
        }
        return;
    }
    if (load_eqn_msc("NCOL", s1)) {
        if ((notAlreadySet.NCOL || force) ||
            ((mask != NULL) && (mask->NCOL == 1))) {
            auto_ncol = atoi(s2);
            notAlreadySet.NCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("DSMIN", s1)) {
        if ((notAlreadySet.DSMIN || force) ||
            ((mask != NULL) && (mask->DSMIN == 1))) {
            auto_dsmin = atof(s2);
            notAlreadySet.DSMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("DSMAX", s1)) {
        if ((notAlreadySet.DSMAX || force) ||
            ((mask != NULL) && (mask->DSMAX == 1))) {
            auto_dsmax = atof(s2);
            notAlreadySet.DSMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("DS", s1)) {
        if ((notAlreadySet.DS || force) ||
            ((mask != NULL) && (mask->DS == 1))) {
            auto_ds = atof(s2);
            notAlreadySet.DS = 0;
        }

        return;
    }
    if (load_eqn_msc("PARMIN", s1)) {
        if ((notAlreadySet.XMAX || force) ||
            ((mask != NULL) && (mask->XMAX == 1))) {
            auto_rl0 = atof(s2);
            notAlreadySet.XMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("PARMAX", s1)) {
        if ((notAlreadySet.PARMAX || force) ||
            ((mask != NULL) && (mask->PARMAX == 1))) {
            auto_rl1 = atof(s2);
            notAlreadySet.PARMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("NORMMIN", s1)) {
        if ((notAlreadySet.NORMMIN || force) ||
            ((mask != NULL) && (mask->NORMMIN == 1))) {
            auto_a0 = atof(s2);
            notAlreadySet.NORMMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("NORMMAX", s1)) {
        if ((notAlreadySet.NORMMAX || force) ||
            ((mask != NULL) && (mask->NORMMAX == 1))) {
            auto_a1 = atof(s2);
            notAlreadySet.NORMMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("EPSL", s1)) {
        if ((notAlreadySet.EPSL || force) ||
            ((mask != NULL) && (mask->EPSL == 1))) {
            auto_epsl = atof(s2);
            notAlreadySet.EPSL = 0;
        }
        return;
    }

    if (load_eqn_msc("EPSU", s1)) {
        if ((notAlreadySet.EPSU || force) ||
            ((mask != NULL) && (mask->EPSU == 1))) {
            auto_epsu = atof(s2);
            notAlreadySet.EPSU = 0;
        }
        return;
    }
    if (load_eqn_msc("EPSS", s1)) {
        if ((notAlreadySet.EPSS || force) ||
            ((mask != NULL) && (mask->EPSS == 1))) {
            auto_epss = atof(s2);
            notAlreadySet.EPSS = 0;
        }
        return;
    }
    if (load_eqn_msc("RUNNOW", s1)) {
        if ((notAlreadySet.RUNNOW || force) ||
            ((mask != NULL) && (mask->RUNNOW == 1))) {
            RunImmediately = atoi(s2);
            notAlreadySet.RUNNOW = 0;
        }
        return;
    }

    if (load_eqn_msc("SEC", s1)) {
        if ((notAlreadySet.SEC || force) ||
            ((mask != NULL) && (mask->SEC == 1))) {
            SEc = atoi(s2);
            notAlreadySet.SEC = 0;
        }
        return;
    }
    if (load_eqn_msc("UEC", s1)) {
        if ((notAlreadySet.UEC || force) ||
            ((mask != NULL) && (mask->UEC == 1))) {
            UEc = atoi(s2);
            notAlreadySet.UEC = 0;
        }
        return;
    }
    if (load_eqn_msc("SPC", s1)) {
        if ((notAlreadySet.SPC || force) ||
            ((mask != NULL) && (mask->SPC == 1))) {
            SPc = atoi(s2);
            notAlreadySet.SPC = 0;
        }
        return;
    }
    if (load_eqn_msc("UPC", s1)) {
        if ((notAlreadySet.UPC || force) ||
            ((mask != NULL) && (mask->UPC == 1))) {
            UPc = atoi(s2);
            notAlreadySet.UPC = 0;
        }
        return;
    }

    if (load_eqn_msc("AUTOEVAL", s1)) {
        if ((notAlreadySet.AUTOEVAL || force) ||
            ((mask != NULL) && (mask->AUTOEVAL == 1))) {
            f = atoi(s2);
            tabular_set_auto_eval_flags(f);
            notAlreadySet.AUTOEVAL = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOXMAX", s1)) {
        if ((notAlreadySet.AUTOXMAX || force) ||
            ((mask != NULL) && (mask->AUTOXMAX == 1))) {
            auto_xmax = atof(s2);
            notAlreadySet.AUTOXMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOYMAX", s1)) {
        if ((notAlreadySet.AUTOYMAX || force) ||
            ((mask != NULL) && (mask->AUTOYMAX == 1))) {
            auto_ymax = atof(s2);
            notAlreadySet.AUTOYMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOXMIN", s1)) {
        if ((notAlreadySet.AUTOXMIN || force) ||
            ((mask != NULL) && (mask->AUTOXMIN == 1))) {
            auto_xmin = atof(s2);
            notAlreadySet.AUTOXMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOYMIN", s1)) {
        if ((notAlreadySet.AUTOYMIN || force) ||
            ((mask != NULL) && (mask->AUTOYMIN == 1))) {
            auto_ymin = atof(s2);
            notAlreadySet.AUTOYMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOVAR", s1)) {
        if ((notAlreadySet.AUTOVAR || force) ||
            ((mask != NULL) && (mask->AUTOVAR == 1))) {
            browse_find_variable(s2, &i);
            if (i > 0) {
                auto_var = i - 1;
            }
            notAlreadySet.AUTOVAR = 0;
        }
        return;
    }

    // postscript options

    if (load_eqn_msc("PS_FONT", s1)) {
        if ((notAlreadySet.PS_FONT || force) ||
            ((mask != NULL) && (mask->PS_FONT == 1))) {
            strcpy(PS_FONT, s2);
            notAlreadySet.PS_FONT = 0;
        }
        return;
    }

    if (load_eqn_msc("PS_LW", s1)) {
        if ((notAlreadySet.PS_LW || force) ||
            ((mask != NULL) && (mask->PS_LW == 1))) {
            PS_LW = atof(s2);
            notAlreadySet.PS_LW = 0;
        }
        return;
    }

    if (load_eqn_msc("PS_FSIZE", s1)) {
        if ((notAlreadySet.PS_FSIZE || force) ||
            ((mask != NULL) && (mask->PS_FSIZE == 1))) {
            PS_FONTSIZE = atoi(s2);
            notAlreadySet.PS_FSIZE = 0;
        }
        return;
    }

    if (load_eqn_msc("PS_COLOR", s1)) {
        if ((notAlreadySet.PS_COLOR || force) ||
            ((mask != NULL) && (mask->PS_COLOR == 1))) {
            PSColorFlag = atoi(s2);
            PS_Color = PSColorFlag;
            notAlreadySet.PS_COLOR = 0;
        }
        return;
    }
    if (load_eqn_msc("TUTORIAL", s1)) {
        if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
            ggets_plintf("TUTORIAL option must be 0 or 1.\n");
            exit(-1);
        }
        if ((notAlreadySet.TUTORIAL || force) ||
            ((mask != NULL) && (mask->TUTORIAL == 1))) {
            do_tutorial = atoi(s2);
            notAlreadySet.TUTORIAL = 0;
        }
        return;
    }
    if (load_eqn_msc("S1", s1)) {
        if ((notAlreadySet.SLIDER1 || force) ||
            ((mask != NULL) && (mask->SLIDER1 == 1))) {
            strncpy(SLIDER1VAR, s2, 20);
            SLIDER1VAR[19] = '\0';
            notAlreadySet.SLIDER1 = 0;
        }
        return;
    }

    if (load_eqn_msc("S2", s1)) {
        if ((notAlreadySet.SLIDER2 || force) ||
            ((mask != NULL) && (mask->SLIDER2 == 1))) {
            strncpy(SLIDER2VAR, s2, 20);
            SLIDER2VAR[19] = '\0';
            notAlreadySet.SLIDER2 = 0;
        }
        return;
    }
    if (load_eqn_msc("S3", s1)) {
        if ((notAlreadySet.SLIDER3 || force) ||
            ((mask != NULL) && (mask->SLIDER3 == 1))) {
            strncpy(SLIDER3VAR, s2, 20);
            SLIDER3VAR[19] = '\0';
            notAlreadySet.SLIDER3 = 0;
        }
        return;
    }
    if (load_eqn_msc("SLO1", s1)) {
        if ((notAlreadySet.SLIDER1LO || force) ||
            ((mask != NULL) && (mask->SLIDER1LO == 1))) {
            SLIDER1LO = atof(s2);
            notAlreadySet.SLIDER1LO = 0;
        }
        return;
    }

    if (load_eqn_msc("SLO2", s1)) {
        if ((notAlreadySet.SLIDER2LO || force) ||
            ((mask != NULL) && (mask->SLIDER2LO == 1))) {
            SLIDER2LO = atof(s2);
            notAlreadySet.SLIDER2LO = 0;
        }
        return;
    }
    if (load_eqn_msc("SLO3", s1)) {
        if ((notAlreadySet.SLIDER3LO || force) ||
            ((mask != NULL) && (mask->SLIDER3LO == 1))) {
            SLIDER3LO = atof(s2);
            notAlreadySet.SLIDER3LO = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI1", s1)) {
        if ((notAlreadySet.SLIDER1HI || force) ||
            ((mask != NULL) && (mask->SLIDER1HI == 1))) {
            SLIDER1HI = atof(s2);
            notAlreadySet.SLIDER1HI = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI2", s1)) {
        if ((notAlreadySet.SLIDER2HI || force) ||
            ((mask != NULL) && (mask->SLIDER2HI == 1))) {
            SLIDER2HI = atof(s2);
            notAlreadySet.SLIDER2HI = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI3", s1)) {
        if ((notAlreadySet.SLIDER3HI || force) ||
            ((mask != NULL) && (mask->SLIDER3HI == 1))) {
            SLIDER3HI = atof(s2);
            notAlreadySet.SLIDER3HI = 0;
        }
        return;
    }

    /* postprocessing options
       This is rally only relevant for batch jobs as it
       writes files then
    */

    if (load_eqn_msc("POSTPROCESS", s1)) {
        if ((notAlreadySet.POSTPROCESS || force) ||
            ((mask != NULL) && (mask->POSTPROCESS == 1))) {
            post_process = atoi(s2);
            notAlreadySet.POSTPROCESS = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTLO", s1)) {
        if ((notAlreadySet.HISTLO || force) ||
            ((mask != NULL) && (mask->HISTLO == 1))) {
            hist_inf.xlo = atof(s2);
            notAlreadySet.HISTLO = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTHI", s1)) {
        if ((notAlreadySet.HISTHI || force) ||
            ((mask != NULL) && (mask->HISTHI == 1))) {
            hist_inf.xhi = atof(s2);
            notAlreadySet.HISTHI = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTBINS", s1)) {
        if ((notAlreadySet.HISTBINS || force) ||
            ((mask != NULL) && (mask->HISTBINS == 1))) {
            hist_inf.nbins = atoi(s2);
            notAlreadySet.HISTBINS = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTCOL", s1)) {
        if ((notAlreadySet.HISTCOL || force) ||
            ((mask != NULL) && (mask->HISTCOL == 1))) {
            browse_find_variable(s2, &i);
            if (i > (-1)) {
                hist_inf.col = i;
            }
            notAlreadySet.HISTCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTLO2", s1)) {
        if ((notAlreadySet.HISTLO2 || force) ||
            ((mask != NULL) && (mask->HISTLO2 == 1))) {
            hist_inf.ylo = atof(s2);
            notAlreadySet.HISTLO2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTHI2", s1)) {
        if ((notAlreadySet.HISTHI2 || force) ||
            ((mask != NULL) && (mask->HISTHI2 == 1))) {
            hist_inf.yhi = atof(s2);
            notAlreadySet.HISTHI2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTBINS2", s1)) {
        if ((notAlreadySet.HISTBINS2 || force) ||
            ((mask != NULL) && (mask->HISTBINS2 == 1))) {
            hist_inf.nbins2 = atoi(s2);
            notAlreadySet.HISTBINS2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTCOL2", s1)) {
        if ((notAlreadySet.HISTCOL2 || force) ||
            ((mask != NULL) && (mask->HISTCOL2 == 1))) {
            browse_find_variable(s2, &i);
            if (i > (-1)) {
                hist_inf.col2 = i;
            }
            notAlreadySet.HISTCOL2 = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECCOL", s1)) {
        if ((notAlreadySet.SPECCOL || force) ||
            ((mask != NULL) && (mask->SPECCOL == 1))) {
            browse_find_variable(s2, &i);
            if (i > (-1)) {
                spec_col = i;
            }
            notAlreadySet.SPECCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECCOL2", s1)) {
        if ((notAlreadySet.SPECCOL2 || force) ||
            ((mask != NULL) && (mask->SPECCOL2 == 1))) {
            browse_find_variable(s2, &i);
            if (i > (-1)) {
                spec_col2 = i;
            }
            notAlreadySet.SPECCOL2 = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECWIDTH", s1)) {
        if ((notAlreadySet.SPECWIDTH || force) ||
            ((mask != NULL) && (mask->SPECWIDTH == 1))) {
            spec_wid = atoi(s2);
            notAlreadySet.SPECWIDTH = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECWIN", s1)) {
        if ((notAlreadySet.SPECWIN || force) ||
            ((mask != NULL) && (mask->SPECWIN == 1))) {
            spec_win = atoi(s2);
            notAlreadySet.SPECWIN = 0;
        }
        return;
    }

    if (load_eqn_msc("DFGRID", s1)) {
        if ((notAlreadySet.DFGRID || force) ||
            ((mask != NULL) && (mask->DFGRID == 1))) {
            DF_GRID = atoi(s2);
            notAlreadySet.DFGRID = 0;
        }
        return;
    }
    if (load_eqn_msc("DFDRAW", s1)) {
        if ((notAlreadySet.DFBATCH || force) ||
            ((mask != NULL) && (mask->DFBATCH == 1))) {
            DFBatch = atoi(s2);
            notAlreadySet.DFBATCH = 0;
        }
        return;
    }
    if (load_eqn_msc("NCDRAW", s1)) {
        if ((notAlreadySet.NCBATCH || force) ||
            ((mask != NULL) && (mask->NCBATCH == 1))) {
            NCBatch = atoi(s2);
            notAlreadySet.NCBATCH = 0;
        }
        return;
    }

    // colorize customizing !!
    if (load_eqn_msc("COLORVIA", s1)) {
        if ((notAlreadySet.COLORVIA || force) ||
            ((mask != NULL) && (mask->COLORVIA == 1))) {
            strcpy(ColorVia, s2);
        }
        notAlreadySet.COLORVIA = 0;
        return;
    }
    if (load_eqn_msc("COLORIZE", s1)) {
        if ((notAlreadySet.COLORIZE || force) ||
            ((mask != NULL) && (mask->COLORIZE == 1))) {
            ColorizeFlag = atoi(s2);
        }
        notAlreadySet.COLORIZE = 0;
        return;
    }
    if (load_eqn_msc("COLORLO", s1)) {
        if ((notAlreadySet.COLORLO || force) ||
            ((mask != NULL) && (mask->COLORLO == 1))) {
            ColorViaLo = atof(s2);
        }
        notAlreadySet.COLORLO = 0;
        return;
    }
    if (load_eqn_msc("COLORHI", s1)) {
        if ((notAlreadySet.COLORHI || force) ||
            ((mask != NULL) && (mask->COLORHI == 1))) {
            ColorViaHi = atof(s2);
        }
        notAlreadySet.COLORHI = 0;
        return;
    }

    ggets_plintf("!! Option %s not recognized\n", s1);
}
