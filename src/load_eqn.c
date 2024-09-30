#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include "xmalloc.h"

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
int32 run_immediately = 0;

int32 IX_PLT[10];
int32 IY_PLT[10];
int32 IZ_PLT[10];
int32 NPltV;
int32 multi_win = 0;
double x_lo[10];
double y_lo[10];
double x_hi[10];
double y_hi[10];
int32 start_line_type = 1;
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
int32 stor_flag;
int32 in_flag;
int32 max_stor;
double x_3d[2];
double y_3d[2];
double z_3d[2];
int32 ix_plt;
int32 iy_plt;
int32 iz_plt;
int32 axes;
int32 TIMPLOT;
int32 PLOT_3D;
double my_xlo;
double my_ylo;
double my_xhi;
double my_yhi;
double torus_period = 6.2831853071795864770;
int32 TORUS = 0;
int32 n_equations;
char options[100];

/*   Numerical stuff ....   */

double delta_t;
double TEND;
double T0;
double TRANS;
double NULL_ERR;
double evec_err;
double newt_err;
double bound;
double delay;
double TOLER;
double atoler;
double h_min;
double h_max;
double bvp_eps;
double bvp_tol;

double POIPLN;

int32 euler_max_iter;
double euler_tol;
int32 nmesh;
int32 njmp;
int32 METHOD;
int32 color_flag;
int32 NC_ITER;
int32 evec_iter;
int32 bpv_maxit;
int32 bvp_flag;

int32 POIMAP;
int32 POIVAR;
int32 POISGN;
int32 SOS;
int32 FFT;
int32 NULL_HERE;
int32 POIEXT;
int32 hist;
int32 h_var;
int32 hist_ind;
int32 forever;

/*  control of range stuff  */

int32 end_sing;
int32 shoot;
int32 par_fol;

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
    lunch_io_double(&torus_period, fp, f, "Torus period");
    if (TORUS) {
        for (int32 i = 0; i < n_equations; i++) {
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
    if (got_file == 1 && (std == 0) && (dp = (struct dirent *)opendir(this_file)) != NULL) {
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

    if (strlen(user_draw_win_color) == 0) {
        sprintf(user_draw_win_color, "#%s", "FFFFFF");
    }

    if (user_gradients < 0) {
        user_gradients = 1;
    }
    return;
}

void
load_eqn_set_all_vals(void) {
    FILE *fp;

    if (not_already_set.TIMEPLOT) {
        TIMPLOT = 1;
        not_already_set.TIMEPLOT = 0;
    }
    if (not_already_set.forever) {
        forever = 0;
        not_already_set.forever = 0;
    }
    if (not_already_set.bvp_tol) {
        bvp_tol = 1.e-5;
        not_already_set.bvp_tol = 0;
    }
    if (not_already_set.bvp_eps) {
        bvp_eps = 1.e-5;
        not_already_set.bvp_eps = 0;
    }
    if (not_already_set.bpv_maxit) {
        bpv_maxit = 20;
        not_already_set.bpv_maxit = 0;
    }
    if (not_already_set.bvp_flag) {
        bvp_flag = 0;
        not_already_set.bvp_flag = 0;
    }
    if (not_already_set.nmesh) {
        nmesh = 40;
        not_already_set.nmesh = 0;
    }
    if (not_already_set.NOUT) {
        njmp = 1;
        not_already_set.NOUT = 0;
    }
    if (not_already_set.SOS) {
        SOS = 0;
        not_already_set.SOS = 0;
    }
    if (not_already_set.FFT) {
        FFT = 0;
        not_already_set.FFT = 0;
    }
    if (not_already_set.hist) {
        hist = 0;
        not_already_set.hist = 0;
    }
    if (not_already_set.plt_fmt_flag) {
        plt_fmt_flag = 0;
        not_already_set.plt_fmt_flag = 0;
    }
    if (not_already_set.axes) {
        axes = 0;
        not_already_set.axes = 0;
    }
    if (not_already_set.TOLER) {
        TOLER = 0.001;
        not_already_set.TOLER = 0;
    }
    if (not_already_set.atoler) {
        atoler = 0.001;
        not_already_set.atoler = 0;
    }
    if (not_already_set.euler_max_iter) {
        euler_max_iter = 10;
        not_already_set.euler_max_iter = 0;
    }
    if (not_already_set.euler_tol) {
        euler_tol = 1.e-7;
        not_already_set.euler_tol = 0;
    }
    if (not_already_set.delay) {
        delay = 0.0;
        not_already_set.delay = 0;
    }
    if (not_already_set.DTMIN) {
        h_min = 1e-12;
        not_already_set.DTMIN = 0;
    }
    if (not_already_set.evec_iter) {
        evec_iter = 100;
        not_already_set.evec_iter = 0;
    }
    if (not_already_set.evec_err) {
        evec_err = .001;
        not_already_set.evec_err = 0;
    }
    if (not_already_set.NULL_ERR) {
        NULL_ERR = .001;
        not_already_set.NULL_ERR = 0;
    }
    if (not_already_set.newt_err) {
        newt_err = .001;
        not_already_set.newt_err = 0;
    }
    if (not_already_set.NULL_HERE) {
        NULL_HERE = 0;
        not_already_set.NULL_HERE = 0;
    }
    del_stab_flag = DFNORMAL;
    if (not_already_set.DTMAX) {
        h_max = 1.000;
        not_already_set.DTMAX = 0;
    }
    if (not_already_set.POIMAP) {
        POIMAP = 0;
        not_already_set.POIMAP = 0;
    }
    if (not_already_set.POIVAR) {
        POIVAR = 1;
        not_already_set.POIVAR = 0;
    }
    if (not_already_set.POIEXT) {
        POIEXT = 0;
        not_already_set.POIEXT = 0;
    }
    if (not_already_set.POISGN) {
        POISGN = 1;
        not_already_set.POISGN = 0;
    }
    if (not_already_set.POIPLN) {
        POIPLN = 0.0;
        not_already_set.POIPLN = 0;
    }

    storind = 0;
    mov_ind = 0;

    stor_flag = 0;

    in_flag = 0;
    solver = odesol_rung_kut;
    PLOT_3D = 0;
    if (not_already_set.METHOD) {
        METHOD = 3;
        not_already_set.METHOD = 0;
    }
    if (not_already_set.XLO) {
        my_xlo = 0.0;
        x_3d[0] = my_xlo;
        not_already_set.XLO = 0;
        not_already_set.XMIN = 0;
    }
    if (not_already_set.XHI) {
        my_xhi = 20.0;
        x_3d[1] = my_xhi;
        not_already_set.XHI = 0;
        not_already_set.XMAX = 0;
    }
    if (not_already_set.YLO) {
        my_ylo = -1;
        y_3d[0] = my_ylo;
        not_already_set.YLO = 0;
        not_already_set.YMIN = 0;
    }
    if (not_already_set.YHI) {
        my_yhi = 1;
        y_3d[0] = my_yhi;
        not_already_set.YHI = 0;
        not_already_set.YMAX = 0;
    }

    if (not_already_set.bound) {
        bound = 100;
        not_already_set.bound = 0;
    }
    if (not_already_set.max_stor) {
        max_stor = 5000;
        not_already_set.max_stor = 0;
    }

    if (not_already_set.T0) {
        T0 = 0.0;
        not_already_set.T0 = 0;
    }
    if (not_already_set.TRANS) {
        TRANS = 0.0;
        not_already_set.TRANS = 0;
    }
    if (not_already_set.DT) {
        delta_t = .05;
        not_already_set.DT = 0;
    }

    if (not_already_set.XMIN) {
        x_3d[0] = -12;
        not_already_set.XMIN = 0;
        not_already_set.XLO = 0;
    }
    if (not_already_set.XMAX) {
        x_3d[1] = 12;
        not_already_set.XMAX = 0;
        not_already_set.XHI = 0;
    }
    if (not_already_set.YMIN) {
        y_3d[0] = -12;
        not_already_set.YMIN = 0;
        not_already_set.YLO = 0;
    }
    if (not_already_set.YMAX) {
        y_3d[1] = 12;
        not_already_set.YMAX = 0;
        not_already_set.YHI = 0;
    }
    if (not_already_set.ZMIN) {
        z_3d[0] = -12;
        not_already_set.ZMIN = 0;
    }
    if (not_already_set.ZMAX) {
        z_3d[1] = 12;
        not_already_set.ZMAX = 0;
    }

    if (not_already_set.TEND) {
        TEND = 20.00;
        not_already_set.TEND = 0;
    }
    if (not_already_set.ix_plt) {
        ix_plt = 0;
        not_already_set.ix_plt = 0;
    }
    if (not_already_set.iy_plt) {
        iy_plt = 1;
        not_already_set.iy_plt = 0;
    }
    if (not_already_set.iz_plt) {
        iz_plt = 1;
        not_already_set.iz_plt = 0;
    }

    if (not_already_set.NPLOT) {
        if (n_equations > 2) {
            if (not_already_set.iz_plt) {
                iz_plt = 2;
            }
        }
        NPltV = 1;
        for (int32 i = 0; i < 10; i++) {
            IX_PLT[i] = ix_plt;
            IY_PLT[i] = iy_plt;
            IZ_PLT[i] = iz_plt;
            x_lo[i] = 0;
            y_lo[i] = -1;
            x_hi[i] = 20;
            y_hi[i] = 1;
        }
        not_already_set.NPLOT = 0;
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

    if (iz_plt > n_equations) {
        iz_plt = n_equations;
    }
    if (iy_plt > n_equations) {
        iy_plt = n_equations;
    }
    if (ix_plt == 0 || iy_plt == 0) {
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
    if (my_xlo >= my_xhi) {
        my_xlo = -2.0;
        my_xhi = 2.0;
    }
    if (my_ylo >= my_yhi) {
        my_ylo = -2.0;
        my_yhi = 2.0;
    }
    if (axes < 5) {
        x_3d[0] = my_xlo;
        y_3d[0] = my_ylo;
        x_3d[1] = my_xhi;
        y_3d[1] = my_yhi;
    }
    storage_init_stor(max_stor, n_equations + 1);
    if (axes >= 5) {
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
    if (not_already_set.BIG_FONT_NAME) {
        strcpy(font_name_big, ptr);
        not_already_set.BIG_FONT_NAME = 0;
    }

    fgets(bob, 80, fp);
    ptr = form_ode_get_first(bob, " ");
    if (not_already_set.SMALL_FONT_NAME) {
        strcpy(font_name_small, ptr);
        not_already_set.SMALL_FONT_NAME = 0;
    }

    if (not_already_set.paper_white) {
        load_eqn_fil_int(fp, &paper_white);
        not_already_set.paper_white = 0;
    }
    if (not_already_set.ix_plt) {
        load_eqn_fil_int(fp, &ix_plt);
        not_already_set.ix_plt = 0;
    }
    if (not_already_set.iy_plt) {
        load_eqn_fil_int(fp, &iy_plt);
        not_already_set.iy_plt = 0;
    }
    if (not_already_set.iz_plt) {
        load_eqn_fil_int(fp, &iz_plt);
        not_already_set.iz_plt = 0;
    }
    if (not_already_set.axes) {
        load_eqn_fil_int(fp, &axes);
        not_already_set.paper_white = 0;
    }
    if (not_already_set.NOUT) {
        load_eqn_fil_int(fp, &njmp);
        not_already_set.NOUT = 0;
    }
    if (not_already_set.nmesh) {
        load_eqn_fil_int(fp, &nmesh);
        not_already_set.nmesh = 0;
    }
    if (not_already_set.METHOD) {
        load_eqn_fil_int(fp, &METHOD);
        not_already_set.METHOD = 0;
    }

    if (not_already_set.TIMEPLOT) {
        load_eqn_fil_int(fp, &TIMPLOT);
        not_already_set.TIMEPLOT = 0;
    }
    if (not_already_set.max_stor) {
        load_eqn_fil_int(fp, &max_stor);
        not_already_set.max_stor = 0;
    }
    if (not_already_set.TEND) {
        load_eqn_fil_flt(fp, &TEND);
        not_already_set.TEND = 0;
    }
    if (not_already_set.DT) {
        load_eqn_fil_flt(fp, &delta_t);
        not_already_set.DT = 0;
    }
    if (not_already_set.T0) {
        load_eqn_fil_flt(fp, &T0);
        not_already_set.T0 = 0;
    }
    if (not_already_set.TRANS) {
        load_eqn_fil_flt(fp, &TRANS);
        not_already_set.TRANS = 0;
    }
    if (not_already_set.bound) {
        load_eqn_fil_flt(fp, &bound);
        not_already_set.bound = 0;
    }
    if (not_already_set.DTMIN) {
        load_eqn_fil_flt(fp, &h_min);
        not_already_set.DTMIN = 0;
    }
    if (not_already_set.DTMAX) {
        load_eqn_fil_flt(fp, &h_max);
        not_already_set.DTMIN = 0;
    }
    if (not_already_set.TOLER) {
        load_eqn_fil_flt(fp, &TOLER);
        not_already_set.TOLER = 0;
    }
    if (not_already_set.delay) {
        load_eqn_fil_flt(fp, &delay);
        not_already_set.delay = 0;
    }
    if (not_already_set.XLO) {
        load_eqn_fil_flt(fp, &my_xlo);
        not_already_set.XLO = 0;
    }
    if (not_already_set.XHI) {
        load_eqn_fil_flt(fp, &my_xhi);
        not_already_set.XHI = 0;
    }
    if (not_already_set.YLO) {
        load_eqn_fil_flt(fp, &my_ylo);
        not_already_set.YLO = 0;
    }
    if (not_already_set.YHI) {
        load_eqn_fil_flt(fp, &my_yhi);
        not_already_set.YHI = 0;
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
        ggets_plintf(" %s not added -- too many must be less than %d \n", name, MAX_INTERN_SET);
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
    ggets_plintf(" added %s doing %s \n", intern_set[j].name, intern_set[j].does);
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
                        if (strlen(user_draw_win_color)!=0)
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
load_eqn_set_internopts_xpprc_and_cli(void) {
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
    but prevents overwriting across cli, .xpprc etc.
    */
    tempNAS = xmalloc(sizeof(*tempNAS));
    *tempNAS = not_already_set;

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
        if (override_quiet == 0)  // Will be 1 if -quiet was specified on the command line.
        {
            xpp_verbose = (atoi(s2) == 0);
        }
        return;
    }
    if (load_eqn_msc("LOGFILE", s1)) {
        if (override_logfile == 0)  // Will be 1 if -logfile was specified on the command line.
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
        if ((not_already_set.BIG_FONT_NAME || force) ||
            ((mask != NULL) && (mask->BIG_FONT_NAME == 1))) {
            strcpy(font_name_big, s2);
            not_already_set.BIG_FONT_NAME = 0;
        }
        return;
    }
    if ((load_eqn_msc("SMALLFONT", s1)) || (load_eqn_msc("SMALL", s1))) {
        if ((not_already_set.SMALL_FONT_NAME || force) ||
            ((mask != NULL) && (mask->SMALL_FONT_NAME == 1))) {
            strcpy(font_name_small, s2);
            not_already_set.SMALL_FONT_NAME = 0;
        }
        return;
    }
    if (load_eqn_msc("FORECOLOR", s1)) {
        if ((not_already_set.user_black || force) || ((mask != NULL) && (mask->user_black == 1))) {
            sprintf(user_black, "#%s", s2);
            not_already_set.user_black = 0;
        }
        return;
    }
    if (load_eqn_msc("BACKCOLOR", s1)) {
        if ((not_already_set.user_white || force) || ((mask != NULL) && (mask->user_white == 1))) {
            sprintf(user_white, "#%s", s2);
            not_already_set.user_white = 0;
        }
        return;
    }
    if (load_eqn_msc("MWCOLOR", s1)) {
        if ((not_already_set.user_main_win_color || force) ||
            ((mask != NULL) && (mask->user_main_win_color == 1))) {
            sprintf(user_main_win_color, "#%s", s2);
            not_already_set.user_main_win_color = 0;
        }
        return;
    }
    if (load_eqn_msc("DWCOLOR", s1)) {
        if ((not_already_set.user_draw_win_color || force) ||
            ((mask != NULL) && (mask->user_draw_win_color == 1))) {
            sprintf(user_draw_win_color, "#%s", s2);
            not_already_set.user_draw_win_color = 0;
        }
        return;
    }
    if (load_eqn_msc("GRADS", s1)) {
        if ((not_already_set.user_gradients || force) ||
            ((mask != NULL) && (mask->user_gradients == 1))) {
            if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
                ggets_plintf("GRADS option must be 0 or 1.\n");
                exit(-1);
            }
            user_gradients = atoi(s2);
            not_already_set.user_gradients = 0;
        }
        return;
    }

    if (load_eqn_msc("PLOTFMT", s1)) {
        if ((not_already_set.PLOTFORMAT || force) || ((mask != NULL) && (mask->PLOTFORMAT == 1))) {
            strcpy(plot_format, s2);
            not_already_set.PLOTFORMAT = 0;
        }
        return;
    }

    if (load_eqn_msc("BACKIMAGE", s1)) {
        if ((not_already_set.user_bg_bitmap || force) ||
            ((mask != NULL) && (mask->user_bg_bitmap == 1))) {
            strcpy(user_bg_bitmap, s2);
            not_already_set.user_bg_bitmap = 0;
        }
        return;
    }
    if (load_eqn_msc("WIDTH", s1)) {
        if ((not_already_set.user_min_width || force) ||
            ((mask != NULL) && (mask->user_min_width == 1))) {
            user_min_width = atoi(s2);
            not_already_set.user_min_width = 0;
        }
        return;
    }
    if (load_eqn_msc("HEIGHT", s1)) {
        if ((not_already_set.user_min_height || force) ||
            ((mask != NULL) && (mask->user_min_height == 1))) {
            user_min_height = atoi(s2);
            not_already_set.user_min_height = 0;
        }
        return;
    }
    if (load_eqn_msc("YNC", s1)) {
        if ((not_already_set.YNullColor || force) || ((mask != NULL) && (mask->YNullColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                YNullColor = i;
            }
            not_already_set.YNullColor = 0;
        }
        return;
    }
    if (load_eqn_msc("XNC", s1)) {
        if ((not_already_set.XNullColor || force) || ((mask != NULL) && (mask->XNullColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                XNullColor = i;
                not_already_set.XNullColor = 0;
            }
        }
        return;
    }

    if (load_eqn_msc("SMC", s1)) {
        if ((not_already_set.stable_manifold_color || force) ||
            ((mask != NULL) && (mask->stable_manifold_color == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                stable_manifold_color = i;
                not_already_set.stable_manifold_color = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("UMC", s1)) {
        if ((not_already_set.UnstableManifoldColor || force) ||
            ((mask != NULL) && (mask->UnstableManifoldColor == 1))) {
            i = atoi(s2);
            if (i > -1 && i < 11) {
                UnstableManifoldColor = i;
                not_already_set.UnstableManifoldColor = 0;
            }
        }
        return;
    }

    if (load_eqn_msc("LT", s1)) {
        if ((not_already_set.start_line_type || force) ||
            ((mask != NULL) && (mask->start_line_type == 1))) {
            i = atoi(s2);
            if (i < 2 && i > -6) {
                start_line_type = i;
                graphics_reset_all_line_type();
                not_already_set.start_line_type = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("SEED", s1)) {
        if ((not_already_set.rand_seed || force) || ((mask != NULL) && (mask->rand_seed == 1))) {
            i = atoi(s2);
            if (i >= 0) {
                rand_seed = i;
                markov_nsrand48(rand_seed);
                not_already_set.rand_seed = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("BACK", s1)) {
        if ((not_already_set.paper_white || force) || ((mask != NULL) && (mask->paper_white == 1))) {
            if (s2[0] == 'w' || s2[0] == 'W') {
                paper_white = 1;
            } else {
                paper_white = 0;
            }
            not_already_set.paper_white = 0;
        }
        return;
    }
    if (load_eqn_msc("COLORMAP", s1)) {
        if ((not_already_set.COLORMAP || force) || ((mask != NULL) && (mask->COLORMAP == 1))) {
            i = atoi(s2);
            if (i < 7) {
                custom_color = i;
            }
            not_already_set.COLORMAP = 0;
        }
        return;
    }
    if (load_eqn_msc("NPLOT", s1)) {
        if ((not_already_set.NPLOT || force) || ((mask != NULL) && (mask->NPLOT == 1))) {
            NPltV = atoi(s2);
            not_already_set.NPLOT = 0;
        }
        return;
    }

    if (load_eqn_msc("DLL_LIB", s1)) {
        if ((not_already_set.DLL_LIB || force) || ((mask != NULL) && (mask->DLL_LIB == 1))) {
            sprintf(dll_lib, "%s", s2);
            dll_flag += 1;
            not_already_set.DLL_LIB = 0;
        }
        return;
    }
    if (load_eqn_msc("DLL_FUN", s1)) {
        if ((not_already_set.DLL_FUN || force) || ((mask != NULL) && (mask->DLL_FUN == 1))) {
            sprintf(dll_fun, "%s", s2);
            dll_flag += 2;
            not_already_set.DLL_FUN = 0;
        }
        return;
    }
    // can now initialize several plots
    if (load_eqn_msc("SIMPLOT", s1)) {
        simul_plot_flag = 1;
        return;
    }
    if (load_eqn_msc("MULTIWIN", s1)) {
        multi_win = 1;
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
            browser_find_variable(s2, &i);
            if (i > -1) {
                IX_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(yy, s1)) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                IY_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(zz, s1)) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                IZ_PLT[j] = i;
            }
            return;
        }
        if (load_eqn_msc(xxh, s1)) {
            x_hi[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(xxl, s1)) {
            x_lo[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(yyh, s1)) {
            y_hi[j] = atof(s2);
            return;
        }
        if (load_eqn_msc(yyl, s1)) {
            y_lo[j] = atof(s2);
            return;
        }
    }
    if (load_eqn_msc("XP", s1)) {
        if ((not_already_set.XP || force) || ((mask != NULL) && (mask->XP == 1))) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                ix_plt = i;
            }
            not_already_set.XP = 0;
            not_already_set.ix_plt = 0;
        }
        return;
    }
    if (load_eqn_msc("YP", s1)) {
        if ((not_already_set.YP || force) || ((mask != NULL) && (mask->YP == 1))) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                iy_plt = i;
            }
            not_already_set.YP = 0;
            not_already_set.iy_plt = 0;
        }
        return;
    }
    if (load_eqn_msc("ZP", s1)) {
        if ((not_already_set.ZP || force) || ((mask != NULL) && (mask->ZP == 1))) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                iz_plt = i;
            }

            not_already_set.ZP = 0;
            not_already_set.iz_plt = 0;
        }
        return;
    }
    if (load_eqn_msc("axes", s1)) {
        if ((not_already_set.axes || force) || ((mask != NULL) && (mask->axes == 1))) {
            if (s2[0] == '3') {
                axes = 5;
            } else {
                axes = 0;
            }

            not_already_set.axes = 0;
        }
        return;
    }

    if (load_eqn_msc("njmp", s1)) {
        if ((not_already_set.NOUT || force) || ((mask != NULL) && (mask->NOUT == 1))) {
            njmp = atoi(s2);
            not_already_set.NOUT = 0;
        }
        return;
    }
    if (load_eqn_msc("NOUT", s1)) {
        if ((not_already_set.NOUT || force) || ((mask != NULL) && (mask->NOUT == 1))) {
            njmp = atoi(s2);
            not_already_set.NOUT = 0;
        }
        return;
    }
    if (load_eqn_msc("nmesh", s1)) {
        if ((not_already_set.nmesh || force) || ((mask != NULL) && (mask->nmesh == 1))) {
            nmesh = atoi(s2);
            not_already_set.nmesh = 0;
        }
        return;
    }
    if (load_eqn_msc("METH", s1)) {
        if ((not_already_set.METHOD || force) || ((mask != NULL) && (mask->METHOD == 1))) {
            for (i = 0; i < 15; i++) {
                if (s2[0] == mkey[i] || s2[0] == Mkey[i]) {
                    METHOD = i;
                }
            }

            not_already_set.METHOD = 0;
        }
        return;
    }
    if (load_eqn_msc("VMAXPTS", s1)) {
        if ((not_already_set.VMAXPTS || force) || ((mask != NULL) && (mask->VMAXPTS == 1))) {
            max_points = atoi(s2);
            not_already_set.VMAXPTS = 0;
        }
        return;
    }
    if (load_eqn_msc("max_stor", s1)) {
        if ((not_already_set.max_stor || force) || ((mask != NULL) && (mask->max_stor == 1))) {
            max_stor = atoi(s2);
            not_already_set.max_stor = 0;
        }
        return;
    }
    if (load_eqn_msc("TOR_PER", s1)) {
        if ((not_already_set.TOR_PER || force) || ((mask != NULL) && (mask->TOR_PER == 1))) {
            torus_period = atof(s2);
            TORUS = 1;
            not_already_set.TOR_PER = 0;
        }
        return;
    }
    if (load_eqn_msc("JAC_EPS", s1)) {
        if ((not_already_set.JAC_EPS || force) || ((mask != NULL) && (mask->JAC_EPS == 1))) {
            newt_err = atof(s2);
            not_already_set.JAC_EPS = 0;
        }
        return;
    }
    if (load_eqn_msc("NEWT_TOL", s1)) {
        if ((not_already_set.NEWT_TOL || force) || ((mask != NULL) && (mask->NEWT_TOL == 1))) {
            evec_err = atof(s2);
            not_already_set.NEWT_TOL = 0;
        }
        return;
    }
    if (load_eqn_msc("NEWT_ITER", s1)) {
        if ((not_already_set.NEWT_ITER || force) || ((mask != NULL) && (mask->NEWT_ITER == 1))) {
            evec_iter = atoi(s2);
            not_already_set.NEWT_ITER = 0;
        }
        return;
    }
    if (load_eqn_msc("FOLD", s1)) {
        if ((not_already_set.FOLD || force) || ((mask != NULL) && (mask->FOLD == 1))) {
            browser_find_variable(s2, &i);
            if (i > 0) {
                itor[i - 1] = 1;
                TORUS = 1;
            }
        }
        return;
    }
    if (load_eqn_msc("TOTAL", s1)) {
        if ((not_already_set.TEND || force) || ((mask != NULL) && (mask->TEND == 1))) {
            TEND = atof(s2);
            not_already_set.TEND = 0;
        }
        return;
    }
    if (load_eqn_msc("DTMIN", s1)) {
        if ((not_already_set.DTMIN || force) || ((mask != NULL) && (mask->DTMIN == 1))) {
            h_min = atof(s2);
            not_already_set.DTMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("DTMAX", s1)) {
        if ((not_already_set.DTMAX || force) || ((mask != NULL) && (mask->DTMAX == 1))) {
            h_max = atof(s2);
            not_already_set.DTMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("DT", s1)) {
        if ((not_already_set.DT || force) || ((mask != NULL) && (mask->DT == 1))) {
            delta_t = atof(s2);
            not_already_set.DT = 0;
        }
        return;
    }
    if (load_eqn_msc("T0", s1)) {
        if ((not_already_set.T0 || force) || ((mask != NULL) && (mask->T0 == 1))) {
            T0 = atof(s2);
            not_already_set.T0 = 0;
        }
        return;
    }
    if (load_eqn_msc("TRANS", s1)) {
        if ((not_already_set.TRANS || force) || ((mask != NULL) && (mask->TRANS == 1))) {
            TRANS = atof(s2);
            not_already_set.TRANS = 0;
        }
        return;
    }
    if (load_eqn_msc("bound", s1)) {
        if ((not_already_set.bound || force) || ((mask != NULL) && (mask->bound == 1))) {
            bound = atof(s2);
            not_already_set.bound = 0;
        }
        return;
    }
    if (load_eqn_msc("ATOL", s1)) {
        if ((not_already_set.atoler || force) || ((mask != NULL) && (mask->atoler == 1))) {
            atoler = atof(s2);
            not_already_set.atoler = 0;
        }
        return;
    }
    if (load_eqn_msc("TOL", s1)) {
        if ((not_already_set.TOLER || force) || ((mask != NULL) && (mask->TOLER == 1))) {
            TOLER = atof(s2);
            not_already_set.TOLER = 0;
        }
        return;
    }

    if (load_eqn_msc("delay", s1)) {
        if ((not_already_set.delay || force) || ((mask != NULL) && (mask->delay == 1))) {
            delay = atof(s2);
            not_already_set.delay = 0;
        }
        return;
    }
    if (load_eqn_msc("BANDUP", s1)) {
        if ((not_already_set.BANDUP || force) || ((mask != NULL) && (mask->BANDUP == 1))) {
            cv_bandflag = 1;
            cv_bandupper = atoi(s2);
            not_already_set.BANDUP = 0;
        }
        return;
    }
    if (load_eqn_msc("BANDLO", s1)) {
        if ((not_already_set.BANDLO || force) || ((mask != NULL) && (mask->BANDLO == 1))) {
            cv_bandflag = 1;
            cv_bandlower = atoi(s2);
            not_already_set.BANDLO = 0;
        }
        return;
    }

    if (load_eqn_msc("PHI", s1)) {
        if ((not_already_set.PHI || force) || ((mask != NULL) && (mask->PHI == 1))) {
            PHI0 = atof(s2);
            not_already_set.PHI = 0;
        }
        return;
    }
    if (load_eqn_msc("THETA", s1)) {
        if ((not_already_set.THETA || force) || ((mask != NULL) && (mask->THETA == 1))) {
            THETA0 = atof(s2);
            not_already_set.THETA = 0;
        }
        return;
    }
    if (load_eqn_msc("XLO", s1)) {
        if ((not_already_set.XLO || force) || ((mask != NULL) && (mask->XLO == 1))) {
            my_xlo = atof(s2);
            not_already_set.XLO = 0;
        }
        return;
    }
    if (load_eqn_msc("YLO", s1)) {
        if ((not_already_set.YLO || force) || ((mask != NULL) && (mask->YLO == 1))) {
            my_ylo = atof(s2);
            not_already_set.YLO = 0;
        }
        return;
    }

    if (load_eqn_msc("XHI", s1)) {
        if ((not_already_set.XHI || force) || ((mask != NULL) && (mask->XHI == 1))) {
            my_xhi = atof(s2);
            not_already_set.XHI = 0;
        }
        return;
    }
    if (load_eqn_msc("YHI", s1)) {
        if ((not_already_set.YHI || force) || ((mask != NULL) && (mask->YHI == 1))) {
            my_yhi = atof(s2);
            not_already_set.YHI = 0;
        }
        return;
    }
    if (load_eqn_msc("XMAX", s1)) {
        if ((not_already_set.XMAX || force) || ((mask != NULL) && (mask->XMAX == 1))) {
            x_3d[1] = atof(s2);
            not_already_set.XMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("YMAX", s1)) {
        if ((not_already_set.YMAX || force) || ((mask != NULL) && (mask->YMAX == 1))) {
            y_3d[1] = atof(s2);
            not_already_set.YMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("ZMAX", s1)) {
        if ((not_already_set.ZMAX || force) || ((mask != NULL) && (mask->ZMAX == 1))) {
            z_3d[1] = atof(s2);
            not_already_set.ZMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("XMIN", s1)) {
        if ((not_already_set.XMIN || force) || ((mask != NULL) && (mask->XMIN == 1))) {
            x_3d[0] = atof(s2);
            not_already_set.XMIN = 0;
            if ((not_already_set.XLO || force) || ((mask != NULL) && (mask->XLO == 1))) {
                my_xlo = atof(s2);
                not_already_set.XLO = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("YMIN", s1)) {
        if ((not_already_set.YMIN || force) || ((mask != NULL) && (mask->YMIN == 1))) {
            y_3d[0] = atof(s2);
            not_already_set.YMIN = 0;
            if ((not_already_set.YLO || force) || ((mask != NULL) && (mask->YLO == 1))) {
                my_ylo = atof(s2);
                not_already_set.YLO = 0;
            }
        }
        return;
    }
    if (load_eqn_msc("ZMIN", s1)) {
        if ((not_already_set.ZMIN || force) || ((mask != NULL) && (mask->ZMIN == 1))) {
            z_3d[0] = atof(s2);
            not_already_set.ZMIN = 0;
        }
        return;
    }

    if (load_eqn_msc("POIMAP", s1)) {
        if ((not_already_set.POIMAP || force) || ((mask != NULL) && (mask->POIMAP == 1))) {
            if (s2[0] == 'm' || s2[0] == 'M') {
                POIMAP = 2;
            }
            if (s2[0] == 's' || s2[0] == 'S') {
                POIMAP = 1;
            }
            if (s2[0] == 'p' || s2[0] == 'P') {
                POIMAP = 3;
            }
            not_already_set.POIMAP = 0;
        }
        return;
    }

    if (load_eqn_msc("POIVAR", s1)) {
        if ((not_already_set.POIVAR || force) || ((mask != NULL) && (mask->POIVAR == 1))) {
            browser_find_variable(s2, &i);
            if (i > -1) {
                POIVAR = i;
            }

            not_already_set.POIVAR = 0;
        }
        return;
    }
    if (load_eqn_msc("OUTPUT", s1)) {
        if ((not_already_set.OUTPUT || force) || ((mask != NULL) && (mask->OUTPUT == 1))) {
            strcpy(batch_out, s2);
            not_already_set.OUTPUT = 0;
        }
        return;
    }

    if (load_eqn_msc("POISGN", s1)) {
        if ((not_already_set.POISGN || force) || ((mask != NULL) && (mask->POISGN == 1))) {
            POISGN = atoi(s2);
            not_already_set.POISGN = 0;
        }
        return;
    }

    if (load_eqn_msc("POISTOP", s1)) {
        if ((not_already_set.POISTOP || force) || ((mask != NULL) && (mask->POISTOP == 1))) {
            SOS = atoi(s2);
            not_already_set.POISTOP = 0;
        }
        return;
    }
    if (load_eqn_msc("STOCH", s1)) {
        if ((not_already_set.STOCH || force) || ((mask != NULL) && (mask->STOCH == 1))) {
            stoch_flag = atoi(s2);
            not_already_set.STOCH = 0;
        }
        return;
    }
    if (load_eqn_msc("POIPLN", s1)) {
        if ((not_already_set.POIPLN || force) || ((mask != NULL) && (mask->POIPLN == 1))) {
            POIPLN = atof(s2);
            not_already_set.POIPLN = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGEOVER", s1)) {
        if ((not_already_set.RANGEOVER || force) || ((mask != NULL) && (mask->RANGEOVER == 1))) {
            strcpy(range.item, s2);
            not_already_set.RANGEOVER = 0;
        }

        return;
    }
    if (load_eqn_msc("RANGESTEP", s1)) {
        if ((not_already_set.RANGESTEP || force) || ((mask != NULL) && (mask->RANGESTEP == 1))) {
            range.steps = atoi(s2);
            not_already_set.RANGESTEP = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGELOW", s1)) {
        if ((not_already_set.RANGELOW || force) || ((mask != NULL) && (mask->RANGELOW == 1))) {
            range.plow = atof(s2);
            not_already_set.RANGELOW = 0;
        }

        return;
    }

    if (load_eqn_msc("RANGEHIGH", s1)) {
        if ((not_already_set.RANGEHIGH || force) || ((mask != NULL) && (mask->RANGEHIGH == 1))) {
            range.phigh = atof(s2);
            not_already_set.RANGEHIGH = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGERESET", s1)) {
        if ((not_already_set.RANGERESET || force) || ((mask != NULL) && (mask->RANGERESET == 1))) {
            if (s2[0] == 'y' || s2[0] == 'Y') {
                range.reset = 1;
            } else {
                range.reset = 0;
            }
            not_already_set.RANGERESET = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGEOLDIC", s1)) {
        if ((not_already_set.RANGEOLDIC || force) || ((mask != NULL) && (mask->RANGEOLDIC == 1))) {
            if (s2[0] == 'y' || s2[0] == 'Y') {
                range.oldic = 1;
            } else {
                range.oldic = 0;
            }

            not_already_set.RANGEOLDIC = 0;
        }
        return;
    }

    if (load_eqn_msc("RANGE", s1)) {
        if ((not_already_set.RANGE || force) || ((mask != NULL) && (mask->RANGE == 1))) {
            batch_range = atoi(s2);
            not_already_set.RANGE = 0;
        }
        return;
    }

    if (load_eqn_msc("NTST", s1)) {
        if ((not_already_set.NTST || force) || ((mask != NULL) && (mask->NTST == 1))) {
            auto_ntst = atoi(s2);
            not_already_set.NTST = 0;
        }
        return;
    }
    if (load_eqn_msc("NMAX", s1)) {
        if ((not_already_set.NMAX || force) || ((mask != NULL) && (mask->NMAX == 1))) {
            auto_nmx = atoi(s2);
            not_already_set.NMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("NPR", s1)) {
        if ((not_already_set.NPR || force) || ((mask != NULL) && (mask->NPR == 1))) {
            auto_npr = atoi(s2);
            not_already_set.NPR = 0;
        }
        return;
    }
    if (load_eqn_msc("NCOL", s1)) {
        if ((not_already_set.NCOL || force) || ((mask != NULL) && (mask->NCOL == 1))) {
            auto_ncol = atoi(s2);
            not_already_set.NCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("DSMIN", s1)) {
        if ((not_already_set.DSMIN || force) || ((mask != NULL) && (mask->DSMIN == 1))) {
            auto_dsmin = atof(s2);
            not_already_set.DSMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("DSMAX", s1)) {
        if ((not_already_set.DSMAX || force) || ((mask != NULL) && (mask->DSMAX == 1))) {
            auto_dsmax = atof(s2);
            not_already_set.DSMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("DS", s1)) {
        if ((not_already_set.DS || force) || ((mask != NULL) && (mask->DS == 1))) {
            auto_ds = atof(s2);
            not_already_set.DS = 0;
        }

        return;
    }
    if (load_eqn_msc("PARMIN", s1)) {
        if ((not_already_set.XMAX || force) || ((mask != NULL) && (mask->XMAX == 1))) {
            auto_rl0 = atof(s2);
            not_already_set.XMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("PARMAX", s1)) {
        if ((not_already_set.PARMAX || force) || ((mask != NULL) && (mask->PARMAX == 1))) {
            auto_rl1 = atof(s2);
            not_already_set.PARMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("NORMMIN", s1)) {
        if ((not_already_set.NORMMIN || force) || ((mask != NULL) && (mask->NORMMIN == 1))) {
            auto_a0 = atof(s2);
            not_already_set.NORMMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("NORMMAX", s1)) {
        if ((not_already_set.NORMMAX || force) || ((mask != NULL) && (mask->NORMMAX == 1))) {
            auto_a1 = atof(s2);
            not_already_set.NORMMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("EPSL", s1)) {
        if ((not_already_set.EPSL || force) || ((mask != NULL) && (mask->EPSL == 1))) {
            auto_epsl = atof(s2);
            not_already_set.EPSL = 0;
        }
        return;
    }

    if (load_eqn_msc("EPSU", s1)) {
        if ((not_already_set.EPSU || force) || ((mask != NULL) && (mask->EPSU == 1))) {
            auto_epsu = atof(s2);
            not_already_set.EPSU = 0;
        }
        return;
    }
    if (load_eqn_msc("EPSS", s1)) {
        if ((not_already_set.EPSS || force) || ((mask != NULL) && (mask->EPSS == 1))) {
            auto_epss = atof(s2);
            not_already_set.EPSS = 0;
        }
        return;
    }
    if (load_eqn_msc("RUNNOW", s1)) {
        if ((not_already_set.RUNNOW || force) || ((mask != NULL) && (mask->RUNNOW == 1))) {
            run_immediately = atoi(s2);
            not_already_set.RUNNOW = 0;
        }
        return;
    }

    if (load_eqn_msc("SEC", s1)) {
        if ((not_already_set.SEC || force) || ((mask != NULL) && (mask->SEC == 1))) {
            SEc = atoi(s2);
            not_already_set.SEC = 0;
        }
        return;
    }
    if (load_eqn_msc("UEC", s1)) {
        if ((not_already_set.UEC || force) || ((mask != NULL) && (mask->UEC == 1))) {
            UEc = atoi(s2);
            not_already_set.UEC = 0;
        }
        return;
    }
    if (load_eqn_msc("SPC", s1)) {
        if ((not_already_set.SPC || force) || ((mask != NULL) && (mask->SPC == 1))) {
            SPc = atoi(s2);
            not_already_set.SPC = 0;
        }
        return;
    }
    if (load_eqn_msc("UPC", s1)) {
        if ((not_already_set.UPC || force) || ((mask != NULL) && (mask->UPC == 1))) {
            UPc = atoi(s2);
            not_already_set.UPC = 0;
        }
        return;
    }

    if (load_eqn_msc("AUTOEVAL", s1)) {
        if ((not_already_set.AUTOEVAL || force) || ((mask != NULL) && (mask->AUTOEVAL == 1))) {
            f = atoi(s2);
            tabular_set_auto_eval_flags(f);
            not_already_set.AUTOEVAL = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOXMAX", s1)) {
        if ((not_already_set.AUTOXMAX || force) || ((mask != NULL) && (mask->AUTOXMAX == 1))) {
            auto_xmax = atof(s2);
            not_already_set.AUTOXMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOYMAX", s1)) {
        if ((not_already_set.AUTOYMAX || force) || ((mask != NULL) && (mask->AUTOYMAX == 1))) {
            auto_ymax = atof(s2);
            not_already_set.AUTOYMAX = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOXMIN", s1)) {
        if ((not_already_set.AUTOXMIN || force) || ((mask != NULL) && (mask->AUTOXMIN == 1))) {
            auto_xmin = atof(s2);
            not_already_set.AUTOXMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOYMIN", s1)) {
        if ((not_already_set.AUTOYMIN || force) || ((mask != NULL) && (mask->AUTOYMIN == 1))) {
            auto_ymin = atof(s2);
            not_already_set.AUTOYMIN = 0;
        }
        return;
    }
    if (load_eqn_msc("AUTOVAR", s1)) {
        if ((not_already_set.AUTOVAR || force) || ((mask != NULL) && (mask->AUTOVAR == 1))) {
            browser_find_variable(s2, &i);
            if (i > 0) {
                auto_var = i - 1;
            }
            not_already_set.AUTOVAR = 0;
        }
        return;
    }

    // postscript options

    if (load_eqn_msc("ps_font", s1)) {
        if ((not_already_set.ps_font || force) || ((mask != NULL) && (mask->ps_font == 1))) {
            strcpy(ps_font, s2);
            not_already_set.ps_font = 0;
        }
        return;
    }

    if (load_eqn_msc("ps_lw", s1)) {
        if ((not_already_set.ps_lw || force) || ((mask != NULL) && (mask->ps_lw == 1))) {
            ps_lw = atof(s2);
            not_already_set.ps_lw = 0;
        }
        return;
    }

    if (load_eqn_msc("PS_FSIZE", s1)) {
        if ((not_already_set.PS_FSIZE || force) || ((mask != NULL) && (mask->PS_FSIZE == 1))) {
            ps_fontsize = atoi(s2);
            not_already_set.PS_FSIZE = 0;
        }
        return;
    }

    if (load_eqn_msc("PS_COLOR", s1)) {
        if ((not_already_set.PS_COLOR || force) || ((mask != NULL) && (mask->PS_COLOR == 1))) {
            ps_color_flag = atoi(s2);
            ps_color = ps_color_flag;
            not_already_set.PS_COLOR = 0;
        }
        return;
    }
    if (load_eqn_msc("TUTORIAL", s1)) {
        if (!(load_eqn_msc(s2, "0") || load_eqn_msc(s2, "1"))) {
            ggets_plintf("TUTORIAL option must be 0 or 1.\n");
            exit(-1);
        }
        if ((not_already_set.TUTORIAL || force) || ((mask != NULL) && (mask->TUTORIAL == 1))) {
            do_tutorial = atoi(s2);
            not_already_set.TUTORIAL = 0;
        }
        return;
    }
    if (load_eqn_msc("S1", s1)) {
        if ((not_already_set.slider1 || force) || ((mask != NULL) && (mask->slider1 == 1))) {
            strncpy(slider1var, s2, 20);
            slider1var[19] = '\0';
            not_already_set.slider1 = 0;
        }
        return;
    }

    if (load_eqn_msc("S2", s1)) {
        if ((not_already_set.slider2 || force) || ((mask != NULL) && (mask->slider2 == 1))) {
            strncpy(slider2var, s2, 20);
            slider2var[19] = '\0';
            not_already_set.slider2 = 0;
        }
        return;
    }
    if (load_eqn_msc("S3", s1)) {
        if ((not_already_set.slider3 || force) || ((mask != NULL) && (mask->slider3 == 1))) {
            strncpy(slider3var, s2, 20);
            slider3var[19] = '\0';
            not_already_set.slider3 = 0;
        }
        return;
    }
    if (load_eqn_msc("SLO1", s1)) {
        if ((not_already_set.slider1lo || force) || ((mask != NULL) && (mask->slider1lo == 1))) {
            slider1lo = atof(s2);
            not_already_set.slider1lo = 0;
        }
        return;
    }

    if (load_eqn_msc("SLO2", s1)) {
        if ((not_already_set.slider2lo || force) || ((mask != NULL) && (mask->slider2lo == 1))) {
            slider2lo = atof(s2);
            not_already_set.slider2lo = 0;
        }
        return;
    }
    if (load_eqn_msc("SLO3", s1)) {
        if ((not_already_set.slider3lo || force) || ((mask != NULL) && (mask->slider3lo == 1))) {
            slider3lo = atof(s2);
            not_already_set.slider3lo = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI1", s1)) {
        if ((not_already_set.slider1hi || force) || ((mask != NULL) && (mask->slider1hi == 1))) {
            slider1hi = atof(s2);
            not_already_set.slider1hi = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI2", s1)) {
        if ((not_already_set.slider2hi || force) || ((mask != NULL) && (mask->slider2hi == 1))) {
            slider2hi = atof(s2);
            not_already_set.slider2hi = 0;
        }
        return;
    }
    if (load_eqn_msc("SHI3", s1)) {
        if ((not_already_set.slider3hi || force) || ((mask != NULL) && (mask->slider3hi == 1))) {
            slider3hi = atof(s2);
            not_already_set.slider3hi = 0;
        }
        return;
    }

    /* postprocessing options
       This is rally only relevant for batch jobs as it
       writes files then
    */

    if (load_eqn_msc("POSTPROCESS", s1)) {
        if ((not_already_set.POSTPROCESS || force) || ((mask != NULL) && (mask->POSTPROCESS == 1))) {
            post_process = atoi(s2);
            not_already_set.POSTPROCESS = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTLO", s1)) {
        if ((not_already_set.HISTLO || force) || ((mask != NULL) && (mask->HISTLO == 1))) {
            hist_inf.xlo = atof(s2);
            not_already_set.HISTLO = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTHI", s1)) {
        if ((not_already_set.HISTHI || force) || ((mask != NULL) && (mask->HISTHI == 1))) {
            hist_inf.xhi = atof(s2);
            not_already_set.HISTHI = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTBINS", s1)) {
        if ((not_already_set.HISTBINS || force) || ((mask != NULL) && (mask->HISTBINS == 1))) {
            hist_inf.nbins = atoi(s2);
            not_already_set.HISTBINS = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTCOL", s1)) {
        if ((not_already_set.HISTCOL || force) || ((mask != NULL) && (mask->HISTCOL == 1))) {
            browser_find_variable(s2, &i);
            if (i > (-1)) {
                hist_inf.col = i;
            }
            not_already_set.HISTCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTLO2", s1)) {
        if ((not_already_set.HISTLO2 || force) || ((mask != NULL) && (mask->HISTLO2 == 1))) {
            hist_inf.ylo = atof(s2);
            not_already_set.HISTLO2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTHI2", s1)) {
        if ((not_already_set.HISTHI2 || force) || ((mask != NULL) && (mask->HISTHI2 == 1))) {
            hist_inf.yhi = atof(s2);
            not_already_set.HISTHI2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTBINS2", s1)) {
        if ((not_already_set.HISTBINS2 || force) || ((mask != NULL) && (mask->HISTBINS2 == 1))) {
            hist_inf.nbins2 = atoi(s2);
            not_already_set.HISTBINS2 = 0;
        }
        return;
    }

    if (load_eqn_msc("HISTCOL2", s1)) {
        if ((not_already_set.HISTCOL2 || force) || ((mask != NULL) && (mask->HISTCOL2 == 1))) {
            browser_find_variable(s2, &i);
            if (i > (-1)) {
                hist_inf.col2 = i;
            }
            not_already_set.HISTCOL2 = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECCOL", s1)) {
        if ((not_already_set.SPECCOL || force) || ((mask != NULL) && (mask->SPECCOL == 1))) {
            browser_find_variable(s2, &i);
            if (i > (-1)) {
                spec_col = i;
            }
            not_already_set.SPECCOL = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECCOL2", s1)) {
        if ((not_already_set.SPECCOL2 || force) || ((mask != NULL) && (mask->SPECCOL2 == 1))) {
            browser_find_variable(s2, &i);
            if (i > (-1)) {
                spec_col2 = i;
            }
            not_already_set.SPECCOL2 = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECWIDTH", s1)) {
        if ((not_already_set.SPECWIDTH || force) || ((mask != NULL) && (mask->SPECWIDTH == 1))) {
            spec_wid = atoi(s2);
            not_already_set.SPECWIDTH = 0;
        }
        return;
    }

    if (load_eqn_msc("SPECWIN", s1)) {
        if ((not_already_set.SPECWIN || force) || ((mask != NULL) && (mask->SPECWIN == 1))) {
            spec_win = atoi(s2);
            not_already_set.SPECWIN = 0;
        }
        return;
    }

    if (load_eqn_msc("DFGRID", s1)) {
        if ((not_already_set.DFGRID || force) || ((mask != NULL) && (mask->DFGRID == 1))) {
            df_grid = atoi(s2);
            not_already_set.DFGRID = 0;
        }
        return;
    }
    if (load_eqn_msc("DFDRAW", s1)) {
        if ((not_already_set.DFBATCH || force) || ((mask != NULL) && (mask->DFBATCH == 1))) {
            df_batch = atoi(s2);
            not_already_set.DFBATCH = 0;
        }
        return;
    }
    if (load_eqn_msc("NCDRAW", s1)) {
        if ((not_already_set.NCBATCH || force) || ((mask != NULL) && (mask->NCBATCH == 1))) {
            NCBatch = atoi(s2);
            not_already_set.NCBATCH = 0;
        }
        return;
    }

    // colorize customizing !!
    if (load_eqn_msc("COLORVIA", s1)) {
        if ((not_already_set.COLORVIA || force) || ((mask != NULL) && (mask->COLORVIA == 1))) {
            strcpy(color_via, s2);
        }
        not_already_set.COLORVIA = 0;
        return;
    }
    if (load_eqn_msc("COLORIZE", s1)) {
        if ((not_already_set.COLORIZE || force) || ((mask != NULL) && (mask->COLORIZE == 1))) {
            colorize_flag = atoi(s2);
        }
        not_already_set.COLORIZE = 0;
        return;
    }
    if (load_eqn_msc("COLORLO", s1)) {
        if ((not_already_set.COLORLO || force) || ((mask != NULL) && (mask->COLORLO == 1))) {
            color_via_lo = atof(s2);
        }
        not_already_set.COLORLO = 0;
        return;
    }
    if (load_eqn_msc("COLORHI", s1)) {
        if ((not_already_set.COLORHI || force) || ((mask != NULL) && (mask->COLORHI == 1))) {
            color_via_hi = atof(s2);
        }
        not_already_set.COLORHI = 0;
        return;
    }

    ggets_plintf("!! Option %s not recognized\n", s1);
}
