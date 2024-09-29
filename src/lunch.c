#include "functions.h"
#include "parserslow.h"
#include "integers.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <time.h>
#include "struct.h"

#define READEM 1
#define VOLTERRA 6
#define PARAMBOX 1

static int32 set_type = 0;

static void lunch_io_graph(int32 f, FILE *fp);
static void lunch_io_exprs(int32 f, FILE *fp);
static void lunch_io_parameters(int32 f, FILE *fp);
static void lunch_io_numerics(int32 f, FILE *fp);
static void lunch_dump_eqn(FILE *fp);
static void lunch_do_info(FILE *fp);

void
lunch_file_inf(void) {
    int32 ok;
    FILE *fp;
    char filename[XPP_MAX_NAME + 5];
    sprintf(filename, "%s.pars", this_file);
    ggets_ping();
    if (!init_conds_file_selector("Save info", filename, "*.pars*")) {
        return;
    }
    browser_open_write_file(&fp, filename, &ok);
    if (!ok) {
        return;
    }
    init_conds_redraw_params();
    lunch_do_info(fp);
    fclose(fp);
    return;
}

void
lunch_ps_write_pars(FILE *fp) {
    int32 div;
    int32 rem;
    double z;

    fprintf(fp, "\n %%%% %s \n %%%% Parameters ...\n", this_file);
    div = NUPAR / 4;
    rem = NUPAR % 4;
    for (int32 j = 0; j < div; j++) {
        for (int32 i = 0; i < 4; i++) {
            get_val(upar_names[i + 4*j], &z);
            fprintf(fp, "%%%% %s=%.16g   ", upar_names[i + 4*j], z);
        }
        fprintf(fp, "\n");
    }
    for (int32 i = 0; i < rem; i++) {
        get_val(upar_names[i + 4*div], &z);
        fprintf(fp, "%%%% %s=%.16g   ", upar_names[i + 4*div], z);
    }

    fprintf(fp, "\n");
    return;
}

void
lunch_do_info(FILE *fp) {
    static char *method[] = {
        "Discrete", "Euler",    "Mod. Euler", "Runge-Kutta", "Adams",
        "Gear",     "Volterra", "BackEul",    "QualRK",      "Stiff",
        "CVode",    "DoPri5",   "DoPri8(3)",  "Rosenbrock",  "Symplectic"};
    int32 div;
    int32 rem;
    double z;
    char bob[200];
    char fstr[15];
    fprintf(fp, "File: %s \n\n Equations... \n", this_file);
    for (int32 i = 0; i < NEQ; i++) {
        if (i < NODE && METHOD > 0) {
            strcpy(fstr, "d%s/dT=%s\n");
        }
        if (i < NODE && METHOD == 0) {
            strcpy(fstr, "%s(n+1)=%s\n");
        }
        if (i >= NODE) {
            strcpy(fstr, "%s=%s\n");
        }
        fprintf(fp, fstr, uvar_names[i], ode_names[i]);
    }

    if (fix_var > 0) {
        fprintf(fp, "\nwhere ...\n");
        for (int32 i = 0; i < fix_var; i++) {
            fprintf(fp, "%s = %s \n", fixinfo[i].name, fixinfo[i].value);
        }
    }
    if (NFUN > 0) {
        fprintf(fp, "\nUser-defined functions:\n");
        edit_rhs_user_fun_info(fp);
    }

    fprintf(fp, "\n\n Numerical parameters ...\n");

    fprintf(fp, "NJMP=%d  NMESH=%d METHOD=%s evec_iter=%d \n", NJMP, NMESH,
            method[METHOD], evec_iter);
    fprintf(fp, "bvp_eps=%g,bvp_tol=%g,bpv_maxit=%d \n", bvp_eps, bvp_tol,
            bpv_maxit);
    fprintf(fp, "DT=%g T0=%g TRANS=%g TEND=%g bound=%g delay=%g MaxPts=%d\n",
            delta_t, T0, TRANS, TEND, bound, delay, MaxPoints);
    fprintf(fp, "evec_err=%g, NEWT_ERR=%g h_min=%g h_max=%g TOLER=%g \n",
            evec_err, NEWT_ERR, h_min, h_max, TOLER);
    if (POIVAR == 0) {
        strcpy(bob, "T");
    } else {
        strcpy(bob, uvar_names[POIVAR - 1]);
    }
    fprintf(fp, "POIMAP=%d POIVAR=%s POIPLN=%g POISGN=%d \n", POIMAP, bob,
            POIPLN, POISGN);

    fprintf(fp, "\n\n Delay strings ...\n");

    for (int32 i = 0; i < NODE; i++) {
        fprintf(fp, "%s\n", delay_string[i]);
    }
    fprintf(fp, "\n\n BCs ...\n");

    for (int32 i = 0; i < NODE; i++) {
        fprintf(fp, "0=%s\n", my_bc[i].string);
    }
    fprintf(fp, "\n\n ICs ...\n");

    for (int32 i = 0; i < NODE + NMarkov; i++) {
        fprintf(fp, "%s=%.16g\n", uvar_names[i], last_ic[i]);
    }
    fprintf(fp, "\n\n Parameters ...\n");
    div = NUPAR / 4;
    rem = NUPAR % 4;
    for (int32 j = 0; j < div; j++) {
        for (int32 i = 0; i < 4; i++) {
            get_val(upar_names[i + 4*j], &z);
            fprintf(fp, "%s=%.16g   ", upar_names[i + 4*j], z);
        }
        fprintf(fp, "\n");
    }
    for (int32 i = 0; i < rem; i++) {
        get_val(upar_names[i + 4*div], &z);
        fprintf(fp, "%s=%.16g   ", upar_names[i + 4*div], z);
    }

    fprintf(fp, "\n");
    return;
}

int32
lunch_read(FILE *fp) {
    int32 f = READEM;
    int32 ne;
    int32 np;
    int32 temp;
    char bob[256];

    fgets(bob, 255, fp);
    if (bob[0] == '#') {
        set_type = 1;
        lunch_io_int(&ne, fp, f, " ");
    } else {
        ne = atoi(bob);
        set_type = 0;
    }
    lunch_io_int(&np, fp, f, " ");
    if (ne != NEQ || np != NUPAR) {
        ggets_plintf("Set file has incompatible parameters\n");
        return 0;
    }
    lunch_io_numerics(f, fp);
    if (METHOD == VOLTERRA) {
        lunch_io_int(&temp, fp, f, " ");
        volterra_allocate(temp, 1);
        MyStart = 1;
    }
    numerics_chk_delay();
    lunch_io_exprs(f, fp);
    lunch_io_graph(f, fp);
    if (set_type == 1) {
        adjoints_dump_transpose_info(fp, f);
        adjoints_dump_h_stuff(fp, f);
        array_plot_dump(fp, f);
        load_eqn_dump_torus(fp, f);
        integrate_dump_range(fp, f);
    }

    return 1;
}

void
do_lunch(int32 f) {
    // f=1 to read and 0 to write
    int32 ne;
    int32 np;
    int32 ok;
    int32 temp;
    char bob[256];
    FILE *fp;
    time_t ttt;
    char filename[XPP_MAX_NAME + 4];
    sprintf(filename, "%s.set", this_file);

    if (f == READEM) {
        ggets_ping();
        if (!init_conds_file_selector("Load SET File", filename, "*.set")) {
            return;
        }

        fp = fopen(filename, "r");
        if (fp == NULL) {
            ggets_err_msg("Cannot open file");
            return;
        }
        fgets(bob, 255, fp);
        if (bob[0] == '#') {
            set_type = 1;
            lunch_io_int(&ne, fp, f, " ");
        } else {
            ne = atoi(bob);
            set_type = 0;
        }
        lunch_io_int(&np, fp, f, " ");
        if (ne != NEQ || np != NUPAR) {
            ggets_err_msg("Incompatible parameters");
            fclose(fp);
            return;
        }
        lunch_io_numerics(f, fp);
        if (METHOD == VOLTERRA) {
            lunch_io_int(&temp, fp, f, " ");
            volterra_allocate(temp, 1);
            MyStart = 1;
        }
        numerics_chk_delay();
        lunch_io_exprs(f, fp);
        lunch_io_graph(f, fp);
        if (set_type == 1) {
            adjoints_dump_transpose_info(fp, f);
            adjoints_dump_h_stuff(fp, f);
            array_plot_dump(fp, f);
            load_eqn_dump_torus(fp, f);
            integrate_dump_range(fp, f);
        }
        fclose(fp);
        return;
    }
    if (!init_conds_file_selector("Save SET File", filename, "*.set")) {
        return;
    }
    browser_open_write_file(&fp, filename, &ok);
    if (!ok) {
        return;
    }
    init_conds_redraw_params();
    ttt = time(0);
    fprintf(fp, "## Set file for %s on %s", this_file, ctime(&ttt));
    lunch_io_int(&NEQ, fp, f, "Number of equations and auxiliaries");
    lunch_io_int(&NUPAR, fp, f, "Number of parameters");
    lunch_io_numerics(f, fp);
    if (METHOD == VOLTERRA) {
        lunch_io_int(&MaxPoints, fp, f, "Max points for volterra");
    }
    lunch_io_exprs(f, fp);
    lunch_io_graph(f, fp);
    adjoints_dump_transpose_info(fp, f);
    adjoints_dump_h_stuff(fp, f);
    array_plot_dump(fp, f);
    load_eqn_dump_torus(fp, f);
    integrate_dump_range(fp, f);
    lunch_dump_eqn(fp);
    fclose(fp);
    return;
}

void
lunch_dump_eqn(FILE *fp) {
    char fstr[15];
    fprintf(fp, "RHS etc ...\n");
    for (int32 i = 0; i < NEQ; i++) {
        if (i < NODE && METHOD > 0) {
            strcpy(fstr, "d%s/dT=%s\n");
        }
        if (i < NODE && METHOD == 0) {
            strcpy(fstr, "%s(n+1)=%s\n");
        }
        if (i >= NODE) {
            strcpy(fstr, "%s=%s\n");
        }
        fprintf(fp, fstr, uvar_names[i], ode_names[i]);
    }

    if (fix_var > 0) {
        fprintf(fp, "\nwhere ...\n");
        for (int32 i = 0; i < fix_var; i++) {
            fprintf(fp, "%s = %s \n", fixinfo[i].name, fixinfo[i].value);
        }
    }
    if (NFUN > 0) {
        fprintf(fp, "\nUser-defined functions:\n");
        edit_rhs_user_fun_info(fp);
    }
    return;
}

void
lunch_io_numerics(int32 f, FILE *fp) {
    char *method[] = {"Discrete",  "Euler", "Mod. Euler", "Runge-Kutta",
                      "Adams",     "Gear",  "Volterra",   "BackEul",
                      "Qual RK",   "Stiff", "CVode",      "DorPrin5",
                      "DorPri8(3)"};
    char *pmap[] = {"Poincare None", "Poincare Section", "Poincare Max",
                    "Period"};
    char temp[256];
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Numerical stuff\n");
    }
    lunch_io_int(&NJMP, fp, f, " nout");
    lunch_io_int(&NMESH, fp, f, " nullcline mesh");
    lunch_io_int(&METHOD, fp, f, method[METHOD]);
    if (f == READEM) {
        numerics_do_meth();
        storage_alloc_meth();
    }
    lunch_io_double(&TEND, fp, f, "total");
    lunch_io_double(&delta_t, fp, f, "DeltaT");
    lunch_io_double(&T0, fp, f, "T0");
    lunch_io_double(&TRANS, fp, f, "Transient");
    lunch_io_double(&bound, fp, f, "Bound");
    lunch_io_double(&h_min, fp, f, "DtMin");
    lunch_io_double(&h_max, fp, f, "DtMax");
    lunch_io_double(&TOLER, fp, f, "Tolerance");
    // fix stuff concerning the tolerance
    if (f == READEM) {
        if (set_type == 1) {
            lunch_io_double(&atoler, fp, f, "Abs. Tolerance");
        } else {
            atoler = TOLER*10;
        }
    } else {
        lunch_io_double(&atoler, fp, f, "Abs. Tolerance");
    }

    lunch_io_double(&delay, fp, f, "Max Delay");
    lunch_io_int(&evec_iter, fp, f, "Eigenvector iterates");
    lunch_io_double(&evec_err, fp, f, "Eigenvector tolerance");
    lunch_io_double(&NEWT_ERR, fp, f, "Newton tolerance");
    lunch_io_double(&POIPLN, fp, f, "Poincare plane");
    lunch_io_double(&bvp_tol, fp, f, "Boundary value tolerance");
    lunch_io_double(&bvp_eps, fp, f, "Boundary value epsilon");
    lunch_io_int(&bpv_maxit, fp, f, "Boundary value iterates");
    lunch_io_int(&POIMAP, fp, f, pmap[POIMAP]);

    lunch_io_int(&POIVAR, fp, f, "Poincare variable");
    lunch_io_int(&POISGN, fp, f, "Poincare sign");
    lunch_io_int(&SOS, fp, f, "Stop on Section");
    lunch_io_int(&delay_flag, fp, f, "Delay flag");
    lunch_io_double(&MyTime, fp, f, "Current time");
    lunch_io_double(&LastTime, fp, f, "Last Time");
    lunch_io_int(&MyStart, fp, f, "MyStart");
    lunch_io_int(&in_flag, fp, f, "in_flag");
    return;
}

void
lunch_io_parameter_file(char *fn, int32 flag) {
    char fnx[256];
    char c;
    int32 j = 0;
    int32 np;
    FILE *fp;
    time_t ttt;
    for (usize i = 6; i < strlen(fn); i++) {
        c = fn[i];
        if (c != ' ') {
            fnx[j] = c;
            j++;
        }
    }
    fnx[j] = 0;
    if (flag == READEM) {
        fp = fopen(fnx, "r");
        if (fp == NULL) {
            ggets_err_msg("Cannot open file");
            return;
        }
        lunch_io_int(&np, fp, flag, " ");
        if (np != NUPAR) {
            printf("%d", np);
            printf("%d", NUPAR);
            ggets_err_msg("Incompatible parameters");
            fclose(fp);
            return;
        }
        lunch_io_parameters(flag, fp);
        fclose(fp);
        init_conds_redo_stuff();

        return;
    }
    fp = fopen(fnx, "w");
    if (fp == NULL) {
        ggets_err_msg("Cannot open file");
        return;
    }
    lunch_io_int(&NUPAR, fp, flag, "Number params");
    lunch_io_parameters(flag, fp);
    ttt = time(0);
    fprintf(fp, "\n\nFile:%s\n%s", this_file, ctime(&ttt));
    fclose(fp);
    return;
}

void
lunch_io_ic_file(char *fn, int32 flag) {
    char fnx[256];
    char c;
    int32 j = 0;
    int32 chk = 0;
    FILE *fp;
    char msg[256];

    for (usize i = 0; i < strlen(fn); i++) {
        c = fn[i];
        if (c != ' ') {
            fnx[j] = c;
            j++;
        }
    }
    fnx[j] = 0;
    if (flag == READEM) {
        fp = fopen(fnx, "r");
        if (fp == NULL) {
            ggets_err_msg("Cannot open file");
            return;
        }
        for (int32 i = 0; i < NODE; i++) {
            chk = fscanf(fp, "%lg", &last_ic[i]);
            if (chk != 1) {
                sprintf(
                    msg,
                    "Expected %d initial conditions but only found %d in %s.",
                    NODE, i, fn);
                ggets_err_msg(msg);
                return;
            }
        }

        while (chk != EOF) {
            chk = fscanf(fp, "%lg", &last_ic[NODE]);
            if (chk != EOF) {
                sprintf(msg, "Found more than %d initial conditions in %s.",
                        NODE, fn);
                ggets_err_msg(msg);
                return;
            }
        }
        fclose(fp);
    }
}

void
lunch_io_parameters(int32 f, FILE *fp) {
    int32 index;
    char junk[256];
    double z;
    for (int32 i = 0; i < NUPAR; i++) {
        if (f != READEM) {
            get_val(upar_names[i], &z);
            lunch_io_double(&z, fp, f, upar_names[i]);
        } else {
            lunch_io_double(&z, fp, f, " ");
            set_val(upar_names[i], z);

            if (!xpp_batch) {
                index = init_conds_find_user_name(PARAMBOX, upar_names[i]);
                if (index >= 0) {
                    sprintf(junk, "%.16g", z);
                    init_conds_set_edit_params(&ParamBox, index, junk);
                    init_conds_draw_one_box(ParamBox, index);
                }
            }
        }
    }

    if (!xpp_batch) {
        init_conds_reset_sliders();
    }
    return;
}

void
lunch_io_exprs(int32 f, FILE *fp) {
    char temp[256];
    double z;

    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Delays\n");
    }
    for (int32 i = 0; i < NODE; i++) {
        lunch_io_string(delay_string[i], 100, fp, f);
    }
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Bndry conds\n");
    }
    for (int32 i = 0; i < NODE; i++) {
        lunch_io_string(my_bc[i].string, 100, fp, f);
    }
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Old ICs\n");
    }
    for (int32 i = 0; i < NODE + NMarkov; i++) {
        lunch_io_double(&last_ic[i], fp, f, uvar_names[i]);
    }
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Ending  ICs\n");
    }
    for (int32 i = 0; i < NODE + NMarkov; i++) {
        lunch_io_double(&MyData[i], fp, f, uvar_names[i]);
    }
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Parameters\n");
    }
    for (int32 i = 0; i < NUPAR; i++) {
        if (f != READEM) {
            get_val(upar_names[i], &z);
            lunch_io_double(&z, fp, f, upar_names[i]);
        } else {
            lunch_io_double(&z, fp, f, " ");
            set_val(upar_names[i], z);
        }
    }

    if (f == READEM && Xup) {
        init_conds_redraw_bcs();
        init_conds_redraw_ics();
        init_conds_redraw_delays();
        init_conds_redraw_params();
    }
    return;
}

void
lunch_io_graph(int32 f, FILE *fp) {
    char temp[256];
    if (f == READEM && set_type == 1) {
        fgets(temp, 255, fp);  // skip a line
    }
    if (f != READEM) {
        fprintf(fp, "# Graphics\n");
    }
    for (int32 j = 0; j < 3; j++) {
        for (int32 k = 0; k < 3; k++) {
            lunch_io_double(&(MyGraph->rm[k][j]), fp, f, "rm");
        }
    }
    for (int32 j = 0; j < MAXPERPLOT; j++) {
        lunch_io_int(&(MyGraph->xv[j]), fp, f, " ");
        lunch_io_int(&(MyGraph->yv[j]), fp, f, " ");
        lunch_io_int(&(MyGraph->zv[j]), fp, f, " ");
        lunch_io_int(&(MyGraph->line[j]), fp, f, " ");
        lunch_io_int(&(MyGraph->color[j]), fp, f, " ");
    }

    lunch_io_double(&(MyGraph->ZPlane), fp, f, " ");
    lunch_io_double(&(MyGraph->ZView), fp, f, " ");
    lunch_io_int(&(MyGraph->PerspFlag), fp, f, " ");
    lunch_io_int(&(MyGraph->ThreeDFlag), fp, f, "3DFlag");
    lunch_io_int(&(MyGraph->TimeFlag), fp, f, "Timeflag");
    lunch_io_int(&(MyGraph->ColorFlag), fp, f, "Colorflag");
    lunch_io_int(&(MyGraph->grtype), fp, f, "Type");
    lunch_io_double(&(MyGraph->color_scale), fp, f, "color scale");
    lunch_io_double(&(MyGraph->min_scale), fp, f, " minscale");

    lunch_io_double(&(MyGraph->xmax), fp, f, " xmax");
    lunch_io_double(&(MyGraph->xmin), fp, f, " xmin");
    lunch_io_double(&(MyGraph->ymax), fp, f, " ymax");
    lunch_io_double(&(MyGraph->ymin), fp, f, " ymin");
    lunch_io_double(&(MyGraph->zmax), fp, f, " zmax");
    lunch_io_double(&(MyGraph->zmin), fp, f, " zmin");
    lunch_io_double(&(MyGraph->xbar), fp, f, " ");
    lunch_io_double(&(MyGraph->dx), fp, f, " ");
    lunch_io_double(&(MyGraph->ybar), fp, f, " ");
    lunch_io_double(&(MyGraph->dy), fp, f, " ");
    lunch_io_double(&(MyGraph->zbar), fp, f, " ");
    lunch_io_double(&(MyGraph->dz), fp, f, " ");

    lunch_io_double(&(MyGraph->Theta), fp, f, " Theta");
    lunch_io_double(&(MyGraph->Phi), fp, f, " Phi");
    lunch_io_int(&(MyGraph->xshft), fp, f, " xshft");
    lunch_io_int(&(MyGraph->yshft), fp, f, " yshft");
    lunch_io_int(&(MyGraph->zshft), fp, f, " zshft");
    lunch_io_double(&(MyGraph->xlo), fp, f, " xlo");
    lunch_io_double(&(MyGraph->ylo), fp, f, " ylo");
    lunch_io_double(&(MyGraph->oldxlo), fp, f, " ");
    lunch_io_double(&(MyGraph->oldylo), fp, f, " ");
    lunch_io_double(&(MyGraph->xhi), fp, f, " xhi");
    lunch_io_double(&(MyGraph->yhi), fp, f, " yhi");
    lunch_io_double(&(MyGraph->oldxhi), fp, f, " ");
    lunch_io_double(&(MyGraph->oldyhi), fp, f, " ");
    if (f == READEM && Xup) {
        graf_par_redraw_the_graph();
    }
    return;
}

void
lunch_io_int(int32 *i, FILE *fp, int32 f, char *ss) {
    char bob[256];
    if (f == READEM) {
        fgets(bob, 255, fp);
        *i = atoi(bob);
    } else {
        fprintf(fp, "%d   %s\n", *i, ss);
    }
    return;
}

void
lunch_io_double(double *z, FILE *fp, int32 f, char *ss) {
    char bob[256];
    if (f == READEM) {
        fgets(bob, 255, fp);
        *z = atof(bob);
    } else {
        fprintf(fp, "%.16g  %s\n", *z, ss);
    }
    return;
}

void
lunch_io_string(char *s, int32 len, FILE *fp, int32 f) {
    usize i;
    if (f == READEM) {
        fgets(s, len, fp);
        i = 0;
        while (i < strlen(s)) {
            if (s[i] == '\n') {
                s[i] = 0;
            }
            i++;
        }
    } else {
        fprintf(fp, "%s\n", s);
    }
}
