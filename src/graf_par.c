
#include "graf_par.h"

#include "integrate.h"

#include "menudrive.h"
#include "aniparse.h"
#include "arrayplot.h"
#include "color.h"
#include "init_conds.h"
#include "rubber.h"

#include "ggets.h"
#include "graphics.h"
#include "menu.h"
#include "pop_list.h"
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <math.h>
#include "xpplim.h"
#include "struct.h"
#include "browse.h"
#include "mykeydef.h"
#include "many_pops.h"
#include "kinescope.h"
#include "nullcline.h"
#include "axes2.h"
#include "my_ps.h"
#include "my_svg.h"
#include "load_eqn.h"
#include "integers.h"
#include <libgen.h>

NCLINE nclines[MAXNCLINE];
extern CURVE frz[MAXFRZ];
extern GRAPH *MyGraph;
extern Display *display;
extern Window main_win, draw_win, info_pop;
extern int32 DCURY;
extern int32 storind;
extern int32 PS_FONTSIZE;
extern int32 PS_Port;
/*extern char PS_FONT[100];*/
extern char PS_FONT[XPP_MAX_NAME];
extern double PS_LW;
extern BROWSER my_browser;
extern double x_3d[2], y_3d[2], z_3d[2];
/*Default is now color*/
int32 PS_Color = 1;

extern char PlotFormat[100];

#define SPER 3
#define UPER 4
#define SEQ 1
#define UEQ 2

#define lsSEQ 0
#define lsUEQ 1
#define lsSPER 8
#define lsUPER 9

MOV3D mov3d = {"theta", "N", 45, 45, 7};

BD my_bd;

extern int32 DLeft, DRight, DTop, DBottom, VTic, HTic, VChar, HChar;

extern double T0, TEND;
extern float **storage;

double FreezeKeyX, FreezeKeyY;
int32 FreezeKeyFlag, AutoFreezeFlag = 0;
int32 CurrentCurve = 0;
extern char this_file[XPP_MAX_NAME];
extern char this_internset[XPP_MAX_NAME];

extern int32 PltFmtFlag;
extern char uvar_names[MAXODE][12];

extern char *info_message, *no_hint[], *wind_hint[], *view_hint[], *frz_hint[];
extern char *graf_hint[], *cmap_hint[];

int32 colorline[] = {0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 0};
char *color_names[] = {"WHITE",        "RED",    "REDORANGE",   "ORANGE",
                       "YELLOWORANGE", "YELLOW", "YELLOWGREEN", "GREEN",
                       "BLUEGREEN",    "BLUE",   "PURPLE",      "BLACK"};

void
change_view_com(int32 com) {

    if (com == 2) {
        make_my_aplot("Array!");
        edit_aplot();
        return;
    }
    if (com == 3) {
        new_vcr();
        return;
    }

    MyGraph->grtype = 5 * com;
    if (MyGraph->grtype < 5)
        get_2d_view(CurrentCurve);
    else
        get_3d_view(CurrentCurve);
    check_flags();
    redraw_the_graph();
}

void
ind_to_sym(int32 ind, char *str) {
    if (ind == 0)
        strcpy(str, "T");
    else
        strcpy(str, uvar_names[ind - 1]);
    return;
}

void
check_flags(void) {
    if (MyGraph->grtype > 4)
        MyGraph->ThreeDFlag = 1;
    else
        MyGraph->ThreeDFlag = 0;
    if ((MyGraph->xv[0] == 0) || (MyGraph->yv[0] == 0) ||
        ((MyGraph->zv[0] == 0) && (MyGraph->ThreeDFlag == 1)))
        MyGraph->TimeFlag = 1;
    else
        MyGraph->TimeFlag = 0;
    return;
}

void
get_2d_view(int32 ind) {
    static char *n[] = {"*0X-axis", "*0Y-axis", "Xmin",   "Ymin",
                        "Xmax",     "Ymax",     "Xlabel", "Ylabel"};
    char values[8][MAX_LEN_SBOX];
    int32 status, i;
    int32 i1 = MyGraph->xv[ind], i2 = MyGraph->yv[ind];
    char n1[15], n2[15];
    ind_to_sym(i1, n1);
    ind_to_sym(i2, n2);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%s", n2);
    snprintf(values[2], sizeof(values[2]), "%g", MyGraph->xmin);
    snprintf(values[3], sizeof(values[3]), "%g", MyGraph->ymin);
    snprintf(values[4], sizeof(values[4]), "%g", MyGraph->xmax);
    snprintf(values[5], sizeof(values[5]), "%g", MyGraph->ymax);
    snprintf(values[6], sizeof(values[6]), "%s", MyGraph->xlabel);
    snprintf(values[7], sizeof(values[7]), "%s", MyGraph->ylabel);
    MyGraph->ThreeDFlag = 0;
    status = do_string_box(8, 4, 2, "2D View", n, values, 31);
    if (status != 0) {
        /*  get variable names  */
        find_variable(values[0], &i);
        if (i > -1)
            MyGraph->xv[ind] = i;
        find_variable(values[1], &i);
        if (i > -1)
            MyGraph->yv[ind] = i;

        MyGraph->xmin = atof(values[2]);
        MyGraph->ymin = atof(values[3]);
        MyGraph->xmax = atof(values[4]);
        MyGraph->ymax = atof(values[5]);
        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        snprintf(MyGraph->xlabel, sizeof(MyGraph->xlabel), "%s", values[6]);
        snprintf(MyGraph->ylabel, sizeof(MyGraph->ylabel), "%s", values[7]);
        check_windows();
        /*	      plintf(" x=%d y=%d xlo=%f ylo=%f xhi=%f yhi=%f \n",
                             MyGraph->xv[ind],MyGraph->yv[ind],MyGraph->xlo,
                             MyGraph->ylo,MyGraph->xhi,MyGraph->yhi);
        */
    }
    return;
}

void
axes_opts(void) {
    static char *n[] = {"X-origin",    "Y-origin",   "Z-origin",  "X-org(1=on)",
                        "Y-org(1=on)", "Z-org(1=on", "PSFontSize"};
    char values[7][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%g", MyGraph->xorg);
    snprintf(values[1], sizeof(values[1]), "%g", MyGraph->yorg);
    snprintf(values[2], sizeof(values[2]), "%g", MyGraph->zorg);
    snprintf(values[3], sizeof(values[3]), "%d", MyGraph->xorgflag);
    snprintf(values[4], sizeof(values[4]), "%d", MyGraph->yorgflag);
    snprintf(values[5], sizeof(values[5]), "%d", MyGraph->zorgflag);
    snprintf(values[6], sizeof(values[6]), "%d", PS_FONTSIZE);
    status = do_string_box(7, 7, 1, "Axes options", n, values, 25);
    if (status != 0) {
        MyGraph->xorg = atof(values[0]);
        MyGraph->yorg = atof(values[1]);
        MyGraph->zorg = atof(values[2]);
        MyGraph->xorgflag = atoi(values[3]);
        MyGraph->yorgflag = atoi(values[4]);
        MyGraph->zorgflag = atoi(values[5]);
        PS_FONTSIZE = atoi(values[6]);
        redraw_the_graph();
    }
    return;
}

void
get_3d_view(int32 ind) {
    static char *n[] = {"*0X-axis", "*0Y-axis", "*0Z-axis", "Xmin",
                        "Xmax",     "Ymin",     "Ymax",     "Zmin",
                        "Zmax",     "XLo",      "XHi",      "YLo",
                        "YHi",      "Xlabel",   "Ylabel",   "Zlabel"};
    char values[16][MAX_LEN_SBOX];
    int32 status, i, i1 = MyGraph->xv[ind], i2 = MyGraph->yv[ind],
                     i3 = MyGraph->zv[ind];
    char n1[15], n2[15], n3[15];
    ind_to_sym(i1, n1);
    ind_to_sym(i2, n2);
    ind_to_sym(i3, n3);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%s", n2);
    snprintf(values[2], sizeof(values[2]), "%s", n3);
    snprintf(values[3], sizeof(values[3]), "%g", MyGraph->xmin);
    snprintf(values[5], sizeof(values[5]), "%g", MyGraph->ymin);
    snprintf(values[7], sizeof(values[7]), "%g", MyGraph->zmin);
    snprintf(values[4], sizeof(values[4]), "%g", MyGraph->xmax);
    snprintf(values[6], sizeof(values[6]), "%g", MyGraph->ymax);
    snprintf(values[8], sizeof(values[8]), "%g", MyGraph->zmax);
    snprintf(values[9], sizeof(values[9]), "%g", MyGraph->xlo);
    snprintf(values[11], sizeof(values[11]), "%g", MyGraph->ylo);
    snprintf(values[10], sizeof(values[10]), "%g", MyGraph->xhi);
    snprintf(values[12], sizeof(values[12]), "%g", MyGraph->yhi);
    snprintf(values[13], sizeof(values[13]), "%s", MyGraph->xlabel);
    snprintf(values[14], sizeof(values[14]), "%s", MyGraph->ylabel);
    snprintf(values[15], sizeof(values[15]), "%s", MyGraph->zlabel);
    MyGraph->ThreeDFlag = 1;
    status = do_string_box(16, 6, 3, "3D View", n, values, 31);
    if (status != 0) {
        /*  get variable names  */
        find_variable(values[0], &i);
        if (i > -1)
            MyGraph->xv[ind] = i;
        find_variable(values[1], &i);
        if (i > -1)
            MyGraph->yv[ind] = i;
        find_variable(values[2], &i);
        if (i > -1)
            MyGraph->zv[ind] = i;
        snprintf(MyGraph->xlabel, sizeof(MyGraph->xlabel), "%s", values[13]);
        snprintf(MyGraph->ylabel, sizeof(MyGraph->ylabel), "%s", values[14]);
        snprintf(MyGraph->zlabel, sizeof(MyGraph->zlabel), "%s", values[15]);

        MyGraph->xmin = atof(values[3]);
        MyGraph->ymin = atof(values[5]);
        MyGraph->zmin = atof(values[7]);
        MyGraph->xmax = atof(values[4]);
        MyGraph->ymax = atof(values[6]);
        MyGraph->zmax = atof(values[8]);
        MyGraph->xlo = atof(values[9]);
        MyGraph->ylo = atof(values[11]);
        MyGraph->xhi = atof(values[10]);
        MyGraph->yhi = atof(values[12]);
        check_windows();
        /*      plintf("%f %f %f %f %f %f \n %f %f %f %f",
                     MyGraph->xmin,MyGraph->xmax,
                     MyGraph->ymin,MyGraph->ymax,
                     MyGraph->zmin,MyGraph->zmax,
                     MyGraph->xlo,MyGraph->xhi,
                     MyGraph->ylo,MyGraph->yhi);
    */
    }
    return;
}

void
check_val(double *x1, double *x2, double *xb, double *xd) {
    double temp;

    /*
      see get_max for details
    */

    if (*x1 == *x2) {
        temp = .05 * lmax(fabs(*x1), 1.0);
        *x1 = *x1 - temp;
        *x2 = *x2 + temp;
    }
    if (*x1 > *x2) {
        temp = *x2;
        *x2 = *x1;
        *x1 = temp;
    }
    *xb = .5 * (*x1 + *x2);
    *xd = 2.0 / (*x2 - *x1);
    return;
}

void
get_max(int32 index, double *vmin, double *vmax) {
    float x0, x1, z;
    double temp;
    int32 i;
    x0 = my_browser.data[index][0];
    x1 = x0;
    for (i = 0; i < my_browser.maxrow; i++) {
        z = my_browser.data[index][i];
        if (z < x0)
            x0 = z;
        if (z > x1)
            x1 = z;
    }
    *vmin = (double)x0;
    *vmax = (double)x1;
    if (fabs(*vmin - *vmax) < REAL_SMALL) {
        temp = .05 * lmax(fabs(*vmin), 1.0);
        *vmin = *vmin - temp;
        *vmax = *vmax + temp;
    }
    return;
}

void
corner_cube(double *xlo, double *xhi, double *ylo, double *yhi) {
    float x, y;
    float x1, x2, y1, y2;
    threedproj(-1., -1., -1., &x, &y);
    x1 = x;
    x2 = x;
    y1 = y;
    y2 = y;
    threedproj(-1., -1., 1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(-1., 1., -1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(-1., 1., 1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(1., -1., -1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(1., -1., 1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(1., 1., 1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    threedproj(1., 1., -1., &x, &y);
    if (x < x1)
        x1 = x;
    if (x > x2)
        x2 = x;
    if (y < y1)
        y1 = y;
    if (y > y2)
        y2 = y;
    *xlo = x1;
    *ylo = y1;
    *xhi = x2;
    *yhi = y2;
    return;
}

void
default_window(void) {
    if (MyGraph->ThreeDFlag) {
        MyGraph->xmax = x_3d[1];
        MyGraph->ymax = y_3d[1];
        MyGraph->zmax = z_3d[1];
        MyGraph->xmin = x_3d[0];
        MyGraph->ymin = y_3d[0];
        MyGraph->zmin = z_3d[0];

        corner_cube(&(MyGraph->xlo), &(MyGraph->xhi), &(MyGraph->ylo),
                    &(MyGraph->yhi));
        check_windows();
    } else {
        MyGraph->xmax = x_3d[1];
        MyGraph->ymax = y_3d[1];

        MyGraph->xmin = x_3d[0];
        MyGraph->ymin = y_3d[0];

        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        check_windows();
    }

    redraw_the_graph();
}

void
fit_window(void) {
    double Mx = -1.e25, My = -1.e25, Mz = -1.e25, mx = -Mx, my = -My, mz = -Mz;
    int32 i, n = MyGraph->nvars;
    if (storind < 2)
        return;
    if (MyGraph->ThreeDFlag) {
        for (i = 0; i < n; i++) {

            get_max(MyGraph->xv[i], &(MyGraph->xmin), &(MyGraph->xmax));
            Mx = lmax(MyGraph->xmax, Mx);
            mx = -lmax(-MyGraph->xmin, -mx);

            get_max(MyGraph->yv[i], &(MyGraph->ymin), &(MyGraph->ymax));
            My = lmax(MyGraph->ymax, My);
            my = -lmax(-MyGraph->ymin, -my);

            get_max(MyGraph->zv[i], &(MyGraph->zmin), &(MyGraph->zmax));
            Mz = lmax(MyGraph->zmax, Mz);
            mz = -lmax(-MyGraph->zmin, -mz);
        }
        MyGraph->xmax = Mx;
        MyGraph->ymax = My;
        MyGraph->zmax = Mz;
        MyGraph->xmin = mx;
        MyGraph->ymin = my;
        MyGraph->zmin = mz;

        corner_cube(&(MyGraph->xlo), &(MyGraph->xhi), &(MyGraph->ylo),
                    &(MyGraph->yhi));
        check_windows();
    } else {
        for (i = 0; i < n; i++) {
            get_max(MyGraph->xv[i], &(MyGraph->xmin), &(MyGraph->xmax));
            Mx = lmax(MyGraph->xmax, Mx);
            mx = -lmax(-MyGraph->xmin, -mx);

            get_max(MyGraph->yv[i], &(MyGraph->ymin), &(MyGraph->ymax));
            My = lmax(MyGraph->ymax, My);
            my = -lmax(-MyGraph->ymin, -my);
        }
        MyGraph->xmax = Mx;
        MyGraph->ymax = My;

        MyGraph->xmin = mx;
        MyGraph->ymin = my;

        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        check_windows();
    }
    redraw_the_graph();
}

void
check_windows(void) {
    double zip, zap;
    check_val(&MyGraph->xmin, &MyGraph->xmax, &MyGraph->xbar, &MyGraph->dx);
    check_val(&MyGraph->ymin, &MyGraph->ymax, &MyGraph->ybar, &MyGraph->dy);
    check_val(&MyGraph->zmin, &MyGraph->zmax, &MyGraph->zbar, &MyGraph->dz);
    check_val(&MyGraph->xlo, &MyGraph->xhi, &zip, &zap);
    check_val(&MyGraph->ylo, &MyGraph->yhi, &zip, &zap);
    return;
}

void
user_window(void) {
    static char *n[] = {"X Lo", "X Hi", "Y Lo", "Y Hi"};
    char values[4][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%g", MyGraph->xlo);
    snprintf(values[2], sizeof(values[2]), "%g", MyGraph->ylo);
    snprintf(values[1], sizeof(values[1]), "%g", MyGraph->xhi);
    snprintf(values[3], sizeof(values[3]), "%g", MyGraph->yhi);
    status = do_string_box(4, 2, 2, "Window", n, values, 28);
    if (status != 0) {

        MyGraph->xlo = atof(values[0]);
        MyGraph->ylo = atof(values[2]);
        MyGraph->xhi = atof(values[1]);
        MyGraph->yhi = atof(values[3]);
        if (MyGraph->grtype < 5) {
            MyGraph->xmin = MyGraph->xlo;
            MyGraph->xmax = MyGraph->xhi;
            MyGraph->ymin = MyGraph->ylo;
            MyGraph->ymax = MyGraph->yhi;
        }
        check_windows();
    }
    redraw_the_graph();
}

void
xi_vs_t(void) {
    /*  a short cut   */
    char name[20], value[20];
    int32 i = MyGraph->yv[0];

    ind_to_sym(i, value);
    snprintf(name, sizeof(name), "Plot vs t: ");
    new_string(name, value);
    find_variable(value, &i);

    if (i > -1) {
        MyGraph->yv[0] = i;
        MyGraph->grtype = 0;
        MyGraph->xv[0] = 0;
        if (storind >= 2) {
            get_max(MyGraph->xv[0], &(MyGraph->xmin), &(MyGraph->xmax));
            get_max(MyGraph->yv[0], &(MyGraph->ymin), &(MyGraph->ymax));

        } else {
            MyGraph->xmin = T0;
            MyGraph->xmax = TEND;
        }
        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        check_windows();
        check_flags();
        set_normal_scale();
        redraw_the_graph();
    }
    return;
}

void
redraw_the_graph(void) {
    blank_screen(draw_win);
    set_normal_scale();
    do_axes();
    hi_lite(draw_win);
    restore(0, my_browser.maxrow);
    draw_label(draw_win);
    draw_freeze(draw_win);
    redraw_dfield();
    if (MyGraph->Nullrestore)
        restore_nullclines();
    return;
}

void
movie_rot(double start, double increment, int32 nclip, int32 angle) {
    int32 i;
    double thetaold = MyGraph->Theta, phiold = MyGraph->Phi;
    reset_film();
    for (i = 0; i <= nclip; i++) {

        if (angle == 0)
            make_rot(start + i*increment, phiold);
        else
            make_rot(thetaold, start + i*increment);
        redraw_the_graph();
        film_clip();
    }
    MyGraph->Theta = thetaold;
    MyGraph->Phi = phiold;
    return;
}

void
test_rot(void) {
    int32 done = 0;
    int32 kp;
    XEvent ev;
    double theta = MyGraph->Theta, phi = MyGraph->Phi;
    redraw_cube(theta, phi);
    while (done == 0) {
        XNextEvent(display, &ev);
        if (ev.type == KeyPress) {
            kp = get_key_press(&ev);
            switch (kp) {
            case UP:
                phi = phi + 1;
                redraw_cube(theta, phi);
                break;
            case DOWN:
                phi = phi - 1;
                redraw_cube(theta, phi);
                break;
            case LEFT:
                theta = theta + 1;
                redraw_cube(theta, phi);
                break;
            case RIGHT:
                theta = theta - 1;
                redraw_cube(theta, phi);
                break;
            case FINE:
                done = 1;
                break;
            case ESC:
                done = -1;
                break;
            }
        }
    }
    if (done == 1) {
        MyGraph->Phi = phi;
        MyGraph->Theta = theta;
    }
    redraw_the_graph();
}

void
get_3d_par_com(void) {

    static char *n[] = {"Persp (1=On)",
                        "ZPlane",
                        "ZView",
                        "Theta",
                        "Phi",
                        "Movie(Y/N)",
                        "Vary (theta/phi)",
                        "Start angle",
                        "Increment",
                        "Number increments"};
    char values[10][MAX_LEN_SBOX];
    int32 status;

    int32 nclip = 8, angle = 0;
    double start, increment = 45;
    if (MyGraph->grtype < 5)
        return;

    snprintf(values[0], sizeof(values[0]), "%d", MyGraph->PerspFlag);
    snprintf(values[1], sizeof(values[1]), "%g", MyGraph->ZPlane);
    snprintf(values[2], sizeof(values[2]), "%g", MyGraph->ZView);
    snprintf(values[3], sizeof(values[3]), "%g", MyGraph->Theta);
    snprintf(values[4], sizeof(values[4]), "%g", MyGraph->Phi);
    snprintf(values[5], sizeof(values[5]), "%s", mov3d.yes);
    snprintf(values[6], sizeof(values[6]), "%s", mov3d.angle);
    snprintf(values[7], sizeof(values[7]), "%g", mov3d.start);
    snprintf(values[8], sizeof(values[8]), "%g", mov3d.incr);
    snprintf(values[9], sizeof(values[9]), "%d", mov3d.nclip);

    status = do_string_box(10, 5, 2, "3D Parameters", n, values, 28);
    if (status != 0) {
        MyGraph->PerspFlag = atoi(values[0]);
        MyGraph->ZPlane = atof(values[1]);
        MyGraph->ZView = atof(values[2]);
        MyGraph->Theta = atof(values[3]);
        MyGraph->Phi = atof(values[4]);
        if (values[5][0] == 'y' || values[5][0] == 'Y') {
            strcpy(mov3d.yes, values[5]);
            strcpy(mov3d.angle, values[6]);
            start = atof(values[7]);
            increment = atof(values[8]);
            nclip = atoi(values[9]);
            mov3d.start = start;
            mov3d.incr = increment;
            mov3d.nclip = nclip;
            angle = 0;
            if (mov3d.angle[0] == 'p' || mov3d.angle[0] == 'P')
                angle = 1;
            /*     XRaiseWindow(display,MyGraph->w); */
            movie_rot(start, increment, nclip, angle);
        }

        make_rot(MyGraph->Theta, MyGraph->Phi);
        /*  Redraw the picture   */
        redraw_the_graph();
    }
    return;
}

void
get_3d_par_noper(void) {

    static char *n[] = {
        "Theta",       "Phi",       "Movie(Y/N)",       "Vary (theta/phi)",
        "Start angle", "Increment", "Number increments"};
    char values[10][MAX_LEN_SBOX];
    int32 status;

    int32 nclip = 8, angle = 0;
    double start, increment = 45;
    if (MyGraph->grtype < 5)
        return;

    snprintf(values[0], sizeof(values[0]), "%g", MyGraph->Theta);
    snprintf(values[1], sizeof(values[1]), "%g", MyGraph->Phi);
    snprintf(values[2], sizeof(values[2]), "%s", mov3d.yes);
    snprintf(values[3], sizeof(values[3]), "%s", mov3d.angle);
    snprintf(values[4], sizeof(values[4]), "%g", mov3d.start);
    snprintf(values[5], sizeof(values[5]), "%g", mov3d.incr);
    snprintf(values[6], sizeof(values[6]), "%d", mov3d.nclip);

    status = do_string_box(7, 7, 1, "3D Parameters", n, values, 28);
    if (status != 0) {
        /* MyGraph->PerspFlag=atoi(values[0]);
                   MyGraph->ZPlane=atof(values[1]);
                   MyGraph->ZView=atof(values[2]); */
        MyGraph->Theta = atof(values[0]);
        MyGraph->Phi = atof(values[1]);
        if (values[2][0] == 'y' || values[2][0] == 'Y') {
            strcpy(mov3d.yes, values[2]);
            strcpy(mov3d.angle, values[3]);
            start = atof(values[4]);
            increment = atof(values[5]);
            nclip = atoi(values[6]);
            mov3d.start = start;
            mov3d.incr = increment;
            mov3d.nclip = nclip;
            angle = 0;
            if (mov3d.angle[0] == 'p' || mov3d.angle[0] == 'P')
                angle = 1;
            /*     XRaiseWindow(display,MyGraph->w); */
            movie_rot(start, increment, nclip, angle);
        }

        make_rot(MyGraph->Theta, MyGraph->Phi);
        /*  Redraw the picture   */
        redraw_the_graph();
    }
    return;
}

void
update_view(float xlo, float xhi, float ylo, float yhi) {
    MyGraph->xlo = xlo;
    MyGraph->ylo = ylo;
    MyGraph->xhi = xhi;
    MyGraph->yhi = yhi;
    if (MyGraph->grtype < 5) {
        MyGraph->xmin = MyGraph->xlo;
        MyGraph->xmax = MyGraph->xhi;
        MyGraph->ymin = MyGraph->ylo;
        MyGraph->ymax = MyGraph->yhi;
    }
    check_windows();

    redraw_the_graph();
}

void
scroll_window(void) {
    XEvent ev;
    int32 i = 0, j = 0;
    int32 state = 0;
    float x, y, x0, y0;
    float xlo = MyGraph->xlo;
    float ylo = MyGraph->ylo;
    float xhi = MyGraph->xhi;
    float yhi = MyGraph->yhi;
    float dx = 0;
    float dy = 0;
    int32 alldone = 0;
    XSelectInput(display, draw_win,
                 KeyPressMask | ButtonPressMask | ButtonReleaseMask |
                     PointerMotionMask | ButtonMotionMask | ExposureMask);
    while (!alldone) {
        XNextEvent(display, &ev);
        switch (ev.type) {
        case KeyPress:
            alldone = 1;
            break;
        case Expose:
            do_expose(ev);
            break;
        case ButtonPress:
            if (state == 0) {
                i = ev.xkey.x;
                j = ev.xkey.y;
                scale_to_real(i, j, &x0, &y0);
                state = 1;
            }
            break;
        case MotionNotify:
            if (state == 1) {
                i = ev.xmotion.x;
                j = ev.xmotion.y;
                scale_to_real(i, j, &x, &y);
                dx = -(x - x0) / 2;
                dy = -(y - y0) / 2;

                update_view(xlo + dx, xhi + dx, ylo + dy, yhi + dy);
            }
            break;
        case ButtonRelease:
            state = 0;
            xlo = xlo + dx;
            xhi = xhi + dx;
            ylo = ylo + dy;
            yhi = yhi + dy;
            break;
        }
    }
    return;
}

void
window_zoom_com(int32 c) {
    int32 i1, i2, j1, j2;
    switch (c) {
    case 0:
        user_window();
        break;
    case 1:

        /*  XSelectInput(display,w,
    KeyPressMask|ButtonPressMask|ButtonReleaseMask|
             PointerMotionMask|ButtonMotionMask|ExposureMask);
                while(1)
                {
                     XNextEvent(display,&ev);
                     switch(ev.type){
                } */
        if (rubber(&i1, &j1, &i2, &j2, draw_win, RUBBOX) == 0)
            break;
        zoom_in(i1, j1, i2, j2);

        break;
    case 2:
        if (rubber(&i1, &j1, &i2, &j2, draw_win, RUBBOX) == 0)
            break;
        zoom_out(i1, j1, i2, j2);
        break;
    case 3:
        fit_window();
        break;
    case 4:
        default_window();
        break;
    case 5:
        scroll_window();
        break;
    }
    set_normal_scale();
    return;
}

void
zoom_in(int32 i1, int32 j1, int32 i2, int32 j2) {
    float x1, y1, x2, y2;
    float dx = MyGraph->xhi - MyGraph->xlo;
    float dy = MyGraph->yhi - MyGraph->ylo;
    scale_to_real(i1, j1, &x1, &y1);
    scale_to_real(i2, j2, &x2, &y2);
    if (x1 == x2 || y1 == y2) {

        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        dx = dx / 2;
        dy = dy / 2;

        /*Shrink by thirds and center (track) about the point clicked*/
        MyGraph->xlo = x1 - dx / 2;
        MyGraph->xhi = x1 + dx / 2;

        MyGraph->ylo = y1 - dy / 2;
        MyGraph->yhi = y1 + dy / 2;

    } else {
        MyGraph->xlo = x1;
        MyGraph->ylo = y1;
        MyGraph->xhi = x2;
        MyGraph->yhi = y2;
    }
    if (MyGraph->grtype < 5) {
        MyGraph->xmin = MyGraph->xlo;
        MyGraph->xmax = MyGraph->xhi;
        MyGraph->ymin = MyGraph->ylo;
        MyGraph->ymax = MyGraph->yhi;
    }
    check_windows();
    redraw_the_graph();
    draw_help();
    return;
}

void
zoom_out(int32 i1, int32 j1, int32 i2, int32 j2) {

    float x1, y1, x2, y2;
    float bx, mux, by, muy;
    float dx = MyGraph->xhi - MyGraph->xlo;
    float dy = MyGraph->yhi - MyGraph->ylo;
    scale_to_real(i1, j1, &x1, &y1);
    scale_to_real(i2, j2, &x2, &y2);

    /*
    if(x1==x2||y1==y2)return;
    */
    /*
    plintf("%f %f %f %f \n ",x1,y1,x2,y2);
    plintf("%f %f %f %f
    \n",MyGraph->xlo,MyGraph->ylo,MyGraph->xhi,MyGraph->yhi);
   */
    if (x1 == x2 || y1 == y2) {

        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        /*Grow by thirds and center (track) about the point clicked*/
        dx = dx * 2;
        dy = dy * 2;

        MyGraph->xlo = x1 - dx / 2;
        MyGraph->xhi = x1 + dx / 2;

        MyGraph->ylo = y1 - dy / 2;
        MyGraph->yhi = y1 + dy / 2;
    } else {
        if (x1 > x2) {
            bx = x1;
            x1 = x2;
            x2 = bx;
        }
        if (y1 > y2) {
            by = y1;
            y1 = y2;
            y2 = by;
        }

        bx = dx*dx / (x2 - x1);
        mux = (x1 - MyGraph->xlo) / dx;
        MyGraph->xlo = MyGraph->xlo - bx*mux;
        MyGraph->xhi = MyGraph->xlo + bx;

        by = dy*dy / (y2 - y1);
        muy = (y1 - MyGraph->ylo) / dy;
        MyGraph->ylo = MyGraph->ylo - by*muy;
        MyGraph->yhi = MyGraph->ylo + by;
    }
    if (MyGraph->grtype < 5) {
        MyGraph->xmin = MyGraph->xlo;
        MyGraph->xmax = MyGraph->xhi;
        MyGraph->ymin = MyGraph->ylo;
        MyGraph->ymax = MyGraph->yhi;
    }
    check_windows();
    redraw_the_graph();
    draw_help();
    return;
}

void
graph_all(int32 *list, int32 n, int32 type) {
    int32 i;
    if (type == 0) {
        for (i = 0; i < n; i++) {
            MyGraph->xv[i] = 0;
            MyGraph->yv[i] = list[i];
            MyGraph->line[i] = MyGraph->line[0];
            MyGraph->color[i] = i;
        }
        MyGraph->nvars = n;
        MyGraph->grtype = 0;
        MyGraph->ThreeDFlag = 0;
    }
    if (type == 1) {
        MyGraph->nvars = 1;
        MyGraph->xv[0] = list[0];
        MyGraph->yv[0] = list[1];
        MyGraph->grtype = 0;
        MyGraph->ThreeDFlag = 0;
        if (n == 3) {
            MyGraph->zv[0] = list[2];
            MyGraph->grtype = 5;
            MyGraph->ThreeDFlag = 1;
        }
    }
    check_flags();
    fit_window();
    return;
}

int32
find_color(int32 in) {
    int32 i;
    for (i = 0; i <= 10; i++)
        if (in == colorline[i])
            return i;
    return 0;
}

int32
alter_curve(char *title, int32 in_it, int32 n) {
    static char *nn[] = {"*0X-axis", "*0Y-axis", "*0Z-axis", "*4Color",
                         "Line type"};
    char values[5][MAX_LEN_SBOX];
    int32 status, i;
    int32 i1 = MyGraph->xv[in_it], i2 = MyGraph->yv[in_it],
          i3 = MyGraph->zv[in_it];
    char n1[15], n2[15], n3[15];

    ind_to_sym(i1, n1);
    ind_to_sym(i2, n2);
    ind_to_sym(i3, n3);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%s", n2);
    snprintf(values[2], sizeof(values[2]), "%s", n3);
    snprintf(values[3], sizeof(values[3]), "%d", MyGraph->color[in_it]);
    snprintf(values[4], sizeof(values[4]), "%d", MyGraph->line[in_it]);
    status = do_string_box(5, 5, 1, title, nn, values, 25);
    if (status != 0) {
        find_variable(values[0], &i);
        if (i > -1)
            MyGraph->xv[n] = i;
        find_variable(values[1], &i);
        if (i > -1)
            MyGraph->yv[n] = i;
        find_variable(values[2], &i);
        if (i > -1)
            MyGraph->zv[n] = i;

        MyGraph->line[n] = atoi(values[4]);
        i = atoi(values[3]);
        if (i < 0 || i > 10)
            i = 0;
        MyGraph->color[n] = i;

        return 1;
    }
    return 0;
}

void
edit_curve(void) {
    char bob[20];
    int32 crv = 0;
    snprintf(bob, sizeof(bob), "Edit 0-%d :", MyGraph->nvars - 1);
    ping();
    new_int(bob, &crv);
    if (crv >= 0 && crv < MyGraph->nvars) {
        snprintf(bob, sizeof(bob), "Edit curve %d", crv);
        alter_curve(bob, crv, crv);
    }
    return;
}

void
new_curve(void) {
    if (alter_curve("New Curve", 0, MyGraph->nvars))
        MyGraph->nvars = MyGraph->nvars + 1;
    return;
}

void
create_ps(void) {
    /*char filename[256];*/
    char filename[XPP_MAX_NAME];
    static char *nn[] = {"BW-0/Color-1", "Land(0)/Port(1)", "Axes fontsize",
                         "Font", "Linewidth"};
    int32 status;
    char values[5][MAX_LEN_SBOX];
    snprintf(values[0], sizeof(values[0]), "%d", PS_Color);
    snprintf(values[1], sizeof(values[1]), "%d", PS_Port);
    snprintf(values[2], sizeof(values[2]), "%d", PS_FONTSIZE);
    snprintf(values[3], sizeof(values[3]), "%s", PS_FONT);
    snprintf(values[4], sizeof(values[4]), "%g", PS_LW);
    status = do_string_box(5, 5, 1, "Postscript parameters", nn, values, 25);
    if (status != 0) {
        PS_Color = atoi(values[0]);
        PS_Port = atoi(values[1]);
        PS_FONTSIZE = atoi(values[2]);
        PS_LW = atof(values[4]);
        snprintf(PS_FONT, sizeof(PS_FONT), "%s", values[3]);
        snprintf(filename, sizeof(filename), "%s.ps", this_file);
        ping();

        if (!file_selector("Print postscript", filename, "*.ps"))
            return;
        if (ps_init(filename, PS_Color)) {
            ps_restore();
            ping();
        }
    }
    return;
}

void
dump_ps(int32 i) {
    char filename[XPP_MAX_NAME];
    if (i < 0) {
        snprintf(filename, sizeof(filename), "%s%s.%s", this_file,
                 this_internset, PlotFormat);
    } else {
        snprintf(filename, sizeof(filename), "%s%s_%04d.%s", this_file,
                 this_internset, i, PlotFormat);
    }

    if (strcmp(PlotFormat, "ps") == 0) {
        if (ps_init(filename, PS_Color)) {
            ps_restore();
        }
    } else if (strcmp(PlotFormat, "svg") == 0) {
        if (svg_init(filename, PS_Color)) {
            svg_restore();
        }
    }
    return;
}

void
create_svg(void) {

    char filename[XPP_MAX_NAME];
    strcpy(filename, this_file);
    filename[strlen(filename) - 4] = '\0';
    strcat(filename, ".svg");
    /*snprintf(filename, sizeof(filename),"%s.svg",tmp);*/
    if (!file_selector("Print svg", filename, "*.svg"))
        return;
    if (svg_init(filename, PS_Color)) {
        svg_restore();
        ping();
    }
    return;
}

/*
ps_test()
{
 double xlo=MyGraph->xlo,xhi=MyGraph->xhi,ylo=MyGraph->ylo,yhi=MyGraph->yhi;
 text_abs((float)xlo,(float)ylo,"lolo");
 text_abs((float)xlo,(float)yhi,"lohi");
 text_abs((float)xhi,(float)ylo,"hilo");
 text_abs((float)xhi,(float)yhi,"hihi");
 ps_end();
}

*/

void
change_cmap_com(int32 i) {
    NewColormap(i);
    return;
}

void
freeze_com(int32 c) {

    switch (c) {
    case 0:
        freeze_crv(0);
        break;
    case 1:
        delete_frz();
        break;
    case 2:
        edit_frz();
        break;
    case 3:
        kill_frz();
        break;
        /*case 4:
        key_frz();
        break; */
    case 5:
        frz_bd();
        break;
    case 6:
        free_bd();
        break;
    case 7:
        AutoFreezeFlag = 1 - AutoFreezeFlag;
        break;
    }
    return;
}

void
set_key(int32 x, int32 y) {
    float xp, yp;
    scale_to_real(x, y, &xp, &yp);
    FreezeKeyX = xp;
    FreezeKeyY = yp;
    FreezeKeyFlag = 1;
    return;
}

void
draw_freeze_key(void) {
    int32 ix, iy;
    int32 i, y0;
    int32 ix2;
    int32 dy = 2 * HChar;
    if (FreezeKeyFlag == SCRNFMT)
        return;
    if (PltFmtFlag == PSFMT)
        dy = -dy;
    scale_to_screen((float)FreezeKeyX, (float)FreezeKeyY, &ix, &iy);
    ix2 = ix + 4 * HChar;
    y0 = iy;
    for (i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && frz[i].w == draw_win && strlen(frz[i].key) > 0) {
            set_linestyle(abs(frz[i].color));
            line(ix, y0, ix2, y0);
            set_linestyle(0);
            put_text(ix2 + HChar, y0, frz[i].key);
            y0 += dy;
        }
    }
    return;
}

void
key_frz_com(int32 c) {
    int32 x, y;
    switch (c) {
    case 0:
        FreezeKeyFlag = 0;
        break;
    case 1:
        MessageBox("Position with mouse");
        if (get_mouse_xy(&x, &y, draw_win)) {
            set_key(x, y);
            draw_freeze_key();
        }
        KillMessageBox();
    }
    return;
}

void
edit_frz(void) {
    int32 i;
    i = get_frz_index(draw_win);
    if (i < 0)
        return;
    edit_frz_crv(i);
    return;
}

void
delete_frz_crv(int32 i) {
    if (frz[i].use == 0)
        return;
    frz[i].use = 0;
    frz[i].name[0] = 0;
    frz[i].key[0] = 0;
    free(frz[i].xv);
    free(frz[i].yv);
    if (frz[i].type > 0)
        free(frz[i].zv);
    return;
}

void
delete_frz(void) {
    int32 i;
    i = get_frz_index(draw_win);
    if (i < 0)
        return;
    delete_frz_crv(i);
    return;
}

void
kill_frz(void) {
    int32 i;
    for (i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && frz[i].w == draw_win)
            delete_frz_crv(i);
    }
    return;
}

int32
freeze_crv(int32 ind) {
    int32 i;
    i = create_crv(ind);
    if (i < 0)
        return -1;
    edit_frz_crv(i);
    return 1;
}

void
auto_freeze_it(void) {
    if (AutoFreezeFlag == 0)
        return;
    create_crv(0);
    return;
}

int32
create_crv(int32 ind) {
    int32 i, type, j;
    int32 ix, iy, iz;

    for (i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 0) {
            ix = MyGraph->xv[ind];
            iy = MyGraph->yv[ind];
            iz = MyGraph->zv[ind];
            if (my_browser.maxrow <= 2) {
                err_msg("No Curve to freeze");
                return -1;
            }
            frz[i].xv = malloc(sizeof(float) * my_browser.maxrow);
            frz[i].yv = malloc(sizeof(float) * my_browser.maxrow);
            if ((type = MyGraph->grtype) > 0)
                frz[i].zv = malloc(sizeof(float) * my_browser.maxrow);
            if ((type > 0 && frz[i].zv == NULL) ||
                (type == 0 && frz[i].yv == NULL)) {
                err_msg("Cant allocate storage for curve");
                return -1;
            }
            frz[i].use = 1;
            frz[i].len = my_browser.maxrow;
            for (j = 0; j < my_browser.maxrow; j++) {
                frz[i].xv[j] = my_browser.data[ix][j];
                frz[i].yv[j] = my_browser.data[iy][j];
                if (type > 0)
                    frz[i].zv[j] = my_browser.data[iz][j];
            }
            frz[i].type = type;
            frz[i].w = draw_win;
            snprintf(frz[i].name, sizeof(frz[i].name), "crv%c", 'a' + i);
            strncpy(frz[i].key, frz[i].name, sizeof(frz[i].key));
            return i;
        }
    }
    err_msg("All curves used");
    return -1;
}

void
edit_frz_crv(int32 i) {
    static char *nn[] = {"*4Color", "Key", "Name"};
    char values[3][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", frz[i].color);
    snprintf(values[1], sizeof(values[1]), "%s", frz[i].key);
    snprintf(values[2], sizeof(values[2]), "%s", frz[i].name);
    status = do_string_box(3, 3, 1, "Edit Freeze", nn, values, 25);
    if (status != 0) {
        frz[i].color = atoi(values[0]);
        snprintf(frz[i].key, sizeof(frz[i].key), "%s", values[1]);
        snprintf(frz[i].name, sizeof(frz[i].name), "%s", values[2]);
    }
    return;
}

void
draw_frozen_cline(int32 index, Window w) {
    if (nclines[index].use == 0 || nclines[index].w != w)
        return;
    return;
}

void
draw_freeze(Window w) {
    int32 i, j, type = MyGraph->grtype, lt = 0;
    float oldxpl, oldypl, oldzpl = 0.0, xpl, ypl, zpl = 0.0;
    float *xv, *yv, *zv;
    for (i = 0; i < MAXNCLINE; i++)
        draw_frozen_cline(i, w);
    for (i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && frz[i].w == w && frz[i].type == type) {
            if (frz[i].color < 0) {
                set_linestyle(-frz[i].color);
                lt = 1;
            } else
                set_linestyle(frz[i].color);
            xv = frz[i].xv;
            yv = frz[i].yv;
            zv = frz[i].zv;
            oldxpl = xv[0];
            oldypl = yv[0];
            if (type > 0)
                oldzpl = zv[0];
            for (j = 0; j < frz[i].len; j++) {
                xpl = xv[j];
                ypl = yv[j];
                if (type > 0)
                    zpl = zv[j];
                if (lt == 0) {
                    if (type == 0)
                        line_abs(oldxpl, oldypl, xpl, ypl);
                    else
                        line_3d(oldxpl, oldypl, oldzpl, xpl, ypl, zpl);
                } else {
                    if (type == 0)
                        point_abs(xpl, ypl);
                    else
                        point_3d(xpl, ypl, zpl);
                }
                oldxpl = xpl;
                oldypl = ypl;
                if (type > 0)
                    oldzpl = zpl;
            }
        }
    }
    draw_freeze_key();
    draw_bd(w);
    return;
}

/*  Bifurcation curve importing */

void
init_bd(void) {
    my_bd.nbifcrv = 0;
    return;
}

void
draw_bd(Window w) {
    int32 i, j, len;
    float oldxpl, oldypl, xpl, ypl, *x, *y;
    if (w == my_bd.w && my_bd.nbifcrv > 0) {
        for (i = 0; i < my_bd.nbifcrv; i++) {
            set_linestyle(my_bd.color[i]);
            len = my_bd.npts[i];
            x = my_bd.x[i];
            y = my_bd.y[i];
            xpl = x[0];
            ypl = y[0];
            for (j = 0; j < len; j++) {
                oldxpl = xpl;
                oldypl = ypl;
                xpl = x[j];
                ypl = y[j];
                line_abs(oldxpl, oldypl, xpl, ypl);
            }
        }
    }
    return;
}

void
free_bd(void) {
    int32 i;
    if (my_bd.nbifcrv > 0) {
        for (i = 0; i < my_bd.nbifcrv; i++) {
            free(my_bd.x[i]);
            free(my_bd.y[i]);
        }
        my_bd.nbifcrv = 0;
    }
    return;
}

void
add_bd_crv(float *x, float *y, int32 len, int32 type, int32 ncrv) {
    int32 i;
    if (ncrv >= MAXBIFCRV)
        return;
    my_bd.x[ncrv] = malloc(sizeof(float) * len);
    my_bd.y[ncrv] = malloc(sizeof(float) * len);
    for (i = 0; i < len; i++) {
        my_bd.x[ncrv][i] = x[i];
        my_bd.y[ncrv][i] = y[i];
    }
    my_bd.npts[ncrv] = len;
    i = lsSEQ;
    if (type == UPER)
        i = lsUPER;
    if (type == SPER)
        i = lsSPER;
    if (type == UEQ)
        i = lsUEQ;
    my_bd.color[ncrv] = i;
    return;
}

void
frz_bd(void) {
    FILE *fp;
    /*char filename[256];*/
    char filename[XPP_MAX_NAME];
    snprintf(filename, sizeof(filename), "diagram.dat");
    ping();
    if (!file_selector("Import Diagram", filename, "*.dat"))
        return;
    /* if(new_string("Diagram to import: ",filename)==0)return; */
    if ((fp = fopen(filename, "r")) == NULL) {
        err_msg("Couldn't open file");
        return;
    }
    read_bd(fp);
}

void
read_bd(FILE *fp) {
    int32 oldtype, type, oldbr, br, ncrv = 0, len, f2;
    float x[8000], ylo[8000], yhi[8000];
    len = 0;
    fscanf(fp, "%g %g %g %d %d %d", &x[len], &ylo[len], &yhi[len], &oldtype,
           &oldbr, &f2);
    len++;
    while (!feof(fp)) {
        fscanf(fp, "%g %g %g %d %d %d", &x[len], &ylo[len], &yhi[len], &type,
               &br, &f2);
        if (type == oldtype && br == oldbr)
            len++;
        else {
            /* if(oldbr==br)len++; */ /* extend to point of instability */
            add_bd_crv(x, ylo, len, oldtype, ncrv);
            ncrv++;
            if (oldtype == UPER || oldtype == SPER) {
                add_bd_crv(x, yhi, len, oldtype, ncrv);
                ncrv++;
            }
            if (oldbr == br)
                len--;
            x[0] = x[len];
            ylo[0] = ylo[len];
            yhi[0] = yhi[len];

            len = 1;
        }
        oldbr = br;
        oldtype = type;
    }
    /*  save this last one */
    if (len > 1) {
        add_bd_crv(x, ylo, len, oldtype, ncrv);
        ncrv++;
        if (oldtype == UPER || oldtype == SPER) {
            add_bd_crv(x, yhi, len, oldtype, ncrv);
            ncrv++;
        }
    }
    plintf(" got %d bifurcation curves\n", ncrv);
    fclose(fp);
    my_bd.nbifcrv = ncrv;
    my_bd.w = draw_win;
    return;
}

int32
get_frz_index(Window w) {
    char *n[MAXFRZ];
    char key[MAXFRZ], ch;

    int32 i;
    int32 count = 0;
    Window temp = main_win;
    for (i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && w == frz[i].w) {
            n[count] = malloc(20);
            sprintf(n[count], "%s", frz[i].name);
            key[count] = 'a' + i;

            count++;
        }
    }
    if (count == 0)
        return -1;
    key[count] = 0;
    ch = (char)pop_up_list(&temp, "Curves", n, key, count, 12, 0, 10,
                           8 * DCURY + 8, no_hint, info_pop, info_message);
    for (i = 0; i < count; i++)
        free(n[i]);
    return (int32)(ch - 'a');
}

void
export_graf_data(void) {
    FILE *fp;
    /*char filename[256];*/
    char filename[XPP_MAX_NAME];
    snprintf(filename, sizeof(filename), "curve.dat");
    ping();
    if (!file_selector("Export graph data", filename, "*.dat"))
        return;
    /* if(new_string("Data filename:",filename)==0)return; */
    if ((fp = fopen(filename, "w")) == NULL) {
        err_msg("Couldn't open file");
        return;
    }
    export_data(fp);
    fclose(fp);
    return;
}

void
add_a_curve_com(int32 c) {

    switch (c) {
    case 0:
        if (MyGraph->nvars >= MAXPERPLOT) {
            err_msg("Too many plots!");
            return;
        }
        new_curve();
        break;
    case 1:
        if (MyGraph->nvars > 1)
            MyGraph->nvars = MyGraph->nvars - 1;
        break;
    case 2:
        MyGraph->nvars = 1;
        break;
    case 3:
        edit_curve();
        break;
    case 4:
        create_ps();
        break;
    case 5:
        create_svg();
        break;
        /* case 6: freeze();
           break; */
    case 7:
        axes_opts();
        break;
    case 8:
        export_graf_data();
        break;
        /*  case 9: change_cmap();
            break; */
    }
    check_flags();
    redraw_the_graph();
}
