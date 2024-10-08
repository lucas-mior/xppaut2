#include "functions.h"
#include "xmalloc.h"

#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <math.h>
#include "struct.h"
#include "integers.h"
#include <libgen.h>

static NullCline nclines[MAXNCLINE];
int32 ps_color = 1;

#define SPER 3
#define UPER 4
#define UEQ 2

#define LS_SEQ 0
#define LS_UEQ 1
#define LS_SPER 8
#define LS_UPER 9

static struct Mov3D {
    char angle[20];
    char yes[3];
    double start;
    double incr;
    int32 nclip;
} mov3d = {"theta", "N", 45, 45, 7};

static struct BD {
    double *x[MAXBIFCRV];
    double *y[MAXBIFCRV];
    int32 color[MAXBIFCRV];
    int32 npts[MAXBIFCRV];
    int32 nbifcrv;
    Window window;
} my_bd;

static double freeze_key_x;
static double freeze_key_y;
static double freeze_key_x;
static double freeze_key_y;
static int32 freeze_key_flag;
int32 auto_freeze_flag = 0;
static int32 current_curve = 0;

int32 colorline[] = {0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 0};
char *color_names[] = {"WHITE",       "RED",   "REDORANGE", "ORANGE", "YELLOWORANGE", "YELLOW",
                       "YELLOWGREEN", "GREEN", "BLUEGREEN", "BLUE",   "PURPLE",       "BLACK"};

static int32 graf_par_get_frz_index(Window window);
static void graf_par_read_bd(FILE *fp);
static void graf_par_add_bd_crv(double *x, double *y, int32 len, int32 type, int32 ncrv);
static void graf_par_draw_frozen_cline(int32 index, Window window);
static void graf_par_edit_frz_crv(int32 i);
static int32 graf_par_create_crv(int32 ind);
static void graf_par_delete_frz(void);
static void graf_par_delete_frz_crv(int32 i);
static void graf_par_edit_frz(void);
static void graf_par_draw_freeze_key(void);
static void graf_par_set_key(int32 x, int32 y);
static int32 graf_par_alter_curve(char *title, int32 in_it, int32 n);
static void graf_par_zoom_out(int32 i1, int32 j1, int32 i2, int32 j2);
static void graf_par_zoom_in(int32 i1, int32 j1, int32 i2, int32 j2);
static void graf_par_movie_rot(double start, double increment, int32 nclip, int32 angle);
static void graf_par_fit_window(void);
static void graf_par_corner_cube(double *xlo, double *xhi, double *ylo, double *yhi);
static void graf_par_check_val(double *x1, double *x2, double *xb, double *xd);
static void graf_par_check_flags(void);
static void graf_par_update_view(double xlo, double xhi, double ylo, double yhi);

void
graf_par_change_view_com(int32 com) {
    if (com == 2) {
        array_plot_make_my("Array!");
        array_plot_edit();
        return;
    }
    if (com == 3) {
        ani_new_vcr();
        return;
    }

    MyGraph->grtype = 5*com;
    if (MyGraph->grtype < 5) {
        int32 ind = current_curve;
        static char *n[] = {"*0X-axis", "*0Y-axis", "Xmin",   "Ymin",
                            "Xmax",     "Ymax",     "Xlabel", "Ylabel"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status;
        int32 i;
        int32 i1 = MyGraph->xv[ind];
        int32 i2 = MyGraph->yv[ind];
        char n1[15];
        char n2[15];
        graf_par_ind_to_sym(i1, n1);
        graf_par_ind_to_sym(i2, n2);
        snprintf(values[0], sizeof(values[0]), "%s", n1);
        snprintf(values[1], sizeof(values[1]), "%s", n2);
        snprintf(values[2], sizeof(values[2]), "%g", MyGraph->xmin);
        snprintf(values[3], sizeof(values[3]), "%g", MyGraph->ymin);
        snprintf(values[4], sizeof(values[4]), "%g", MyGraph->xmax);
        snprintf(values[5], sizeof(values[5]), "%g", MyGraph->ymax);
        snprintf(values[6], sizeof(values[6]), "%s", MyGraph->xlabel);
        snprintf(values[7], sizeof(values[7]), "%s", MyGraph->ylabel);
        MyGraph->ThreeDFlag = 0;
        status = pop_list_do_string_box(8, 4, 2, "2D View", n, values, 31);
        if (status != 0) {
            //  get variable names
            browser_find_variable(values[0], &i);
            if (i > -1) {
                MyGraph->xv[ind] = i;
            }
            browser_find_variable(values[1], &i);
            if (i > -1) {
                MyGraph->yv[ind] = i;
            }

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
            graf_par_check_windows();
        }
    } else {
        int32 ind = current_curve;
        static char *n[] = {"*0X-axis", "*0Y-axis", "*0Z-axis", "Xmin",  "Xmax", "Ymin",
                            "Ymax",     "Zmin",     "Zmax",     "XLo",   "XHi",  "YLo",
                            "YHi",      "Xlabel",   "Ylabel",   "Zlabel"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status, i, i1 = MyGraph->xv[ind], i2 = MyGraph->yv[ind], i3 = MyGraph->zv[ind];
        char n1[15];
        char n2[15];
        char n3[15];
        graf_par_ind_to_sym(i1, n1);
        graf_par_ind_to_sym(i2, n2);
        graf_par_ind_to_sym(i3, n3);
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
        status = pop_list_do_string_box(16, 6, 3, "3D View", n, values, 31);
        if (status != 0) {
            //  get variable names
            browser_find_variable(values[0], &i);
            if (i > -1) {
                MyGraph->xv[ind] = i;
            }
            browser_find_variable(values[1], &i);
            if (i > -1) {
                MyGraph->yv[ind] = i;
            }
            browser_find_variable(values[2], &i);
            if (i > -1) {
                MyGraph->zv[ind] = i;
            }
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
            graf_par_check_windows();
        }
    }
    graf_par_check_flags();
    graf_par_redraw_the_graph();
}

void
graf_par_ind_to_sym(int32 ind, char *str) {
    if (ind == 0) {
        strcpy(str, "T");
    } else {
        strcpy(str, uvar_names[ind - 1]);
    }
    return;
}

void
graf_par_check_flags(void) {
    if (MyGraph->grtype > 4) {
        MyGraph->ThreeDFlag = 1;
    } else {
        MyGraph->ThreeDFlag = 0;
    }
    if ((MyGraph->xv[0] == 0) || (MyGraph->yv[0] == 0) ||
        ((MyGraph->zv[0] == 0) && (MyGraph->ThreeDFlag == 1))) {
        MyGraph->TimeFlag = 1;
    } else {
        MyGraph->TimeFlag = 0;
    }
    return;
}

void
graf_par_check_val(double *x1, double *x2, double *xb, double *xd) {
    double temp;

    // see get_max for details

    if (*x1 == *x2) {
        temp = .05*LMAX(fabs(*x1), 1.0);
        *x1 = *x1 - temp;
        *x2 = *x2 + temp;
    }
    if (*x1 > *x2) {
        temp = *x2;
        *x2 = *x1;
        *x1 = temp;
    }
    *xb = .5*(*x1 + *x2);
    *xd = 2.0 / (*x2 - *x1);
    return;
}

void
graf_par_get_max(int32 index, double *vmin, double *vmax) {
    double x0;
    double x1;
    double z;
    double temp;
    x0 = browser_my.data[index][0];
    x1 = x0;
    for (int32 i = 0; i < browser_my.maxrow; i++) {
        z = browser_my.data[index][i];
        if (z < x0) {
            x0 = z;
        }
        if (z > x1) {
            x1 = z;
        }
    }
    *vmin = (double)x0;
    *vmax = (double)x1;
    if (fabs(*vmin - *vmax) < REAL_SMALL) {
        temp = .05*LMAX(fabs(*vmin), 1.0);
        *vmin = *vmin - temp;
        *vmax = *vmax + temp;
    }
    return;
}

void
graf_par_corner_cube(double *xlo, double *xhi, double *ylo, double *yhi) {
    double x;
    double y;
    double x1;
    double x2;
    double y1;
    double y2;
    graphics_threedproj(-1., -1., -1., &x, &y);
    x1 = x;
    x2 = x;
    y1 = y;
    y2 = y;
    graphics_threedproj(-1., -1., 1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(-1., 1., -1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(-1., 1., 1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(1., -1., -1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(1., -1., 1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(1., 1., 1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    graphics_threedproj(1., 1., -1., &x, &y);
    if (x < x1) {
        x1 = x;
    }
    if (x > x2) {
        x2 = x;
    }
    if (y < y1) {
        y1 = y;
    }
    if (y > y2) {
        y2 = y;
    }
    *xlo = x1;
    *ylo = y1;
    *xhi = x2;
    *yhi = y2;
    return;
}

void
graf_par_default_window(void) {
    if (MyGraph->ThreeDFlag) {
        MyGraph->xmax = x_3d[1];
        MyGraph->ymax = y_3d[1];
        MyGraph->zmax = z_3d[1];
        MyGraph->xmin = x_3d[0];
        MyGraph->ymin = y_3d[0];
        MyGraph->zmin = z_3d[0];

        graf_par_corner_cube(&(MyGraph->xlo), &(MyGraph->xhi), &(MyGraph->ylo), &(MyGraph->yhi));
        graf_par_check_windows();
    } else {
        MyGraph->xmax = x_3d[1];
        MyGraph->ymax = y_3d[1];

        MyGraph->xmin = x_3d[0];
        MyGraph->ymin = y_3d[0];

        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        graf_par_check_windows();
    }

    graf_par_redraw_the_graph();
}

void
graf_par_fit_window(void) {
    double Mx = -1.e25;
    double My = -1.e25;
    double Mz = -1.e25;
    double mx = -Mx;
    double my = -My;
    double mz = -Mz;
    int32 n = MyGraph->nvars;

    if (storind < 2) {
        return;
    }
    if (MyGraph->ThreeDFlag) {
        for (int32 i = 0; i < n; i++) {
            graf_par_get_max(MyGraph->xv[i], &(MyGraph->xmin), &(MyGraph->xmax));
            Mx = LMAX(MyGraph->xmax, Mx);
            mx = -LMAX(-MyGraph->xmin, -mx);

            graf_par_get_max(MyGraph->yv[i], &(MyGraph->ymin), &(MyGraph->ymax));
            My = LMAX(MyGraph->ymax, My);
            my = -LMAX(-MyGraph->ymin, -my);

            graf_par_get_max(MyGraph->zv[i], &(MyGraph->zmin), &(MyGraph->zmax));
            Mz = LMAX(MyGraph->zmax, Mz);
            mz = -LMAX(-MyGraph->zmin, -mz);
        }
        MyGraph->xmax = Mx;
        MyGraph->ymax = My;
        MyGraph->zmax = Mz;
        MyGraph->xmin = mx;
        MyGraph->ymin = my;
        MyGraph->zmin = mz;

        graf_par_corner_cube(&(MyGraph->xlo), &(MyGraph->xhi), &(MyGraph->ylo), &(MyGraph->yhi));
        graf_par_check_windows();
    } else {
        for (int32 i = 0; i < n; i++) {
            graf_par_get_max(MyGraph->xv[i], &(MyGraph->xmin), &(MyGraph->xmax));
            Mx = LMAX(MyGraph->xmax, Mx);
            mx = -LMAX(-MyGraph->xmin, -mx);

            graf_par_get_max(MyGraph->yv[i], &(MyGraph->ymin), &(MyGraph->ymax));
            My = LMAX(MyGraph->ymax, My);
            my = -LMAX(-MyGraph->ymin, -my);
        }
        MyGraph->xmax = Mx;
        MyGraph->ymax = My;

        MyGraph->xmin = mx;
        MyGraph->ymin = my;

        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        graf_par_check_windows();
    }
    graf_par_redraw_the_graph();
}

void
graf_par_check_windows(void) {
    double zip;
    double zap;
    graf_par_check_val(&MyGraph->xmin, &MyGraph->xmax, &MyGraph->xbar, &MyGraph->dx);
    graf_par_check_val(&MyGraph->ymin, &MyGraph->ymax, &MyGraph->ybar, &MyGraph->dy);
    graf_par_check_val(&MyGraph->zmin, &MyGraph->zmax, &MyGraph->zbar, &MyGraph->dz);
    graf_par_check_val(&MyGraph->xlo, &MyGraph->xhi, &zip, &zap);
    graf_par_check_val(&MyGraph->ylo, &MyGraph->yhi, &zip, &zap);
    return;
}

void
graf_par_xi_vs_t(void) {
    //  a short cut
    char name[20];
    char value[20];
    int32 i = MyGraph->yv[0];

    graf_par_ind_to_sym(i, value);
    snprintf(name, sizeof(name), "Plot vs t: ");
    ggets_new_string(name, value);
    browser_find_variable(value, &i);

    if (i > -1) {
        MyGraph->yv[0] = i;
        MyGraph->grtype = 0;
        MyGraph->xv[0] = 0;
        if (storind >= 2) {
            graf_par_get_max(MyGraph->xv[0], &(MyGraph->xmin), &(MyGraph->xmax));
            graf_par_get_max(MyGraph->yv[0], &(MyGraph->ymin), &(MyGraph->ymax));

        } else {
            MyGraph->xmin = T0;
            MyGraph->xmax = TEND;
        }
        MyGraph->xlo = MyGraph->xmin;
        MyGraph->ylo = MyGraph->ymin;
        MyGraph->xhi = MyGraph->xmax;
        MyGraph->yhi = MyGraph->ymax;
        graf_par_check_windows();
        graf_par_check_flags();
        graphics_set_normal_scale();
        graf_par_redraw_the_graph();
    }
    return;
}

void
graf_par_redraw_the_graph(void) {
    ggets_blank_screen(draw_win);
    graphics_set_normal_scale();
    axes_do();
    many_pops_hi_lite(draw_win);
    integrate_restore(0, browser_my.maxrow);
    many_pops_draw_label(draw_win);
    graf_par_draw_freeze(draw_win);
    nullcline_redraw_dfield();
    if (MyGraph->Nullrestore) {
        restore_nullclines();
    }
    return;
}

void
graf_par_movie_rot(double start, double increment, int32 nclip, int32 angle) {
    double thetaold = MyGraph->Theta;
    double phiold = MyGraph->Phi;
    kinescope_reset_film();
    for (int32 i = 0; i <= nclip; i++) {
        if (angle == 0) {
            graphics_make_rot(start + i*increment, phiold);
        } else {
            graphics_make_rot(thetaold, start + i*increment);
        }
        graf_par_redraw_the_graph();
        kinescope_film_clip();
    }
    MyGraph->Theta = thetaold;
    MyGraph->Phi = phiold;
    return;
}

void
graf_par_get_3d_com(void) {
    static char *n[] = {
        "Persp (1=On)", "ZPlane",           "ZView",       "Theta",     "Phi",
        "Movie(Y/N)",   "Vary (theta/phi)", "Start angle", "Increment", "Number increments"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;

    int32 nclip = 8;
    int32 angle = 0;
    double start;
    double increment = 45;
    if (MyGraph->grtype < 5) {
        return;
    }

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

    status = pop_list_do_string_box(10, 5, 2, "3D Parameters", n, values, 28);
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
            if (mov3d.angle[0] == 'p' || mov3d.angle[0] == 'P') {
                angle = 1;
            }
            graf_par_movie_rot(start, increment, nclip, angle);
        }

        graphics_make_rot(MyGraph->Theta, MyGraph->Phi);
        //  Redraw the picture
        graf_par_redraw_the_graph();
    }
    return;
}

static void
graf_par_update_view(double xlo, double xhi, double ylo, double yhi) {
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
    graf_par_check_windows();

    graf_par_redraw_the_graph();
}

void
graf_par_window_zoom_com(int32 c) {
    int32 i1;
    int32 i2;
    int32 j1;
    int32 j2;
    switch (c) {
    case 0: {
        // graf par user window
        static char *n[] = {"X Lo", "X Hi", "Y Lo", "Y Hi"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status;
        snprintf(values[0], sizeof(values[0]), "%g", MyGraph->xlo);
        snprintf(values[2], sizeof(values[2]), "%g", MyGraph->ylo);
        snprintf(values[1], sizeof(values[1]), "%g", MyGraph->xhi);
        snprintf(values[3], sizeof(values[3]), "%g", MyGraph->yhi);
        status = pop_list_do_string_box(4, 2, 2, "Window", n, values, 28);
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
            graf_par_check_windows();
        }
        graf_par_redraw_the_graph();
        break;
    }
    case 1:
        if (rubber(&i1, &j1, &i2, &j2, draw_win, RUBBOX) == 0) {
            break;
        }
        graf_par_zoom_in(i1, j1, i2, j2);

        break;
    case 2:
        if (rubber(&i1, &j1, &i2, &j2, draw_win, RUBBOX) == 0) {
            break;
        }
        graf_par_zoom_out(i1, j1, i2, j2);
        break;
    case 3:
        graf_par_fit_window();
        break;
    case 4:
        graf_par_default_window();
        break;
    case 5: {
        // graf par scroll window
        XEvent event;
        int32 i = 0;
        int32 j = 0;
        int32 state = 0;
        double x;
        double y;
        double x0 = 0;
        double y0 = 0;
        double xlo = MyGraph->xlo;
        double ylo = MyGraph->ylo;
        double xhi = MyGraph->xhi;
        double yhi = MyGraph->yhi;
        double dx = 0;
        double dy = 0;
        int32 alldone = 0;
        XSelectInput(display, draw_win,
                     KeyPressMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask |
                         ButtonMotionMask | ExposureMask);
        while (!alldone) {
            XNextEvent(display, &event);
            switch (event.type) {
            case KeyPress:
                alldone = 1;
                break;
            case Expose:
                many_pops_do_expose(event);
                break;
            case ButtonPress:
                if (state == 0) {
                    i = event.xkey.x;
                    j = event.xkey.y;
                    graphics_scale_to_real(i, j, &x0, &y0);
                    state = 1;
                }
                break;
            case MotionNotify:
                if (state == 1) {
                    i = event.xmotion.x;
                    j = event.xmotion.y;
                    graphics_scale_to_real(i, j, &x, &y);
                    dx = -(x - x0) / 2;
                    dy = -(y - y0) / 2;

                    graf_par_update_view(xlo + dx, xhi + dx, ylo + dy, yhi + dy);
                }
                break;
            case ButtonRelease:
                state = 0;
                xlo = xlo + dx;
                xhi = xhi + dx;
                ylo = ylo + dy;
                yhi = yhi + dy;
                break;
            default:
                break;
            }
        }
        break;
    }
    default:
        break;
    }
    graphics_set_normal_scale();
    return;
}

void
graf_par_zoom_in(int32 i1, int32 j1, int32 i2, int32 j2) {
    double x1;
    double y1;
    double x2;
    double y2;
    double dx = MyGraph->xhi - MyGraph->xlo;
    double dy = MyGraph->yhi - MyGraph->ylo;
    graphics_scale_to_real(i1, j1, &x1, &y1);
    graphics_scale_to_real(i2, j2, &x2, &y2);
    if (x1 == x2 || y1 == y2) {
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        dx = dx / 2;
        dy = dy / 2;

        // Shrink by thirds and center (track) about the point clicked
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
    graf_par_check_windows();
    graf_par_redraw_the_graph();
    menu_draw_help();
    return;
}

void
graf_par_zoom_out(int32 i1, int32 j1, int32 i2, int32 j2) {
    double x1;
    double y1;
    double x2;
    double y2;
    double bx;
    double mux;
    double by;
    double muy;
    double dx = MyGraph->xhi - MyGraph->xlo;
    double dy = MyGraph->yhi - MyGraph->ylo;
    graphics_scale_to_real(i1, j1, &x1, &y1);
    graphics_scale_to_real(i2, j2, &x2, &y2);

    /*
    if(x1==x2||y1==y2)return;
    */
    /*
    ggets_plintf("%f %f %f %f \n ",x1,y1,x2,y2);
    ggets_plintf("%f %f %f %f
    \n",MyGraph->xlo,MyGraph->ylo,MyGraph->xhi,MyGraph->yhi);
   */
    if (x1 == x2 || y1 == y2) {
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        // Grow by thirds and center (track) about the point clicked
        dx = dx*2;
        dy = dy*2;

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
    graf_par_check_windows();
    graf_par_redraw_the_graph();
    menu_draw_help();
    return;
}

void
graf_par_graph_all(int32 *list, int32 n, int32 type) {
    if (type == 0) {
        for (int32 i = 0; i < n; i++) {
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
    graf_par_check_flags();
    graf_par_fit_window();
    return;
}

int32
graf_par_alter_curve(char *title, int32 in_it, int32 n) {
    static char *nn[] = {"*0X-axis", "*0Y-axis", "*0Z-axis", "*4Color", "Line type"};
    char values[LENGTH(nn)][MAX_LEN_SBOX];
    int32 status;
    int32 i;
    int32 i1 = MyGraph->xv[in_it], i2 = MyGraph->yv[in_it], i3 = MyGraph->zv[in_it];
    char n1[15];
    char n2[15];
    char n3[15];

    graf_par_ind_to_sym(i1, n1);
    graf_par_ind_to_sym(i2, n2);
    graf_par_ind_to_sym(i3, n3);
    snprintf(values[0], sizeof(values[0]), "%s", n1);
    snprintf(values[1], sizeof(values[1]), "%s", n2);
    snprintf(values[2], sizeof(values[2]), "%s", n3);
    snprintf(values[3], sizeof(values[3]), "%d", MyGraph->color[in_it]);
    snprintf(values[4], sizeof(values[4]), "%d", MyGraph->line[in_it]);
    status = pop_list_do_string_box(5, 5, 1, title, nn, values, 25);
    if (status != 0) {
        browser_find_variable(values[0], &i);
        if (i > -1) {
            MyGraph->xv[n] = i;
        }
        browser_find_variable(values[1], &i);
        if (i > -1) {
            MyGraph->yv[n] = i;
        }
        browser_find_variable(values[2], &i);
        if (i > -1) {
            MyGraph->zv[n] = i;
        }

        MyGraph->line[n] = atoi(values[4]);
        i = atoi(values[3]);
        if (i < 0 || i > 10) {
            i = 0;
        }
        MyGraph->color[n] = i;

        return 1;
    }
    return 0;
}

void
graf_par_dump_ps(int32 i) {
#define FILENAME_SIZE (sizeof(this_file) + sizeof(this_internset) + sizeof(plot_format) + 10)
    char filename[FILENAME_SIZE];
#undef FILENAME_SIZE
    if (i < 0) {
        snprintf(filename, sizeof(filename), "%s%s.%s", this_file, this_internset, plot_format);
    } else {
        snprintf(filename, sizeof(filename), "%s%s_%04d.%s", this_file, this_internset, i,
                 plot_format);
    }

    if (strcmp(plot_format, "ps") == 0) {
        if (ps_init(filename, ps_color)) {
            many_pops_ps_restore();
        }
    } else if (strcmp(plot_format, "svg") == 0) {
        if (svg_init(filename)) {
            many_pops_svg_restore();
        }
    }
    return;
}

void
graf_par_change_cmap_com(int32 i) {
    color_new_map(i);
    return;
}

void
graf_par_freeze_com(int32 c) {
    switch (c) {
    case 0: {
        // freeze crv
        int32 crv;
        crv = graf_par_create_crv(0);
        if (crv < 0) {
            break;
        }
        graf_par_edit_frz_crv(crv);
        break;
    }
    case 1:
        graf_par_delete_frz();
        break;
    case 2:
        graf_par_edit_frz();
        break;
    case 3:
        // kill frz
        for (int32 i = 0; i < MAXFRZ; i++) {
            if (frz[i].use == 1 && frz[i].window == draw_win) {
                graf_par_delete_frz_crv(i);
            }
        }
        break;
    case 5: {
        // frz bd
        FILE *fp;
        char filename[XPP_MAX_NAME];
        snprintf(filename, sizeof(filename), "diagram.dat");
        ggets_ping();
        if (!init_conds_file_selector("Import Diagram", filename, "*.dat")) {
            return;
        }
        if ((fp = fopen(filename, "r")) == NULL) {
            ggets_err_msg("Couldn't open file");
            return;
        }
        graf_par_read_bd(fp);
        break;
    }
    case 6:
        // free bd
        if (my_bd.nbifcrv > 0) {
            for (int32 i = 0; i < my_bd.nbifcrv; i++) {
                free(my_bd.x[i]);
                free(my_bd.y[i]);
            }
            my_bd.nbifcrv = 0;
        }
        break;
    case 7:
        auto_freeze_flag = 1 - auto_freeze_flag;
        break;
    default:
        break;
    }
    return;
}

void
graf_par_set_key(int32 x, int32 y) {
    double xp;
    double yp;
    graphics_scale_to_real(x, y, &xp, &yp);
    freeze_key_x = xp;
    freeze_key_y = yp;
    freeze_key_flag = 1;
    return;
}

void
graf_par_draw_freeze_key(void) {
    int32 ix;
    int32 iy;
    int32 y0;
    int32 ix2;
    int32 dy = 2*h_char;
    if (freeze_key_flag == SCRNFMT) {
        return;
    }
    if (plt_fmt_flag == PSFMT) {
        dy = -dy;
    }
    graphics_scale_to_screen((double)freeze_key_x, (double)freeze_key_y, &ix, &iy);
    ix2 = ix + 4*h_char;
    y0 = iy;
    for (int32 i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && frz[i].window == draw_win && strlen(frz[i].key) > 0) {
            graphics_set_linestyle(ABS(frz[i].color));
            graphics_line(ix, y0, ix2, y0);
            graphics_set_linestyle(0);
            graphics_put_text(ix2 + h_char, y0, frz[i].key);
            y0 += dy;
        }
    }
    return;
}

void
graf_par_key_frz_com(int32 c) {
    int32 x;
    int32 y;
    switch (c) {
    case 0:
        freeze_key_flag = 0;
        break;
    case 1:
        menudrive_message_box("Position with mouse");
        if (ggets_mouse_xy(&x, &y, draw_win)) {
            graf_par_set_key(x, y);
            graf_par_draw_freeze_key();
        }
        menudrive_message_box_kill();
        break;
    default:
        break;
    }
    return;
}

void
graf_par_edit_frz(void) {
    int32 i;
    i = graf_par_get_frz_index(draw_win);
    if (i < 0) {
        return;
    }
    graf_par_edit_frz_crv(i);
    return;
}

void
graf_par_delete_frz_crv(int32 i) {
    if (frz[i].use == 0) {
        return;
    }
    frz[i].use = 0;
    frz[i].name[0] = 0;
    frz[i].key[0] = 0;
    free(frz[i].xv);
    free(frz[i].yv);
    if (frz[i].type > 0) {
        free(frz[i].zv);
    }
    return;
}

void
graf_par_delete_frz(void) {
    int32 i;
    i = graf_par_get_frz_index(draw_win);
    if (i < 0) {
        return;
    }
    graf_par_delete_frz_crv(i);
    return;
}

void
graf_par_auto_freeze_it(void) {
    if (auto_freeze_flag == 0) {
        return;
    }
    graf_par_create_crv(0);
    return;
}

int32
graf_par_create_crv(int32 ind) {
    int32 type;
    int32 ix;
    int32 iy;
    int32 iz;

    for (int32 i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 0) {
            ix = MyGraph->xv[ind];
            iy = MyGraph->yv[ind];
            iz = MyGraph->zv[ind];
            if (browser_my.maxrow <= 2) {
                ggets_err_msg("No Curve to freeze");
                return -1;
            }
            frz[i].xv = xmalloc(sizeof(*(frz[i].xv))*(usize)browser_my.maxrow);
            frz[i].yv = xmalloc(sizeof(*(frz[i].yv))*(usize)browser_my.maxrow);
            if ((type = MyGraph->grtype) > 0) {
                frz[i].zv = xmalloc(sizeof(*(frz[i].zv))*(usize)browser_my.maxrow);
            }
            if ((type > 0 && frz[i].zv == NULL) || (type == 0 && frz[i].yv == NULL)) {
                ggets_err_msg("Cant allocate storage for curve");
                return -1;
            }
            frz[i].use = 1;
            frz[i].len = browser_my.maxrow;
            for (int32 j = 0; j < browser_my.maxrow; j++) {
                frz[i].xv[j] = browser_my.data[ix][j];
                frz[i].yv[j] = browser_my.data[iy][j];
                if (type > 0) {
                    frz[i].zv[j] = browser_my.data[iz][j];
                }
            }
            frz[i].type = (int16)type;
            frz[i].window = draw_win;
            snprintf(frz[i].name, sizeof(frz[i].name), "crv%c", 'a' + i);
            strncpy(frz[i].key, frz[i].name, sizeof(frz[i].key));
            return i;
        }
    }
    ggets_err_msg("All curves used");
    return -1;
}

void
graf_par_edit_frz_crv(int32 i) {
    static char *nn[] = {"*4Color", "Key", "Name"};
    char values[LENGTH(nn)][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", frz[i].color);
    snprintf(values[1], sizeof(values[1]), "%s", frz[i].key);
    snprintf(values[2], sizeof(values[2]), "%s", frz[i].name);
    status = pop_list_do_string_box(3, 3, 1, "Edit Freeze", nn, values, 25);
    if (status != 0) {
        frz[i].color = atoi(values[0]);
        strncpy(frz[i].key, values[1], sizeof(frz[i].key));
        strncpy(frz[i].name, values[2], sizeof(frz[i].name));
    }
    return;
}

void
graf_par_draw_frozen_cline(int32 index, Window window) {
    if (nclines[index].use == 0 || nclines[index].window != window) {
        return;
    }
    return;
}

void
graf_par_draw_freeze(Window window) {
    int32 type = MyGraph->grtype;
    int32 lt = 0;
    double oldxpl;
    double oldypl;
    double oldzpl = 0.0;
    double xpl;
    double ypl;
    double zpl = 0.0;
    double *xv, *yv, *zv;
    for (int32 i = 0; i < MAXNCLINE; i++) {
        graf_par_draw_frozen_cline(i, window);
    }
    for (int32 i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && frz[i].window == window && frz[i].type == type) {
            if (frz[i].color < 0) {
                graphics_set_linestyle(-frz[i].color);
                lt = 1;
            } else {
                graphics_set_linestyle(frz[i].color);
            }
            xv = frz[i].xv;
            yv = frz[i].yv;
            zv = frz[i].zv;
            oldxpl = xv[0];
            oldypl = yv[0];
            if (type > 0) {
                oldzpl = zv[0];
            }
            for (int32 j = 0; j < frz[i].len; j++) {
                xpl = xv[j];
                ypl = yv[j];
                if (type > 0) {
                    zpl = zv[j];
                }
                if (lt == 0) {
                    if (type == 0) {
                        graphics_line_abs(oldxpl, oldypl, xpl, ypl);
                    } else {
                        graphics_line_3d(oldxpl, oldypl, oldzpl, xpl, ypl, zpl);
                    }
                } else {
                    if (type == 0) {
                        graphics_point_abs(xpl, ypl);
                    } else {
                        graphics_point_3d(xpl, ypl, zpl);
                    }
                }
                oldxpl = xpl;
                oldypl = ypl;
                if (type > 0) {
                    oldzpl = zpl;
                }
            }
        }
    }
    graf_par_draw_freeze_key();
    {
        // draw bd
        int32 len;
        double *x, *y;
        if (window == my_bd.window && my_bd.nbifcrv > 0) {
            for (int32 i2 = 0; i2 < my_bd.nbifcrv; i2++) {
                graphics_set_linestyle(my_bd.color[i2]);
                len = my_bd.npts[i2];
                x = my_bd.x[i2];
                y = my_bd.y[i2];
                xpl = x[0];
                ypl = y[0];
                for (int32 j2 = 0; j2 < len; j2++) {
                    oldxpl = xpl;
                    oldypl = ypl;
                    xpl = x[j2];
                    ypl = y[j2];
                    graphics_line_abs(oldxpl, oldypl, xpl, ypl);
                }
            }
        }
    }
    return;
}

/*  Bifurcation curve importing */

void
graf_par_init_bd(void) {
    my_bd.nbifcrv = 0;
    return;
}

void
graf_par_add_bd_crv(double *x, double *y, int32 len, int32 type, int32 ncrv) {
    int32 i;
    if (ncrv >= MAXBIFCRV) {
        return;
    }
    my_bd.x[ncrv] = xmalloc(sizeof(*(my_bd.x[ncrv]))*(usize)len);
    my_bd.y[ncrv] = xmalloc(sizeof(*(my_bd.y[ncrv]))*(usize)len);
    for (i = 0; i < len; i++) {
        my_bd.x[ncrv][i] = x[i];
        my_bd.y[ncrv][i] = y[i];
    }
    my_bd.npts[ncrv] = len;
    i = LS_SEQ;
    if (type == UPER) {
        i = LS_UPER;
    }
    if (type == SPER) {
        i = LS_SPER;
    }
    if (type == UEQ) {
        i = LS_UEQ;
    }
    my_bd.color[ncrv] = i;
    return;
}

void
graf_par_read_bd(FILE *fp) {
    int32 oldtype;
    int32 type;
    int32 oldbr;
    int32 br;
    int32 ncrv = 0;
    int32 len;
    int32 f2;
    double x[8000];
    double ylo[8000];
    double yhi[8000];
    len = 0;
    fscanf(fp, "%lf %lf %lf %d %d %d", &x[len], &ylo[len], &yhi[len], &oldtype, &oldbr, &f2);
    len++;
    while (!feof(fp)) {
        fscanf(fp, "%lf %lf %lf %d %d %d", &x[len], &ylo[len], &yhi[len], &type, &br, &f2);
        if (type == oldtype && br == oldbr) {
            len++;
        } else {
            /* if(oldbr==br)len++; */  // extend to point of instability
            graf_par_add_bd_crv(x, ylo, len, oldtype, ncrv);
            ncrv++;
            if (oldtype == UPER || oldtype == SPER) {
                graf_par_add_bd_crv(x, yhi, len, oldtype, ncrv);
                ncrv++;
            }
            if (oldbr == br) {
                len--;
            }
            x[0] = x[len];
            ylo[0] = ylo[len];
            yhi[0] = yhi[len];

            len = 1;
        }
        oldbr = br;
        oldtype = type;
    }
    //  save this last one
    if (len > 1) {
        graf_par_add_bd_crv(x, ylo, len, oldtype, ncrv);
        ncrv++;
        if (oldtype == UPER || oldtype == SPER) {
            graf_par_add_bd_crv(x, yhi, len, oldtype, ncrv);
            ncrv++;
        }
    }
    ggets_plintf(" got %d bifurcation curves\n", ncrv);
    fclose(fp);
    my_bd.nbifcrv = ncrv;
    my_bd.window = draw_win;
    return;
}

int32
graf_par_get_frz_index(Window window) {
    char *n[MAXFRZ];
    char key[MAXFRZ];
    char ch;

    int32 count = 0;
    Window temp = main_win;
    for (int32 i = 0; i < MAXFRZ; i++) {
        if (frz[i].use == 1 && window == frz[i].window) {
            n[count] = xmalloc(20);
            sprintf(n[count], "%s", frz[i].name);
            key[count] = 'a' + (char)i;

            count++;
        }
    }
    if (count == 0) {
        return -1;
    }
    key[count] = 0;
    ch = (char)pop_list_popup_list_new(&temp, "Curves", n, key, count, 12, 0, 10, 8*dcur_y + 8,
                                       no_hint, info_pop, info_message);
    for (int32 i = 0; i < count; i++) {
        free(n[i]);
    }
    return (int32)(ch - 'a');
}

void
graf_par_add_a_curve_com(int32 c) {
    switch (c) {
    case 0:
        if (MyGraph->nvars >= MAXPERPLOT) {
            ggets_err_msg("Too many plots!");
            return;
        }
        // new curve
        if (graf_par_alter_curve("New Curve", 0, MyGraph->nvars)) {
            MyGraph->nvars = MyGraph->nvars + 1;
        }
        break;
    case 1:
        if (MyGraph->nvars > 1) {
            MyGraph->nvars = MyGraph->nvars - 1;
        }
        break;
    case 2:
        MyGraph->nvars = 1;
        break;
    case 3: {
        // edit curve
        char bob[21];
        int32 crv = 0;
        snprintf(bob, sizeof(bob), "Edit 0-%d :", MyGraph->nvars - 1);
        ggets_ping();
        ggets_new_int(bob, &crv);
        if (crv >= 0 && crv < MyGraph->nvars) {
            snprintf(bob, sizeof(bob), "Edit curve %d", crv);
            graf_par_alter_curve(bob, crv, crv);
        }
        break;
    }
    case 4: {
        // create ps
        char filename[XPP_MAX_NAME + 3];
        static char *nn[] = {"BW-0/Color-1", "Land(0)/Port(1)", "Axes fontsize", "Font",
                             "Linewidth"};
        int32 status;
        char values[LENGTH(nn)][MAX_LEN_SBOX];
        snprintf(values[0], sizeof(values[0]), "%d", ps_color);
        snprintf(values[1], sizeof(values[1]), "%d", ps_port);
        snprintf(values[2], sizeof(values[2]), "%d", ps_fontsize);
        strncpy(values[3], ps_font, sizeof(values[3]));
        snprintf(values[4], sizeof(values[4]), "%g", ps_lw);
        status = pop_list_do_string_box(5, 5, 1, "Postscript parameters", nn, values, 25);
        if (status != 0) {
            ps_color = atoi(values[0]);
            ps_port = atoi(values[1]);
            ps_fontsize = atoi(values[2]);
            ps_lw = atof(values[4]);
            snprintf(ps_font, sizeof(ps_font), "%s", values[3]);
            snprintf(filename, sizeof(filename), "%s.ps", this_file);
            ggets_ping();

            if (!init_conds_file_selector("Print postscript", filename, "*.ps")) {
                return;
            }
            if (ps_init(filename, ps_color)) {
                many_pops_ps_restore();
                ggets_ping();
            }
        }
        break;
    }
    case 5: {
        // create svg
        char filename[XPP_MAX_NAME];
        strcpy(filename, this_file);
        filename[strlen(filename) - 4] = '\0';
        strcat(filename, ".svg");
        if (!init_conds_file_selector("Print svg", filename, "*.svg")) {
            return;
        }
        if (svg_init(filename)) {
            many_pops_svg_restore();
            ggets_ping();
        }
        break;
    }
        /* case 6: freeze();
           break; */
    case 7: {
        // graf par axes opts
        static char *n[] = {"X-origin",    "Y-origin",   "Z-origin",  "X-org(1=on)",
                            "Y-org(1=on)", "Z-org(1=on", "PSFontSize"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status;
        snprintf(values[0], sizeof(values[0]), "%g", MyGraph->xorg);
        snprintf(values[1], sizeof(values[1]), "%g", MyGraph->yorg);
        snprintf(values[2], sizeof(values[2]), "%g", MyGraph->zorg);
        snprintf(values[3], sizeof(values[3]), "%d", MyGraph->xorgflag);
        snprintf(values[4], sizeof(values[4]), "%d", MyGraph->yorgflag);
        snprintf(values[5], sizeof(values[5]), "%d", MyGraph->zorgflag);
        snprintf(values[6], sizeof(values[6]), "%d", ps_fontsize);
        status = pop_list_do_string_box(7, 7, 1, "Axes options", n, values, 25);
        if (status != 0) {
            MyGraph->xorg = atof(values[0]);
            MyGraph->yorg = atof(values[1]);
            MyGraph->zorg = atof(values[2]);
            MyGraph->xorgflag = atoi(values[3]);
            MyGraph->yorgflag = atoi(values[4]);
            MyGraph->zorgflag = atoi(values[5]);
            ps_fontsize = atoi(values[6]);
            graf_par_redraw_the_graph();
        }
        break;
    }
    case 8: {
        // export graf data
        FILE *fp;
        char filename[XPP_MAX_NAME];
        snprintf(filename, sizeof(filename), "curve.dat");
        ggets_ping();
        if (!init_conds_file_selector("Export graph data", filename, "*.dat")) {
            return;
        }
        // if(ggets_new_string("Data filename:",filename)==0)return;
        if ((fp = fopen(filename, "w")) == NULL) {
            ggets_err_msg("Couldn't open file");
            return;
        }
        integrate_export_data(fp);
        fclose(fp);
        break;
    }
        /*  case 9: change_cmap();
            break; */
    default:
        break;
    }
    graf_par_check_flags();
    graf_par_redraw_the_graph();
}
