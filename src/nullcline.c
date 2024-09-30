#include "functions.h"
#include "integers.h"
#include "xmalloc.h"
#include <stdbool.h>

#include "parserslow.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "xpplim.h"
#include "struct.h"
#include <stdio.h>

#define MAX_NULL 10000

typedef struct Point {
    double x;
    double y;
    double z;
} Point;

static int32 NCSuppress = 0;
static int32 df_supress = 0;
int32 df_batch = 0;
int32 NCBatch = 0;

static int32 NullStyle = 0;  // 1 is with little vertical/horizontal lines

int32 XNullColor = 2;
int32 YNullColor = 7;
static int32 num_x_n;
static int32 num_y_n;
static int32 num_index;
static int32 null_ix;
static int32 null_iy;
static int32 WHICH_CRV;
static double *X_n;
static double *Y_n;
static double *saver;
static double *NTop;
static double *NBot;
int32 df_grid = 16;
int32 df_flag = 0;
static int32 df_ix = -1;
static int32 df_iy = -1;
static int32 dfield_type = 0;

int32 doing_dfield = 0;

char color_via[15] = "speed";
double color_via_lo = 0;
double color_via_hi = 1;
int32 colorize_flag = 0;

static RangeInfo ncrange;

static NullClines *ncperm;
static int32 n_nstore = 0;
static int32 ncline_cnt;

static int32 nullcline_interpolate(Point p1, Point p2, double z, double *x, double *y);
static void nullcline_do_cline(int32 ngrid, double x1, double y1, double x2, double y2);
static void nullcline_quad_contour(Point p1, Point p2, Point p3, Point p4);
static double nullcline_fnull(double x, double y);
static void nullcline_store(double x1, double y1, double x2, double y2);
static void nullcline_restor(double *v, int32 n, int32 d);
static void nullcline_dump(FILE *fp, double *x, int32 nx, double *y, int32 ny);
static void nullcline_redraw_froz(int32 flag);
static void nullcline_save_frozen(char *fn);

void
nullcline_froz_cline_stuff_com(int32 i) {
    int32 delay2 = 200;

    if (n_nstore == 0) {
        // nullcline start
        n_nstore = 1;
        ncperm = xmalloc(sizeof(*ncperm));
        ncperm->p = NULL;
        ncperm->n = NULL;
        ncperm->nmx = 0;
        ncperm->nmy = 0;
        ncperm->n_ix = -5;
        ncperm->n_iy = -5;
        ncrange.xlo = 0;
        ncrange.xhi = 1;
        ncrange.nstep = 10;
        snprintf(ncrange.rv, sizeof(ncrange.rv), " ");
    }
    switch (i) {
    case 0:
        if (NULL_HERE == 0) {
            return;
        }
        nullcline_add_froz(X_n, num_x_n, null_ix, Y_n, num_y_n, null_iy);
        break;
    case 1: {
        // nullcline_clear_froz
        NullClines *z;
        NullClines *znew;
        z = ncperm;
        while (z->n != NULL) {
            z = z->n;
        }
        //  this is the bottom but there is nothing here that has been stored

        znew = z->p;
        if (znew == NULL) {
            break;
        }
        free(z);
        z = znew;
        // now we are deleting everything
        while (z->p != NULL) {
            znew = z->p;
            z->n = NULL;
            z->p = NULL;
            free(z->xn);
            free(z->yn);
            free(z);
            z = znew;
        }
        if (ncperm->nmx > 0) {
            free(ncperm->xn);
            ncperm->nmx = 0;
        }
        if (ncperm->nmy > 0) {
            free(ncperm->yn);
            ncperm->nmy = 0;
        }
        ncperm->n = NULL;
        n_nstore = 1;
        ncline_cnt = 0;
        break;
    }
    case 3:
        ggets_new_int("delay2 (msec)", &delay2);
        if (delay2 <= 0) {
            delay2 = 0;
        }
        nullcline_redraw_froz(delay2);
        break;
    case 2: {
        // nullcline do range
        static char *n[] = {"*2Range parameter", "Steps", "Low", "High"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status;
        int32 i2;
        double z;
        double dz;
        double zold;
        double xmin;
        double xmax;
        double y_tp;
        double y_bot;
        int32 col1 = XNullColor;
        int32 col2 = YNullColor;
        int32 course = NMESH;
        /* if(paper_white){
          col1=1;
          col2=9;
          } */
        snprintf(values[0], sizeof(values[0]), "%s", ncrange.rv);
        snprintf(values[1], sizeof(values[1]), "%d", ncrange.nstep);
        snprintf(values[2], sizeof(values[2]), "%g", ncrange.xlo);
        snprintf(values[3], sizeof(values[3]), "%g", ncrange.xhi);
        status = pop_list_do_string_box(4, 4, 1, "Range Clines", n, values, 45);
        if (status == 0) {
            break;
        }

        strcpy(ncrange.rv, values[0]);
        ncrange.nstep = atoi(values[1]);
        ncrange.xlo = atof(values[2]);
        ncrange.xhi = atof(values[3]);
        if (ncrange.nstep <= 0) {
            break;
        }
        dz = (ncrange.xhi - ncrange.xlo) / (double)ncrange.nstep;
        if (dz <= 0.0) {
            break;
        }
        get_val(ncrange.rv, &zold);

        for (i2 = NODE; i2 < NODE + NMarkov; i2++) {
            set_ivar(i2 + 1 + fix_var, last_ic[i2]);
        }
        xmin = (double)MyGraph->xmin;
        xmax = (double)MyGraph->xmax;
        y_tp = (double)MyGraph->ymax;
        y_bot = (double)MyGraph->ymin;
        null_ix = MyGraph->xv[0];
        null_iy = MyGraph->yv[0];

        for (i2 = 0; i2 <= ncrange.nstep; i2++) {
            z = (double)i2*dz + ncrange.xlo;
            set_val(ncrange.rv, z);
            if (NULL_HERE == 0) {
                if ((X_n = xmalloc(4*MAX_NULL*sizeof(double))) != NULL &&
                    (Y_n = xmalloc(4*MAX_NULL*sizeof(double))) != NULL) {

                    NULL_HERE = 1;
                }
                NTop = xmalloc((usize)(course + 1)*sizeof(*NTop));
                NBot = xmalloc((usize)(course + 1)*sizeof(*NBot));
                if (NTop == NULL || NBot == NULL) {
                    NULL_HERE = 0;
                }
            } else {
                free(NTop);
                free(NBot);
                NTop = xmalloc((usize)(course + 1)*sizeof(*NTop));
                NBot = xmalloc((usize)(course + 1)*sizeof(*NBot));
                if (NTop == NULL || NBot == NULL) {
                    NULL_HERE = 0;
                    return;
                }
            }

            WHICH_CRV = null_ix;
            graphics_set_linestyle(col1);
            new_nullcline(course, xmin, y_bot, xmax, y_tp, X_n, &num_x_n);

            WHICH_CRV = null_iy;
            graphics_set_linestyle(col2);
            new_nullcline(course, xmin, y_bot, xmax, y_tp, Y_n, &num_y_n);
            nullcline_add_froz(X_n, num_x_n, null_ix, Y_n, num_y_n, null_iy);
        }
        set_val(ncrange.rv, zold);
        break;
    }
    default:
        break;
    }
    return;
}

void
nullcline_silent_dfields(void) {
    if (df_batch == 5 || df_batch == 4) {
        df_supress = 1;
        graphics_init_ps();
        nullcline_do_batch_dfield();
        df_supress = 0;
    }
    return;
}

void
silent_nullclines(void) {
    FILE *fp;
    if (NCBatch != 2) {
        return;
    }
    NCSuppress = 1;
    nullcline_new_clines_com(0);
    fp = fopen("nullclines.dat", "w");
    if (fp == NULL) {
        ggets_plintf("Cannot open nullcline file\n");
        return;
    }
    nullcline_dump(fp, X_n, num_x_n, Y_n, num_y_n);
    fclose(fp);
    NCSuppress = 0;
    return;
}

int32
get_nullcline_floats(double **v, int32 *n, int32 who, int32 type) {
    // type=0,1
    NullClines *z;
    if (who < 0) {
        if (type == 0) {
            *v = X_n;
            *n = num_x_n;
        } else {
            *v = Y_n;
            *n = num_y_n;
        }
        if (v == NULL) {
            return 1;
        }
        return 0;
    }
    if (who > ncline_cnt || n_nstore == 0) {
        return 1;
    }
    z = ncperm;
    for (int32 i = 0; i < who; i++) {
        z = z->n;
    }
    if (z == NULL) {
        return 1;
    }
    if (type == 0) {
        *v = z->xn;
        *n = z->nmx;
    } else {
        *v = z->yn;
        *n = z->nmy;
    }
    if (v == NULL) {
        return 1;
    }
    return 0;
}

void
nullcline_save_frozen(char *fn) {
    NullClines *z;
    FILE *fp;
    char fnx[256];
    char ch;
    int32 i = 1;
    if (n_nstore == 0) {
        return;
    }
    ch = (char)menudrive_two_choice("YES", "NO", "Save Frozen Clines?", "yn");
    if (ch == 'n') {
        return;
    }
    z = ncperm;
    while (true) {
        if (z == NULL || (z->nmx == 0 && z->nmy == 0)) {
            return;
        }
        snprintf(fnx, sizeof(fnx), "%s.%d", fn, i);
        fp = fopen(fnx, "w");
        if (fp == NULL) {
            ggets_err_msg("Cant open file!");
            return;
        }
        nullcline_dump(fp, z->xn, z->nmx, z->yn, z->nmy);
        fclose(fp);
        i++;
        z = z->n;
        if (z == NULL) {
            break;
        }
    }
    return;
}

void
nullcline_redraw_froz(int32 flag) {
    NullClines *z;
    int32 col1 = XNullColor;
    int32 col2 = YNullColor;
    /* if(paper_white){
      col1=1;
      col2=9;
      } */
    if (n_nstore == 0) {
        return;
    }
    z = ncperm;
    while (true) {
        if (z == NULL || (z->nmx == 0 && z->nmy == 0)) {
            return;
        }

        if (MyGraph->xv[0] == z->n_ix && MyGraph->yv[0] == z->n_iy && MyGraph->ThreeDFlag == 0) {
            if (flag > 0) {
                browser_wait_a_sec(flag);
                main_clr_scrn();
            }
            graphics_set_linestyle(col1);
            nullcline_restor(z->xn, z->nmx, 1);
            graphics_set_linestyle(col2);
            nullcline_restor(z->yn, z->nmy, 2);
            if (flag > 0) {
                menudrive_flush_display();
            }
        }
        z = z->n;
        if (z == NULL) {
            break;
        }
    }
    return;
}

void
nullcline_add_froz(double *xn, int32 nmx, int32 n_ix, double *yn, int32 nmy, int32 n_iy) {
    NullClines *z;
    NullClines *znew;
    z = ncperm;
    // move to end
    while (z->n != NULL) {
        z = (z->n);
    }
    z->xn = xmalloc(4*(usize)nmx*sizeof(*(z->xn)));
    for (int32 i = 0; i < 4*nmx; i++) {
        z->xn[i] = xn[i];
    }
    z->yn = xmalloc(4*(usize)nmy*sizeof(*(z->yn)));
    for (int32 i = 0; i < 4*nmy; i++) {
        z->yn[i] = yn[i];
    }
    z->nmx = nmx;
    z->nmy = nmy;
    z->n_ix = n_ix;
    z->n_iy = n_iy;
    z->n = xmalloc(sizeof(*(z->n)));
    znew = z->n;
    znew->n = NULL;
    znew->p = z;
    znew->nmx = 0;
    znew->nmy = 0;
    znew->n_ix = -5;
    znew->n_iy = -5;
    ncline_cnt++;
    return;
}

void
nullcline_get_max_dfield(double *y, double *ydot, double u0, double v0, double du, double dv,
                         int32 n, int32 inx, int32 iny, double *mdf) {
    double amp;
    double dxp;
    double dyp;
    *mdf = 0.0;
    for (int32 i = 0; i <= n; i++) {
        y[inx] = u0 + du*i;
        for (int32 j = 0; j <= n; j++) {
            y[iny] = v0 + dv*j;
            rhs_function(0.0, y, ydot, NODE);
            main_rhs_extra(y, 0.0, NODE, n_equations);
            graphics_scale_dxdy(ydot[inx], ydot[iny], &dxp, &dyp);
            amp = hypot(dxp, dyp);
            if (amp > *mdf) {
                *mdf = amp;
            }
        }
    }
    return;
}
/*  all the nifty 2D stuff here    */

void
nullcline_do_batch_nclines(void) {
    if (!xpp_batch) {
        return;
    }
    if (!NCBatch) {
        return;
    }
    if (NCBatch == 1) {
        nullcline_new_clines_com(0);
        return;
    }
    return;
}

void
nullcline_set_colorization_stuff(void) {
    numerics_user_set_color_par(colorize_flag, color_via, color_via_lo, color_via_hi);
    return;
}

void
nullcline_do_batch_dfield(void) {
    if (!xpp_batch) {
        return;
    }
    switch (df_batch) {
    case 0:
        return;
    case 1:
        df_flag = 1;
        dfield_type = 1;
        df_ix = MyGraph->xv[0];
        df_iy = MyGraph->yv[0];
        nullcline_redraw_dfield();
        return;
    case 2:
        df_flag = 1;
        dfield_type = 0;
        df_ix = MyGraph->xv[0];
        df_iy = MyGraph->yv[0];
        nullcline_redraw_dfield();
        return;
    case 3:
        df_flag = 2;
        dfield_type = 0;
        df_ix = MyGraph->xv[0];
        df_iy = MyGraph->yv[0];
        nullcline_redraw_dfield();
        return;
    case 4:
        df_flag = 1;
        dfield_type = 1;
        df_ix = MyGraph->xv[0];
        df_iy = MyGraph->yv[0];
        nullcline_redraw_dfield();
        return;
    case 5:
        df_flag = 1;
        dfield_type = 0;
        df_ix = MyGraph->xv[0];
        df_iy = MyGraph->yv[0];
        nullcline_redraw_dfield();
        return;
    default:
        break;
    }
    return;
}

void
nullcline_redraw_dfield(void) {
    int32 inx = MyGraph->xv[0] - 1;
    int32 iny = MyGraph->yv[0] - 1;
    double y[MAX_ODE];
    double ydot[MAX_ODE];
    double xv1;
    double xv2;
    double v1[MAX_ODE];
    double v2[MAX_ODE];
    FILE *fp = NULL;

    double amp;
    double mdf;

    double du;
    double dv;
    double u0;
    double v0;
    double dxp;
    double dyp;
    double dz;
    double dup;
    double dvp;

    int32 grid = df_grid;
    if (df_flag == 0 || MyGraph->TimeFlag || MyGraph->xv[0] == MyGraph->yv[0] ||
        MyGraph->ThreeDFlag || df_ix != MyGraph->xv[0] || df_iy != MyGraph->yv[0]) {
        return;
    }
    if (df_supress == 1) {
        fp = fopen("dirfields.dat", "w");
        if (fp == NULL) {
            return;
        }
    }

    du = (MyGraph->xhi - MyGraph->xlo) / (double)grid;
    dv = (MyGraph->yhi - MyGraph->ylo) / (double)grid;

    dup = (double)(d_right - d_left) / (double)grid;
    dvp = (double)(d_top - d_buttom) / (double)grid;
    dz = hypot(dup, dvp)*(.25 + .75*dfield_type);
    u0 = MyGraph->xlo;
    v0 = MyGraph->ylo;
    if (!df_supress) {
        graphics_set_linestyle(MyGraph->color[0]);
    }
    integrate_get_ic(2, y);
    nullcline_get_max_dfield(y, ydot, u0, v0, du, dv, grid, inx, iny, &mdf);
    if (plt_fmt_flag == SVGFMT) {
        doing_dfield = 1;
        fprintf(svgfile, "<g>\n");
    }
    for (int32 i = 0; i <= grid; i++) {
        y[inx] = u0 + du*i;
        for (int32 j = 0; j <= grid; j++) {
            y[iny] = v0 + dv*j;
            rhs_function(0.0, y, ydot, NODE);
            main_rhs_extra(y, 0.0, NODE, n_equations);
            if (MyGraph->ColorFlag || df_flag == 2) {
                v1[0] = 0.0;
                v2[0] = 0.0;
                for (int32 k = 0; k < n_equations; k++) {
                    v1[k + 1] = (double)y[k];
                    v2[k + 1] = v1[k + 1] + (double)ydot[k];
                }
                if (!df_supress) {
                    integrate_comp_color(v1, v2, NODE, 1.0);
                }
            }
            if (df_flag == 1 || df_flag == 4) {
                graphics_scale_dxdy(ydot[inx], ydot[iny], &dxp, &dyp);
                if (dfield_type == 1) {
                    ydot[inx] /= mdf;
                    ydot[iny] /= mdf;
                } else {
                    amp = hypot(dxp, dyp);
                    if (amp != 0.0) {
                        ydot[inx] /= amp;
                        ydot[iny] /= amp;
                    }
                }
                xv1 = y[inx] + ydot[inx]*dz;
                xv2 = y[iny] + ydot[iny]*dz;
                if (!df_supress) {
                    graphics_bead_abs((double)xv1, (double)xv2);
                    graphics_line_abs((double)y[inx], (double)y[iny], (double)xv1, (double)xv2);
                } else {
                    if (fp) {
                        fprintf(fp, "%g %g %g %g\n", y[inx], y[iny], xv1, xv2);
                    }
                }
            }
            if (df_flag == 2 && j > 0 && i < grid) {
                graphics_frect_abs((double)y[inx], (double)y[iny], (double)du, (double)dv);
            }
        }
    }

    if (plt_fmt_flag == SVGFMT) {
        doing_dfield = 0;
        fprintf(svgfile, "</g>\n");
    }
    if (df_supress == 1) {
        fclose(fp);
    }
    df_supress = 0;
    return;
}

void
nullcline_direct_field_com(int32 c) {
    int32 start;
    int32 inx = MyGraph->xv[0] - 1;
    int32 iny = MyGraph->yv[0] - 1;
    double y[MAX_ODE];
    double ydot[MAX_ODE];
    double xv1;
    double xv2;
    double dtold = delta_t;
    double v1[MAX_ODE];
    double v2[MAX_ODE];

    double amp;
    double mdf;
    double t;
    double du;
    double dv;
    double u0;
    double v0;
    double dxp;
    double dyp;
    double dz;
    double dup;
    double dvp;
    double oldtrans = TRANS;

    int32 grid = df_grid;

    if (MyGraph->TimeFlag || MyGraph->xv[0] == MyGraph->yv[0] || MyGraph->ThreeDFlag) {
        return;
    }

    if (c == 2) {
        df_flag = 0;
        return;
    }
    if (c == 0) {
        dfield_type = 1;
    }
    if (c == 4) {
        dfield_type = 0;
    }
    ggets_new_int("Grid:", &grid);
    if (grid <= 1) {
        return;
    }
    df_grid = grid;
    du = (MyGraph->xhi - MyGraph->xlo) / (double)grid;
    dv = (MyGraph->yhi - MyGraph->ylo) / (double)grid;

    dup = (double)(d_right - d_left) / (double)grid;
    dvp = (double)(d_top - d_buttom) / (double)grid;
    dz = hypot(dup, dvp)*(.25 + .75*dfield_type);
    u0 = MyGraph->xlo;
    v0 = MyGraph->ylo;
    graphics_set_linestyle(MyGraph->color[0]);
    if (c != 1) {
        df_flag = 1;
        if (c == 3) {
            df_flag = 2;
            du = (MyGraph->xhi - MyGraph->xlo) / (double)(grid + 1);
            dv = (MyGraph->yhi - MyGraph->ylo) / (double)(grid + 1);
        }
        df_ix = inx + 1;
        df_iy = iny + 1;
        integrate_get_ic(2, y);
        nullcline_get_max_dfield(y, ydot, u0, v0, du, dv, grid, inx, iny, &mdf);
        if (plt_fmt_flag == SVGFMT) {
            doing_dfield = 1;
            fprintf(svgfile, "<g>\n");
        }

        for (int32 i = 0; i <= grid; i++) {
            y[inx] = u0 + du*i;
            for (int32 j = 0; j <= grid; j++) {
                y[iny] = v0 + dv*j;
                rhs_function(0.0, y, ydot, NODE);
                main_rhs_extra(y, 0.0, NODE, n_equations);
                if (MyGraph->ColorFlag || df_flag == 2) {
                    v1[0] = 0.0;
                    v2[0] = 0.0;
                    for (int32 k = 0; k < n_equations; k++) {
                        v1[k + 1] = (double)y[k];
                        v2[k + 1] = v1[k + 1] + (double)ydot[k];
                    }
                    integrate_comp_color(v1, v2, NODE, 1.0);
                }
                if (df_flag == 1) {
                    graphics_scale_dxdy(ydot[inx], ydot[iny], &dxp, &dyp);
                    if (dfield_type == 0) {
                        amp = hypot(dxp, dyp);
                        if (amp != 0.0) {
                            ydot[inx] /= amp;
                            ydot[iny] /= amp;
                        }
                    } else {
                        ydot[inx] /= mdf;
                        ydot[iny] /= mdf;
                    }
                    xv1 = y[inx] + ydot[inx]*dz;
                    xv2 = y[iny] + ydot[iny]*dz;
                    graphics_bead_abs((double)xv1, (double)xv2);
                    graphics_line_abs((double)y[inx], (double)y[iny], (double)xv1, (double)xv2);
                }
                if (df_flag == 2 && j > 0 && i < grid) {
                    graphics_frect_abs((double)y[inx], (double)y[iny], (double)du, (double)dv);
                }
            }
        }
        TRANS = oldtrans;
        if (plt_fmt_flag == SVGFMT) {
            doing_dfield = 0;
            fprintf(svgfile, "</g>\n");
        }
        return;
    }
    stor_flag = 0;

    suppress_bounds = 1;
    for (int32 k = 0; k < 2; k++) {
        for (int32 i = 0; i <= grid; i++) {
            for (int32 j = 0; j <= grid; j++) {
                integrate_get_ic(2, y);
                y[inx] = u0 + du*i;
                y[iny] = v0 + dv*j;
                t = 0.0;
                start = 1;
                /*if(integrate(&t,y,TEND,delta_t,1,njmp,&start)==1){
                  TRANS=oldtrans;
                  delta_t=dtold;
                  return;
                  stor_flag=1;
                  } */
                integrate(&t, y, TEND, delta_t, 1, njmp, &start);
            }
        }
        delta_t = -delta_t;
    }
    suppress_bounds = 0;
    delta_t = dtold;
    if (plt_fmt_flag == SVGFMT) {
        doing_dfield = 0;
        fprintf(svgfile, "</g>\n");
    }
    return;
}

/* animated nullclines stuff
 * added Aug 31 97
 * just redraws them
 * It will allow you to either freeze a range of them
 * or just one at a time

 * freeze   -   store the current set
 * range_freeze - compute over some range of parameters
 * clear - delete all but the current set
 * animate - replay all frozen ones (not current set )
 */

void
restore_nullclines(void) {
    int32 col1 = XNullColor;
    int32 col2 = YNullColor;
    /* if(paper_white){
      col1=1;
      col2=9;
      } */
    if (NULL_HERE == 0) {
        return;
    }
    if (MyGraph->xv[0] == null_ix && MyGraph->yv[0] == null_iy && MyGraph->ThreeDFlag == 0) {
        graphics_set_linestyle(col1);
        nullcline_restor(X_n, num_x_n, 1);
        graphics_set_linestyle(col2);
        nullcline_restor(Y_n, num_y_n, 2);
    }
    nullcline_redraw_froz(0);
    return;
}

void
nullcline_dump(  // gnuplot format
    FILE *fp, double *x, int32 nx, double *y, int32 ny) {
    fprintf(fp, "# X-nullcline\n");
    for (int32 i = 0; i < nx - 1; i++) {
        fprintf(fp, "%g %g 1 \n", x[4*i], x[4*i + 1]);
        fprintf(fp, "%g %g 1 \n", x[4*i + 2], x[4*i + 3]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n# Y-nullcline\n");
    for (int32 i = 0; i < ny - 1; i++) {
        fprintf(fp, "%g %g 2 \n", y[4*i], y[4*i + 1]);
        fprintf(fp, "%g %g 2 \n", y[4*i + 2], y[4*i + 3]);
        fprintf(fp, "\n");
    }
    return;
}

void
nullcline_restor(  // d=1 for x and 2 for y
    double *v, int32 n, int32 d) {
    int32 i4;
    double xm;
    double ym;
    int32 x1;
    int32 y1;
    if (plt_fmt_flag == SVGFMT) {
        fprintf(svgfile, "<g>\n");
    }

    for (int32 i = 0; i < n; i++) {
        i4 = 4*i;
        graphics_line_abs(v[i4], v[i4 + 1], v[i4 + 2], v[i4 + 3]);
        if (NullStyle == 1) {
            xm = .5*(v[i4] + v[i4 + 2]);
            ym = .5*(v[i4 + 1] + v[i4 + 3]);
            graphics_scale_to_screen(xm, ym, &x1, &y1);
            switch (d) {
            case 1:
                graphics_line(x1, y1 - 4, x1, y1 + 4);

                break;
            case 2:
                graphics_line(x1 - 4, y1, x1 + 4, y1);
                break;
            default:
                break;
            }
        }
    }

    if (plt_fmt_flag == SVGFMT) {
        fprintf(svgfile, "</g>\n");
    }
    return;
}

void
nullcline_create_new_cline(void) {
    if (NULL_HERE) {
        nullcline_new_clines_com(0);
    }
    return;
}

void
nullcline_new_clines_com(int32 c) {
    int32 course = NMESH;
    double xmin;
    double xmax;
    double y_tp;
    double y_bot;
    int32 col1 = XNullColor;
    int32 col2 = YNullColor;

    if (MyGraph->ThreeDFlag || MyGraph->TimeFlag || MyGraph->xv[0] == MyGraph->yv[0]) {
        return;
    }

    if (c == 1) {
        restore_nullclines();
        return;
    }
    if (c == 2) {
        MyGraph->Nullrestore = 1;
        return;
    }
    if (c == 3) {
        MyGraph->Nullrestore = 0;
        return;
    }
    if (c == 4) {
        menudrive_froz_cline_stuff();
        return;
    }
    if (c == 5) {
        // save_the_nullclines
        FILE *fp;
        char filename[256];
        if (NULL_HERE == 0) {
            return;
        }
        snprintf(filename, sizeof(filename), "nc.dat");
        ggets_ping();
        if (!init_conds_file_selector("Save nullclines", filename, "*.dat")) {
            return;
        }
        fp = fopen(filename, "w");
        if (fp == NULL) {
            ggets_err_msg("Cant open file!");
            return;
        }
        nullcline_dump(fp, X_n, num_x_n, Y_n, num_y_n);
        fclose(fp);
        nullcline_save_frozen(filename);
        return;
    }
    if (c == 0) {
        for (int32 i = NODE; i < NODE + NMarkov; i++) {
            set_ivar(i + 1 + fix_var, last_ic[i]);
        }
        xmin = (double)MyGraph->xmin;
        xmax = (double)MyGraph->xmax;
        y_tp = (double)MyGraph->ymax;
        y_bot = (double)MyGraph->ymin;
        null_ix = MyGraph->xv[0];
        null_iy = MyGraph->yv[0];
        if (NULL_HERE == 0) {
            if ((X_n = xmalloc(4*MAX_NULL*sizeof(double))) != NULL &&
                (Y_n = xmalloc(4*MAX_NULL*sizeof(double))) != NULL) {

                NULL_HERE = 1;
            }
            NTop = xmalloc((usize)(course + 1)*sizeof(*NTop));
            NBot = xmalloc((usize)(course + 1)*sizeof(*NBot));
            if (NTop == NULL || NBot == NULL) {
                NULL_HERE = 0;
            }
        } else {
            free(NTop);
            free(NBot);
            NTop = xmalloc((usize)(course + 1)*sizeof(*NTop));
            NBot = xmalloc((usize)(course + 1)*sizeof(*NBot));
            if (NTop == NULL || NBot == NULL) {
                NULL_HERE = 0;
                return;
            }
        }

        WHICH_CRV = null_ix;
        if (!NCSuppress) {
            graphics_set_linestyle(col1);
        }
        new_nullcline(course, xmin, y_bot, xmax, y_tp, X_n, &num_x_n);
        ggets_ping();

        WHICH_CRV = null_iy;
        if (!NCSuppress) {
            graphics_set_linestyle(col2);
        }
        new_nullcline(course, xmin, y_bot, xmax, y_tp, Y_n, &num_y_n);
        ggets_ping();
    }
    return;
}

void
new_nullcline(int32 course, double xlo, double ylo, double xhi, double yhi, double *stor,
              int32 *npts) {
    num_index = 0;
    saver = stor;
    nullcline_do_cline(course, xlo, ylo, xhi, yhi);
    *npts = num_index;
    return;
}

void
nullcline_store(double x1, double y1, double x2, double y2) {
    int32 i;
    if (num_index >= MAX_NULL) {
        return;
    }
    i = 4*num_index;
    saver[i] = x1;
    saver[i + 1] = y1;
    saver[i + 2] = x2;
    saver[i + 3] = y2;
    num_index++;
    return;
}

double
nullcline_fnull(double x, double y) {
    double y1[MAX_ODE];
    double ydot[MAX_ODE];
    for (int32 i = 0; i < NODE; i++) {
        y1[i] = last_ic[i];
    }

    y1[null_ix - 1] = (double)x;
    y1[null_iy - 1] = (double)y;
    rhs_function(0.0, y1, ydot, NODE);
    return (double)ydot[WHICH_CRV - 1];
}

int32
nullcline_interpolate(Point p1, Point p2, double z, double *x, double *y) {
    double scale;
    if (p1.z == p2.z) {
        return 0;
    }
    scale = (z - p1.z) / (p2.z - p1.z);
    *x = p1.x + scale*(p2.x - p1.x);
    *y = p1.y + scale*(p2.y - p1.y);
    return 1;
}

void
nullcline_quad_contour(Point p1, Point p2, Point p3, Point p4) {
    double x[4];
    double y[4];
    int32 count = 0;
    if (p1.z*p2.z <= 0.0) {
        if (nullcline_interpolate(p1, p2, 0.0, &x[count], &y[count])) {
            count++;
        }
    }
    if (p2.z*p3.z <= 0.0) {
        if (nullcline_interpolate(p3, p2, 0.0, &x[count], &y[count])) {
            count++;
        }
    }
    if (p3.z*p4.z <= 0.0) {
        if (nullcline_interpolate(p3, p4, 0.0, &x[count], &y[count])) {
            count++;
        }
    }
    if (p1.z*p4.z <= 0.0) {
        if (nullcline_interpolate(p1, p4, 0.0, &x[count], &y[count])) {
            count++;
        }
    }

    if (count == 2) {
        if (!NCSuppress) {
            graphics_line_abs(x[0], y[0], x[1], y[1]);
        }
        nullcline_store(x[0], y[0], x[1], y[1]);
    }
    return;
}

void
nullcline_do_cline(int32 ngrid, double x1, double y1, double x2, double y2) {
    double dx = (x2 - x1) / (double)ngrid;
    double dy = (y2 - y1) / (double)ngrid;
    double x;
    double y;
    Point p[5];
    int32 nx = ngrid + 1;
    int32 ny = ngrid + 1;

    y = y2;
    for (int32 i = 0; i < nx; i++) {
        x = x1 + i*dx;
        NBot[i] = nullcline_fnull(x, y);
    }

    for (int32 j = 1; j < ny; j++) {
        y = y2 - j*dy;
        NTop[0] = NBot[0];
        NBot[0] = nullcline_fnull(x1, y);
        for (int32 i = 1; i < nx; i++) {
            x = x1 + i*dx;
            NTop[i] = NBot[i];
            NBot[i] = nullcline_fnull(x, y);
            p[0].x = x - dx;
            p[0].y = y + dy;
            p[0].z = NTop[i - 1];
            p[1].x = x;
            p[1].y = y + dy;
            p[1].z = NTop[i];
            p[3].x = x - dx;
            p[3].y = y;
            p[3].z = NBot[i - 1];
            p[2].x = x;
            p[2].y = y;
            p[2].z = NBot[i];
            /*      Uncomment for triangle contour
                 p[4].x=.25*(p[0].x+p[1].x+p[2].x+p[3].x);
                p[4].y=.25*(p[0].y+p[1].y+p[2].y+p[3].y);
                p[4].z=.25*(p[0].z+p[1].z+p[2].z+p[3].z);

                triangle_contour(p[0],p[1],p[4]);
                triangle_contour(p[1],p[4],p[2]);
                triangle_contour(p[4],p[3],p[2]);
                triangle_contour(p[0],p[4],p[3]); */
            //   Uncomment for quad contour
            nullcline_quad_contour(p[0], p[1], p[2], p[3]);
        }
    }
}
