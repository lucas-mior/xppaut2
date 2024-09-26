#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "integers.h"
#include "struct.h"
#include "xpplim.h"

#define SIGNIF (0.01) /* less than one hundredth of a tic mark */
#define CheckZero(x, tic) (fabs(x) < ((tic)*SIGNIF) ? 0.0 : (x))

int32 axes2_doing = 0;
int32 axes2_doing_box = 0;

static void axes2_get_title_str(char *s1, char *s2, char *s3);
static void axes2_make_title(char *str);
static double axes2_dbl_raise(double x, int32 y);
static double axes2_make_tics(double tmin, double tmax);
static void axes2_find_max_min_tic(double *tmin, double *tmax, double tic);
static void axes2_draw_ytics(char *s1, double start, double incr, double end);
static void axes2_draw_xtics(char *s2, double start, double incr, double end);

void
axes2_get_title_str(char *s1, char *s2, char *s3) {
    int32 i;
    if ((i = MyGraph->xv[0]) == 0)
        strcpy(s1, "T");
    else
        strcpy(s1, uvar_names[i - 1]);

    if ((i = MyGraph->yv[0]) == 0)
        strcpy(s2, "T");
    else
        strcpy(s2, uvar_names[i - 1]);

    if ((i = MyGraph->zv[0]) == 0)
        strcpy(s3, "T");
    else
        strcpy(s3, uvar_names[i - 1]);
    return;
}

void
axes2_make_title(char *str) {
    int32 i;
    char name1[20];
    char name2[20];
    char name3[20];
    if ((i = MyGraph->xv[0]) == 0)
        strcpy(name1, "T");
    else
        strcpy(name1, uvar_names[i - 1]);

    if ((i = MyGraph->yv[0]) == 0)
        strcpy(name2, "T");
    else
        strcpy(name2, uvar_names[i - 1]);

    if ((i = MyGraph->zv[0]) == 0)
        strcpy(name3, "T");
    else
        strcpy(name3, uvar_names[i - 1]);

    if (MyGraph->grtype >= 5)
        sprintf(str, "%s vs %s vs %s", name3, name2, name1);
    else
        sprintf(str, "%s vs %s", name2, name1);
    return;
}

double
axes2_dbl_raise(double x, int32 y) {
    double val;

    val = 1.0;
    for (int32 i = 0; i < abs(y); i++)
        val *= x;
    if (y < 0)
        return 1.0 / val;
    return val;
}

double
axes2_make_tics(double tmin, double tmax) {
    register double xr, xnorm, tics, tic, l10;

    xr = fabs(tmin - tmax);

    l10 = log10(xr);
    xnorm =
        pow(10.0, l10 - (double)((l10 >= 0.0) ? (int32)l10 : ((int32)l10 - 1)));
    if (xnorm <= 2)
        tics = 0.2;
    else if (xnorm <= 5)
        tics = 0.5;
    else
        tics = 1.0;
    tic = tics *
          axes2_dbl_raise(10.0, (l10 >= 0.0) ? (int32)l10 : ((int32)l10 - 1));
    return tic;
}

void
axes2_find_max_min_tic(double *tmin, double *tmax, double tic) {
    double t1 = *tmin;
    t1 = tic*floor(*tmin / tic);
    if (t1 < *tmin)
        t1 += tic;
    *tmin = t1;
    t1 = tic*ceil(*tmax / tic);
    if (t1 > *tmax)
        t1 -= tic;
    *tmax = t1;
    return;
}

void
axes2_redraw_cube_pt(double theta, double phi) {
    char bob[50];
    graphics_set_linestyle(0);
    graphics_make_rot(theta, phi);
    main_clr_scrn();

    sprintf(bob, "theta=%g phi=%g", theta, phi);
    many_pops_canvas_xy(bob);
    return;
}

void
axes2_do(void) {
    char s1[20], s2[20], s3[20];
    axes2_get_title_str(s1, s2, s3);
    graphics_set_linestyle(0);
    if (Xup) {
        /* axes2 re title */
        char bob[40];
        axes2_make_title(bob);
        many_pops_title_text(bob);
        many_pops_small_gr();
    }

    switch (MyGraph->grtype) {
    case 0:
        axes2_box(MyGraph->xlo, MyGraph->xhi, MyGraph->ylo, MyGraph->yhi,
                  MyGraph->xlabel, MyGraph->ylabel, 1);
        break;
    case 5: {
        double tx, ty, tz;
        double x1, y1, z1, x2, y2, z2, dt = .03;
        double x0 = MyGraph->xorg, y0 = MyGraph->yorg, z0 = MyGraph->zorg;
        char bob[20];

        double xmin = MyGraph->xmin, xmax = MyGraph->xmax, ymin = MyGraph->ymin;
        double ymax = MyGraph->ymax, zmin = MyGraph->zmin, zmax = MyGraph->zmax;
        double x4 = xmin, y4 = ymin, z4 = zmin, x5 = xmax, y5 = ymax, z5 = zmax;
        double x3, y3, z3, x6, y6, z6;

        axes2_doing = 1;

        tx = axes2_make_tics(xmin, xmax);
        ty = axes2_make_tics(ymin, ymax);
        tz = axes2_make_tics(zmin, zmax);
        axes2_find_max_min_tic(&xmin, &xmax, tx);
        axes2_find_max_min_tic(&zmin, &zmax, tz);
        axes2_find_max_min_tic(&ymin, &ymax, ty);
        graphics_scale3d((double)xmin, (double)ymin, (double)zmin, &x1, &y1,
                         &z1);
        graphics_scale3d((double)xmax, (double)ymax, (double)zmax, &x2, &y2,
                         &z2);

        graphics_scale3d(x4, y4, z4, &x3, &y3, &z3);
        graphics_scale3d(x5, y5, z5, &x6, &y6, &z6);
        graphics_set_linestyle(-2);
        graphics_line3d(-1., -1., -1., 1., -1., -1.);
        graphics_line3d(1., -1., -1., 1., 1., -1.);
        graphics_line3d(1., 1., -1., -1., 1., -1.);
        graphics_line3d(-1., 1., -1., -1., -1., -1.);
        graphics_line3d(-1., -1., 1., 1., -1., 1.);
        graphics_line3d(1., -1., 1., 1., 1., 1.);
        graphics_line3d(1., 1., 1., -1., 1., 1.);
        graphics_line3d(-1., 1., 1., -1., -1., 1.);
        graphics_line3d(1., 1., 1., 1., 1., -1.);
        graphics_line3d(-1., 1., 1., -1., 1., -1.);
        graphics_line3d(-1., -1., 1., -1., -1., -1.);
        graphics_line3d(1., -1., 1., 1., -1., -1.);

        graphics_line3dn(-1. - dt, -1., z2, -1. + dt, -1., z2);
        graphics_line3dn(-1. - dt, -1., z1, -1. + dt, -1., z1);
        graphics_line3dn(x2, -1. - dt, -1.0, x2, -1.0 + dt, -1.0);
        graphics_line3dn(x1, -1. - dt, -1.0, x1, -1.0 + dt, -1.0);
        graphics_line3dn(1.0 - dt, y1, -1.0, 1.0 + dt, y1, -1.0);
        graphics_line3dn(1.0 - dt, y2, -1.0, 1.0 + dt, y2, -1.0);

        graphics_set_linestyle(-1);

        if (MyGraph->zorgflag)
            graphics_line_3d(x0, y0, z4, x0, y0, z5);
        if (MyGraph->yorgflag)
            graphics_line_3d(x0, y4, z0, x0, y5, z0);
        if (MyGraph->xorgflag)
            graphics_line_3d(x4, y0, z0, x5, y0, z0);

        dt = .06;
        TextJustify = 2;
        sprintf(bob, "%g", xmin);
        graphics_text3d(x1, -1 - 2.*dt, -1.0, bob);
        sprintf(bob, "%g", xmax);
        graphics_text3d(x2, -1 - 2.*dt, -1.0, bob);
        graphics_text3d(0.0, -1 - dt, -1.0, MyGraph->xlabel);
        TextJustify = 0;
        sprintf(bob, "%g", ymin);
        graphics_text3d(1 + dt, y1, -1.0, bob);
        sprintf(bob, "%g", ymax);
        graphics_text3d(1 + dt, y2, -1.0, bob);
        graphics_text3d(1 + dt, 0.0, -1.0, MyGraph->ylabel);
        TextJustify = 2;
        sprintf(bob, "%g", zmin);
        graphics_text3d(-1. - dt, -1 - dt, z1, bob);
        sprintf(bob, "%g", zmax);
        graphics_text3d(-1. - dt, -1 - dt, z2, bob);
        graphics_text3d(-1. - dt, -1. - dt, 0.0, MyGraph->zlabel);
        TextJustify = 0;

        axes2_doing = 0;
        break;
    }
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    if (Xup)
        many_pops_small_base();
    return;
}

void
axes2_redraw_cube(double theta, double phi) {
    char bob[50];
    graphics_set_linestyle(0);
    graphics_make_rot(theta, phi);
    ggets_blank_screen(draw_win);

    /* axes2 draw unit cube */
    graphics_line3d(-1., -1., -1., 1., -1., -1.);
    graphics_line3d(1., -1., -1., 1., 1., -1.);
    graphics_line3d(1., 1., -1., -1., 1., -1.);
    graphics_line3d(-1., 1., -1., -1., -1., -1.);
    graphics_line3d(-1., -1., 1., 1., -1., 1.);
    graphics_line3d(1., -1., 1., 1., 1., 1.);
    graphics_line3d(1., 1., 1., -1., 1., 1.);
    graphics_line3d(-1., 1., 1., -1., -1., 1.);
    graphics_line3d(1., 1., 1., 1., 1., -1.);
    graphics_line3d(-1., 1., 1., -1., 1., -1.);
    graphics_line3d(-1., -1., 1., -1., -1., -1.);
    graphics_line3d(1., -1., 1., 1., -1., -1.);

    sprintf(bob, "theta=%g phi=%g", theta, phi);
    many_pops_canvas_xy(bob);
    return;
}

void
axes2_box(double x_min, double x_max, double y_min, double y_max, char *sx,
          char *sy, int32 flag) {
    double ytic;
    double xtic;

    int32 xaxis_y;
    int32 yaxis_x;

    int32 ybot = DBottom, ytop = DTop;
    int32 xleft = DLeft, xright = DRight;

    axes2_doing = 1;

    if (ybot > ytop) {
        ytop = ybot;
        ybot = DTop;
    }

    ytic = axes2_make_tics(y_min, y_max);
    xtic = axes2_make_tics(x_min, x_max);
    graphics_scale_to_screen((double)MyGraph->xorg, (double)MyGraph->yorg,
                             &yaxis_x, &xaxis_y);
    graphics_set_linestyle(-1);
    if (MyGraph->xorgflag && flag)
        if (xaxis_y >= ybot && xaxis_y <= ytop)
            graphics_line(xleft, xaxis_y, xright, xaxis_y);
    if (MyGraph->yorgflag && flag)
        if (yaxis_x >= xleft && yaxis_x <= xright)
            graphics_line(yaxis_x, ybot, yaxis_x, ytop);
    graphics_set_linestyle(-2);
    axes2_doing_box = 1;
    graphics_line(xleft, ybot, xright, ybot);
    graphics_line(xright, ybot, xright, ytop);
    axes2_doing_box = 0;
    graphics_line(xright, ytop, xleft, ytop);
    graphics_line(xleft, ytop, xleft, ybot);
    axes2_draw_ytics(sy, ytic*floor(y_min / ytic), ytic,
                     ytic*ceil(y_max / ytic));
    axes2_draw_xtics(sx, xtic*floor(x_min / xtic), xtic,
                     xtic*ceil(x_max / xtic));
    TextJustify = 0;
    graphics_set_linestyle(0);

    axes2_doing = 0;
    return;
}

void
axes2_draw_ytics(char *s1, double start, double incr, double end) {
    double ticvalue;
    double place;
    double y_min = YMin, y_max = YMax, x_min = XMin;
    char bob[100];
    int32 xt, yt, s = 1;
    TextJustify = 2; /* Right justification  */
    for (ticvalue = start; ticvalue <= end; ticvalue += incr) {
        place = CheckZero(ticvalue, incr);
        if (ticvalue < y_min || ticvalue > y_max)
            continue;
        sprintf(bob, "%g", place);
        graphics_scale_to_screen((double)x_min, (double)place, &xt, &yt);
        axes2_doing_box = 0;
        graphics_line(DLeft, yt, DLeft + HTic, yt);
        axes2_doing_box = 1;
        graphics_line(DRight, yt, DRight - HTic, yt);
        axes2_doing_box = 0;
        graphics_put_text(DLeft - (int32)(1.25*HChar), yt, bob);
    }
    graphics_scale_to_screen((double)x_min, (double)y_max, &xt, &yt);
    if (DTop < DBottom)
        s = -1;
    if (PltFmtFlag == SVGFMT) {
        fprintf(svgfile,
                "\n      <text class=\"xppyaxislabelv\" text-anchor=\"middle\" "
                "x=\"%d\"  y=\"%d\"\n",
                0, 0);
        fprintf(svgfile,
                "      transform=\"rotate(-90,75,180) translate(75,180)\"\n");
        fprintf(svgfile, "      >%s</text>\n", s1);

        fprintf(svgfile,
                "\n      <text class=\"xppyaxislabelh\" text-anchor=\"end\" "
                "x=\"%d\"  y=\"%d\"\n",
                DLeft - HChar, yt + 2*s*VChar);
        fprintf(svgfile, "      >%s</text>\n", s1);
    } else {
        graphics_put_text(DLeft - HChar, yt + 2*s*VChar, s1);
    }
    return;
}

void
axes2_draw_xtics(char *s2, double start, double incr, double end) {
    double ticvalue;
    double place;
    double y_min = YMin, x_min = XMin, x_max = XMax;

    char bob[100];
    int32 xt, yt = 0;
    int32 s = 1;
    if (DTop < DBottom)
        s = -1;
    TextJustify = 1; /* Center justification  */
    for (ticvalue = start; ticvalue <= end; ticvalue += incr) {
        place = CheckZero(ticvalue, incr);
        if (ticvalue < x_min || ticvalue > x_max)
            continue;
        sprintf(bob, "%g", place);
        graphics_scale_to_screen((double)place, y_min, &xt, &yt);
        axes2_doing_box = 0;
        graphics_line(xt, DBottom, xt, DBottom + s*VTic);
        axes2_doing_box = 1;
        graphics_line(xt, DTop, xt, DTop - s*VTic);
        axes2_doing_box = 0;
        graphics_put_text(xt, yt - (int32)(1.25*VChar*s), bob);
    }
    graphics_put_text((DLeft + DRight) / 2, yt - (int32)(2.5*VChar*s), s2);
    return;
}
