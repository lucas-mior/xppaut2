#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "functions.h"
#include "integers.h"

#define GREYSCALE -1

static struct DevScale {
    double xmin, xmax, ymin, ymax;
    double xscale, yscale, xoff, yoff;
    double tx, ty, angle, slant; /* text attributes   */
    double linecol, letx, lety;
    int32 linewid;
} ps_scale;

static FILE *my_plot_file;

static void ps_replot(double **z, int32 col0, int32 row0, int32 nskip,
                      int32 ncskip, int32 maxrow, int32 maxcol, int32 nacross,
                      int32 ndown, double zmin, double zmax, int32 type);
static void array_print_ps_begin(double xlo, double ylo, double xhi, double yhi,
                                 double sx, double sy);
static void array_print_ps_convert(double x, double y, double *xs, double *ys);
static void ps_col_scale(double y0, double x0, double dy, double dx, int32 n,
                         double zlo, double zhi, int32 type);
static void ps_boxit(double tlo, double thi, double jlo, double jhi, double zlo,
                     double zhi, char *sx, char *sy, char *sb, int32 type);
static void array_print_ps_close(void);
static void array_print_ps_setline(double fill, int32 thick);
static void ps_text2(char *str, double xr, double yr, int32 icent);
static void array_print_ps_set_text(double angle, double slant, double x_size,
                                    double y_size);
static void array_print_ps_rect(double x, double y, double wid, double len);
static void array_print_ps_bar(double x, double y, double wid, double len,
                               double fill, int32 flag);
static void ps_rgb_bar(double x, double y, double wid, double len, double fill,
                       int32 flag, int32 rgb);
static void ps_hsb_bar(double x, double y, double wid, double len, double fill,
                       int32 flag);

int32
array_print(char *filename, char *xtitle, char *ytitle, char *bottom,
            int32 nacross, int32 ndown, int32 col0, int32 row0, int32 nskip,
            int32 ncskip, int32 maxrow, int32 maxcol, double **data,
            double zmin, double zmax, double tlo, double thi, int32 type) {
    double xx = (double)ndown;
    double yy = (double)(nacross / ncskip);
    my_plot_file = fopen(filename, "w");
    if (my_plot_file == NULL)
        return -1;
    array_print_ps_begin(0.0, 0.0, xx, yy, 10., 7.);
    ps_replot(data, col0, row0, nskip, ncskip, maxrow, maxcol, nacross, ndown,
              zmin, zmax, type);
    ps_boxit(tlo, thi, 0.0, yy, zmin, zmax, xtitle, ytitle, bottom, type);
    array_print_ps_close();
    return 0;
}

void
ps_replot(double **z, int32 col0, int32 row0, int32 nskip, int32 ncskip,
          int32 maxrow, int32 maxcol, int32 nacross, int32 ndown, double zmin,
          double zmax, int32 type) {
    int32 i, j, ib, jb;

    double fill, x, y;
    double dx = (ps_scale.xmax - ps_scale.xmin);
    double dy = (ps_scale.ymax - ps_scale.ymin);
    double xhi = .95*dx, yhi = .85*dy;
    double delx;
    double dely;
    delx = .8*dx / (double)ndown;
    dely = .8*dy / (double)(nacross / ncskip);
    for (i = 0; i < nacross / ncskip; i++) {
        ib = col0 + i*ncskip;
        if (ib > maxcol)
            return;
        for (j = 0; j < ndown; j++) {
            jb = row0 + j*nskip;
            if (jb < maxrow && jb >= 0) {
                fill = (z[ib][jb] - zmin) / (zmax - zmin);
                if (fill < 0.0)
                    fill = 0.0;
                if (fill > 1.0)
                    fill = 1.0;
                fill = 1 - fill;
                x = xhi - delx - j*delx;
                y = yhi - dely - i*dely;
                if (type == GREYSCALE)
                    array_print_ps_bar(x, y, delx, dely, fill, 0);
                else
                    ps_rgb_bar(x, y, delx, dely, fill, 0, type);
            }
        }
    }
    return;
}

void
array_print_ps_begin(double xlo, double ylo, double xhi, double yhi, double sx,
                     double sy) {
    double x0, y0, x1, y1;
    ps_scale.xmin = xlo;
    ps_scale.ymin = ylo;
    ps_scale.ymax = yhi;
    ps_scale.xmax = xhi;
    ps_scale.xoff = 300;
    ps_scale.yoff = 300;
    ps_scale.angle = -90.;
    ps_scale.xscale = 1800.*sx*.2 / (xhi - xlo);
    ps_scale.yscale = 1800.*sy*.2 / (yhi - ylo);

    array_print_ps_set_text(-90., 0.0, 18.0, 18.0);
    ps_scale.letx = ps_scale.tx / ps_scale.xscale;
    ps_scale.lety = ps_scale.ty / ps_scale.yscale;
    array_print_ps_convert(xlo, ylo, &x0, &y0);
    array_print_ps_convert(xhi, yhi, &x1, &y1);
    fprintf(my_plot_file, "%s\n", "%!");
    fprintf(my_plot_file, "%s %g %g %g %g\n", "%%BoundingBox: ", .2*x0,
            .2*y0, .2*x1, .2*y1);
    fprintf(my_plot_file, "20 dict begin\n");
    fprintf(my_plot_file, "gsave\n");
    fprintf(my_plot_file, "/m {moveto} def\n");
    fprintf(my_plot_file, "/l {lineto} def\n");
    fprintf(my_plot_file, "/Cshow { currentpoint stroke moveto\n");
    fprintf(my_plot_file,
            "  dup stringwidth pop -2 div vshift rmoveto show } def\n");
    fprintf(my_plot_file, "/Lshow { currentpoint stroke moveto\n");
    fprintf(my_plot_file, "  0 vshift rmoveto show } def\n");
    fprintf(my_plot_file, "/Rshow { currentpoint stroke moveto\n");
    fprintf(my_plot_file,
            "  dup stringwidth pop neg vshift rmoveto show } def\n");
    fprintf(my_plot_file, "/C {setrgbcolor} def\n");
    fprintf(my_plot_file, "/G {setgray} def\n");
    fprintf(my_plot_file, "/S {stroke} def\n");
    fprintf(my_plot_file, "/HSB {sethsbcolor} def\n");
    fprintf(my_plot_file, "/RGB {setrgbcolor} def\n");
    fprintf(my_plot_file, "/FS {fill stroke} def\n");
    fprintf(my_plot_file, "630 -20 translate\n");
    fprintf(my_plot_file, "90 rotate\n");
    fprintf(my_plot_file, ".2 .2 scale\n");
    fprintf(my_plot_file, "/basefont /Times-Roman findfont def\n");
    array_print_ps_setline(0.0, 4);
    return;
}

void
array_print_ps_convert(double x, double y, double *xs, double *ys) {
    *xs = (x - ps_scale.xmin)*ps_scale.xscale + ps_scale.xoff;
    *ys = (y - ps_scale.ymin)*ps_scale.yscale + ps_scale.yoff;
    return;
}

void
ps_col_scale(double y0, double x0, double dy, double dx, int32 n, double zlo,
             double zhi, int32 type) {
    int32 i;
    char s[100];

    double dz = 1. / (double)(n - 1);

    for (i = 0; i < n; i++) {
        if (type == GREYSCALE)
            array_print_ps_bar(x0, y0 - (i + 1)*dy, dx, dy,
                               1 - (double)i*dz, 0);
        else
            ps_rgb_bar(x0, y0 - (i + 1)*dy, dx, dy, 1. - (double)i*dz, 0,
                       type);
    }
    fprintf(my_plot_file, "0 G\n");
    sprintf(s, "%g", zlo);
    ps_text2(s, x0 + .5*dx, y0 + .01*dx, 2);
    sprintf(s, "%g", zhi);
    ps_text2(s, x0 + .5*dx, y0 - n*dy - dy / 2, 0);
    return;
}

void
ps_boxit(double tlo, double thi, double jlo, double jhi, double zlo, double zhi,
         char *sx, char *sy, char *sb, int32 type) {
    char str[100];
    int32 i = ps_scale.linewid;
    double z = ps_scale.linecol;
    double dx = ps_scale.xmax - ps_scale.xmin;
    double dy = ps_scale.ymax - ps_scale.ymin;
    double xlo = .15*dx, ylo = .05*dy, xhi = .95*dx, yhi = .85*dy;

    array_print_ps_setline(0.0, 10);
    array_print_ps_rect(xlo, ylo, .8*dx, .8*dy);
    array_print_ps_setline(z, i);

    ps_text2(sx, xhi + .01*dx, .5*(yhi + ylo), 1);
    ps_text2(sy, .5*(xhi + xlo), yhi + .01*dy, 2);
    sprintf(str, "%g", tlo);
    ps_text2(str, xhi - .01*dx, yhi + .01*dy, 2);
    sprintf(str, "%g", thi);
    ps_text2(str, xlo, yhi + .01*dy, 2);
    sprintf(str, "%g", jlo);
    ps_text2(str, xhi + .01*dx, yhi, 0);
    sprintf(str, "%g", jhi);
    ps_text2(str, xhi + .01*dx, ylo + .01, 2);
    ps_col_scale(yhi - .15*dy, xlo - .1*dx, .025*dy, .05*dx, 20, zlo,
                 zhi, type);
    ps_text2(sb, xlo - .035*dx, .5*(yhi + ylo), 1);
    return;
}

void
array_print_ps_close(void) {
    fprintf(my_plot_file, "showpage\n");
    fprintf(my_plot_file, "grestore\n");
    fprintf(my_plot_file, "end\n");
    fclose(my_plot_file);
    return;
}

void
array_print_ps_setline(double fill, int32 thick) {
    fprintf(my_plot_file, "%f G\n %d setlinewidth \n", fill, thick);
    ps_scale.linewid = thick;
    ps_scale.linecol = fill;
    return;
}

void
ps_text2(char *str, double xr, double yr, int32 icent /* ignores for now  */
) {
    double slant = .0174532*ps_scale.slant;
    double x;
    double y;
    double sizex = ps_scale.tx, sizey = ps_scale.ty, rot = ps_scale.angle;
    double a = sizex*cos(slant), b = sizey*sin(slant),
           c = -sizex*sin(slant), d = sizey*cos(slant);
    array_print_ps_convert(xr, yr, &x, &y);
    fprintf(my_plot_file, "%d %d m\n", (int32)x, (int32)y);
    fprintf(my_plot_file, "gsave \n %f rotate \n", rot);
    fprintf(my_plot_file,
            "basefont [%.4f %.4f %.4f %.4f 0 0] makefont setfont\n", a, b, c,
            d);
    switch (icent) {
    case 0:
        fprintf(my_plot_file, "( %s ) show \n grestore\n", str);
        break;
    case 1: /* centered */
        fprintf(my_plot_file,
                "(%s) dup stringwidth pop -2 div 0 rmoveto show \n grestore\n",
                str);
        break;
    case 2: /* left edge */
        fprintf(my_plot_file,
                "(%s) dup stringwidth pop neg 0 rmoveto show \n grestore\n",
                str);
        break;
    case 3: /* right edge */
        fprintf(my_plot_file,
                "(%s) dup stringwidth pop  0 rmoveto show \n grestore\n", str);
        break;
    default:
        break;
    }
    return;
}

void
array_print_ps_set_text(double angle, double slant, double x_size,
                        double y_size) {
    ps_scale.tx = x_size*5.0;
    ps_scale.ty = y_size*5.0;
    ps_scale.angle = angle;
    ps_scale.slant = slant;
    return;
}

void
array_print_ps_rect(double x, double y, double wid, double len) {
    double x1, y1, x2, y2;
    array_print_ps_convert(x, y, &x1, &y1);
    array_print_ps_convert(x + wid, y + len, &x2, &y2);
    fprintf(my_plot_file,
            "%d %d m \n %d %d l \n %d %d l \n %d %d l \n %d %d l \n S \n",
            (int32)x1, (int32)y1, (int32)x2, (int32)y1, (int32)x2, (int32)y2,
            (int32)x1, (int32)y2, (int32)x1, (int32)y1);
    return;
}

void
array_print_ps_bar(double x, double y, double wid, double len, double fill,
                   int32 flag) {
    double x1, y1, x2, y2;
    fprintf(my_plot_file, "%f G\n", fill);
    array_print_ps_convert(x, y, &x1, &y1);
    array_print_ps_convert(x + wid, y + len, &x2, &y2);
    fprintf(my_plot_file, "%d %d m \n %d %d l \n %d %d l \n %d %d l \n FS\n",
            (int32)x1, (int32)y1, (int32)x2, (int32)y1, (int32)x2, (int32)y2,
            (int32)x1, (int32)y2);

    if (flag) {
        fprintf(my_plot_file, "0 G\n");
        array_print_ps_rect(x, y, wid, len);
    }
    return;
}

void
ps_rgb_bar(double x, double y, double wid, double len, double fill, int32 flag,
           int32 rgb) {
    double x1, y1, x2, y2;
    double r = 0.0, g = 0.0, b = 0.0;
    if (rgb == 2) {
        ps_hsb_bar(x, y, wid, len, fill, flag);
        return;
    }
    if (fill < 0.0)
        fill = 0.0;
    if (fill > 1.0)
        fill = 1.0;
    switch (rgb) {
    case 0:
        fill = 1. - fill;
        b = (double)sqrt((double)(1.0 - fill*fill));
        r = (double)sqrt((double)(fill*(2.0 - fill)));
        break;
    case 1:
        if (fill > .4999)
            r = 0.0;
        else
            r = (double)sqrt((double)(1. - 4*fill*fill));
        g = (double)2*sqrt((double)fill*(1. - fill));

        if (fill < .5001)
            b = 0.0;
        else
            b = (double)sqrt((double)(4*(fill - .5)*(1.5 - fill)));
        break;
    default:
        fprintf(stderr, "ps_rgb_bar: rgb is neither 0 nor 1.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(my_plot_file, "%f %f %f RGB\n", r, g, b);
    array_print_ps_convert(x, y, &x1, &y1);
    array_print_ps_convert(x + wid, y + len, &x2, &y2);
    /*   fprintf(my_plot_file,"%f %f m \n %f %f l \n %f %f l \n %f %f l \n
       FS\n", x1,y1,x2,y1,x2,y2,x1,y2); */
    fprintf(my_plot_file, "%d %d m \n %d %d l \n %d %d l \n %d %d l \n FS\n",
            (int32)x1, (int32)y1, (int32)x2, (int32)y1, (int32)x2, (int32)y2,
            (int32)x1, (int32)y2);
    if (flag) {
        fprintf(my_plot_file, "0 G\n");
        array_print_ps_rect(x, y, wid, len);
    }
    return;
}

void
ps_hsb_bar(double x, double y, double wid, double len, double fill,
           int32 flag) {
    double x1, y1, x2, y2;
    fprintf(my_plot_file, "%f 1.0 1.0 HSB\n", fill);
    array_print_ps_convert(x, y, &x1, &y1);
    array_print_ps_convert(x + wid, y + len, &x2, &y2);
    /* fprintf(my_plot_file,"%f %f m \n %f %f l \n %f %f l \n %f %f l \n FS\n",
             x1,y1,x2,y1,x2,y2,x1,y2); */
    fprintf(my_plot_file, "%d %d m \n %d %d l \n %d %d l \n %d %d l \n FS\n",
            (int32)x1, (int32)y1, (int32)x2, (int32)y1, (int32)x2, (int32)y2,
            (int32)x1, (int32)y2);
    if (flag) {
        fprintf(my_plot_file, "0 G\n");
        array_print_ps_rect(x, y, wid, len);
    }
}
