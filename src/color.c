#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "functions.h"
#include "integers.h"

#define RED 20
#define REDORANGE 21
#define ORANGE 22
#define YELLOWORANGE 23
#define YELLOW 24
#define YELLOWGREEN 25
#define GREEN 26
#define BLUEGREEN 27
#define BLUE 28
#define PURPLE 29

#define C_NORM 0
#define C_PERIODIC 1
#define C_HOT 2
#define C_COOL 3
#define C_REDBLUE 4
#define C_GRAY 5
#define C_CUBHLX 6

static int32 color_mode = 1;
int32 color_min;
int32 color_total;
int32 COLOR;
int32 color_max;
int32 custom_color = 0;

#define MAX_COLORS 256
#define COL_TOTAL 150

static XColor color[MAX_COLORS];

static int32 color_bfun(double y, int32 per);
static int32 color_gfun(double y, int32 per);
static int32 color_rfun(double y, int32 per);
static void color_make_cmaps(int32 *r, int32 *g, int32 *b, int32 n, int32 type);

void
color_set_s(int32 col) {
    if (col < 0)
        XSetForeground(display, small_gc, GrBack);
    if (col == 0)
        XSetForeground(display, small_gc, GrFore);
    else {
        if (COLOR)
            XSetForeground(display, small_gc, (ulong)color_map(col));
        else
            XSetForeground(display, small_gc, GrFore);
    }
    return;
}

void
color_set(int32 col) {
    if (col < 0)
        XSetForeground(display, gc_graph, GrBack);
    if (col == 0)
        XSetForeground(display, gc_graph, GrFore);
    else {

        if (COLOR)
            XSetForeground(display, gc_graph, (ulong)color_map(col));
        else
            XSetForeground(display, gc_graph, GrFore);
    }
    return;
}

/* this makes alot of nice color maps */
void
color_make_cmaps(int32 *r, int32 *g, int32 *b, int32 n, int32 type) {
    double x;
    int32 i1;
    int32 i2;
    int32 i3;
    double pii = 3.1415926;
    /* for CUBHLX  */
    double start = .5;
    double rots = -1.5;
    double hue = 1.2;
    double gamma = 1.;
    double angle;
    double amp;
    double rr;
    double gg;
    double bb;

    switch (type) {
    case C_NORM:
        for (int32 i = 0; i < n; i++) {
            x = (double)i / ((double)n);
            r[i] = color_rfun(1 - x, 0) << 8;
            g[i] = color_gfun(1 - x, 0) << 8;
            b[i] = color_bfun(1 - x, 0) << 8;
        }
        break;
    case C_PERIODIC:
        for (int32 i = 0; i < n; i++) {
            x = (double)i / ((double)n);
            r[i] = color_rfun(x, 1) << 8;
            g[i] = color_gfun(x, 1) << 8;
            b[i] = color_bfun(x, 1) << 8;
        }
        break;
    case C_HOT:
        i1 = (int32)(.375*(double)n);
        i2 = 2*i1;
        i3 = n - i2;

        for (int32 i = 0; i < i1; i++) {
            x = 256*255*(double)i / ((double)i1);

            r[i] = (int32)x;
            g[i] = 0;
            b[i] = 0;
            g[i + i1] = (int32)x;
            b[i + i1] = 0;
        }

        for (int32 i = i1; i < n; i++)
            r[i] = 256*255;
        for (int32 i = i2; i < n; i++) {
            x = 256*255*(double)(i - i2) / ((double)i3);

            g[i] = 256*255;
            b[i] = (int32)x;
        }
        break;
    case C_COOL:
        for (int32 i = 0; i < n; i++) {
            x = (double)i / ((double)n);
            r[i] = (int32)(256*255*x);
            b[i] = (int32)(256*255*(1 - x));
            g[i] = 256*255;
        }
        break;
    case C_REDBLUE:
        for (int32 i = 0; i < n; i++) {
            x = (double)i / ((double)n);
            r[i] = (int32)(256*255*x);
            b[i] = (int32)(256*255*(1 - x));
            g[i] = 0;
        }
        break;

    case C_GRAY:
        for (int32 i = 0; i < n; i++) {
            r[i] = i*256*255 / n;
            b[i] = i*256*255 / n;
            g[i] = i*256*255 / n;
        }
        break;
        /* https://www.mrao.cam.ac.uk/~dag/CUBEHELIX/ */
    case C_CUBHLX:
        for (int32 i = 0; i < n; i++) {
            x = (double)i / ((double)n);
            angle = 2*pii*(start / 3.0 + 1 + rots*x);
            x = pow(x, gamma);
            amp = hue*x*(1 - x) / 2.0;
            rr = x + amp*(-.14861*cos(angle) + 1.78277*sin(angle));
            gg = x + amp*(-.29227*cos(angle) - .90649*sin(angle));
            bb = x + amp*(1.97294*cos(angle));

            if (rr < 0.0)
                rr = 0.0;
            if (rr > 1.0)
                rr = 1.0;
            if (gg < 0.0)
                gg = 0.0;
            if (gg > 1.0)
                gg = 1.0;
            if (bb < 0.0)
                bb = 0.0;
            if (bb > 1.0)
                bb = 1.0;
            r[i] = 256*(int32)(255.*rr);
            b[i] = 256*(int32)(255.*bb);
            g[i] = 256*(int32)(255.*gg);
        }
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

int32
color_rfun(double y, int32 per) {
    double x;
    x = y;
    if ((y > .666666) && (per == 1))
        x = 1. - y;

    if (x > .33333333333)
        return 0;
    return (int32)(3.*255*sqrt((.333334 - x)*(x + .33334)));
}

int32
color_gfun(double y, int32 per) {
    (void)per;
    if (y > .666666)
        return 0;
    return (int32)(3.*255*sqrt((.6666667 - y)*(y)));
}

int32
color_bfun(double y, int32 per) {
    (void)per;
    if (y < .333334)
        return 0;
    return (int32)(2.79*255*sqrt((1.05 - y)*(y - .333333333)));
}

void
color_new_map(int32 type) {
    /*  printf(" My color map = %d\n",type); */
    if (TrueColorFlag == 0) {
        ggets_err_msg("New colormaps not supported without TrueColor");
        return;
    }
    custom_color = type;
    color_map_make();
    return;
}

void
color_get_ps(int32 i, double *r, double *g, double *b) {
    double z = 1. / (65535);
    *r = z*(double)color[i].red;
    *g = z*(double)color[i].green;
    *b = z*(double)color[i].blue;
    return;
}

void
color_get_svg(int32 i, int32 *r, int32 *g, int32 *b) {
    *r = color[i].red / 255;
    *g = color[i].green / 255;
    *b = color[i].blue / 255;
    return;
}

void
color_map_make(void) {
    Colormap cmap;
    int32 clo = 20;

    int32 r[256];
    int32 g[256];
    int32 b[256];

    cmap = (Colormap)NULL;
    color_min = 30;
    color_max = MAX_COLORS - 1;
    color_total = color_max - color_min + 1;
    if (color_total > COL_TOTAL)
        color_total = COL_TOTAL;
    color_max = color_min + color_total;
    if (Xup)
        cmap = DefaultColormap(display, screen);
    for (int32 i = 0; i < clo; i++) {
        color[i].pixel = (ulong)i;
    }
    for (int32 i = 20; i < 30; i++) {
        color[i].red = 0;
        color[i].blue = 0;
        color[i].green = 0;
        color[i].flags = DoRed | DoGreen | DoBlue;
    }

    color[RED].red = 255;
    color[BLUE].blue = 255;
    color[GREEN].green = 225;
    color[YELLOWGREEN].red = 200;
    color[YELLOWGREEN].blue = 75;
    color[YELLOWGREEN].green = 235;
    color[REDORANGE].red = 240;
    color[REDORANGE].green = 100;
    color[ORANGE].red = 255;
    color[ORANGE].green = 165;
    color[YELLOWORANGE].red = 255;
    color[YELLOWORANGE].green = 205;
    color[YELLOW].red = 200;
    color[YELLOW].green = 200;
    color[BLUEGREEN].blue = 200;
    color[BLUEGREEN].green = 200;
    color[PURPLE].red = 160;
    color[PURPLE].green = 32;
    color[PURPLE].blue = 240;
    for (int32 i = 20; i < 30; i++) {
        color[i].red = (ushort)(color[i].red << 8);
        color[i].blue = (ushort)(color[i].blue << 8);
        color[i].green = (ushort)(color[i].green << 8);
        color[i].flags = DoRed | DoGreen | DoBlue;
        if (Xup) {
            XAllocColor(display, cmap, &color[i]);
        }
    }

    color_make_cmaps(r, g, b, color_total + 1, custom_color);
    for (int32 i = color_min; i <= color_max; i++) {
        color[i].red = (ushort)r[i - color_min];
        color[i].green = (ushort)g[i - color_min];
        color[i].blue = (ushort)b[i - color_min];

        color[i].flags = DoRed | DoGreen | DoBlue;
        if (Xup) {
            XAllocColor(display, cmap, &color[i]);
        }
    }
    return;
}

uint32
color_map(int32 i) {
    if (i == -1)
        return GrBack;

    if (i == 0)
        return GrFore;
    if (color_mode) {
        if (i < 0)
            i = 0;
        if (i >= color_max)
            i = color_max;
        return (uint32)color[i].pixel;
    } else {
        return (uint32)i;
    }
}
