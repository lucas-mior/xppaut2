#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "struct.h"

#define MAXPERPLOT 10
#define DEGTORAD .0174532
#define EP1 1.000001
#define SYMSIZE .00175

double THETA0 = 45;
double PHI0 = 45;

int32 PS_Port = 0;
int32 PointRadius = 0;

static char dashes[10][5] = {{0},       {1, 6, 0}, {0},       {4, 2, 0},
                             {1, 3, 0}, {4, 4, 0}, {1, 5, 0}, {4, 4, 4, 1, 0},
                             {4, 2, 0}, {1, 3, 0}};

XFontStruct *symfonts[5];
XFontStruct *romfonts[5];
int32 avsymfonts[5];
int32 avromfonts[5];

/*  This is an improved graphics driver for XPP
    It requires only a few commands
    All positions are integers

    graphics_point(x,y)        Draws point to (x,y) with pointtype PointStyle
    graphics_line(x1,y1,x2,y2) Draws line with linetype LineStyle
    graphics_put_text(x1,y,text)  Draws text with TextAngle, Justify
    init_device()     Sets up the default for tics, plotting area,
                      and anything else

    close_device()  closes files, etc

    The device is assumed to go from (0,0) to (XDMax,YDMax)

    DLeft,DRight,DTop,DBottom are the actual graph areas
    VTic Htic are the actual sizes of the tics in device pixels
    VChar HChar are the height and width used for character spacing

    linetypes   -2  thick plain lines
                -1  thin plain lines

*/

int32 DLeft;
int32 DRight;
int32 DTop;
int32 DBottom;
int32 VTic;
int32 HTic;
int32 VChar;
int32 HChar;
int32 XDMax;
int32 YDMax;
double XMin;
double YMin;
double XMax;
double YMax;
int32 PointType = -1;
int32 TextJustify;
int32 TextAngle;

static void graphics_draw_symbol(double x, double y, double size,
                                 int32 my_symb);
static void graphics_line_nabs(double x1_out, double y1_out, double x2_out,
                               double y2_out);
static void graphics_pers_line(double x, double y, double z, double xp,
                               double yp, double zp);
static void graphics_init(int32 i);
static void graphics_line_x11(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
static void graphics_rect_x11(int32 x, int32 y, int32 w, int32 h);
static void graphics_bead_x11(int32 x, int32 y);
static void graphics_point_x11(int32 xp, int32 yp);
static void graphics_bead(int32 x1, int32 y1);
static void graphics_rot_3dvec(double x, double y, double z, double *xp,
                               double *yp, double *zp);

void
graphics_get_scale(double *x1, double *y1, double *x2, double *y2) {
    *x1 = XMin;
    *y1 = YMin;
    *x2 = XMax;
    *y2 = YMax;
    return;
}

void
graphics_set_scale(double x1, double y1, double x2, double y2) {
    XMin = x1;
    YMin = y1;
    XMax = x2;
    YMax = y2;
    return;
}

void
graphics_get_draw_area(void) {
    graphics_get_draw_area_flag(1);
    return;
}

void
graphics_get_draw_area_flag(int32 flag) {
    int32 x;
    int32 y;
    uint32 w;
    uint32 h;
    uint32 bw;
    uint32 de;
    Window root;

    if (flag == 1) {
        XGetGeometry(display, draw_win, &root, &x, &y, &w, &h, &bw, &de);
        MyGraph->x11Wid = (int32)w;
        MyGraph->x11Hgt = (int32)h;
    } else {
        w = (uint32)MyGraph->x11Wid;
        h = (uint32)MyGraph->x11Hgt;
    }

    XDMax = (int32)w;
    YDMax = (int32)h;
    VTic = MAX(h / 100, 1);
    HTic = MAX(w / 150, 1);
    VChar = DCURYs;
    HChar = DCURXs;

    DLeft = 12*HChar;
    DRight = XDMax - 3*HChar - HTic;
    DBottom = YDMax - 1 - VChar*7 / 2;
    DTop = VChar*5 / 2 + 1;
    h = (uint32)(DBottom - DTop);
    w = (uint32)(DRight - DLeft);
    MyGraph->Width = (int32)w;
    MyGraph->Height = (int32)h;
    MyGraph->x0 = DLeft;
    MyGraph->y0 = DTop;
    graphics_set_normal_scale();
    return;
}

void
graphics_change_current_linestyle(int32 new, int32 *old) {
    *old = MyGraph->color[0];
    MyGraph->color[0] = new;
    return;
}

void
graphics_set_normal_scale(void) {
    XMin = MyGraph->xlo;
    YMin = MyGraph->ylo;
    XMax = MyGraph->xhi;
    YMax = MyGraph->yhi;
    return;
}

void
graphics_point(int32 x, int32 y) {
    if (PltFmtFlag == PSFMT) {
        ps_point(x, y);
    } else if (PltFmtFlag == SVGFMT) {
        svg_point(x, y);
    } else {
        graphics_point_x11(x, y);
    }
    return;
}

void
graphics_line(int32 x1, int32 y1, int32 x2, int32 y2) {
    if (PltFmtFlag == PSFMT) {
        ps_line(x1, y1, x2, y2);
    } else if (PltFmtFlag == SVGFMT) {
        svg_line(x1, y1, x2, y2);
    } else {
        graphics_line_x11(x1, y1, x2, y2);
    }
}
/* draw a little filled circle */

void
graphics_bead(int32 x1, int32 y1) {
    if (PltFmtFlag == PSFMT) {
        return;
    } else if (PltFmtFlag == SVGFMT) {
        svg_bead();
    } else {
        graphics_bead_x11(x1, y1);
    }
    return;
}

void
graphics_frect(int32 x1, int32 y1, int32 w, int32 h) {
    if (PltFmtFlag == PSFMT) {
        ps_frect(x1, y1, w, h);
    } else if (PltFmtFlag == SVGFMT) {
        svg_frect(x1, y1, w, h);
    } else {
        graphics_rect_x11(x1, y1, w, h);
    }
    return;
}

void
graphics_put_text(int32 x, int32 y, char *str) {
    if (PltFmtFlag == PSFMT) {
        ps_text(x, y, str);
    } else if (PltFmtFlag == SVGFMT) {
        svg_text(x, y, str);
    } else {
        graphics_put_text_x11(x, y, str);
    }
    return;
}

void
graphics_init_x11(void) {
    graphics_get_draw_area();
    return;
}

void
graphics_init_ps(void) {
    if (!PS_Port) {
        XDMax = 7200;
        YDMax = 5040;
        VTic = 63;
        HTic = 63;
        VChar = 140;
        HChar = 84;
        DLeft = 12*HChar;
        DRight = XDMax - 3*HChar - HTic;
        DTop = YDMax - 1 - VChar*7 / 2;
        DBottom = VChar*5 / 2 + 1;
    } else {
        YDMax = 7200;
        XDMax = 5040;
        VTic = 63;
        HTic = 63;
        VChar = 140;
        HChar = 84;
        DLeft = 12*HChar;
        DRight = XDMax - 3*HChar - HTic;
        DTop = YDMax - 1 - VChar*7 / 2;
        DBottom = VChar*5 / 2 + 1;
    }
    return;
}

void
graphics_init_svg(void) {
    XDMax = 640;
    YDMax = 400;
    VTic = 9;
    HTic = 9;
    VChar = 20;
    HChar = 12;
    DLeft = 12*HChar;
    DRight = XDMax - 3*HChar - HTic;
    DBottom = YDMax - 1 - VChar*7 / 2;
    DTop = VChar*5 / 2 + 1;
    return;
}

void
graphics_point_x11(int32 xp, int32 yp) {
    int32 r = PointRadius;
    int32 r2 = (int32)(r / 1.41421356 + 0.5);
    int32 wh = 2*r2;
    if (PointRadius == 0) {
        XDrawPoint(display, draw_win, gc_graph, xp, yp);
    } else {
        XFillArc(display, draw_win, gc_graph, xp - r2, yp - r2, (uint)wh,
                 (uint)wh, 0, 360*64);
    }
    return;
}

void
graphics_set_linestyle(int32 ls) {
    if (PltFmtFlag == PSFMT) {
        ps_linetype(ls);
    } else if (PltFmtFlag == SVGFMT) {
        svg_linetype(ls);
    } else {
        /* graphics set line style x11 */
        int32 type = 0;
        if (ls == -2) { /*  Border  */
            color_set(0);
            XSetLineAttributes(display, gc_graph, 2, LineSolid, CapButt,
                               JoinBevel);
            return;
        }
        /*width=0;
         */
        if (ls == -1) {
            color_set(0);
            XSetDashes(display, gc_graph, 0, dashes[1], (int)strlen(dashes[1]));
            XSetLineAttributes(display, gc_graph, 0, LineOnOffDash, CapButt,
                               JoinBevel);
            return;
        }
        if (!COLOR) { /* Mono  */
            ls = (ls % 8) + 2;
            if (ls == 2) {
                type = LineSolid;
            } else {
                type = LineOnOffDash;
                XSetDashes(display, gc_graph, 0, dashes[ls],
                           (int)strlen(dashes[ls]));
            }
            color_set(0);
            XSetLineAttributes(display, gc_graph, 0, type, CapButt, JoinBevel);
            return;
        }
        /* color system  */
        ls = ls % 11;
        XSetLineAttributes(display, gc_graph, 0, LineSolid, CapButt, JoinBevel);
        color_set(colorline[ls]);
    }
    return;
}

void
graphics_bead_x11(int32 x, int32 y) {
    XFillArc(display, draw_win, gc_graph, x - 2, y - 2, 4, 4, 0, 360*64);
    return;
}

void
graphics_rect_x11(int32 x, int32 y, int32 w, int32 h) {
    XFillRectangle(display, draw_win, gc_graph, x, y, (uint)w, (uint)h);
    return;
}

void
graphics_draw_many_lines(void) {
    int32 NLINE = 500000;
    for (int32 i = 0; i < NLINE; i++) {
        XDrawLine(display, draw_win, gc_graph, rand() % 200, rand() % 200,
                  rand() % 200, rand() % 200);
    }
    printf("Done\n");
    return;
}

void
graphics_line_x11(int32 xp1, int32 yp1, int32 xp2, int32 yp2) {
    XDrawLine(display, draw_win, gc_graph, xp1, yp1, xp2, yp2);
    return;
}

void
graphics_put_text_x11(int32 x, int32 y, char *str) {
    int32 sw = (int32)strlen(str)*DCURXs;
    switch (TextJustify) {
    case 0:
        sw = 0;
        break;
    case 1:
        sw = -sw / 2;
        break;
    case 2:
        sw = -sw;
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    XSetForeground(display, small_gc, GrFore);
    XDrawString(display, draw_win, small_gc, x + sw, y + DCURYs / 3, str,
                (int32)strlen(str));
    XSetForeground(display, small_gc, GrBack);
    return;
}

void
graphics_special_put_text_x11(int32 x, int32 y, char *str, int32 size) {
    int32 i = 0;
    int32 j = 0;
    int32 cx = x;
    int32 cy = y;
    int32 cf = 0;
    int32 cs;
    int32 n = (int32)strlen(str);
    int32 dx = 0;
    char tmp[256];
    char c;
    int32 sub;
    int32 sup;
    cs = size;
    if (avromfonts[size] == 1) {
        sup = romfonts[size]->ascent;
        sub = sup / 2;
    } else {
        sup = font_small->ascent;
        sub = sup / 2;
    }
    while (i < n) {
        c = str[i];
        if (c == '\\') {
            i++;
            c = str[i];
            tmp[j] = 0; /* end the current buffer */

            graphics_fancy_put_text_x11(cx, cy, tmp, cs,
                                        cf); /* render the current buffer */
            if (cf == 0) {
                if (avromfonts[cs] == 1) {
                    dx = XTextWidth(romfonts[cs], tmp, (int)strlen(tmp));
                } else {
                    dx = XTextWidth(font_small, tmp, (int)strlen(tmp));
                }
            }
            if (cf == 1) {
                if (avsymfonts[cs] == 1) {
                    dx = XTextWidth(symfonts[cs], tmp, (int)strlen(tmp));
                } else {
                    dx = XTextWidth(font_small, tmp, (int)strlen(tmp));
                }
            }
            cx += dx;
            j = 0;
            if (c == '0') {
                cf = 0;
            }
            if (c == 'n') {
                cy = y;
                cs = size;
            }
            if (c == 's') {
                cy = cy + sub;
                if (size > 0) {
                    cs = size - 1;
                }
            }
            if (c == 'S') {
                if (size > 0) {
                    cs = size - 1;
                }
                cy = cy - sup;
            }
            if (c == '1') {
                cf = 1;
            }

            i++;
        } else {
            tmp[j] = c;
            j++;
            i++;
        }
    }
    tmp[j] = 0;
    graphics_fancy_put_text_x11(cx, cy, tmp, cs, cf);
    return;
}

void
graphics_fancy_put_text_x11(int32 x, int32 y, char *str, int32 size,
                            int32 font) {
    /*int32 yoff;
     */
    if (strlen(str) == 0) {
        return;
    }
    switch (font) {
    case 1:
        if (avsymfonts[size] == 1) {
            XSetFont(display, font_gc, symfonts[size]->fid);
        } else {
            XSetFont(display, font_gc, font_small->fid);
        }
        XSetForeground(display, font_gc, GrFore);
        XDrawString(display, draw_win, font_gc, x, y, str, (int)strlen(str));
        XSetForeground(display, font_gc, GrBack);
        break;
    default:
        if (avromfonts[size] == 1) {
            XSetFont(display, font_gc, romfonts[size]->fid);

        } else {
            XSetFont(display, font_gc, font_small->fid);
        }
        XSetForeground(display, font_gc, GrFore);
        XDrawString(display, draw_win, font_gc, x, y, str, (int)strlen(str));
        XSetForeground(display, font_gc, GrBack);
        break;
    }
    return;
}

void
graphics_scale_dxdy(double x, double y, double *i, double *j) {
    double dx = (DRight - DLeft) / (XMax - XMin);
    double dy = (DTop - DBottom) / (YMax - YMin);
    *i = x*dx;
    *j = y*dy;
    return;
}

void
graphics_scale_to_screen(/* not really the screen!  */
                         double x, double y, int32 *i, int32 *j) {
    double dx = (DRight - DLeft) / (XMax - XMin);
    double dy = (DTop - DBottom) / (YMax - YMin);
    *i = (int32)((x - XMin)*dx) + DLeft;
    *j = (int32)((y - YMin)*dy) + DBottom;
    return;
}

void
graphics_scale_to_real(/* Not needed except for X */
                       int32 i, int32 j, double *x, double *y) {
    int32 i1;
    int32 j1;
    double x1;
    double y1;
    graphics_get_draw_area();
    i1 = i - DLeft;
    j1 = j - DBottom;
    x1 = (double)i1;
    y1 = (double)j1;
    *x = (MyGraph->xhi - MyGraph->xlo)*x1 / ((double)(DRight - DLeft)) +
         MyGraph->xlo;
    *y = (MyGraph->yhi - MyGraph->ylo)*y1 / ((double)(DTop - DBottom)) +
         MyGraph->ylo;
    return;
}

void
graphics_reset_all_line_type(void) {
    for (int32 j = 0; j < MAXPOP; j++) {
        for (int32 k = 0; k < MAXPERPLOT; k++) {
            graph[j].line[k] = START_LINE_TYPE;
        }
    }
    return;
}

void
graphics_init_all(void) {
    for (int32 i = 0; i < MAXPOP; i++) {
        graphics_init(i);
    }
    MyGraph = &graph[0];
    graphics_set_normal_scale();
    return;
}

void
graphics_set_extra(void) {
    if (NPltV < 2) {
        return;
    }
    if (NPltV > 8) {
        NPltV = 8;
    }
    if (MultiWin == 0) {
        MyGraph->nvars = NPltV;
        for (int32 i = 1; i < NPltV; i++) {
            MyGraph->xv[i] = IX_PLT[i + 1];
            MyGraph->yv[i] = IY_PLT[i + 1];
            MyGraph->zv[i] = IZ_PLT[i + 1];
            MyGraph->color[i] = i;
        }
        return;
    }
    if (Xup) {
        for (int32 i = 1; i < NPltV; i++) {
            many_pops_create_a_pop();
            graph[i].xv[0] = IX_PLT[i + 1];
            graph[i].yv[0] = IY_PLT[i + 1];
            graph[i].zv[0] = IZ_PLT[i + 1]; /* irrelevant probably */
            graph[i].grtype = 0;            /* force 2D */
            graph[i].xlo = X_LO[i + 1];
            graph[i].xhi = X_HI[i + 1];
            graph[i].ylo = Y_LO[i + 1];
            graph[i].yhi = Y_HI[i + 1];
            /*  printf(" %g %g %g %g
             * \n",X_LO[i+1],X_HI[i+1],Y_LO[i+1],Y_HI[i+1]); */
        }
        many_pops_set_active_windows();
        many_pops_make_active(0, 1);
    }
    return;
}

void
graphics_reset_graph(void) {
    if (AXES >= 5) {
        PLOT_3D = 1;
    } else {
        PLOT_3D = 0;
    }
    MyGraph->xv[0] = IXPLT;
    MyGraph->yv[0] = IYPLT;
    MyGraph->zv[0] = IZPLT;
    MyGraph->xmax = x_3d[1];
    MyGraph->ymax = y_3d[1];
    MyGraph->zmax = z_3d[1];
    MyGraph->xbar = .5*(x_3d[1] + x_3d[0]);
    MyGraph->ybar = .5*(y_3d[1] + y_3d[0]);
    MyGraph->zbar = .5*(z_3d[1] + z_3d[0]);
    MyGraph->dx = 2. / (x_3d[1] - x_3d[0]);
    MyGraph->dy = 2. / (y_3d[1] - y_3d[0]);
    MyGraph->dz = 2. / (z_3d[1] - z_3d[0]);
    MyGraph->xmin = x_3d[0];
    MyGraph->ymin = y_3d[0];
    MyGraph->zmin = z_3d[0];
    MyGraph->xlo = MY_XLO;
    MyGraph->ylo = MY_YLO;
    MyGraph->xhi = MY_XHI;
    MyGraph->yhi = MY_YHI;
    MyGraph->grtype = AXES;
    graf_par_check_windows();
    graphics_set_normal_scale();
    graf_par_redraw_the_graph();
}

void
graphics_get_graph(void) {
    x_3d[0] = MyGraph->xmin;
    x_3d[1] = MyGraph->xmax;
    y_3d[0] = MyGraph->ymin;
    y_3d[1] = MyGraph->ymax;
    z_3d[0] = MyGraph->zmin;
    z_3d[1] = MyGraph->zmax;
    MY_XLO = MyGraph->xlo;
    MY_YLO = MyGraph->ylo;
    MY_XHI = MyGraph->xhi;
    MY_YHI = MyGraph->yhi;
    IXPLT = MyGraph->xv[0];
    IYPLT = MyGraph->yv[0];
    IZPLT = MyGraph->zv[0];
    PLOT_3D = MyGraph->ThreeDFlag;
    if (PLOT_3D) {
        AXES = 5;
    } else {
        AXES = 0;
    }
    AXES = MyGraph->grtype;
    return;
}

void
graphics_init(int32 i) {
    if (AXES <= 3) {
        AXES = 0;
    }
    for (int32 j = 0; j < 3; j++) {
        for (int32 k = 0; k < 3; k++) {
            if (k == j) {
                graph[i].rm[k][j] = 1.0;
            } else {
                graph[i].rm[k][j] = 0.0;
            }
        }
    }
    graph[i].nvars = 1;
    for (int32 j = 0; j < MAXPERPLOT; j++) {
        graph[i].xv[j] = IXPLT;
        graph[i].yv[j] = IYPLT;
        graph[i].zv[j] = IZPLT;
        graph[i].line[j] = START_LINE_TYPE;
        graph[i].color[j] = 0;
    }

    /*sprintf(graph[i].xlabel,"");
    sprintf(graph[i].ylabel,"");
    sprintf(graph[i].zlabel,"");
    */
    graph[i].xlabel[0] = '\0';
    graph[i].ylabel[0] = '\0';
    graph[i].zlabel[0] = '\0';

    graph[i].Use = 0;
    graph[i].state = 0;
    graph[i].Restore = 1;
    graph[i].Nullrestore = 0;
    graph[i].ZPlane = -1000.0;
    graph[i].ZView = 1000.0;
    graph[i].PerspFlag = 0;
    graph[i].ThreeDFlag = PLOT_3D;
    graph[i].TimeFlag = TIMPLOT;
    graph[i].ColorFlag = 0;
    graph[i].grtype = AXES;
    graph[i].color_scale = 1.0;
    graph[i].min_scale = 0.0;
    strcpy(graph[i].gr_info, "");
    graph[i].xmax = x_3d[1];
    graph[i].ymax = y_3d[1];
    graph[i].zmax = z_3d[1];
    graph[i].xbar = .5*(x_3d[1] + x_3d[0]);
    graph[i].ybar = .5*(y_3d[1] + y_3d[0]);
    graph[i].zbar = .5*(z_3d[1] + z_3d[0]);
    graph[i].dx = 2. / (x_3d[1] - x_3d[0]);
    graph[i].dy = 2. / (y_3d[1] - y_3d[0]);
    graph[i].dz = 2. / (z_3d[1] - z_3d[0]);
    graph[i].xmin = x_3d[0];
    graph[i].ymin = y_3d[0];
    graph[i].zmin = z_3d[0];
    graph[i].xorg = 0.0;
    graph[i].yorg = 0.0;
    graph[i].yorg = 0.0;
    graph[i].xorgflag = 1;
    graph[i].yorgflag = 1;
    graph[i].zorgflag = 1;
    graph[i].Theta = THETA0;
    graph[i].Phi = PHI0;
    graph[i].xshft = 0;
    graph[i].yshft = 0;
    graph[i].zshft = 0;
    graph[i].xlo = MY_XLO;
    graph[i].ylo = MY_YLO;
    graph[i].oldxlo = MY_XLO;
    graph[i].oldylo = MY_YLO;
    graph[i].xhi = MY_XHI;
    graph[i].yhi = MY_YHI;
    graph[i].oldxhi = MY_XHI;
    graph[i].oldyhi = MY_YHI;
    MyGraph = &graph[i];
    graphics_make_rot(THETA0, PHI0);
    return;
}

void
graphics_copy_graph(/*  Graph[i]=Graph[l]  */
                    int32 i, int32 l) {
    graph[i].Use = graph[l].Use;
    graph[i].Restore = graph[l].Restore;
    graph[i].Nullrestore = graph[l].Nullrestore;
    for (int32 j = 0; j < 3; j++) {
        for (int32 k = 0; k < 3; k++) {
            graph[i].rm[k][j] = graph[l].rm[k][j];
        }
    }
    graph[i].nvars = graph[l].nvars;
    for (int32 j = 0; j < MAXPERPLOT; j++) {
        graph[i].xv[j] = graph[l].xv[j];
        graph[i].yv[j] = graph[l].yv[j];
        graph[i].zv[j] = graph[l].zv[j];
        graph[i].line[j] = graph[l].line[j];
        graph[i].color[j] = graph[l].color[j];
    }

    graph[i].ZPlane = graph[l].ZPlane;
    graph[i].ZView = graph[l].ZView;
    graph[i].PerspFlag = graph[l].PerspFlag;
    graph[i].ThreeDFlag = graph[l].ThreeDFlag;
    graph[i].TimeFlag = graph[l].TimeFlag;
    graph[i].ColorFlag = graph[l].ColorFlag;
    graph[i].grtype = graph[l].grtype;
    graph[i].color_scale = graph[l].color_scale;
    graph[i].min_scale = graph[l].min_scale;

    graph[i].xmax = graph[l].xmax;
    graph[i].xmin = graph[l].xmin;
    graph[i].ymax = graph[l].ymax;
    graph[i].ymin = graph[l].ymin;
    graph[i].zmax = graph[l].zmax;
    graph[i].zmin = graph[l].zmin;
    graph[i].xbar = graph[l].xbar;
    graph[i].dx = graph[l].dx;
    graph[i].ybar = graph[l].ybar;
    graph[i].dy = graph[l].dy;
    graph[i].zbar = graph[l].zbar;
    graph[i].dz = graph[l].dz;

    graph[i].Theta = graph[l].Theta;
    graph[i].Phi = graph[l].Phi;
    graph[i].xshft = graph[l].xshft;
    graph[i].yshft = graph[l].yshft;
    graph[i].zshft = graph[l].zshft;
    graph[i].xlo = graph[l].xlo;
    graph[i].ylo = graph[l].ylo;
    graph[i].oldxlo = graph[l].oldxlo;
    graph[i].oldylo = graph[l].oldylo;
    graph[i].xhi = graph[l].xhi;
    graph[i].yhi = graph[l].yhi;
    graph[i].oldxhi = graph[l].oldxhi;
    graph[i].oldyhi = graph[l].oldyhi;
    return;
}

void
graphics_make_rot(double theta, double phi) {
    double ct = cos(DEGTORAD*theta);
    double st = sin(DEGTORAD*theta);
    double sp = sin(DEGTORAD*phi);
    double cp = cos(DEGTORAD*phi);
    MyGraph->Theta = theta;
    MyGraph->Phi = phi;
    MyGraph->rm[0][0] = ct;
    MyGraph->rm[0][1] = st;
    MyGraph->rm[0][2] = 0.0;
    MyGraph->rm[1][0] = -cp*st;
    MyGraph->rm[1][1] = cp*ct;
    MyGraph->rm[1][2] = sp;
    MyGraph->rm[2][0] = st*sp;
    MyGraph->rm[2][1] = -sp*ct;
    MyGraph->rm[2][2] = cp;
    return;
}

void
graphics_scale3d(double x, double y, double z, double *xp, double *yp,
                 double *zp) {
    *xp = (x - MyGraph->xbar)*MyGraph->dx;
    *yp = (y - MyGraph->ybar)*MyGraph->dy;
    *zp = (z - MyGraph->zbar)*MyGraph->dz;
    return;
}

int32
graphics_threedproj(double x2p, double y2p, double z2p, double *xp,
                    double *yp) {
    double x1p;
    double y1p;
    double z1p;
    double s;
    /*  if(fabs(x2p)>1||fabs(y2p)>1||fabs(z2p)>1)return 0; */
    graphics_rot_3dvec(x2p, y2p, z2p, &x1p, &y1p, &z1p);

    if (MyGraph->PerspFlag == 0) {
        *xp = x1p;
        *yp = y1p;
        return 1;
    }
    if ((z1p >= (double)(MyGraph->ZView)) ||
        (z1p < (double)(MyGraph->ZPlane))) {
        return 0;
    }
    s = (double)(MyGraph->ZView - MyGraph->ZPlane) /
        ((double)(MyGraph->ZView) - z1p);
    x1p = s*x1p;
    y1p = s*y1p;
    *xp = x1p;
    *yp = y1p;
    return 1;
}

void
graphics_text3d(double x, double y, double z, char *s) {
    double xp;
    double yp;
    if (graphics_threedproj(x, y, z, &xp, &yp)) {
        graphics_text_abs(xp, yp, s);
    }
    return;
}

int32
graphics_threed_proj(double x, double y, double z, double *xp, double *yp) {
    double x1p;
    double y1p;
    double z1p;
    double s;
    double x2p;
    double y2p;
    double z2p;
    graphics_scale3d(x, y, z, &x2p, &y2p, &z2p); /* scale to a cube  */
    /* if(fabs(x2p)>1||fabs(y2p)>1||fabs(z2p)>1)return 0; */
    graphics_rot_3dvec(x2p, y2p, z2p, &x1p, &y1p, &z1p);

    if (MyGraph->PerspFlag == 0) {
        *xp = x1p;
        *yp = y1p;
        return 1;
    }
    if ((z1p >= (double)(MyGraph->ZView)) ||
        (z1p < (double)(MyGraph->ZPlane))) {
        return 0;
    }
    s = (double)(MyGraph->ZView - MyGraph->ZPlane) /
        ((double)(MyGraph->ZView) - z1p);
    x1p = s*x1p;
    y1p = s*y1p;
    *xp = x1p;
    *yp = y1p;
    return 1;
}

void
graphics_point_3d(double x, double y, double z) {
    double xp;
    double yp;
    if (graphics_threed_proj(x, y, z, &xp, &yp)) {
        graphics_point_abs(xp, yp);
    }
    return;
}

void
graphics_line3dn(/* unscaled version  unclipped */
                 double xs1, double ys1, double zs1, double xsp1, double ysp1,
                 double zsp1) {
    double xs;
    double ys;
    double zs;
    double xsp;
    double ysp;
    double zsp;
    graphics_rot_3dvec(xs1, ys1, zs1, &xs, &ys, &zs); /* rotate the line */
    graphics_rot_3dvec(xsp1, ysp1, zsp1, &xsp, &ysp, &zsp);
    if (MyGraph->PerspFlag) {
        graphics_pers_line(xs, ys, zs, xsp, ysp, zsp);
    } else {
        graphics_line_nabs(xs, ys, xsp, ysp);
    }
    return;
}

void
graphics_line3d(/* unscaled version     */
                double x01, double y01, double z01, double x02, double y02,
                double z02) {
    double xs;
    double ys;
    double zs;
    double xs1;
    double ys1;
    double zs1;
    double xsp;
    double ysp;
    double zsp;
    double xsp1;
    double ysp1;
    double zsp1;
    if (!graphics_clip3d(x01, y01, z01, x02, y02, z02, &xs1, &ys1, &zs1, &xsp1,
                         &ysp1, &zsp1)) {
        return;
    }
    graphics_rot_3dvec(xs1, ys1, zs1, &xs, &ys, &zs); /* rotate the line */
    graphics_rot_3dvec(xsp1, ysp1, zsp1, &xsp, &ysp, &zsp);
    if (MyGraph->PerspFlag) {
        graphics_pers_line(xs, ys, zs, xsp, ysp, zsp);
    } else {
        graphics_line_abs(xs, ys, xsp, ysp);
    }
    return;
}

void
graphics_line_3d(double x, double y, double z, double xp, double yp,
                 double zp) {
    double xs;
    double ys;
    double zs;
    double xs1;
    double ys1;
    double zs1;
    double xsp;
    double ysp;
    double zsp;
    double xsp1;
    double ysp1;
    double zsp1;
    double x01;
    double x02;
    double y01;
    double y02;
    double z01;
    double z02;
    graphics_scale3d(x, y, z, &x01, &y01, &z01); /* scale to a cube  */
    graphics_scale3d(xp, yp, zp, &x02, &y02, &z02);
    if (!graphics_clip3d(x01, y01, z01, x02, y02, z02, &xs1, &ys1, &zs1, &xsp1,
                         &ysp1, &zsp1)) {
        return;
    }
    graphics_rot_3dvec(xs1, ys1, zs1, &xs, &ys, &zs); /* rotate the line */
    graphics_rot_3dvec(xsp1, ysp1, zsp1, &xsp, &ysp, &zsp);
    if (MyGraph->PerspFlag) {
        graphics_pers_line(xs, ys, zs, xsp, ysp, zsp);
    } else {
        graphics_line_abs(xs, ys, xsp, ysp);
    }
    return;
}

void
graphics_pers_line(double x, double y, double z, double xp, double yp,
                   double zp) {
    double Zv = (double)MyGraph->ZView;
    double Zp = (double)MyGraph->ZPlane;
    double d = Zv - Zp, s;
    double eps = .005*d;

    if (((zp >= Zv) && (z >= Zv)) || ((zp < Zp) && (z < Zp))) {
        return;
    }
    if (zp > Zv) {
        s = (Zv - eps - z) / (zp - z);
        zp = Zv - eps;
        yp = y + s*(yp - y);
        xp = x + s*(xp - x);
    }
    if (z > Zv) {
        s = (Zv - eps - zp) / (z - zp);
        z = Zv - eps;
        y = yp + s*(y - yp);
        x = xp + s*(x - xp);
    }
    if (zp < Zp) {
        s = (Zp - z) / (zp - z);
        zp = Zp;
        yp = y + s*(yp - y);
        xp = x + s*(xp - x);
    }
    if (z < Zp) {
        s = (Zp - zp) / (z - zp);
        z = Zp;
        y = yp + s*(y - yp);
        x = xp + s*(x - xp);
    }
    s = d / (Zv - zp);
    xp = xp*s;
    yp = yp*s;
    s = d / (Zv - z);
    x = s*x;
    y = s*y;
    graphics_line_abs(x, y, xp, yp);
    return;
}

void
graphics_rot_3dvec(double x, double y, double z, double *xp, double *yp,
                   double *zp) {
    double vt[3];
    double vnew[3];
    vt[0] = x;
    vt[1] = y;
    vt[2] = z;

    for (int32 i = 0; i < 3; i++) {
        vnew[i] = 0.0;
        for (int32 j = 0; j < 3; j++) {
            vnew[i] = vnew[i] + MyGraph->rm[i][j]*vt[j];
        }
    }
    *xp = vnew[0];
    *yp = vnew[1];
    *zp = vnew[2];
    return;
}

void
graphics_point_abs(double x1, double y1) {
    int32 xp;
    int32 yp;

    double x_left = XMin;
    double x_right = XMax;
    double y_top = YMax;
    double y_bottom = YMin;
    if ((x1 > x_right) || (x1 < x_left) || (y1 > y_top) || (y1 < y_bottom)) {
        return;
    }
    graphics_scale_to_screen(x1, y1, &xp, &yp);
    graphics_point(xp, yp);
    return;
}

void
graphics_line_nabs(double x1_out, double y1_out, double x2_out, double y2_out) {
    int32 xp1;
    int32 yp1;
    int32 xp2;
    int32 yp2;

    graphics_scale_to_screen(x1_out, y1_out, &xp1, &yp1);
    graphics_scale_to_screen(x2_out, y2_out, &xp2, &yp2);
    graphics_line(xp1, yp1, xp2, yp2);
    return;
}

void
graphics_bead_abs(double x1, double y1) {
    int32 i1;
    int32 j1;
    double x_left = XMin;
    double x_right = XMax;
    double y_top = YMax;
    double y_bottom = YMin;
    if ((x1 > x_right) || (x1 < x_left) || (y1 > y_top) || (y1 < y_bottom)) {
        return;
    }
    graphics_scale_to_screen(x1, y1, &i1, &j1);
    graphics_bead(i1, j1);
    return;
}

void
graphics_frect_abs(double x1, double y1, double w, double h) {
    int32 i1;
    int32 i2;
    int32 j1;
    int32 j2;
    int32 ih;
    int32 iw;
    double x2 = x1 + w;
    double y2 = y1 + h;
    graphics_scale_to_screen(x1, y1, &i1, &j1);
    graphics_scale_to_screen(x2, y2, &i2, &j2);
    iw = ABS(i2 - i1);
    ih = ABS(j2 - j1);
    graphics_frect(i1, j1, iw + 1, ih + 1);
    return;
}

void
graphics_line_abs(double x1, double y1, double x2, double y2) {
    double x1_out;
    double y1_out;
    double x2_out;
    double y2_out;

    int32 xp1;
    int32 yp1;
    int32 xp2;
    int32 yp2;
    if (graphics_clip(x1, x2, y1, y2, &x1_out, &y1_out, &x2_out, &y2_out)) {
        graphics_scale_to_screen(x1_out, y1_out, &xp1, &yp1);
        graphics_scale_to_screen(x2_out, y2_out, &xp2, &yp2);
        graphics_line(xp1, yp1, xp2, yp2);
    }
    return;
}

void
graphics_text_abs(double x, double y, char *text) {
    int32 xp;
    int32 yp;
    graphics_scale_to_screen(x, y, &xp, &yp);
    graphics_put_text(xp, yp, text);
    return;
}

void
graphics_fillin_text(char *old, char *new) {
    int32 i;
    int32 l = (int32)strlen(old);
    int32 j;
    int32 m;
    int32 ans;
    char name[256];
    char c;
    char c2;
    double z;
    char val[25];
    i = 0;
    j = 0;
    while (true) {
        c = old[i];

        if (c == '\\') {
            c2 = old[i + 1];
            if (c2 != '{') {
                goto na;
            }
            if (c2 == '{') {
                m = 0;
                i = i + 2;
                while (true) {
                    c2 = old[i];
                    if (c2 == '}') {
                        name[m] = 0;
                        ans = calc_do_calc(name, &z);
                        if (ans != -1) {
                            sprintf(val, "%g", z);

                            for (usize k = 0; k < strlen(val); k++) {
                                new[j] = val[k];
                                j++;
                            }

                            break;
                        } else {
                            new[j] = '?';
                            j++;
                        }
                    } else {
                        name[m] = c2;
                        m++;
                    }
                    i++;
                    if (i >= l) { /* oops - end of string */
                        new[j] = '?';
                        new[j + 1] = 0;
                        return;
                    }
                }
            } /* ok - we have found matching and are done */
            goto nc; /* sometimes its just easier to use the !#$$# goto */
        }
    na:
        new[j] = c;
        j++;
    nc: /* normal characters */
        i++;
        if (i >= l) {
            break;
        }
    }
    new[j] = 0;
    return;
}

void
graphics_fancy_text_abs(double x, double y, char *old, int32 size) {
    int32 xp;
    int32 yp;
    char text[256];
    graphics_scale_to_screen(x, y, &xp, &yp);
    graphics_fillin_text(old, text);
    if (PltFmtFlag == PSFMT) {
        ps_special_put_text(xp, yp, text, size);
    } else if (PltFmtFlag == SVGFMT) {
        special_put_text_svg(xp, yp, text, size);
    } else {
        graphics_special_put_text_x11(xp, yp, text, size);
    }
    return;
}

int32
graphics_clip3d(double x1, double y1, double z1, double x2, double y2,
                double z2, double *x1p, double *y1p, double *z1p, double *x2p,
                double *y2p, double *z2p) {
    int32 istack, ix1 = 0, ix2 = 0, iy1 = 0, iy2 = 0, iz1 = 0, iz2 = 0,
                  iflag = 0;

    double wh;
    double wv;
    double wo;
    double xhat;
    double yhat;
    double zhat;
    double del;

    istack = 1;
    *x1p = x1;
    *y1p = y1;
    *z1p = z1;
    *x2p = x2;
    *y2p = y2;
    *z2p = z2;
    if (x1 < -1.) {
        ix1 = -1;
    }
    if (x1 > 1.) {
        ix1 = 1;
    }
    if (x2 < -1.) {
        ix2 = -1;
    }
    if (x2 > 1.) {
        ix2 = 1;
    }
    if (y1 < -1.) {
        iy1 = -1;
    }
    if (y1 > 1.) {
        iy1 = 1;
    }
    if (y2 < -1.) {
        iy2 = -1;
    }
    if (y2 > 1.) {
        iy2 = 1;
    }
    if (z1 < -1.) {
        iz1 = -1;
    }
    if (z1 > 1.) {
        iz1 = 1;
    }
    if (z2 < -1.) {
        iz2 = -1;
    }
    if (z2 > 1.) {
        iz2 = 1;
    }

    if ((ABS(ix1) + ABS(ix2) + abs(iy1) + abs(iy2) + abs(iz1) + abs(iz2)) ==
        0) {
        return 1;
    }

    /*  Both are outside the cube  */

    if ((ix1 == ix2) && (ix1 != 0)) {
        return 0;
    }
    if ((iy1 == iy2) && (iy1 != 0)) {
        return 0;
    }
    if ((iz1 == iz2) && (iz1 != 0)) {
        return 0;
    }

    if (ix1 == 0) {
        goto C2;
    }
    wv = -1;
    if (ix1 > 0) {
        wv = 1;
    }
    *x1p = wv;
    del = (wv - x2) / (x1 - x2);
    yhat = (y1 - y2)*del + y2;
    zhat = (z1 - z2)*del + z2;
    if (fabs(zhat) <= EP1 && fabs(yhat) <= EP1) {
        *y1p = yhat;
        *z1p = zhat;
        iflag = 1;
        goto C3;
    }
    istack = 0;
C2:
    if (iy1 == 0) {
        goto C22;
    }
    wh = -1;
    if (iy1 > 0) {
        wh = 1;
    }
    *y1p = wh;
    del = (wh - y2) / (y1 - y2);
    xhat = (x1 - x2)*del + x2;
    zhat = (z1 - z2)*del + z2;
    if (fabs(zhat) <= EP1 && fabs(xhat) <= EP1) {
        *x1p = xhat;
        *z1p = zhat;
        iflag = 1;
        goto C3;
    }
    istack = 0;
C22:
    if (iz1 == 0) {
        goto C3;
    }
    wo = -1;
    if (iz1 > 0) {
        wo = 1;
    }
    *z1p = wo;

    del = (wo - z2) / (z1 - z2);
    xhat = del*(x1 - x2) + x2;
    yhat = del*(y1 - y2) + y2;

    if (fabs(xhat) <= EP1 && fabs(yhat) <= EP1) {
        *x1p = xhat;
        *y1p = yhat;
        iflag = 1;
    } else {
        istack = 0;
    }
C3:
    istack += iflag;
    if ((ix2 == 0) || (istack == 0)) {
        goto C44;
    }
    wv = -1;
    if (ix2 > 0) {
        wv = 1;
    }
    *x2p = wv;
    del = (wv - x1) / (x2 - x1);
    yhat = (y2 - y1)*del + y1;
    zhat = (z2 - z1)*del + z1;
    if (fabs(yhat) <= EP1 && fabs(zhat) <= EP1) {
        *y2p = yhat;
        *z2p = zhat;
        return 1;
    }
C44:
    if (iy2 == 0 || istack == 0) {
        goto C4;
    }
    wh = -1;
    if (iy2 > 0) {
        wh = 1;
    }
    *y2p = wh;
    del = (wh - y1) / (y2 - y1);
    xhat = (x2 - x1)*del + x1;
    zhat = (z2 - z1)*del + z1;
    if (fabs(xhat) <= EP1 && fabs(zhat) <= EP1) {
        *z2p = zhat;
        *x2p = xhat;
        return 1;
    }
C4:
    if (iz2 == 0 || istack == 0) {
        return iflag;
    }
    wo = -1;
    if (iz2 > 0) {
        wo = 1;
    }
    *z2p = wo;
    del = (wo - z1) / (z2 - z1);
    xhat = (x2 - x1)*del + x1;
    yhat = (y2 - y1)*del + y1;
    if (fabs(xhat) <= EP1 && fabs(yhat) <= EP1) {
        *x2p = xhat;
        *y2p = yhat;
        return 1;
    }
    return iflag;
}

/************************************************************ *
 *  Clipping algorithm                                         *
 *   on input,                                                 *
 *           (x1,y1) and (x2,y2) are endpoints for line        *
 *           (x_left,y_bottom) and (x_right,y_top) is window*output:*value
 *is 1 for drawing, 0 for no drawing          * (x1_out,y1_out),(x2_out,y2_out)
 *are endpoints     * of clipped line                                  *
 ***************************************************************/
int32
graphics_clip(double x1, double x2, double y1, double y2, double *x1_out,
              double *y1_out, double *x2_out, double *y2_out) {
    int32 istack;
    int32 ix1;
    int32 ix2;
    int32 iy1;
    int32 iy2;
    int32 isum;
    int32 iflag;
    double wh;
    double xhat;
    double yhat;
    double wv;
    double x_left = XMin;
    double x_right = XMax;
    double y_top = YMax;
    double y_bottom = YMin;
    istack = 1;
    ix1 = ix2 = iy1 = iy2 = iflag = 0;
    *y1_out = y1;
    *y2_out = y2;
    *x1_out = x1;
    *x2_out = x2;
    if (x1 < x_left) {
        ix1 = -1;
    }
    if (x1 > x_right) {
        ix1 = 1;
    }
    if (x2 < x_left) {
        ix2 = -1;
    }
    if (x2 > x_right) {
        ix2 = 1;
    }
    if (y2 < y_bottom) {
        iy2 = -1;
    }
    if (y2 > y_top) {
        iy2 = 1;
    }
    if (y1 < y_bottom) {
        iy1 = -1;
    }
    if (y1 > y_top) {
        iy1 = 1;
    }
    isum = ABS(ix1) + ABS(ix2) + abs(iy1) + abs(iy2);
    if (isum == 0) {
        return 1; /* both inside window so plottem' */
    }

    if (((ix1 == ix2) && (ix1 != 0)) || ((iy1 == iy2) && (iy1 != 0))) {
        return 0;
    }
    if (ix1 == 0) {
        goto C2;
    }
    wv = x_left;
    if (ix1 > 0) {
        wv = x_right;
    }
    *x1_out = wv;
    yhat = (y1 - y2)*(wv - x2) / (x1 - x2) + y2;
    if ((yhat <= y_top) && (yhat >= y_bottom)) {
        *y1_out = yhat;
        iflag = 1;
        goto C3;
    }
    istack = 0;
C2:
    if (iy1 == 0) {
        goto C3;
    }
    wh = y_bottom;
    if (iy1 > 0) {
        wh = y_top;
    }
    *y1_out = wh;
    xhat = (x1 - x2)*(wh - y2) / (y1 - y2) + x2;
    if ((xhat <= x_right) && (xhat >= x_left)) {
        *x1_out = xhat;
        iflag = 1;
    } else {
        istack = 0;
    }
C3:
    istack += iflag;
    if ((ix2 == 0) || (istack == 0)) {
        goto C4;
    }
    wv = x_left;
    if (ix2 > 0) {
        wv = x_right;
    }
    *x2_out = wv;
    yhat = (y2 - y1)*(wv - x1) / (x2 - x1) + y1;
    if ((yhat <= y_top) && (yhat >= y_bottom)) {
        *y2_out = yhat;
        return 1;
    }
C4:
    if ((iy2 == 0) || (istack == 0)) {
        return iflag;
    }
    wh = y_bottom;
    if (iy2 > 0) {
        wh = y_top;
    }
    *y2_out = wh;
    xhat = (x2 - x1)*(wh - y1) / (y2 - y1) + x1;
    if ((xhat <= x_right) && (xhat >= x_left)) {
        *x2_out = xhat;
        return 1;
    }
    return iflag;
}

void
graphics_eq_symb(double *x, int32 type) {
    double dx = 6.0*(double)(MyGraph->xhi - MyGraph->xlo)*SYMSIZE;
    double dy = 6.0*(double)(MyGraph->yhi - MyGraph->ylo)*SYMSIZE;
    int32 ix = MyGraph->xv[0] - 1, iy = MyGraph->yv[0] - 1,
          iz = MyGraph->zv[0] - 1;
    if (!Xup) {
        return;
    }
    if (MyGraph->TimeFlag) {
        return;
    }
    color_set(0);
    if (MyGraph->ThreeDFlag) {
        dx = 6.0*SYMSIZE / MyGraph->dx;
        dy = 6.0*SYMSIZE / MyGraph->dy;
        graphics_line_3d((double)x[ix] + dx, (double)x[iy], (double)x[iz],
                         (double)x[ix] - dx, (double)x[iy], (double)x[iz]);
        graphics_line_3d((double)x[ix], (double)x[iy] + dy, (double)x[iz],
                         (double)x[ix], (double)x[iy] - dy, (double)x[iz]);
        return;
    }
    graphics_draw_symbol((double)x[ix], (double)x[iy], SYMSIZE, type);
    graphics_point_abs((double)x[ix], (double)x[iy]);
    return;
}

void
graphics_draw_symbol(double x, double y, double size, int32 my_symb) {
    double dx = (double)(MyGraph->xhi - MyGraph->xlo)*size;
    double dy = (double)(MyGraph->yhi - MyGraph->ylo)*size;
    static int32 sym_dir[4][48] = {
        /*          box              */
        {0, -6, -6, 1, 12, 0, 1, 0, 12, 1, -12, 0, 1, 0, -12, 3,
         0, 0,  3,  0, 0,  3, 0, 0, 3,  0, 0,   3, 0, 0, 3,   0,
         0, 3,  0,  0, 3,  0, 0, 3, 0,  0, 3,   0, 0, 3, 0,   0},

        /*          triangle         */
        {0, -6, -6, 1, 12, 0, 1, -6, 12, 1, -6, -12, 3, 0, 0, 3,
         0, 0,  3,  0, 0,  3, 0, 0,  3,  0, 0,  3,   0, 0, 3, 0,
         0, 3,  0,  0, 3,  0, 0, 3,  0,  0, 3,  0,   0, 3, 0, 0},

        /*          cross            */
        {0, -6, 0, 1, 12, 0, 0, -6, -6, 1, 0, 12, 3, 0, 0, 3,
         0, 0,  3, 0, 0,  3, 0, 0,  3,  0, 0, 3,  0, 0, 3, 0,
         0, 3,  0, 0, 3,  0, 0, 3,  0,  0, 3, 0,  0, 3, 0, 0},

        /*          circle           */
        {0,  6,  0, 1,  -1, 3, 1, -2, 2, 1, -3, 1, 1, -3, -1, 1,
         -2, -2, 1, -1, -3, 1, 1, -3, 1, 2, -2, 1, 3, -1, 1,  3,
         1,  1,  2, 2,  1,  1, 3, 3,  0, 0, 3,  0, 0, 3,  0,  0},
    };
    int32 ind = 0;
    int32 pen = 0;
    double x1 = x;
    double y1 = y;
    double x2;
    double y2;

    while (pen != 3) {
        x2 = sym_dir[my_symb][3*ind + 1]*dx + x1;
        y2 = sym_dir[my_symb][3*ind + 2]*dy + y1;
        pen = sym_dir[my_symb][3*ind];
        if (pen != 0) {
            graphics_line_abs(x1, y1, x2, y2);
        }
        x1 = x2;
        y1 = y2;
        ind++;
    }
}
