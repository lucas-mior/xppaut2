/*   routines for plotting arrays as functions of time

     makes a window
     of  N X M pixels
     user specifies   starting variable  x0 and ending variable xn
                      starting time  ending time
                      max var  min var

                                TITLE

                   [Kill]  [Edit]  [Print]  [Style] [Fit] [Range]
   ________________________________________________________
           1 |  |      tic marks              |  | N
            ---------------------------------------
     T0
          -                                                   MAX
          -                                                   ---
          -                                                   | |
                                                              | |
                                                              | |
                                                              | |
                                                              | |
                                                              | |
          -
          -                                                   MIN
     TN     ---------------------------------------

    and it creates a color plot

*/

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.bitmap"
#include "functions.h"
#include "integers.h"
#include "xpplim.h"

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define READEM 1

#define FIRSTCOLOR 30
#define FIX_MIN_SIZE 2

int32 array_plot_range;
static int32 array_plot_range_count = 0;
static char array_plot_range_stem[256] = "rangearray";
static int32 array_plot_still = 1, array_plot_tag = 0;
static struct ArrayPlot {
    Window base, wclose, wedit, wprint, wstyle, wscale, wmax, wmin;
    Window wplot, wredraw, wtime, wgif, wrange, wfit;
    int32 index0, indexn, alive, nacross, ndown, plotdef;
    int32 height, width, ploth, plotw;
    int32 nstart, nskip, ncskip;
    char name[20];
    double tstart, tend, zmin, zmax, dt;
    char xtitle[256];
    char ytitle[256];
    char filename[256];
    char bottom[256];
    int32 type;
} array_plot;

static int32 plot3d_auto_redraw = 0;
static FILE *ap_fp;
static GC array_plot_gc;
static int32 first_aplot_press;

static void set_acolor(int32 col);
static void tag_aplot(char *);
static void array_plot_gif_all(char *, int32);
static void scale_aplot(struct ArrayPlot *ap, double *zmax, double *zmin);
static void wborder(Window window, int32 i, struct ArrayPlot ap);
static void create_array_plot(struct ArrayPlot *ap, char *wname, char *iname);
static void print_aplot(struct ArrayPlot *ap);
static void draw_scale(struct ArrayPlot ap);
static void array_plot_draw(struct ArrayPlot ap);
static void reset_aplot_axes(struct ArrayPlot ap);
static int32 array_plot_edit2(struct ArrayPlot *ap);
static void redraw_aplot(struct ArrayPlot ap);
static void array_plot_display(Window window, struct ArrayPlot ap);
static void array_plot_destroy(void);
static void array_plot_button(Window window);
static void get_root(char *s, char *sroot, int32 *num);
static void gif_aplot(void);

void
array_plot_draw_one(char *bob) {
    char filename[300];

    redraw_aplot(array_plot);
    if (array_plot_tag)
        tag_aplot(bob);
    XFlush(display);
    snprintf(filename, sizeof(filename), "%s.%d.gif", array_plot_range_stem,
             array_plot_range_count);
    array_plot_gif_all(filename, array_plot_still);
    array_plot_range_count++;
    return;
}

static void
set_up_aplot_range(void) {
    static char *n[] = {"Basename", "Still(1/0)", "Tag(0/1)"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    double *x;
    strncpy(values[0], array_plot_range_stem, sizeof(values[0]));
    snprintf(values[1], sizeof(values[1]), "%d", array_plot_still);
    snprintf(values[2], sizeof(values[2]), "%d", array_plot_tag);
    status = do_string_box(3, 3, 1, "Array range saving", n, values, 28);
    if (status != 0) {
        snprintf(array_plot_range_stem, sizeof(array_plot_range_stem), "%s",
                 values[0]);
        array_plot_still = atoi(values[1]);
        array_plot_tag = atoi(values[2]);
        array_plot_range = 1;
        array_plot_range_count = 0;
        x = &MyData[0];
        integrate_do_range(x, 0);
    }
    return;
}

static void
fit_aplot(void) {
    double zmax;
    double zmin;
    scale_aplot(&array_plot, &zmax, &zmin);
    array_plot.zmin = zmin;
    array_plot.zmax = zmax;
    redraw_aplot(array_plot);
    return;
}

void
array_plot_optimize(int32 *plist) {
    int32 i0 = plist[0] - 1;
    int32 i1 = plist[1] - 1;
    int32 nr;
    int32 ns;
    double zmax, zmin;
    int32 nrows = my_browser.maxrow;
    int32 ncol = i1 + 1 - i0;
    if (ncol < 2 || nrows < 2)
        return;
    array_plot_make_my("Array!");

    array_plot.index0 = i0 + 1;
    strcpy(array_plot.name, uvar_names[i0]);
    array_plot.nacross = ncol;
    nr = 201;
    if (nrows < nr)
        nr = nrows;
    array_plot.ndown = nr;
    ns = nrows / nr;
    array_plot.nskip = ns;
    array_plot.ncskip = 1;
    scale_aplot(&array_plot, &zmax, &zmin);
    array_plot.zmin = zmin;
    array_plot.zmax = zmax;
    array_plot.plotdef = 1;
    reset_aplot_axes(array_plot);
    redraw_aplot(array_plot);
    return;
}

void
array_plot_make_my(char *name) {
    if (array_plot.alive == 1)
        return;
    create_array_plot(&array_plot, name, name);
    return;
}

void
scale_aplot(struct ArrayPlot *ap, double *zmax, double *zmin) {
    int32 i, j, ib, jb, row0 = ap->nstart, col0 = ap->index0;
    int32 nrows = my_browser.maxrow;
    double z;
    ib = col0;
    jb = row0;
    *zmax = my_browser.data[ib][jb];
    *zmin = *zmax;
    for (i = 0; i < ap->nacross / ap->ncskip; i++) {
        ib = col0 + i*ap->ncskip;
        if (ib <= my_browser.maxcol) {
            for (j = 0; j < ap->ndown; j++) {
                jb = row0 + ap->nskip*j;
                if (jb < nrows && jb >= 0) {
                    z = my_browser.data[ib][jb];
                    if (z < *zmin)
                        *zmin = z;
                    if (z > *zmax)
                        *zmax = z;
                }
            }
        }
    }
    if (*zmin >= *zmax)
        *zmax = fabs(*zmin) + 1 + *zmin;
    return;
}

void
array_plot_expose(Window window) {
    if (array_plot.alive)
        array_plot_display(window, array_plot);
    return;
}

void
array_plot_do_events(XEvent event) {
    int32 x;
    int32 y;
    if (array_plot.alive == 0)
        return;
    switch (event.type) {
        /* case Expose:
         case MapNotify:
           array_plot_display(ev.xany.window,array_plot);
           break;
         */
    case MotionNotify:
        if (event.xany.window == array_plot.wplot) {
            /*printf("%d\n",ev.xmotion.y-first_aplot_press); */
            array_plot.nstart =
                array_plot.nstart - event.xmotion.y + first_aplot_press;
            if (array_plot.nstart < 0)
                array_plot.nstart = 0;
            redraw_aplot(array_plot);
        }
        break;
    case ConfigureNotify:
        if (event.xconfigure.window != array_plot.base)
            return;
        x = event.xconfigure.width;
        y = event.xconfigure.height;
        array_plot.width = x;
        array_plot.height = y;
        array_plot.ploth = y - 55;
        array_plot.plotw = x - 30 - 10*DCURXs;
        XMoveResizeWindow(display, array_plot.wplot, 20 + 10*DCURXs, 45,
                          (uint)array_plot.plotw, (uint)array_plot.ploth);
        break;
    case EnterNotify:
        wborder(event.xexpose.window, 2, array_plot);
        break;
    case LeaveNotify:
        wborder(event.xexpose.window, 1, array_plot);
        break;
    case ButtonPress:
        if (event.xany.window == array_plot.wplot)
            first_aplot_press = event.xbutton.y;
        /*array_plot_button(ev.xbutton.window,array_plot);*/
        array_plot_button(event.xbutton.window);
        break;
    default:
        break;
    }
    return;
}

void
wborder(Window window, int32 i, struct ArrayPlot ap) {
    /* if(w==ap.wedit||w==ap.wprint||w==ap.wkill||w==ap.wstyle||w==ap.wredraw)
     */
    Window w = window;
    if (w == ap.wedit || w == ap.wprint || w == ap.wclose || w == ap.wredraw ||
        w == ap.wgif || w == ap.wrange || w == ap.wfit)
        XSetWindowBorderWidth(display, w, (uint)i);
    return;
}

void
array_plot_destroy(void) {
    array_plot.alive = 0;
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, array_plot.base);
    XDestroyWindow(display, array_plot.base);
    return;
}

void
array_plot_init_my(void) {
    array_plot.height = 400;
    array_plot.width = 400;
    array_plot.zmin = 0.0;
    array_plot.zmax = 1.0;
    array_plot.alive = 0;
    array_plot.plotdef = 0;
    array_plot.index0 = 1;
    array_plot.indexn = 0;
    array_plot.nacross = 1;
    array_plot.ndown = 50;
    array_plot.nstart = 0;
    array_plot.nskip = 8;
    array_plot.ncskip = 1;
    array_plot.tstart = 0.0;
    array_plot.tend = 20.0;
    strcpy(array_plot.filename, "output.ps");
    strcpy(array_plot.xtitle, "index");
    strcpy(array_plot.ytitle, "time");
    strcpy(array_plot.bottom, "");
    array_plot.type = -1;
    return;
}

void
create_array_plot(struct ArrayPlot *ap, char *wname, char *iname) {
    Window base;
    int32 width;
    int32 height;
    uint32 valuemask = 0;
    XGCValues values;
    XTextProperty winname;
    XTextProperty iconname;
    XSizeHints size_hints;
    /* init_array_plot(ap); */
    width = ap->width;
    height = ap->height;
    base = make_plain_window(RootWindow(display, screen), 0, 0, ap->width,
                             ap->height, 1);
    ap->base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask);
    XStringListToTextProperty(&wname, 1, &winname);
    XStringListToTextProperty(&iname, 1, &iconname);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.min_width = width;
    size_hints.min_height = height;
    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints,
                     NULL, NULL);
    main_fix_window_size(base, width, height, FIX_MIN_SIZE);
    many_pops_make_icon((char *)array_bits, array_width, array_height, base);

    ap->wredraw = browse_button2(base, 0, 0, 0);
    ap->wedit = browse_button2(base, 0, 1, 0);
    ap->wprint = browse_button2(base, 0, 2, 0);
    ap->wclose = browse_button2(base, 0, 3, 0);
    ap->wfit = browse_button2(base, 0, 4, 0);
    ap->wrange = browse_button2(base, 0, 5, 0);
    ap->wgif = browse_button2(base, 1, 0, 0);
    ap->wmax = make_window(base, 10, 45, 10*DCURXs, DCURYs, 1);
    ap->wmin = make_window(base, 10, 51 + DCURYs + color_total, 10*DCURXs,
                           DCURYs, 1);
    ap->wscale = make_window(base, 10 + 4*DCURXs, 48 + DCURYs, 2*DCURXs,
                             color_total, 0);
    ap->wtime = make_window(base, 20 + 10*DCURXs, 30, 20*DCURXs, DCURYs, 0);
    ap->wplot = make_plain_window(base, 20 + 10*DCURXs, 45,
                                  width - 30 - 10*DCURXs, height - 55, 2);
    ap->plotw = width - 30 - 10*DCURXs;
    ap->ploth = height - 55;
    ap->alive = 1;
    array_plot_gc = XCreateGC(display, ap->wplot, valuemask, &values);
    return;
}

void
print_aplot(struct ArrayPlot *ap) {
    double tlo;
    double thi;
    int32 status, errflag;
    static char *n[] = {"Filename", "Top label", "Side label", "Bottom label",
                        "Render(-1,0,1,2)"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 nrows = my_browser.maxrow;
    int32 row0 = ap->nstart;
    int32 col0 = ap->index0;
    int32 jb;
    if (nrows <= 2)
        return;
    if (ap->plotdef == 0 || ap->nacross < 2 || ap->ndown < 2)
        return;
    jb = row0;
    tlo = 0.0;
    thi = 20.0;
    if (jb > 0 && jb < nrows)
        tlo = my_browser.data[0][jb];
    jb = row0 + ap->nskip*(ap->ndown - 1);
    if (jb >= nrows)
        jb = nrows - 1;
    if (jb >= 0)
        thi = my_browser.data[0][jb];
    strncpy(values[0], ap->filename, sizeof(values[0]));
    strncpy(values[1], ap->xtitle, sizeof(values[1]));
    strncpy(values[2], ap->ytitle, sizeof(values[2]));
    strncpy(values[3], ap->bottom, sizeof(values[3]));
    snprintf(values[4], sizeof(values[4]), "%d", ap->type);
    status = do_string_box(5, 5, 1, "Print array_plot", n, values, 40);
    if (status != 0) {
        strcpy(ap->filename, values[0]);
        strcpy(ap->xtitle, values[1]);
        strcpy(ap->ytitle, values[2]);
        strcpy(ap->bottom, values[3]);
        ap->type = atoi(values[4]);
        if (ap->type < -1 || ap->type > 2)
            ap->type = -1;
        errflag =
            array_print(ap->filename, ap->xtitle, ap->ytitle, ap->bottom,
                        ap->nacross, ap->ndown, col0, row0, ap->nskip,
                        ap->ncskip, nrows, my_browser.maxcol, my_browser.data,
                        ap->zmin, ap->zmax, tlo, thi, ap->type);
        if (errflag == -1)
            ggets_err_msg("Couldn't open file");
    }
    return;
}

void
array_plot_button(Window window) {
    if (window == array_plot.wedit)
        array_plot_edit2(&array_plot);
    if (window == array_plot.wfit)
        fit_aplot();
    if (window == array_plot.wrange)
        set_up_aplot_range();
    if (window == array_plot.wredraw)
        redraw_aplot(array_plot);
    if (window == array_plot.wprint)
        print_aplot(&array_plot);
    if (window == array_plot.wclose)
        array_plot_destroy();
    if (window == array_plot.wgif)
        gif_aplot();
    return;
}

void
draw_scale(struct ArrayPlot ap) {
    int32 i;
    int32 y;
    Window window = ap.wscale;
    for (i = 0; i < color_total; i++) {
        y = color_total - i - 1;
        color_set(i + FIRSTCOLOR);
        XDrawLine(display, window, gc_graph, 0, y, 2*DCURXs, y);
    }
    return;
}

void
array_plot_draw(struct ArrayPlot ap) {
    if (plot3d_auto_redraw != 1)
        return;
    redraw_aplot(ap);
}

void
array_plot_edit(void) {
    array_plot_edit2(&array_plot);
    return;
}

void
get_root(char *s, char *sroot, int32 *num) {
    int32 n = (int32)strlen(s);
    int32 i = n - 1, j;

    char me[100];
    *num = 0;
    while (true) {
        if (!isdigit(s[i])) {
            break;
        }
        i--;
        if (i < 0)
            break;
    }
    if (i < 0)
        strcpy(sroot, s);
    else {
        for (j = 0; j <= i; j++)
            sroot[j] = s[j];
        sroot[i + 1] = 0;
    }
    if (i >= 0 && i < n) {
        for (j = i + 1; j < n; j++)
            me[j - i - 1] = s[j];
        me[n - i] = 0;
        *num = atoi(me);
    }
    return;
}

void
reset_aplot_axes(struct ArrayPlot ap) {
    char bob[200];
    char sroot[100];
    int32 num;
    if (ap.alive == 0)
        return;
    get_root(ap.name, sroot, &num);
    snprintf(bob, sizeof(bob), "%s%d..%d", sroot, num, num + ap.nacross - 1);
    XClearWindow(display, ap.wmax);
    XClearWindow(display, ap.wmin);
    array_plot_display(ap.wmax, ap);
    array_plot_display(ap.wmin, ap);
    many_pops_gtitle_text(bob, ap.base);
    return;
}

void
array_plot_dump(FILE *fp, int32 f) {
    char bob[256];
    if (f == READEM)
        fgets(bob, 255, fp);
    else
        fprintf(fp, "# Array plot stuff\n");
    io_string(array_plot.name, 11, fp, f);
    io_int(&array_plot.nacross, fp, f, "NCols");
    io_int(&array_plot.nstart, fp, f, "Row 1");
    io_int(&array_plot.ndown, fp, f, "NRows");
    io_int(&array_plot.nskip, fp, f, "RowSkip");
    io_double(&array_plot.zmin, fp, f, "Zmin");
    io_double(&array_plot.zmax, fp, f, "Zmax");
    return;
}

int32
array_plot_edit2(struct ArrayPlot *ap) {
    int32 i;
    int32 status;
    double zmax, zmin;
    char *n[] = {"*0Column 1", "NCols", "Row 1",         "NRows",  "RowSkip",
                 "Zmin",       "Zmax",  "Autoplot(0/1)", "ColSkip"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    snprintf(values[0], sizeof(values[0]), "%s", ap->name);
    snprintf(values[1], sizeof(values[1]), "%d", ap->nacross);
    snprintf(values[2], sizeof(values[2]), "%d", ap->nstart);
    snprintf(values[3], sizeof(values[3]), "%d", ap->ndown);
    snprintf(values[4], sizeof(values[4]), "%d", ap->nskip);
    snprintf(values[5], sizeof(values[5]), "%g", ap->zmin);
    snprintf(values[6], sizeof(values[6]), "%g", ap->zmax);
    snprintf(values[7], sizeof(values[7]), "%d", plot3d_auto_redraw);
    snprintf(values[8], sizeof(values[8]), "%d", ap->ncskip);
    status = do_string_box(9, 9, 1, "Edit array_plot", n, values, 40);
    if (status != 0) {
        browse_find_variable(values[0], &i);
        if (i > -1) {
            ap->index0 = i;
            strcpy(ap->name, values[0]);
        } else {
            ggets_err_msg("No such columns");
            ap->plotdef = 0;
            return 0;
        }
        zmax = atof(values[6]);
        zmin = atof(values[5]);
        if (zmin < zmax) {
            ap->zmin = zmin;
            ap->zmax = zmax;
        }
        ap->nacross = atoi(values[1]);
        ap->nstart = atoi(values[2]);
        ap->ndown = atoi(values[3]);
        ap->nskip = atoi(values[4]);
        plot3d_auto_redraw = atoi(values[7]);
        ap->plotdef = 1;
        ap->ncskip = atoi(values[8]);
        if (ap->ncskip < 1)
            ap->ncskip = 1;
        reset_aplot_axes(*ap);
    }
    return 1;
}

void
array_plot_close_files(void) {
    if (array_plot_still == 0)
        fclose(ap_fp);
    return;
}

void
array_plot_gif_all(char *filename, int32 still) {
    Pixmap xi;
    int32 x;
    int32 y;
    uint32 h, w, bw, d;
    Window root;
    /* FILE *fp; */
    if (still == 0) {
        if (array_plot_range_count == 0) {
            if ((ap_fp = fopen(filename, "w")) == NULL) {
                ggets_err_msg("Cannot open file ");
                return;
            }
        }
        XGetGeometry(display, array_plot.wplot, &root, &x, &y, &w, &h, &bw, &d);
        xi = XCreatePixmap(display, RootWindow(display, screen), w, h,
                           (uint)DefaultDepth(display, screen));
        XCopyArea(display, array_plot.wplot, xi, array_plot_gc, 0, 0, w, h, 0,
                  0);

        scrngif_add_ani_gif(xi, ap_fp, array_plot_range_count);
        XFreePixmap(display, xi);
        return;
    }

    if (still == 1) {
        if ((ap_fp = fopen(filename, "w")) == NULL) {
            ggets_err_msg("Cannot open file ");
            return;
        }

        XGetGeometry(display, array_plot.wplot, &root, &x, &y, &w, &h, &bw, &d);
        xi = XCreatePixmap(display, RootWindow(display, screen), w, h,
                           (uint)DefaultDepth(display, screen));
        XCopyArea(display, array_plot.wplot, xi, array_plot_gc, 0, 0, w, h, 0,
                  0);
        scrngif_screen_to_gif(xi, ap_fp);
        fclose(ap_fp);
        XFreePixmap(display, xi);
    }
    return;
}

void
gif_aplot(void) {
    char filename[XPP_MAX_NAME + 4];
    snprintf(filename, sizeof(filename), "%s.gif", this_file);
    if (!init_conds_file_selector("GIF plot", filename, "*.gif"))
        return;
    array_plot_gif_all(filename, 1);
    return;
}

void
redraw_aplot(struct ArrayPlot ap) {
    int32 i;
    int32 j;
    Window window = ap.wplot;
    double z, dx, dy, x, y, tlo, thi;
    char bob[100];
    int32 nrows = my_browser.maxrow, colr, cmax = FIRSTCOLOR + color_total;
    int32 row0 = ap.nstart;
    int32 col0 = ap.index0, delx, dely;
    int32 ib, jb, ix, iy;
    if (nrows <= 2)
        return;
    if (ap.plotdef == 0 || ap.nacross < 2 || ap.ndown < 2)
        return;
    XClearWindow(display, ap.wtime);
    XClearWindow(display, window);
    jb = row0;
    tlo = 0.0;
    thi = 20.0;
    if (jb > 0 && jb < nrows)
        tlo = my_browser.data[0][jb];
    jb = row0 + ap.nskip*(ap.ndown - 1);
    if (jb >= nrows)
        jb = nrows - 1;
    if (jb >= 0)
        thi = my_browser.data[0][jb];
    snprintf(bob, sizeof(bob), " %g < t < %g ", tlo, thi);
    XDrawString(display, ap.wtime, small_gc, 0, CURY_OFFs, bob,
                (int32)strlen(bob));
    dx = (double)ap.plotw / (double)(ap.nacross / ap.ncskip);
    dy = (double)ap.ploth / (double)ap.ndown;
    delx = (int32)dx + 1;
    dely = (int32)dy + 1;
    for (i = 0; i < ap.nacross / ap.ncskip; i++) {
        ib = col0 + i*ap.ncskip;
        x = dx*i;
        ix = (int32)x;

        if (ib >= my_browser.maxcol)
            return;
        for (j = 0; j < ap.ndown; j++) {
            jb = row0 + ap.nskip*j;

            if (jb < nrows && jb >= 0) {
                y = j*dy;
                iy = (int32)y;
                /*	  if(j==0)
                          ggets_plintf(" ib=%d ix=%d iy=%d \n",ib,ix,iy); */
                z = (double)color_total*(my_browser.data[ib][jb] - ap.zmin) /
                    (ap.zmax - ap.zmin);
                colr = (int32)z + FIRSTCOLOR;
                if (colr < FIRSTCOLOR)
                    colr = FIRSTCOLOR;
                if (colr > cmax)
                    colr = cmax;
                set_acolor(colr);
                XFillRectangle(display, window, array_plot_gc, ix, iy,
                               (uint)delx, (uint)dely);
            }
        }
    }
    XFlush(display);
    return;
}

void
tag_aplot(char *bob) {
    color_set(0);
    XDrawString(display, array_plot.wplot, small_gc, 0, CURY_OFFs, bob,
                (int32)strlen(bob));
    return;
}

void
set_acolor(int32 col) {
    if (col < 0)
        XSetForeground(display, array_plot_gc, GrBack);
    if (col == 0)
        XSetForeground(display, array_plot_gc, GrFore);
    else {
        if (COLOR)
            XSetForeground(display, array_plot_gc, (uint)color_map(col));
        else
            XSetForeground(display, array_plot_gc, GrFore);
    }
    return;
}

void
array_plot_display(Window window, struct ArrayPlot ap) {
    char bob[200];

    if (window == ap.wplot) {
        array_plot_draw(ap);
        return;
    }
    if (window == ap.wscale) {
        draw_scale(ap);
        return;
    }
    if (window == ap.wmin) {
        snprintf(bob, sizeof(bob), "%g", ap.zmin);
        XDrawString(display, window, small_gc, 0, CURY_OFFs, bob,
                    (int32)strlen(bob));
        return;
    }
    if (window == ap.wmax) {
        snprintf(bob, sizeof(bob), "%g", ap.zmax);
        XDrawString(display, window, small_gc, 0, CURY_OFFs, bob,
                    (int32)strlen(bob));
        return;
    }
    if (window == ap.wedit) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Edit", 4);
        return;
    }
    if (window == ap.wgif) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "GIF", 3);
        return;
    }
    if (window == ap.wredraw) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Redraw", 6);
        return;
    }
    if (window == ap.wfit) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Fit", 3);
        return;
    }
    if (window == ap.wrange) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Range", 5);
        return;
    }
    if (window == ap.wprint) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Print", 5);
        return;
    }

    /*if(w==ap.wkill){
     XDrawString(display,w,small_gc,0,CURY_OFFs,"Kill",4);
     return;
   }
    */

    if (window == ap.wclose) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Close", 5);
        return;
    }
}
