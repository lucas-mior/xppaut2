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

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define READEM 1

#define FIRSTCOLOR 30

int32 array_plot_range;
static int32 array_plot_range_count = 0;
static char array_plot_range_stem[256] = "rangearray";
static int32 array_plot_still = 1;
static int32 array_plot_tag = 0;
static struct ArrayPlot {
    Window base;
    Window wclose;
    Window wedit;
    Window wprint;
    Window wstyle;
    Window wscale;
    Window wmax;
    Window wmin;
    Window wplot;
    Window wredraw;
    Window wtime;
    Window wgif;
    Window wrange;
    Window wfit;
    int32 index0;
    int32 indexn;
    int32 alive;
    int32 nacross;
    int32 ndown;
    int32 plotdef;
    int32 height;
    int32 width;
    int32 ploth;
    int32 plotw;
    int32 nstart;
    int32 nskip;
    int32 ncskip;
    char name[20];
    double tstart;
    double tend;
    double zmin;
    double zmax;
    double dt;
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

static void array_plot_tag2(char *);
static void array_plot_gif_all(char *, int32);
static void array_plot_scale(struct ArrayPlot *ap, double *zmax, double *zmin);
static void array_plot_wborder(Window window, int32 i, struct ArrayPlot ap);
static void create_array_plot(struct ArrayPlot *ap, char *wname, char *iname);
static void array_plot_print(struct ArrayPlot *ap);
static void array_plot_draw_scale(struct ArrayPlot ap);
static void array_plot_draw(struct ArrayPlot ap);
static void array_plot_reset_axes(struct ArrayPlot ap);
static int32 array_plot_edit2(struct ArrayPlot *ap);
static void array_plot_redraw(struct ArrayPlot ap);
static void array_plot_display(Window window, struct ArrayPlot ap);
static void array_plot_button(Window window);
static void array_plot_get_root(char *s, char *sroot, int32 *num);
static void array_plot_gif(void);

void
array_plot_draw_one(char *bob) {
    char filename[300];

    array_plot_redraw(array_plot);
    if (array_plot_tag) {
        array_plot_tag2(bob);
    }
    XFlush(display);
    snprintf(filename, sizeof(filename), "%s.%d.gif", array_plot_range_stem,
             array_plot_range_count);
    array_plot_gif_all(filename, array_plot_still);
    array_plot_range_count++;
    return;
}

void
array_plot_optimize(int32 *plist) {
    int32 i0 = plist[0] - 1;
    int32 i1 = plist[1] - 1;
    int32 nr;
    int32 ns;
    double zmax;
    double zmin;
    int32 nrows = browser_my.maxrow;
    int32 ncol = i1 + 1 - i0;

    if (ncol < 2 || nrows < 2) {
        return;
    }
    array_plot_make_my("Array!");

    array_plot.index0 = i0 + 1;
    strcpy(array_plot.name, uvar_names[i0]);
    array_plot.nacross = ncol;
    nr = 201;
    if (nrows < nr) {
        nr = nrows;
    }
    array_plot.ndown = nr;
    ns = nrows / nr;
    array_plot.nskip = ns;
    array_plot.ncskip = 1;
    array_plot_scale(&array_plot, &zmax, &zmin);
    array_plot.zmin = zmin;
    array_plot.zmax = zmax;
    array_plot.plotdef = 1;
    array_plot_reset_axes(array_plot);
    array_plot_redraw(array_plot);
    return;
}

void
array_plot_make_my(char *name) {
    if (array_plot.alive == 1) {
        return;
    }
    create_array_plot(&array_plot, name, name);
    return;
}

void
array_plot_scale(struct ArrayPlot *ap, double *zmax, double *zmin) {
    int32 ib;
    int32 jb;
    int32 row0 = ap->nstart;
    int32 col0 = ap->index0;
    int32 nrows = browser_my.maxrow;
    double z;
    ib = col0;
    jb = row0;
    *zmax = browser_my.data[ib][jb];
    *zmin = *zmax;
    for (int32 i = 0; i < ap->nacross / ap->ncskip; i++) {
        ib = col0 + i*ap->ncskip;
        if (ib <= browser_my.maxcol) {
            for (int32 j = 0; j < ap->ndown; j++) {
                jb = row0 + ap->nskip*j;
                if (jb < nrows && jb >= 0) {
                    z = browser_my.data[ib][jb];
                    if (z < *zmin) {
                        *zmin = z;
                    }
                    if (z > *zmax) {
                        *zmax = z;
                    }
                }
            }
        }
    }
    if (*zmin >= *zmax) {
        *zmax = fabs(*zmin) + 1 + *zmin;
    }
    return;
}

void
array_plot_expose(Window window) {
    if (array_plot.alive) {
        array_plot_display(window, array_plot);
    }
    return;
}

void
array_plot_do_events(XEvent event) {
    int32 x;
    int32 y;
    if (array_plot.alive == 0) {
        return;
    }
    switch (event.type) {
        /* case Expose:
         case MapNotify:
           array_plot_display(ev.xany.window,array_plot);
           break;
         */
    case MotionNotify:
        if (event.xany.window == array_plot.wplot) {
            array_plot.nstart = array_plot.nstart - event.xmotion.y + first_aplot_press;
            if (array_plot.nstart < 0) {
                array_plot.nstart = 0;
            }
            array_plot_redraw(array_plot);
        }
        break;
    case ConfigureNotify:
        if (event.xconfigure.window != array_plot.base) {
            return;
        }
        x = event.xconfigure.width;
        y = event.xconfigure.height;
        array_plot.width = x;
        array_plot.height = y;
        array_plot.ploth = y - 55;
        array_plot.plotw = x - 30 - 10*dcur_xs;
        XMoveResizeWindow(display, array_plot.wplot, 20 + 10*dcur_xs, 45, (uint)array_plot.plotw,
                          (uint)array_plot.ploth);
        break;
    case EnterNotify:
        array_plot_wborder(event.xexpose.window, 2, array_plot);
        break;
    case LeaveNotify:
        array_plot_wborder(event.xexpose.window, 1, array_plot);
        break;
    case ButtonPress:
        if (event.xany.window == array_plot.wplot) {
            first_aplot_press = event.xbutton.y;
        }
        array_plot_button(event.xbutton.window);
        break;
    default:
        break;
    }
    return;
}

void
array_plot_wborder(Window window, int32 i, struct ArrayPlot ap) {
    /* if(w==ap.wedit||w==ap.wprint||w==ap.wkill||w==ap.wstyle||w==ap.wredraw)
     */
    Window w = window;
    if (w == ap.wedit || w == ap.wprint || w == ap.wclose || w == ap.wredraw || w == ap.wgif ||
        w == ap.wrange || w == ap.wfit) {
        XSetWindowBorderWidth(display, w, (uint)i);
    }
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
    width = ap->width;
    height = ap->height;
    base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0, ap->width, ap->height, 1);
    ap->base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask | StructureNotifyMask);
    XStringListToTextProperty(&wname, 1, &winname);
    XStringListToTextProperty(&iname, 1, &iconname);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.min_width = width;
    size_hints.min_height = height;
    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints, NULL, NULL);
    main_fix_window_size(base, width, height, FIX_MIN_SIZE);
    many_pops_make_icon((char *)array_bits, array_width, array_height, base);

    ap->wredraw = browser_button2(base, 0, 0, 0);
    ap->wedit = browser_button2(base, 0, 1, 0);
    ap->wprint = browser_button2(base, 0, 2, 0);
    ap->wclose = browser_button2(base, 0, 3, 0);
    ap->wfit = browser_button2(base, 0, 4, 0);
    ap->wrange = browser_button2(base, 0, 5, 0);
    ap->wgif = browser_button2(base, 1, 0, 0);
    ap->wmax = pop_list_make_window(base, 10, 45, 10*dcur_xs, dcur_ys, 1);
    ap->wmin = pop_list_make_window(base, 10, 51 + dcur_ys + color_total, 10*dcur_xs, dcur_ys, 1);
    ap->wscale =
        pop_list_make_window(base, 10 + 4*dcur_xs, 48 + dcur_ys, 2*dcur_xs, color_total, 0);
    ap->wtime = pop_list_make_window(base, 20 + 10*dcur_xs, 30, 20*dcur_xs, dcur_ys, 0);
    ap->wplot = pop_list_make_plain_window(base, 20 + 10*dcur_xs, 45, width - 30 - 10*dcur_xs,
                                           height - 55, 2);
    ap->plotw = width - 30 - 10*dcur_xs;
    ap->ploth = height - 55;
    ap->alive = 1;
    array_plot_gc = XCreateGC(display, ap->wplot, valuemask, &values);
    return;
}

void
array_plot_print(struct ArrayPlot *ap) {
    double tlo;
    double thi;
    int32 status;
    int32 errflag;
    static char *n[] = {"Filename", "Top label", "Side label", "Bottom label", "Render(-1,0,1,2)"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 nrows = browser_my.maxrow;
    int32 row0 = ap->nstart;
    int32 col0 = ap->index0;
    int32 jb;
    if (nrows <= 2) {
        return;
    }
    if (ap->plotdef == 0 || ap->nacross < 2 || ap->ndown < 2) {
        return;
    }
    jb = row0;
    tlo = 0.0;
    thi = 20.0;
    if (jb > 0 && jb < nrows) {
        tlo = browser_my.data[0][jb];
    }
    jb = row0 + ap->nskip*(ap->ndown - 1);
    if (jb >= nrows) {
        jb = nrows - 1;
    }
    if (jb >= 0) {
        thi = browser_my.data[0][jb];
    }
    strncpy(values[0], ap->filename, sizeof(values[0]));
    strncpy(values[1], ap->xtitle, sizeof(values[1]));
    strncpy(values[2], ap->ytitle, sizeof(values[2]));
    strncpy(values[3], ap->bottom, sizeof(values[3]));
    snprintf(values[4], sizeof(values[4]), "%d", ap->type);
    status = pop_list_do_string_box(5, 5, 1, "Print array_plot", n, values, 40);
    if (status != 0) {
        strcpy(ap->filename, values[0]);
        strcpy(ap->xtitle, values[1]);
        strcpy(ap->ytitle, values[2]);
        strcpy(ap->bottom, values[3]);
        ap->type = atoi(values[4]);
        if (ap->type < -1 || ap->type > 2) {
            ap->type = -1;
        }
        errflag =
            array_print(ap->filename, ap->xtitle, ap->ytitle, ap->bottom, ap->nacross, ap->ndown,
                        col0, row0, ap->nskip, ap->ncskip, nrows, browser_my.maxcol,
                        browser_my.data, ap->zmin, ap->zmax, tlo, thi, ap->type);
        if (errflag == -1) {
            ggets_err_msg("Couldn't open file");
        }
    }
    return;
}

void
array_plot_button(Window window) {
    if (window == array_plot.wedit) {
        array_plot_edit2(&array_plot);
    }
    if (window == array_plot.wfit) {
        // array plot fit
        double zmax;
        double zmin;
        array_plot_scale(&array_plot, &zmax, &zmin);
        array_plot.zmin = zmin;
        array_plot.zmax = zmax;
        array_plot_redraw(array_plot);
    }
    if (window == array_plot.wrange) {
        // array plot set up range
        static char *n[] = {"Basename", "Still(1/0)", "Tag(0/1)"};
        char values[LENGTH(n)][MAX_LEN_SBOX];
        int32 status;
        double *x;
        strncpy(values[0], array_plot_range_stem, sizeof(values[0]));
        snprintf(values[1], sizeof(values[1]), "%d", array_plot_still);
        snprintf(values[2], sizeof(values[2]), "%d", array_plot_tag);
        status = pop_list_do_string_box(3, 3, 1, "Array range saving", n, values, 28);
        if (status != 0) {
            snprintf(array_plot_range_stem, sizeof(array_plot_range_stem), "%s", values[0]);
            array_plot_still = atoi(values[1]);
            array_plot_tag = atoi(values[2]);
            array_plot_range = 1;
            array_plot_range_count = 0;
            x = &my_data[0];
            integrate_do_range(x, 0);
        }
    }
    if (window == array_plot.wredraw) {
        array_plot_redraw(array_plot);
    }
    if (window == array_plot.wprint) {
        array_plot_print(&array_plot);
    }
    if (window == array_plot.wclose) {
        // array plot destroy
        array_plot.alive = 0;
        browser_wait_a_sec(CLICK_TIME);
        XDestroySubwindows(display, array_plot.base);
        XDestroyWindow(display, array_plot.base);
    }
    if (window == array_plot.wgif) {
        array_plot_gif();
    }
    return;
}

void
array_plot_draw_scale(struct ArrayPlot ap) {
    int32 y;
    Window window = ap.wscale;
    for (int32 i = 0; i < color_total; i++) {
        y = color_total - i - 1;
        color_set(i + FIRSTCOLOR);
        XDrawLine(display, window, gc_graph, 0, y, 2*dcur_xs, y);
    }
    return;
}

void
array_plot_draw(struct ArrayPlot ap) {
    if (plot3d_auto_redraw != 1) {
        return;
    }
    array_plot_redraw(ap);
}

void
array_plot_edit(void) {
    array_plot_edit2(&array_plot);
    return;
}

void
array_plot_get_root(char *s, char *sroot, int32 *num) {
    int32 n = (int32)strlen(s);
    int32 i = n - 1;

    char me[100];
    *num = 0;
    while (true) {
        if (!isdigit(s[i])) {
            break;
        }
        i--;
        if (i < 0) {
            break;
        }
    }
    if (i < 0) {
        strcpy(sroot, s);
    } else {
        for (int32 j = 0; j <= i; j++) {
            sroot[j] = s[j];
        }
        sroot[i + 1] = 0;
    }
    if (i >= 0 && i < n) {
        for (int32 j = i + 1; j < n; j++) {
            me[j - i - 1] = s[j];
        }
        me[n - i] = 0;
        *num = atoi(me);
    }
    return;
}

void
array_plot_reset_axes(struct ArrayPlot ap) {
    char bob[200];
    char sroot[100];
    int32 num;
    if (ap.alive == 0) {
        return;
    }
    array_plot_get_root(ap.name, sroot, &num);
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
    if (f == READEM) {
        fgets(bob, 255, fp);
    } else {
        fprintf(fp, "# Array plot stuff\n");
    }
    lunch_io_string(array_plot.name, 11, fp, f);
    lunch_io_int(&array_plot.nacross, fp, f, "NCols");
    lunch_io_int(&array_plot.nstart, fp, f, "Row 1");
    lunch_io_int(&array_plot.ndown, fp, f, "NRows");
    lunch_io_int(&array_plot.nskip, fp, f, "RowSkip");
    lunch_io_double(&array_plot.zmin, fp, f, "Zmin");
    lunch_io_double(&array_plot.zmax, fp, f, "Zmax");
    return;
}

int32
array_plot_edit2(struct ArrayPlot *ap) {
    int32 i;
    int32 status;
    double zmax;
    double zmin;
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
    status = pop_list_do_string_box(9, 9, 1, "Edit array_plot", n, values, 40);
    if (status != 0) {
        browser_find_variable(values[0], &i);
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
        if (ap->ncskip < 1) {
            ap->ncskip = 1;
        }
        array_plot_reset_axes(*ap);
    }
    return 1;
}

void
array_plot_close_files(void) {
    if (array_plot_still == 0) {
        fclose(ap_fp);
    }
    return;
}

void
array_plot_gif_all(char *filename, int32 still) {
    Pixmap xi;
    int32 x;
    int32 y;
    uint32 h;
    uint32 w;
    uint32 bw;
    uint32 d;
    Window root;
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
        XCopyArea(display, array_plot.wplot, xi, array_plot_gc, 0, 0, w, h, 0, 0);

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
        XCopyArea(display, array_plot.wplot, xi, array_plot_gc, 0, 0, w, h, 0, 0);
        scrngif_screen_to_gif(xi, ap_fp);
        fclose(ap_fp);
        XFreePixmap(display, xi);
    }
    return;
}

void
array_plot_gif(void) {
    char filename[XPP_MAX_NAME + 4];
    snprintf(filename, sizeof(filename), "%s.gif", this_file);
    if (!init_conds_file_selector("GIF plot", filename, "*.gif")) {
        return;
    }
    array_plot_gif_all(filename, 1);
    return;
}

void
array_plot_redraw(struct ArrayPlot ap) {
    Window window = ap.wplot;
    double z;
    double dx;
    double dy;
    double x;
    double y;
    double tlo;
    double thi;
    char bob[100];
    int32 nrows = browser_my.maxrow, colr, cmax = FIRSTCOLOR + color_total;
    int32 row0 = ap.nstart;
    int32 col0 = ap.index0;
    int32 delx;
    int32 dely;
    int32 ib;
    int32 jb;
    int32 ix;
    int32 iy;
    if (nrows <= 2) {
        return;
    }
    if (ap.plotdef == 0 || ap.nacross < 2 || ap.ndown < 2) {
        return;
    }
    XClearWindow(display, ap.wtime);
    XClearWindow(display, window);
    jb = row0;
    tlo = 0.0;
    thi = 20.0;
    if (jb > 0 && jb < nrows) {
        tlo = browser_my.data[0][jb];
    }
    jb = row0 + ap.nskip*(ap.ndown - 1);
    if (jb >= nrows) {
        jb = nrows - 1;
    }
    if (jb >= 0) {
        thi = browser_my.data[0][jb];
    }
    snprintf(bob, sizeof(bob), " %g < t < %g ", tlo, thi);
    XDrawString(display, ap.wtime, small_gc, 0, cury_offs, bob, (int32)strlen(bob));
    dx = (double)ap.plotw / (double)(ap.nacross / ap.ncskip);
    dy = (double)ap.ploth / (double)ap.ndown;
    delx = (int32)dx + 1;
    dely = (int32)dy + 1;
    for (int32 i = 0; i < ap.nacross / ap.ncskip; i++) {
        ib = col0 + i*ap.ncskip;
        x = dx*i;
        ix = (int32)x;

        if (ib >= browser_my.maxcol) {
            return;
        }
        for (int32 j = 0; j < ap.ndown; j++) {
            jb = row0 + ap.nskip*j;

            if (jb < nrows && jb >= 0) {
                y = j*dy;
                iy = (int32)y;
                /*	  if(j==0)
                          ggets_plintf(" ib=%d ix=%d iy=%d \n",ib,ix,iy); */
                z = (double)color_total*(browser_my.data[ib][jb] - ap.zmin) / (ap.zmax - ap.zmin);
                colr = (int32)z + FIRSTCOLOR;
                if (colr < FIRSTCOLOR) {
                    colr = FIRSTCOLOR;
                }
                if (colr > cmax) {
                    colr = cmax;
                }
                if (colr < 0) {
                    XSetForeground(display, array_plot_gc, gr_back);
                }
                if (colr == 0) {
                    XSetForeground(display, array_plot_gc, gr_fore);
                } else {
                    if (COLOR) {
                        XSetForeground(display, array_plot_gc, (uint)color_map(colr));
                    } else {
                        XSetForeground(display, array_plot_gc, gr_fore);
                    }
                }
                XFillRectangle(display, window, array_plot_gc, ix, iy, (uint)delx, (uint)dely);
            }
        }
    }
    XFlush(display);
    return;
}

void
array_plot_tag2(char *bob) {
    color_set(0);
    XDrawString(display, array_plot.wplot, small_gc, 0, cury_offs, bob, (int32)strlen(bob));
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
        array_plot_draw_scale(ap);
        return;
    }
    if (window == ap.wmin) {
        snprintf(bob, sizeof(bob), "%g", ap.zmin);
        XDrawString(display, window, small_gc, 0, cury_offs, bob, (int32)strlen(bob));
        return;
    }
    if (window == ap.wmax) {
        snprintf(bob, sizeof(bob), "%g", ap.zmax);
        XDrawString(display, window, small_gc, 0, cury_offs, bob, (int32)strlen(bob));
        return;
    }
    if (window == ap.wedit) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Edit", 4);
        return;
    }
    if (window == ap.wgif) {
        XDrawString(display, window, small_gc, 0, cury_offs, "GIF", 3);
        return;
    }
    if (window == ap.wredraw) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Redraw", 6);
        return;
    }
    if (window == ap.wfit) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Fit", 3);
        return;
    }
    if (window == ap.wrange) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Range", 5);
        return;
    }
    if (window == ap.wprint) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Print", 5);
        return;
    }

    /*if(w==ap.wkill){
     XDrawString(display,w,small_gc,0,cury_offs,"Kill",4);
     return;
   }
    */

    if (window == ap.wclose) {
        XDrawString(display, window, small_gc, 0, cury_offs, "Close", 5);
        return;
    }
}
