#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/time.h>
#include <sys/time.h>
#include <unistd.h>

#include "browse.bitmap"
#include "functions.h"
#include "integers.h"
#include "mykeydef.h"
#include "parserslow.h"
#include "xpplim.h"

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define xds(a)                                                                 \
    do {                                                                       \
        XDrawString(display, window, small_gc, 5, CURY_OFFs, a, strlen(a));    \
        return;                                                                \
    } while (0)

#define BMAXCOL 20

#define MYMASK                                                                 \
    (ButtonPressMask | ButtonReleaseMask | KeyPressMask | ExposureMask |       \
     StructureNotifyMask | LeaveWindowMask | EnterWindowMask)

#define SIMPMASK                                                               \
    (ButtonPressMask | ButtonReleaseMask | KeyPressMask | ExposureMask |       \
     StructureNotifyMask)

/*  The one and only primitive data browser   */

Browser my_browser;

static double *old_rep;
static int32 REPLACE = 0, R_COL = 0;

static void write_browser_data(FILE *fp, Browser *b);
static int32 browse_check_for_stor(double **data);
static void browse_data_del_col(Browser *b);
static void browse_data_add_col(Browser *b);
static int32 browse_add_stor_col(char *name, char *formula, Browser *b);
static void browse_chk_seq(char *f, int32 *seq, double *a1, double *a2);
static void browse_replace_column(char *var, char *form, double **dat, int32 n);
static void browse_unreplace_column(void);
static void browse_make_d_table(double xlo, double xhi, int32 col, char *filename,
                         Browser b);
static void browse_find_value(int32 col, double val, int32 *row, Browser b);
static void browse_but_on(Browser *b, int32 i, Window window, int32 yn);
static void enter_browser(XEvent event, Browser *b, int32 yn);
static void display_browser(Window window, Browser b);
static void redraw_browser(Browser b);
static void browse_draw_data(Browser b);
static void kill_browser(Browser *b);
static void make_browser(Browser *b, char *wname, char *iname, int32 row,
                         int32 col);
static void expose_browser(XEvent event, Browser b);
static void resize_browser(Window win, Browser *b);
static void browse_button(XEvent event, Browser *b);
static void browse_keypress(XEvent event, int32 *used, Browser *b);
static void browse_data_up(Browser *b);
static void browse_data_down(Browser *b);
static void browse_data_pgup(Browser *b);
static void browse_data_pgdn(Browser *b);
static void browse_data_home(Browser *b);
static void browse_data_end(Browser *b);
static void browse_data_get(Browser *b);
static void browse_data_replace(Browser *b);
static void browse_data_unreplace(Browser *b);
static void browse_data_table(Browser *b);
static void browse_data_find(Browser *b);
static void browse_data_read(Browser *b);
static void browse_data_write(Browser *b);
static void browse_data_left(Browser *b);
static void browse_data_right(Browser *b);
static void browse_data_first(Browser *b);
static void browse_data_last(Browser *b);
static void browse_data_restore(Browser *b);

double **
get_browser_data(void) {
    return my_browser.data;
}

void
set_browser_data(double **data, int32 col0) {
    my_browser.data = data;
    my_browser.col0 = col0;
    return;
}

double *
browse_get_data_col(int32 c) {
    return my_browser.data[c];
}

int32
browse_get_time_now(void) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (int32)now.tv_usec;
}

void
browse_wait_a_sec(int32 msec) {
    struct timeval tim;

    double sec = (double)msec / 1000;
    double t1;
    double t2;
    gettimeofday(&tim, NULL);
    t1 = (double)tim.tv_sec + ((double)tim.tv_usec / 1000000.0);

    while (true) {
        gettimeofday(&tim, NULL);
        t2 = (double)tim.tv_sec + ((double)tim.tv_usec / 1000000.0);

        if ((t2 - t1) > sec)
            break;
    }
    return;
}

int32
get_max_row_browser(void) {
    return my_browser.maxrow;
}

void
write_my_browser_data(FILE *fp) {
    write_browser_data(fp, &my_browser);
    return;
}

void
write_browser_data(FILE *fp, Browser *b) {
    int32 i, j, l;

    for (i = b->istart; i < b->iend; i++) {
        if (N_plist > 0) {
            for (l = 0; l < N_plist; l++) {
                j = plotlist[l];
                fprintf(fp, "%.8g ", b->data[j][i]);
            }
        } else {
            for (j = 0; j < b->maxcol; j++)
                fprintf(fp, "%.8g ", b->data[j][i]);
        }
        fprintf(fp, "\n");
    }
    return;
}

int32
browse_check_for_stor(double **data) {
    if (data != storage) {
        ggets_err_msg("Only data can be in browser");
        return 0;
    } else
        return 1;
}

void
browse_data_del_col(Browser *b) {
    /*  this only works with storage  */
    Window window;
    int32 rev;
    if (browse_check_for_stor(b->data) == 0)
        return;
    XGetInputFocus(display, &window, &rev);
    ggets_err_msg("Sorry - not working very well yet...");
    return;
}

void
browse_data_add_col(Browser *b) {
    Window window;
    int32 rev;
    int32 status;
    char var[20], form[80];
    if (browse_check_for_stor(b->data) == 0)
        return;
    XGetInputFocus(display, &window, &rev);
    strcpy(var, "");
    strcpy(form, "");
    status = dialog_box_get("Add Column", "Name", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status =
            dialog_box_get("Add Column", "Formula:", form, "Add it", "Cancel", 80);
        if (status != 0)
            browse_add_stor_col(var, form, b);
    }
    return;
}

int32
browse_add_stor_col(char *name, char *formula, Browser *b) {
    int32 com[4000], i, j;

    if (add_expr(formula, com, &i)) {
        ggets_err_msg("Bad Formula .... ");
        return 0;
    }
    if ((my_ode[NEQ + FIX_VAR] = xmalloc((usize)(i + 2)*sizeof(int32))) ==
        NULL) {
        ggets_err_msg("Cant allocate formula space");
        return 0;
    }
    if ((storage[NEQ + 1] = xmalloc((usize)MAXSTOR*sizeof(double))) == NULL) {
        ggets_err_msg("Cant allocate space ....");
        free(my_ode[NEQ]);
        return 0;
    }
    if ((ode_names[NEQ] = xmalloc(80)) == NULL) {
        ggets_err_msg("Cannot allocate space ...");
        free(my_ode[NEQ]);
        free(storage[NEQ + 1]);
        return 0;
    }
    strcpy(ode_names[NEQ], formula);
    strupr(ode_names[NEQ]);
    for (j = 0; j <= i; j++)
        my_ode[NEQ + FIX_VAR][j] = com[j];
    strcpy(uvar_names[NEQ], name);
    strupr(uvar_names[NEQ]);
    for (i = 0; i < b->maxrow; i++)
        storage[NEQ + 1][i] = 0.0; /*  zero it all   */
    for (i = 0; i < b->maxrow; i++) {
        for (j = 0; j < NODE + 1; j++)
            set_ivar(j, (double)storage[j][i]);
        for (j = NODE; j < NEQ; j++)
            set_val(uvar_names[j], (double)storage[j + 1][i]);
        storage[NEQ + 1][i] = (double)evaluate(com);
    }
    add_var(uvar_names[NEQ], 0.0); /*  this could be trouble .... */
    NEQ++;
    b->maxcol = NEQ + 1;
    redraw_browser(*b);
    return 1;
}

void
browse_chk_seq(char *f, int32 *seq, double *a1, double *a2) {
    int32 i, j = -1;
    char n1[256], n2[256];
    int32 n = (int32)strlen(f);
    *seq = 0;
    *a1 = 0.0;
    *a2 = 0.0;
    for (i = 0; i < n; i++) {
        if (f[i] == ':') {
            *seq = 1;
            j = i;
        }

        if (f[i] == ';') {
            *seq = 2;
            j = i;
        }
    }
    if (j > -1) {
        for (i = 0; i < j; i++)
            n1[i] = f[i];
        n1[j] = 0;
        for (i = j + 1; i < n; i++)
            n2[i - j - 1] = f[i];
        n2[n - j - 1] = 0;
        *a1 = atof(n1);
        *a2 = atof(n2);
    }
    return;
}

void
browse_replace_column(char *var, char *form, double **dat, int32 n) {
    int32 com[200], i, j;
    int32 intflag = 0;
    int32 dif_var = -1;
    int32 seq = 0;
    double a1 = 0, a2 = 0, da = 0.0;
    double old = 0.0, dt, derv = 0.0;
    double sum = 0.0;
    if (n < 2)
        return;

    dt = NJMP*DELTA_T;
    /* first check for derivative or integral symbol */
    i = 0;
    while (i < (int32)strlen(form)) {
        if (!isspace(form[i]))
            break;
        i++;
    }
    if (form[i] == '&') {
        intflag = 1;
        form[i] = ' ';
    }
    if (form[i] == '@') {
        form[i] = ' ';
        browse_find_variable(form, &dif_var);
        if (dif_var < 0) {
            ggets_err_msg("No such variable");
            return;
        }
    }

    if (dif_var < 0)
        browse_chk_seq(form, &seq, &a1, &a2);
    if (seq == 1) {
        if (a1 == a2)
            seq = 3;
        else
            da = (a2 - a1) / ((double)(n - 1));
    }
    if (seq == 2)
        da = a2;
    if (seq == 3) {
        ggets_err_msg("Illegal sequence");
        return;
    }

    /*  first compile formula ... */

    if (dif_var < 0 && seq == 0) {
        if (add_expr(form, com, &i)) {
            NCON = NCON_START;
            NSYM = NSYM_START;
            ggets_err_msg("Illegal formula...");
            return;
        }
    }
    /* next check to see if column is known ... */

    browse_find_variable(var, &i);
    if (i < 0) {
        ggets_err_msg("No such column...");
        NCON = NCON_START;
        NSYM = NSYM_START;
        return;
    }
    R_COL = i;

    /* Okay the formula is cool so lets allocate and replace  */

    browse_wipe_rep();
    old_rep = xmalloc(sizeof(*old_rep)*(usize)n);
    REPLACE = 1;
    for (i = 0; i < n; i++) {
        old_rep[i] = dat[R_COL][i];
        if (dif_var < 0) {
            if (seq == 0) {
                for (j = 0; j < NODE + 1; j++)
                    set_ivar(j, (double)dat[j][i]);
                for (j = NODE; j < NEQ; j++)
                    set_val(uvar_names[j], (double)dat[j + 1][i]);
                if (intflag) {
                    sum += (double)evaluate(com);
                    dat[R_COL][i] = sum*dt;
                } else
                    dat[R_COL][i] = (double)evaluate(com);
            } else {
                dat[R_COL][i] = (double)(a1 + i*da);
            }
        } else {
            if (i == 0)
                derv = (dat[dif_var][1] - dat[dif_var][0]) / dt;
            if (i == (n - 1))
                derv = (dat[dif_var][i] - old) / dt;
            /* if(i>0&&i<(n-1))derv=(dat[dif_var][i+1]-old)/(2*dt); */
            if (i > 0 && i < (n - 1))
                derv = (dat[dif_var][i + 1] - dat[dif_var][i]) / dt;
            old = dat[dif_var][i];
            dat[R_COL][i] = derv;
        }
    }
    NCON = NCON_START;
    NSYM = NSYM_START;
    return;
}

void
browse_wipe_rep(void) {
    if (!REPLACE)
        return;
    free(old_rep);
    REPLACE = 0;
    return;
}

void
browse_unreplace_column(void) {
    int32 i, n = my_browser.maxrow;
    if (!REPLACE)
        return;
    for (i = 0; i < n; i++)
        my_browser.data[R_COL][i] = old_rep[i];
    browse_wipe_rep();
    return;
}

void
browse_make_d_table(double xlo, double xhi, int32 col, char *filename, Browser b) {
    int32 i, npts, ok;
    FILE *fp;
    browse_open_write_file(&fp, filename, &ok);
    if (!ok)
        return;
    npts = b.iend - b.istart;

    fprintf(fp, "%d\n", npts);
    fprintf(fp, "%g\n%g\n", xlo, xhi);
    for (i = 0; i < npts; i++)
        fprintf(fp, "%10.10g\n", b.data[col][i + b.istart]);
    fclose(fp);
    ggets_ping();
    return;
}

void
browse_find_value(int32 col, double val, int32 *row, Browser b) {
    int32 n = b.maxrow;
    int32 i;
    int32 ihot = 0;
    double err;
    double errm;
    errm = fabs(b.data[col][0] - val);
    for (i = b.row0; i < n; i++) {
        err = fabs(b.data[col][i] - val);
        if (err < errm) {
            ihot = i;
            errm = err;
        }
    }
    *row = ihot;
    return;
}

void
browse_find_variable(char *s, int32 *col) {
    *col = -1;
    if (strcasecmp("T", s) == 0) {
        *col = 0;
        return;
    }
    *col = init_conds_find_user_name(2, s);
    if (*col > -1)
        *col = *col + 1;
    return;
}

void
browse_but_on(Browser *b, int32 i, Window window, int32 yn) {
    uint32 val = 1;
    if (yn)
        val = 2;
    XSetWindowBorderWidth(display, window, val);
    if (yn && TipsFlag && i >= 0) {
        strcpy(b->hinttxt, browse_hint[i]);
        display_browser(b->hint, *b);
    }
    return;
}

void
enter_browser(XEvent event, Browser *b, int32 yn) {
    Window window = event.xexpose.window;
    if (window == b->find)
        browse_but_on(b, 0, window, yn);
    if (window == b->up)
        browse_but_on(b, 1, window, yn);
    if (window == b->down)
        browse_but_on(b, 2, window, yn);
    if (window == b->pgup)
        browse_but_on(b, 3, window, yn);
    if (window == b->pgdn)
        browse_but_on(b, 4, window, yn);
    if (window == b->left)
        browse_but_on(b, 5, window, yn);
    if (window == b->right)
        browse_but_on(b, 6, window, yn);
    if (window == b->home)
        browse_but_on(b, 7, window, yn);
    if (window == b->end)
        browse_but_on(b, 8, window, yn);
    if (window == b->first)
        browse_but_on(b, 9, window, yn);
    if (window == b->last)
        browse_but_on(b, 10, window, yn);
    if (window == b->restore)
        browse_but_on(b, 11, window, yn);
    if (window == b->write)
        browse_but_on(b, 12, window, yn);
    if (window == b->get)
        browse_but_on(b, 13, window, yn);
    if (window == b->repl)
        browse_but_on(b, 14, window, yn);
    if (window == b->unrepl)
        browse_but_on(b, 15, window, yn);
    if (window == b->table)
        browse_but_on(b, 16, window, yn);
    if (window == b->load)
        browse_but_on(b, 17, window, yn);
    if (window == b->time)
        browse_but_on(b, 18, window, yn);
    if (window == b->addcol)
        browse_but_on(b, 19, window, yn);
    if (window == b->delcol)
        browse_but_on(b, 20, window, yn);
    if (window == b->close)
        browse_but_on(b, -1, window, yn);
    return;
}

void
display_browser(Window window, Browser b) {
    int32 i0;
    if (window == b.hint) {
        XClearWindow(display, b.hint);
        XDrawString(display, window, small_gc, 8, CURY_OFFs, b.hinttxt,
                    (int)strlen(b.hinttxt));
        return;
    }

    if (window == b.find)
        xds("Find");
    if (window == b.up)
        xds("Up");
    if (window == b.down)
        xds("Down");
    if (window == b.pgup)
        xds("PgUp");
    if (window == b.pgdn)
        xds("PgDn");
    if (window == b.left)
        xds("Left");
    if (window == b.right)
        xds("Right");
    if (window == b.home)
        xds("Home");
    if (window == b.end)
        xds("End");
    if (window == b.first)
        xds("First");
    if (window == b.last)
        xds("Last");
    if (window == b.restore)
        xds("Restore");
    if (window == b.write)
        xds("Write");
    if (window == b.get)
        xds("Get");
    if (window == b.repl)
        xds("Replace");
    if (window == b.unrepl)
        xds("Unrepl");
    if (window == b.table)
        xds("Table");
    if (window == b.load)
        xds("Load");
    if (window == b.time)
        xds("Time");
    if (window == b.addcol)
        xds("Add col");
    if (window == b.close)
        xds("Close");
    if (window == b.delcol)
        xds("Del col");

    for (int32 i = 0; i < BMAXCOL; i++) {
        if (window == b.label[i]) {
            i0 = i + b.col0 - 1;
            if (i0 < b.maxcol - 1)
                XDrawString(display, window, small_gc, 5, CURY_OFFs,
                            uvar_names[i0], (int)strlen(uvar_names[i0]));
        }
    }
    if (window == b.main)
        browse_draw_data(b);
    return;
}

void
redraw_browser(Browser b) {
    int32 i;
    int32 i0;
    Window window;
    browse_draw_data(b);
    for (i = 0; i < BMAXCOL; i++) {
        window = b.label[i];
        i0 = i + b.col0 - 1;
        if (i0 < (b.maxcol - 1)) {
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, 5, CURY_OFFs, uvar_names[i0],
                        (int)strlen(uvar_names[i0]));
        }
    }
    return;
}

void
refresh_browser(int32 length) {
    my_browser.dataflag = 1;
    my_browser.maxrow = length;
    my_browser.iend = length;
    if (Xup && my_browser.xflag == 1)
        browse_draw_data(my_browser);
    return;
}

void
reset_browser(void) {
    my_browser.maxrow = 0;
    my_browser.dataflag = 0;
    return;
}

void
browse_draw_data(Browser b) {
    int32 i, i0, j, j0;
    int32 x0;
    char string[50];
    int32 dcol = DCURXs*14;
    int32 drow = (DCURYs + 6);
    if (b.dataflag == 0)
        return; /*   no data  */
    XClearWindow(display, b.main);

    /* Do time data first  */

    for (i = 0; i < b.nrow; i++) {
        i0 = i + b.row0;
        if (i0 < b.maxrow) {
            sprintf(string, "%.8g", b.data[0][i0]);
            XDrawString(display, b.main, small_gc, DCURXs / 2 + 5,
                        i*drow + DCURYs, string, (int)strlen(string));
        }
    }

    /* Do data stuff   */
    for (j = 0; j < b.ncol; j++) {
        x0 = (j + 1)*dcol + DCURXs / 2;
        j0 = j + b.col0;
        if (j0 >= b.maxcol)
            return; /* if this one is too big, they all are  */
        for (i = 0; i < b.nrow; i++) {
            i0 = i + b.row0;
            if (i0 < b.maxrow) {
                sprintf(string, "%.7g", b.data[j0][i0]);
                XDrawString(display, b.main, small_gc, x0 + 5,
                            i*drow + DCURYs, string, (int)strlen(string));
            }
        }
    }
    return;
}

void
init_browser(void) {
    my_browser.dataflag = 0;
    my_browser.data = storage;
    my_browser.maxcol = NEQ + 1;
    my_browser.maxrow = 0;
    my_browser.col0 = 1;
    my_browser.row0 = 0;
    my_browser.istart = 0;
    my_browser.iend = 0;
    strcpy(my_browser.hinttxt, "hint");
    return;
}

void
kill_browser(Browser *b) {
    b->xflag = 0;
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, b->base);
    XDestroyWindow(display, b->base);
    return;
}

void
make_new_browser(void) {
    if (my_browser.xflag == 1) {
        XRaiseWindow(display, my_browser.base);
        return;
    }
    make_browser(&my_browser, "Data Viewer", "Data", 20, 5);
    my_browser.xflag = 1;
}

Window
browse_button2(Window root, int32 row, int32 col, int32 iflag) {
    Window window;
    int32 dcol = 12*DCURXs;
    int32 drow = (DCURYs + 6);
    int32 width = 8*DCURXs;
    int32 x;
    int32 y;
    if (iflag == 1)
        dcol = 14*DCURXs;
    x = dcol*col + 4;
    y = drow*row + 4;
    window = make_window(root, x, y, width + 5, DCURYs + 1, 1);
    XSelectInput(display, window, MYMASK);
    return window;
}

Window
browse_button_data(Window root, int32 row, int32 col, char *name, int32 iflag) {
    Window window;
    int32 dcol = 12*DCURXs;
    int32 drow = (DCURYs + 6);
    int32 width = (int32)strlen(name)*DCURXs;

    int32 x;
    int32 y;
    if (iflag == 1)
        dcol = 14*DCURXs;
    x = dcol*col + 4;
    y = drow*row + 4;
    window = make_window(root, x, y, width + 5, DCURYs + 1, 1);
    XSelectInput(display, window, MYMASK);
    return window;
}

void
make_browser(Browser *b, char *wname, char *iname, int32 row, int32 col) {
    int32 i;
    int32 ncol = col;
    int32 width;
    int32 height;
    Window base;
    /* XWMHints wm_hints;
     */
    XTextProperty winname;
    XTextProperty iconname;
    XSizeHints size_hints;
    int32 dcol = DCURXs*17;
    int32 drow = (DCURYs + 6);
    int32 ystart = 8;

    if (ncol < 5)
        ncol = 5;

    height = drow*(row + 6);
    width = ncol*dcol;
    b->nrow = row;
    b->ncol = ncol;
    base =
        make_plain_window(RootWindow(display, screen), 0, 0, width, height, 4);
    b->base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask);
    XStringListToTextProperty(&wname, 1, &winname);
    XStringListToTextProperty(&iname, 1, &iconname);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    /* size_hints.width=width;
     size_hints.height=height; */
    size_hints.min_width = width - 15;
    size_hints.min_height = height;
    /* wm_hints.initial_state=IconicState;
    wm_hints.flags=StateHint;
    */
    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, &iconname, NULL, 0,
                         &size_hints, NULL, &class_hints);
    }
    many_pops_make_icon((char *)browse_bits, browse_width, browse_height, base);
    b->upper = make_window(base, 0, 0, width, ystart + drow*6, 1);
    XSetWindowBackground(display, b->upper, MyMainWinColor);
    b->main =
        make_plain_window(base, 0, ystart + drow*6, width, row*drow, 1);
    XSetWindowBackground(display, b->main, MyDrawWinColor);
    b->find = browse_button2(base, 0, 0, 0);
    b->get = browse_button2(base, 1, 0, 0);
    b->repl = browse_button2(base, 2, 0, 0);
    b->restore = browse_button2(base, 0, 1, 0);
    b->write = browse_button2(base, 1, 1, 0);
    b->load = browse_button2(base, 2, 1, 0);
    b->first = browse_button2(base, 0, 2, 0);
    b->last = browse_button2(base, 1, 2, 0);
    b->unrepl = browse_button2(base, 2, 2, 0);
    b->table = browse_button2(base, 2, 3, 0);
    b->up = browse_button2(base, 0, 3, 0);
    b->down = browse_button2(base, 1, 3, 0);
    b->pgup = browse_button2(base, 0, 4, 0);
    b->pgdn = browse_button2(base, 1, 4, 0);
    b->left = browse_button2(base, 0, 5, 0);
    b->right = browse_button2(base, 1, 5, 0);
    b->home = browse_button2(base, 0, 6, 0);
    b->end = browse_button2(base, 1, 6, 0);
    b->addcol = browse_button2(base, 2, 4, 0);
    b->delcol = browse_button2(base, 2, 5, 0);
    b->close = browse_button2(base, 2, 6, 0);
    b->time = browse_button2(base, 5, 0, 1);
    b->hint = make_window(base, 0, 4*drow, width - 17, drow - 3, 1);
    XSelectInput(display, b->time, SIMPMASK);

    for (i = 0; i < BMAXCOL; i++) {
        b->label[i] = browse_button_data(base, 5, i + 1, "1234567890", 1);
        XSelectInput(display, b->label[i], SIMPMASK);
    }
    if (noicon == 0)
        XIconifyWindow(display, base, screen);
    /*  XMapWindow(display,base);  */
    return;
}

/*   These are the global exporters ...   */

void
expose_my_browser(XEvent event) {
    if (my_browser.xflag == 0)
        return;
    expose_browser(event, my_browser);
    return;
}

void
enter_my_browser(XEvent event, int32 yn) {
    if (my_browser.xflag == 0)
        return;
    enter_browser(event, &my_browser, yn);
    return;
}

void
my_browse_button(XEvent event) {
    if (my_browser.xflag == 0)
        return;
    browse_button(event, &my_browser);
    return;
}

void
my_browse_keypress(XEvent event, int32 *used) {
    if (my_browser.xflag == 0)
        return;
    browse_keypress(event, used, &my_browser);
    return;
}

void
resize_my_browser(Window window) {
    if (my_browser.xflag == 0)
        return;
    resize_browser(window, &my_browser);
}

void
expose_browser(XEvent event, Browser b) {
    if (my_browser.xflag == 0)
        return;
    if (event.type != Expose)
        return;
    display_browser(event.xexpose.window, b);
    return;
}

void
resize_browser(Window window, Browser *b) {
    uint32 w, h, hreal;
    int32 dcol = 17*DCURXs, drow = DCURYs + 6;
    int32 i0;
    int32 newrow;
    int32 newcol;
    if (my_browser.xflag == 0)
        return;
    if (window != b->base)
        return;
    /* w=ev.xconfigure.width;
    h=ev.xconfigure.height; */
    eig_list_get_new_size(window, &w, &h);
    hreal = h;

    /* first make sure the size is is ok  and an integral
       value of the proper width and height
     */
    i0 = (int32)w / dcol;
    if ((w % (uint32)dcol) > 0)
        i0++;
    if (i0 > b->maxcol)
        i0 = b->maxcol;

    w = (uint32)(i0*dcol);
    if (i0 < 5)
        w = 5*(uint32)dcol;
    newcol = i0;
    h = hreal - 8 - 5*(uint32)drow;
    i0 = (int32)h / drow;
    if ((h % (uint32)drow) > 0)
        i0++;
    if (i0 > b->maxrow)
        i0 = b->maxrow;
    h = (uint32)(i0*drow + DCURXs / 2);
    newrow = i0;
    /*  Now resize everything   */
    if (b->ncol == newcol && b->nrow == newrow)
        return;
    b->ncol = newcol;
    b->nrow = newrow;

    XResizeWindow(display, b->base, w - 17, hreal);
    XResizeWindow(display, b->upper, w - 17, (uint)(8 + drow*3));
    XResizeWindow(display, b->main, w - 17, h);

    /* Let the browser know how many rows and columns of data  */
    return;
}

/*  if button is pressed in the browser
    then do the following  */

void
browse_button(XEvent event, Browser *b) {
    XEvent zz;
    int32 done = 1;
    Window w = event.xbutton.window;
    if (my_browser.xflag == 0)
        return;
    if (w == b->up || w == b->down || w == b->pgup || w == b->pgdn ||
        w == b->left || w == b->right) {
        done = 1;
        while (done) {
            if (w == b->up)
                browse_data_up(b);
            if (w == b->down)
                browse_data_down(b);
            if (w == b->pgup)
                browse_data_pgup(b);
            if (w == b->pgdn)
                browse_data_pgdn(b);
            if (w == b->left)
                browse_data_left(b);
            if (w == b->right)
                browse_data_right(b);
            browse_wait_a_sec(100);
            if (XPending(display) > 0) {
                XNextEvent(display, &zz);
                switch (zz.type) {
                case ButtonRelease:
                    done = 0;
                    break;
                default:
                    break;
                }
            }
        }
        return;
    }

    if (w == b->home) {
        browse_data_home(b);
        return;
    }

    if (w == b->end) {
        browse_data_end(b);
        return;
    }

    if (w == b->first) {
        browse_data_first(b);
        return;
    }

    if (w == b->last) {
        browse_data_last(b);
        return;
    }

    if (w == b->restore) {
        browse_data_restore(b);
        return;
    }

    if (w == b->write) {
        browse_data_write(b);
        return;
    }

    if (w == b->get) {
        browse_data_get(b);
        return;
    }

    if (w == b->find) {
        browse_data_find(b);
        return;
    }

    if (w == b->repl) {
        browse_data_replace(b);
        return;
    }

    if (w == b->load) {
        browse_data_read(b);
        return;
    }

    if (w == b->addcol) {
        browse_data_add_col(b);
        return;
    }

    if (w == b->delcol) {
        browse_data_del_col(b);
        return;
    }

    if (w == b->unrepl) {
        browse_data_unreplace(b);
        return;
    }

    if (w == b->table) {
        browse_data_table(b);
        return;
    }

    if (w == b->close) {
        kill_browser(b);
        return;
    }
    return;
}

void
browse_keypress(XEvent event, int32 *used, Browser *b) {
    Window w = event.xkey.window;

    char ks;
    Window w2;
    int32 rev;

    *used = 0;
    if (my_browser.xflag == 0)
        return;
    XGetInputFocus(display, &w2, &rev);

    if (w == b->main || w == b->base || w == b->upper || w2 == b->base) {
        *used = 1;

        ks = (char)ggets_get_key_press(&event);

        /*
         XLookupString(&ev,buf,maxlen,&ks,&comp);

        */

        if (ks == KEY_UP) {
            browse_data_up(b);
            return;
        }

        if (ks == KEY_DOWN) {
            browse_data_down(b);
            return;
        }

        if (ks == KEY_PGUP) {
            browse_data_pgup(b);
            return;
        }

        if (ks == KEY_PGDN) {
            browse_data_pgdn(b);
            return;
        }

        if (ks == KEY_LEFT) {
            browse_data_left(b);
            return;
        }

        if (ks == KEY_RIGHT) {
            browse_data_right(b);
            return;
        }

        if (ks == KEY_HOME) {
            browse_data_home(b);
            return;
        }

        if (ks == KEY_END) {
            browse_data_end(b);
            return;
        }

        if (ks == 's' || ks == 'S') {
            browse_data_first(b);
            return;
        }

        if (ks == 'e' || ks == 'E') {
            browse_data_last(b);
            return;
        }

        if (ks == 'r' || ks == 'R') {
            browse_data_restore(b);
            return;
        }

        if (ks == 'W' || ks == 'w') {
            browse_data_write(b);
            return;
        }

        if (ks == 'g' || ks == 'G') {
            browse_data_get(b);
            return;
        }

        if (ks == 'f' || ks == 'F') {
            browse_data_find(b);
            return;
        }

        if (ks == 'l' || ks == 'L') {
            browse_data_read(b);
            return;
        }

        if (ks == 'u' || ks == 'U') {
            browse_data_unreplace(b);
            return;
        }

        if (ks == 't' || ks == 'T') {
            browse_data_table(b);
            return;
        }

        if (ks == 'p' || ks == 'P') {
            browse_data_replace(b);
            return;
        }

        if (ks == 'a' || ks == 'A') {
            browse_data_add_col(b);
            return;
        }

        if (ks == 'd' || ks == 'D') {
            browse_data_del_col(b);
            return;
        }

        if (ks == KEY_ESC) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
            return;
        }

    } /* end of cases */
    return;
}

void
browse_data_up(Browser *b) {
    if (b->row0 > 0) {
        b->row0--;
        browse_draw_data(*b);
    }
    return;
}

void
browse_data_down(Browser *b) {
    if (b->row0 < (b->maxrow - 1)) {
        b->row0++;
        browse_draw_data(*b);
    }
    return;
}

void
browse_data_pgup(Browser *b) {
    int32 i = b->row0 - b->nrow;
    if (i > 0)
        b->row0 = i;
    else
        b->row0 = 0;
    browse_draw_data(*b);
    return;
}

void
browse_data_pgdn(Browser *b) {
    int32 i = b->row0 + b->nrow;
    if (i < (b->maxrow - 1))
        b->row0 = i;
    else
        b->row0 = b->maxrow - 1;
    browse_draw_data(*b);
    return;
}

void
browse_data_home(Browser *b) {
    b->row0 = 0;
    b->istart = 0;
    b->iend = b->maxrow;
    browse_draw_data(*b);
    return;
}

void
browse_data_end(Browser *b) {
    b->row0 = b->maxrow - 1;
    browse_draw_data(*b);
    return;
}

void
browse_get_data_xyz(double *x, double *y, double *z, int32 i1, int32 i2, int32 i3,
             int32 off) {
    int32 in = my_browser.row0 + off;
    *x = my_browser.data[i1][in];
    *y = my_browser.data[i2][in];
    *z = my_browser.data[i3][in];
    return;
}

void
data_get_my_browser(int32 row) {
    my_browser.row0 = row;
    browse_data_get(&my_browser);
    return;
}

void
browse_data_get(Browser *b) {
    int32 i, in = b->row0;
    set_ivar(0, (double)storage[0][in]);
    for (i = 0; i < NODE; i++) {
        last_ic[i] = (double)storage[i + 1][in];
        set_ivar(i + 1, last_ic[i]);
    }
    for (i = 0; i < NMarkov; i++) {
        last_ic[i + NODE] = (double)storage[i + NODE + 1][in];
        set_ivar(i + 1 + NODE + FIX_VAR, last_ic[i + NODE]);
    }
    for (i = NODE + NMarkov; i < NEQ; i++)
        set_val(uvar_names[i], storage[i + 1][in]);

    init_conds_redraw_ics();
}

void
browse_data_replace(Browser *b) {
    Window window;
    int32 rev;
    int32 status;
    char var[20], form[80];
    XGetInputFocus(display, &window, &rev);
    strcpy(var, uvar_names[0]);
    strcpy(form, uvar_names[0]);
    status = dialog_box_get("Replace", "Variable:", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status =
            dialog_box_get("Replace", "Formula:", form, "Replace", "Cancel", 80);
        if (status != 0)
            browse_replace_column(var, form, b->data, b->maxrow);
        browse_draw_data(*b);
    }

    XSetInputFocus(display, window, rev, CurrentTime);
    return;
}

void
browse_data_unreplace(Browser *b) {
    browse_unreplace_column();
    browse_draw_data(*b);
    return;
}

void
browse_data_table(Browser *b) {
    Window window;
    int32 rev;
    int32 status;

    static char *name[] = {"Variable", "Xlo", "Xhi", "File"};
    char value[LENGTH(name)][MAX_LEN_SBOX];

    double xlo = 0, xhi = 1;
    int32 col;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    snprintf(value[1], sizeof(value[1]), "0.00");
    snprintf(value[2], sizeof(value[2]), "1.00");
    snprintf(value[3], sizeof(value[3]), "%.*s.tab", (int)sizeof(value[0]) - 5,
             value[0]);

    XGetInputFocus(display, &window, &rev);
    status = do_string_box(4, 4, 1, "Tabulate", name, value, 40);
    XSetInputFocus(display, window, rev, CurrentTime);
    if (status == 0)
        return;
    xlo = atof(value[1]);
    xhi = atof(value[2]);
    browse_find_variable(value[0], &col);
    if (col >= 0)
        browse_make_d_table(xlo, xhi, col, value[3], *b);
    return;
}

void
browse_data_find(Browser *b) {
    Window window;
    int32 rev;
    int32 status;

    static char *name[] = {"*0Variable", "Value"};
    char value[LENGTH(name)][MAX_LEN_SBOX];
    int32 col, row = 0;

    double val;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    sprintf(value[1], "0.00");
    XGetInputFocus(display, &window, &rev);
    status = do_string_box(2, 2, 1, "Find Data", name, value, 40);

    XSetInputFocus(display, window, rev, CurrentTime);

    if (status == 0)
        return;
    val = (double)atof(value[1]);
    browse_find_variable(value[0], &col);
    if (col >= 0)
        browse_find_value(col, val, &row, *b);
    if (row >= 0) {
        b->row0 = row;
        browse_draw_data(*b);
    }
    return;
}

void
browse_open_write_file(FILE **fp, char *fil, int32 *ok) {
    char ans;
    *ok = 0;
    *fp = fopen(fil, "r");
    if (*fp != NULL) {
        fclose(*fp);
        ans = (char)menudrive_two_choice("Yes", "No", "File Exists! Overwrite?",
                                         "yn");
        if (ans != 'y')
            return;
    }

    *fp = fopen(fil, "w");
    if (*fp == NULL) {
        pop_list_respond_box("Ok", "Cannot open file");
        *ok = 0;
    } else
        *ok = 1;
    return;
}

void
browse_data_read(Browser *b) {
    int32 status;
    char fil[256];
    char ch;
    FILE *fp;
    int32 k;
    int32 len, count = 0, white = 1;
    double z;

    strcpy(fil, "test.dat");
    /*  XGetInputFocus(display,&w,&rev);
    status=dialog_box_get("Load","Filename:",fil,"Ok","Cancel",40);
    XSetInputFocus(display,w,rev,CurrentTime);
    */
    status = init_conds_file_selector("Load data", fil, "*.dat");
    if (status == 0)
        return;
    fp = fopen(fil, "r");
    if (fp == NULL) {
        pop_list_respond_box("Ok", "Cannot open file");
        return;
    }
    /*  Now we establish the width of the file and read it.
         If there are more columns than available we
         ignore them.

        if there are fewer rows we read whats necessary
        if there are more rows then read until we
        are done or MAX_STOR_ROW.
        This data can be plotted etc like anything else
       */

    do {
        fscanf(fp, "%c", &ch);
        if (!isspace((int32)ch) && (white)) {
            white = 0;
            ++count;
        }
        if (isspace((int32)ch) && (1 - white))
            white = 1;
    } while (ch != '\n');
    rewind(fp);
    len = 0;
    while (!feof(fp)) {
        for (k = 0; k < count; k++) {
            fscanf(fp, "%lf ", &z);
            if (k < b->maxcol)
                b->data[k][len] = z;
        }
        ++len;
        if (len >= MAXSTOR)
            break;
    }
    fclose(fp);
    refresh_browser(len);
    storind = len;
    /*  b->maxrow=len;
    browse_draw_data(*b); */
    return;
}

void
browse_data_write(Browser *b) {
    int32 status;
    char fil[256];
    FILE *fp;
    int32 i;
    int32 j;
    int32 ok;

    strcpy(fil, "test.dat");

    /*
    XGetInputFocus(display,&w,&rev);
     XSetInputFocus(display,command_pop,RevertToParent,CurrentTime);
     strcpy(fil,"test.dat");
     ggets_new_string("Write to:",fil);
    */
    /* status=dialog_box_get("Write","Filename:",fil,"Ok","Cancel",40);

       XSetInputFocus(display,w,rev,CurrentTime); */
    status = init_conds_file_selector("Write data", fil, "*.dat");
    if (status == 0)
        return;
    browse_open_write_file(&fp, fil, &ok);
    if (!ok)
        return;
    for (i = b->istart; i < b->iend; i++) {
        for (j = 0; j < b->maxcol; j++)
            fprintf(fp, "%.8g ", b->data[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return;
}

void
browse_data_left(Browser *b) {
    int32 i = b->col0;
    if (i > 1) {
        b->col0--;
        redraw_browser(*b);
    }
    return;
}

void
browse_data_right(Browser *b) {
    int32 i = b->col0 + b->ncol;
    if (i <= b->maxcol) {
        b->col0++;
        redraw_browser(*b);
    }
    return;
}

void
browse_data_first(Browser *b) {
    b->istart = b->row0;
}

void
browse_data_last(Browser *b) {
    b->iend = b->row0 + 1;
}

void
browse_data_restore(Browser *b) {
    integrate_restore(b->istart, b->iend);
}
