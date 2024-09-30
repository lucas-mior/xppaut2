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

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define XDS(a)                                                                 \
    do {                                                                       \
        XDrawString(display, window, small_gc, 5, cury_offs, a, strlen(a));    \
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

Browser browser_my;

static double *old_rep;
static int32 REPLACE = 0;
static int32 R_COL = 0;

static void write_browser_data(FILE *fp, Browser *b);
static int32 browser_check_for_stor(double **data);
static void browser_data_del_col(Browser *b);
static void browser_data_add_col(Browser *b);
static int32 browser_add_stor_col(char *name, char *formula, Browser *b);
static void browser_chk_seq(char *f, int32 *seq, double *a1, double *a2);
static void browser_replace_column(char *var, char *form, double **dat,
                                   int32 n);
static void browser_make_d_table(double xlo, double xhi, int32 col,
                                 char *filename, Browser b);
static void browser_find_value(int32 col, double val, int32 *row, Browser b);
static void browser_but_on(Browser *b, int32 i, Window window, int32 yn);
static void enter_browser(XEvent event, Browser *b, int32 yn);
static void display_browser(Window window, Browser b);
static void redraw_browser(Browser b);
static void browser_draw_data(Browser b);
static void kill_browser(Browser *b);
static void make_browser(Browser *b, char *wname, char *iname, int32 row,
                         int32 col);
static void expose_browser(XEvent event, Browser b);
static void resize_browser(Window win, Browser *b);
static void browser_button(XEvent event, Browser *b);
static void browser_keypress(XEvent event, int32 *used, Browser *b);
static void browser_data_up(Browser *b);
static void browser_data_down(Browser *b);
static void browser_data_pgup(Browser *b);
static void browser_data_pgdn(Browser *b);
static void browser_data_home(Browser *b);
static void browser_data_end(Browser *b);
static void browser_data_get(Browser *b);
static void browser_data_replace(Browser *b);
static void browser_data_unreplace(Browser *b);
static void browser_data_table(Browser *b);
static void browser_data_find(Browser *b);
static void browser_data_read(Browser *b);
static void browser_data_write(Browser *b);
static void browser_data_left(Browser *b);
static void browser_data_right(Browser *b);
static void browser_data_first(Browser *b);
static void browser_data_last(Browser *b);
static void browser_data_restore(Browser *b);

double **
browser_get_data(void) {
    return browser_my.data;
}

void
browser_set_data(double **data, int32 col0) {
    browser_my.data = data;
    browser_my.col0 = col0;
    return;
}

double *
browser_get_data_col(int32 c) {
    return browser_my.data[c];
}

int32
browser_get_time_now(void) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (int32)now.tv_usec;
}

void
browser_wait_a_sec(int32 msec) {
    struct timeval tim;

    double sec = (double)msec / 1000;
    double t1;
    double t2;
    gettimeofday(&tim, NULL);
    t1 = (double)tim.tv_sec + ((double)tim.tv_usec / 1000000.0);

    while (true) {
        gettimeofday(&tim, NULL);
        t2 = (double)tim.tv_sec + ((double)tim.tv_usec / 1000000.0);

        if ((t2 - t1) > sec) {
            break;
        }
    }
    return;
}

int32
browser_get_max_row(void) {
    return browser_my.maxrow;
}

void
browser_my_write_data(FILE *fp) {
    write_browser_data(fp, &browser_my);
    return;
}

void
write_browser_data(FILE *fp, Browser *b) {
    int32 i;
    int32 j;
    int32 l;

    for (i = b->istart; i < b->iend; i++) {
        if (N_plist > 0) {
            for (l = 0; l < N_plist; l++) {
                j = plotlist[l];
                fprintf(fp, "%.8g ", b->data[j][i]);
            }
        } else {
            for (j = 0; j < b->maxcol; j++) {
                fprintf(fp, "%.8g ", b->data[j][i]);
            }
        }
        fprintf(fp, "\n");
    }
    return;
}

int32
browser_check_for_stor(double **data) {
    if (data != storage) {
        ggets_err_msg("Only data can be in browser");
        return 0;
    } else {
        return 1;
    }
}

void
browser_data_del_col(Browser *b) {
    //  this only works with storage
    Window window;
    int32 rev;

    if (browser_check_for_stor(b->data) == 0) {
        return;
    }
    XGetInputFocus(display, &window, &rev);
    ggets_err_msg("Sorry - not working very well yet...");
    return;
}

void
browser_data_add_col(Browser *b) {
    Window window;
    int32 rev;
    int32 status;
    char var[20];
    char form[80];
    if (browser_check_for_stor(b->data) == 0) {
        return;
    }
    XGetInputFocus(display, &window, &rev);
    strcpy(var, "");
    strcpy(form, "");
    status = dialog_box_get("Add Column", "Name", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status = dialog_box_get("Add Column", "Formula:", form, "Add it",
                                "Cancel", 80);
        if (status != 0) {
            browser_add_stor_col(var, form, b);
        }
    }
    return;
}

int32
browser_add_stor_col(char *name, char *formula, Browser *b) {
    int32 com[4000];
    int32 i;

    if (parserslow_add_expr(formula, com, &i)) {
        ggets_err_msg("Bad Formula .... ");
        return 0;
    }
    if ((my_ode[NEQ + fix_var] = xmalloc((usize)(i + 2)*sizeof(int32))) ==
        NULL) {
        ggets_err_msg("Cant allocate formula space");
        return 0;
    }
    if ((storage[NEQ + 1] = xmalloc((usize)max_stor*sizeof(double))) == NULL) {
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
    for (int32 j = 0; j <= i; j++) {
        my_ode[NEQ + fix_var][j] = com[j];
    }
    strcpy(uvar_names[NEQ], name);
    strupr(uvar_names[NEQ]);
    for (i = 0; i < b->maxrow; i++) {
        storage[NEQ + 1][i] = 0.0;  //  zero it all
    }
    for (i = 0; i < b->maxrow; i++) {
        for (int32 j = 0; j < NODE + 1; j++) {
            set_ivar(j, (double)storage[j][i]);
        }
        for (int32 j = NODE; j < NEQ; j++) {
            set_val(uvar_names[j], (double)storage[j + 1][i]);
        }
        storage[NEQ + 1][i] = (double)evaluate(com);
    }
    parserslow_add_var(uvar_names[NEQ], 0.0);  //  this could be trouble ....
    NEQ++;
    b->maxcol = NEQ + 1;
    redraw_browser(*b);
    return 1;
}

void
browser_chk_seq(char *f, int32 *seq, double *a1, double *a2) {
    int32 j = -1;
    char n1[256];
    char n2[256];
    int32 n = (int32)strlen(f);
    *seq = 0;
    *a1 = 0.0;
    *a2 = 0.0;
    for (int32 i = 0; i < n; i++) {
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
        for (int32 i = 0; i < j; i++) {
            n1[i] = f[i];
        }
        n1[j] = 0;
        for (int32 i = j + 1; i < n; i++) {
            n2[i - j - 1] = f[i];
        }
        n2[n - j - 1] = 0;
        *a1 = atof(n1);
        *a2 = atof(n2);
    }
    return;
}

void
browser_replace_column(char *var, char *form, double **dat, int32 n) {
    int32 com[200];
    int32 i;
    int32 intflag = 0;
    int32 dif_var = -1;
    int32 seq = 0;
    double a1 = 0;
    double a2 = 0;
    double da = 0.0;
    double old = 0.0;
    double dt;
    double derv = 0.0;
    double sum = 0.0;
    if (n < 2) {
        return;
    }

    dt = NJMP*delta_t;
    // first check for derivative or integral symbol
    i = 0;
    while (i < (int32)strlen(form)) {
        if (!isspace(form[i])) {
            break;
        }
        i++;
    }
    if (form[i] == '&') {
        intflag = 1;
        form[i] = ' ';
    }
    if (form[i] == '@') {
        form[i] = ' ';
        browser_find_variable(form, &dif_var);
        if (dif_var < 0) {
            ggets_err_msg("No such variable");
            return;
        }
    }

    if (dif_var < 0) {
        browser_chk_seq(form, &seq, &a1, &a2);
    }
    if (seq == 1) {
        if (a1 == a2) {
            seq = 3;
        } else {
            da = (a2 - a1) / ((double)(n - 1));
        }
    }
    if (seq == 2) {
        da = a2;
    }
    if (seq == 3) {
        ggets_err_msg("Illegal sequence");
        return;
    }

    //  first compile formula ...

    if (dif_var < 0 && seq == 0) {
        if (parserslow_add_expr(form, com, &i)) {
            NCON = NCON_START;
            NSYM = NSYM_START;
            ggets_err_msg("Illegal formula...");
            return;
        }
    }
    // next check to see if column is known ...

    browser_find_variable(var, &i);
    if (i < 0) {
        ggets_err_msg("No such column...");
        NCON = NCON_START;
        NSYM = NSYM_START;
        return;
    }
    R_COL = i;

    // Okay the formula is cool so lets allocate and replace

    browser_wipe_rep();
    old_rep = xmalloc(sizeof(*old_rep)*(usize)n);
    REPLACE = 1;
    for (i = 0; i < n; i++) {
        old_rep[i] = dat[R_COL][i];
        if (dif_var < 0) {
            if (seq == 0) {
                for (int32 j = 0; j < NODE + 1; j++) {
                    set_ivar(j, (double)dat[j][i]);
                }
                for (int32 j = NODE; j < NEQ; j++) {
                    set_val(uvar_names[j], (double)dat[j + 1][i]);
                }
                if (intflag) {
                    sum += (double)evaluate(com);
                    dat[R_COL][i] = sum*dt;
                } else {
                    dat[R_COL][i] = (double)evaluate(com);
                }
            } else {
                dat[R_COL][i] = (double)(a1 + i*da);
            }
        } else {
            if (i == 0) {
                derv = (dat[dif_var][1] - dat[dif_var][0]) / dt;
            }
            if (i == (n - 1)) {
                derv = (dat[dif_var][i] - old) / dt;
            }
            if (i > 0 && i < (n - 1)) {
                derv = (dat[dif_var][i + 1] - dat[dif_var][i]) / dt;
            }
            old = dat[dif_var][i];
            dat[R_COL][i] = derv;
        }
    }
    NCON = NCON_START;
    NSYM = NSYM_START;
    return;
}

void
browser_wipe_rep(void) {
    if (!REPLACE) {
        return;
    }
    free(old_rep);
    REPLACE = 0;
    return;
}

void
browser_make_d_table(double xlo, double xhi, int32 col, char *filename,
                     Browser b) {
    int32 npts;
    int32 ok;
    FILE *fp;
    browser_open_write_file(&fp, filename, &ok);
    if (!ok) {
        return;
    }
    npts = b.iend - b.istart;

    fprintf(fp, "%d\n", npts);
    fprintf(fp, "%g\n%g\n", xlo, xhi);
    for (int32 i = 0; i < npts; i++) {
        fprintf(fp, "%10.10g\n", b.data[col][i + b.istart]);
    }
    fclose(fp);
    ggets_ping();
    return;
}

void
browser_find_value(int32 col, double val, int32 *row, Browser b) {
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
browser_find_variable(char *s, int32 *col) {
    *col = -1;
    if (strcasecmp("T", s) == 0) {
        *col = 0;
        return;
    }
    *col = init_conds_find_user_name(2, s);
    if (*col > -1) {
        *col = *col + 1;
    }
    return;
}

void
browser_but_on(Browser *b, int32 i, Window window, int32 yn) {
    uint32 val = 1;
    if (yn) {
        val = 2;
    }
    XSetWindowBorderWidth(display, window, val);
    if (yn && flag_tips && i >= 0) {
        strcpy(b->hinttxt, browse_hint[i]);
        display_browser(b->hint, *b);
    }
    return;
}

void
enter_browser(XEvent event, Browser *b, int32 yn) {
    Window window = event.xexpose.window;
    if (window == b->find) {
        browser_but_on(b, 0, window, yn);
    }
    if (window == b->up) {
        browser_but_on(b, 1, window, yn);
    }
    if (window == b->down) {
        browser_but_on(b, 2, window, yn);
    }
    if (window == b->pgup) {
        browser_but_on(b, 3, window, yn);
    }
    if (window == b->pgdn) {
        browser_but_on(b, 4, window, yn);
    }
    if (window == b->left) {
        browser_but_on(b, 5, window, yn);
    }
    if (window == b->right) {
        browser_but_on(b, 6, window, yn);
    }
    if (window == b->home) {
        browser_but_on(b, 7, window, yn);
    }
    if (window == b->end) {
        browser_but_on(b, 8, window, yn);
    }
    if (window == b->first) {
        browser_but_on(b, 9, window, yn);
    }
    if (window == b->last) {
        browser_but_on(b, 10, window, yn);
    }
    if (window == b->restore) {
        browser_but_on(b, 11, window, yn);
    }
    if (window == b->write) {
        browser_but_on(b, 12, window, yn);
    }
    if (window == b->get) {
        browser_but_on(b, 13, window, yn);
    }
    if (window == b->repl) {
        browser_but_on(b, 14, window, yn);
    }
    if (window == b->unrepl) {
        browser_but_on(b, 15, window, yn);
    }
    if (window == b->table) {
        browser_but_on(b, 16, window, yn);
    }
    if (window == b->load) {
        browser_but_on(b, 17, window, yn);
    }
    if (window == b->time) {
        browser_but_on(b, 18, window, yn);
    }
    if (window == b->addcol) {
        browser_but_on(b, 19, window, yn);
    }
    if (window == b->delcol) {
        browser_but_on(b, 20, window, yn);
    }
    if (window == b->close) {
        browser_but_on(b, -1, window, yn);
    }
    return;
}

void
display_browser(Window window, Browser b) {
    int32 i0;
    if (window == b.hint) {
        XClearWindow(display, b.hint);
        XDrawString(display, window, small_gc, 8, cury_offs, b.hinttxt,
                    (int)strlen(b.hinttxt));
        return;
    }

    if (window == b.find) {
        XDS("Find");
    }
    if (window == b.up) {
        XDS("Up");
    }
    if (window == b.down) {
        XDS("Down");
    }
    if (window == b.pgup) {
        XDS("PgUp");
    }
    if (window == b.pgdn) {
        XDS("PgDn");
    }
    if (window == b.left) {
        XDS("Left");
    }
    if (window == b.right) {
        XDS("Right");
    }
    if (window == b.home) {
        XDS("Home");
    }
    if (window == b.end) {
        XDS("End");
    }
    if (window == b.first) {
        XDS("First");
    }
    if (window == b.last) {
        XDS("Last");
    }
    if (window == b.restore) {
        XDS("Restore");
    }
    if (window == b.write) {
        XDS("Write");
    }
    if (window == b.get) {
        XDS("Get");
    }
    if (window == b.repl) {
        XDS("Replace");
    }
    if (window == b.unrepl) {
        XDS("Unrepl");
    }
    if (window == b.table) {
        XDS("Table");
    }
    if (window == b.load) {
        XDS("Load");
    }
    if (window == b.time) {
        XDS("Time");
    }
    if (window == b.addcol) {
        XDS("Add col");
    }
    if (window == b.close) {
        XDS("Close");
    }
    if (window == b.delcol) {
        XDS("Del col");
    }

    for (int32 i = 0; i < BMAXCOL; i++) {
        if (window == b.label[i]) {
            i0 = i + b.col0 - 1;
            if (i0 < b.maxcol - 1) {
                XDrawString(display, window, small_gc, 5, cury_offs,
                            uvar_names[i0], (int)strlen(uvar_names[i0]));
            }
        }
    }
    if (window == b.main) {
        browser_draw_data(b);
    }
    return;
}

void
redraw_browser(Browser b) {
    int32 i0;
    Window window;
    browser_draw_data(b);
    for (int32 i = 0; i < BMAXCOL; i++) {
        window = b.label[i];
        i0 = i + b.col0 - 1;
        if (i0 < (b.maxcol - 1)) {
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, 5, cury_offs, uvar_names[i0],
                        (int)strlen(uvar_names[i0]));
        }
    }
    return;
}

void
browser_refresh(int32 length) {
    browser_my.dataflag = 1;
    browser_my.maxrow = length;
    browser_my.iend = length;
    if (Xup && browser_my.xflag == 1) {
        browser_draw_data(browser_my);
    }
    return;
}

void
browser_reset(void) {
    browser_my.maxrow = 0;
    browser_my.dataflag = 0;
    return;
}

void
browser_draw_data(Browser b) {
    int32 i0;
    int32 j0;
    int32 x0;
    char string[50];
    int32 dcol = dcur_xs*14;
    int32 drow = (dcur_ys + 6);
    if (b.dataflag == 0) {
        return;  //   no data
    }
    XClearWindow(display, b.main);

    // Do time data first

    for (int32 i = 0; i < b.nrow; i++) {
        i0 = i + b.row0;
        if (i0 < b.maxrow) {
            sprintf(string, "%.8g", b.data[0][i0]);
            XDrawString(display, b.main, small_gc, dcur_xs / 2 + 5,
                        i*drow + dcur_ys, string, (int)strlen(string));
        }
    }

    // Do data stuff
    for (int32 j = 0; j < b.ncol; j++) {
        x0 = (j + 1)*dcol + dcur_xs / 2;
        j0 = j + b.col0;
        if (j0 >= b.maxcol) {
            return;  // if this one is too big, they all are
        }
        for (int32 i = 0; i < b.nrow; i++) {
            i0 = i + b.row0;
            if (i0 < b.maxrow) {
                sprintf(string, "%.7g", b.data[j0][i0]);
                XDrawString(display, b.main, small_gc, x0 + 5,
                            i*drow + dcur_ys, string, (int)strlen(string));
            }
        }
    }
    return;
}

void
browser_init(void) {
    browser_my.dataflag = 0;
    browser_my.data = storage;
    browser_my.maxcol = NEQ + 1;
    browser_my.maxrow = 0;
    browser_my.col0 = 1;
    browser_my.row0 = 0;
    browser_my.istart = 0;
    browser_my.iend = 0;
    strcpy(browser_my.hinttxt, "hint");
    return;
}

void
kill_browser(Browser *b) {
    b->xflag = 0;
    browser_wait_a_sec(CLICK_TIME);
    XDestroySubwindows(display, b->base);
    XDestroyWindow(display, b->base);
    return;
}

void
browser_make_new(void) {
    if (browser_my.xflag == 1) {
        XRaiseWindow(display, browser_my.base);
        return;
    }
    make_browser(&browser_my, "Data Viewer", "Data", 20, 5);
    browser_my.xflag = 1;
}

Window
browser_button2(Window root, int32 row, int32 col, int32 iflag) {
    Window window;
    int32 dcol = 12*dcur_xs;
    int32 drow = (dcur_ys + 6);
    int32 width = 8*dcur_xs;
    int32 x;
    int32 y;
    if (iflag == 1) {
        dcol = 14*dcur_xs;
    }
    x = dcol*col + 4;
    y = drow*row + 4;
    window = pop_list_make_window(root, x, y, width + 5, dcur_ys + 1, 1);
    XSelectInput(display, window, MYMASK);
    return window;
}

Window
browser_button_data(Window root, int32 row, int32 col, char *name,
                    int32 iflag) {
    Window window;
    int32 dcol = 12*dcur_xs;
    int32 drow = (dcur_ys + 6);
    int32 width = (int32)strlen(name)*dcur_xs;

    int32 x;
    int32 y;
    if (iflag == 1) {
        dcol = 14*dcur_xs;
    }
    x = dcol*col + 4;
    y = drow*row + 4;
    window = pop_list_make_window(root, x, y, width + 5, dcur_ys + 1, 1);
    XSelectInput(display, window, MYMASK);
    return window;
}

void
make_browser(Browser *b, char *wname, char *iname, int32 row, int32 col) {
    int32 ncol = col;
    int32 width;
    int32 height;
    Window base;
    XTextProperty winname;
    XTextProperty iconname;
    XSizeHints size_hints;
    int32 dcol = dcur_xs*17;
    int32 drow = (dcur_ys + 6);
    int32 ystart = 8;

    if (ncol < 5) {
        ncol = 5;
    }

    height = drow*(row + 6);
    width = ncol*dcol;
    b->nrow = row;
    b->ncol = ncol;
    base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0, width,
                                      height, 4);
    b->base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask);
    XStringListToTextProperty(&wname, 1, &winname);
    XStringListToTextProperty(&iname, 1, &iconname);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.min_width = width - 15;
    size_hints.min_height = height;
    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, &iconname, NULL, 0,
                         &size_hints, NULL, &class_hints);
    }
    many_pops_make_icon((char *)browse_bits, browse_width, browse_height, base);
    b->upper = pop_list_make_window(base, 0, 0, width, ystart + drow*6, 1);
    XSetWindowBackground(display, b->upper, my_main_win_color);
    b->main = pop_list_make_plain_window(base, 0, ystart + drow*6, width,
                                         row*drow, 1);
    XSetWindowBackground(display, b->main, my_draw_win_color);
    b->find = browser_button2(base, 0, 0, 0);
    b->get = browser_button2(base, 1, 0, 0);
    b->repl = browser_button2(base, 2, 0, 0);
    b->restore = browser_button2(base, 0, 1, 0);
    b->write = browser_button2(base, 1, 1, 0);
    b->load = browser_button2(base, 2, 1, 0);
    b->first = browser_button2(base, 0, 2, 0);
    b->last = browser_button2(base, 1, 2, 0);
    b->unrepl = browser_button2(base, 2, 2, 0);
    b->table = browser_button2(base, 2, 3, 0);
    b->up = browser_button2(base, 0, 3, 0);
    b->down = browser_button2(base, 1, 3, 0);
    b->pgup = browser_button2(base, 0, 4, 0);
    b->pgdn = browser_button2(base, 1, 4, 0);
    b->left = browser_button2(base, 0, 5, 0);
    b->right = browser_button2(base, 1, 5, 0);
    b->home = browser_button2(base, 0, 6, 0);
    b->end = browser_button2(base, 1, 6, 0);
    b->addcol = browser_button2(base, 2, 4, 0);
    b->delcol = browser_button2(base, 2, 5, 0);
    b->close = browser_button2(base, 2, 6, 0);
    b->time = browser_button2(base, 5, 0, 1);
    b->hint = pop_list_make_window(base, 0, 4*drow, width - 17, drow - 3, 1);
    XSelectInput(display, b->time, SIMPMASK);

    for (int32 i = 0; i < BMAXCOL; i++) {
        b->label[i] = browser_button_data(base, 5, i + 1, "1234567890", 1);
        XSelectInput(display, b->label[i], SIMPMASK);
    }
    if (noicon == 0) {
        XIconifyWindow(display, base, screen);
    }
    return;
}

/*   These are the global exporters ...   */

void
browser_my_expose(XEvent event) {
    if (browser_my.xflag == 0) {
        return;
    }
    expose_browser(event, browser_my);
    return;
}

void
browser_my_enter(XEvent event, int32 yn) {
    if (browser_my.xflag == 0) {
        return;
    }
    enter_browser(event, &browser_my, yn);
    return;
}

void
browser_my_button(XEvent event) {
    if (browser_my.xflag == 0) {
        return;
    }
    browser_button(event, &browser_my);
    return;
}

void
browser_my_keypress(XEvent event, int32 *used) {
    if (browser_my.xflag == 0) {
        return;
    }
    browser_keypress(event, used, &browser_my);
    return;
}

void
browser_my_resize(Window window) {
    if (browser_my.xflag == 0) {
        return;
    }
    resize_browser(window, &browser_my);
}

void
expose_browser(XEvent event, Browser b) {
    if (browser_my.xflag == 0) {
        return;
    }
    if (event.type != Expose) {
        return;
    }
    display_browser(event.xexpose.window, b);
    return;
}

void
resize_browser(Window window, Browser *b) {
    uint32 w;
    uint32 h;
    uint32 hreal;
    int32 dcol = 17*dcur_xs, drow = dcur_ys + 6;
    int32 i0;
    int32 newrow;
    int32 newcol;
    if (browser_my.xflag == 0) {
        return;
    }
    if (window != b->base) {
        return;
    }
    eig_list_get_new_size(window, &w, &h);
    hreal = h;

    /* first make sure the size is is ok  and an integral
       value of the proper width and height
     */
    i0 = (int32)w / dcol;
    if ((w % (uint32)dcol) > 0) {
        i0++;
    }
    if (i0 > b->maxcol) {
        i0 = b->maxcol;
    }

    w = (uint32)(i0*dcol);
    if (i0 < 5) {
        w = 5*(uint32)dcol;
    }
    newcol = i0;
    h = hreal - 8 - 5*(uint32)drow;
    i0 = (int32)h / drow;
    if ((h % (uint32)drow) > 0) {
        i0++;
    }
    if (i0 > b->maxrow) {
        i0 = b->maxrow;
    }
    h = (uint32)(i0*drow + dcur_xs / 2);
    newrow = i0;
    //  Now resize everything
    if (b->ncol == newcol && b->nrow == newrow) {
        return;
    }
    b->ncol = newcol;
    b->nrow = newrow;

    XResizeWindow(display, b->base, w - 17, hreal);
    XResizeWindow(display, b->upper, w - 17, (uint)(8 + drow*3));
    XResizeWindow(display, b->main, w - 17, h);

    // Let the browser know how many rows and columns of data
    return;
}

/*  if button is pressed in the browser
    then do the following  */

void
browser_button(XEvent event, Browser *b) {
    XEvent zz;
    int32 done = 1;
    Window w = event.xbutton.window;
    if (browser_my.xflag == 0) {
        return;
    }
    if (w == b->up || w == b->down || w == b->pgup || w == b->pgdn ||
        w == b->left || w == b->right) {
        done = 1;
        while (done) {
            if (w == b->up) {
                browser_data_up(b);
            }
            if (w == b->down) {
                browser_data_down(b);
            }
            if (w == b->pgup) {
                browser_data_pgup(b);
            }
            if (w == b->pgdn) {
                browser_data_pgdn(b);
            }
            if (w == b->left) {
                browser_data_left(b);
            }
            if (w == b->right) {
                browser_data_right(b);
            }
            browser_wait_a_sec(100);
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
        browser_data_home(b);
        return;
    }

    if (w == b->end) {
        browser_data_end(b);
        return;
    }

    if (w == b->first) {
        browser_data_first(b);
        return;
    }

    if (w == b->last) {
        browser_data_last(b);
        return;
    }

    if (w == b->restore) {
        browser_data_restore(b);
        return;
    }

    if (w == b->write) {
        browser_data_write(b);
        return;
    }

    if (w == b->get) {
        browser_data_get(b);
        return;
    }

    if (w == b->find) {
        browser_data_find(b);
        return;
    }

    if (w == b->repl) {
        browser_data_replace(b);
        return;
    }

    if (w == b->load) {
        browser_data_read(b);
        return;
    }

    if (w == b->addcol) {
        browser_data_add_col(b);
        return;
    }

    if (w == b->delcol) {
        browser_data_del_col(b);
        return;
    }

    if (w == b->unrepl) {
        browser_data_unreplace(b);
        return;
    }

    if (w == b->table) {
        browser_data_table(b);
        return;
    }

    if (w == b->close) {
        kill_browser(b);
        return;
    }
    return;
}

void
browser_keypress(XEvent event, int32 *used, Browser *b) {
    Window w = event.xkey.window;

    char ks;
    Window w2;
    int32 rev;

    *used = 0;
    if (browser_my.xflag == 0) {
        return;
    }
    XGetInputFocus(display, &w2, &rev);

    if (w == b->main || w == b->base || w == b->upper || w2 == b->base) {
        *used = 1;

        ks = (char)ggets_get_key_press(&event);

        /*
         XLookupString(&ev,buf,maxlen,&ks,&comp);

        */

        if (ks == KEY_UP) {
            browser_data_up(b);
            return;
        }

        if (ks == KEY_DOWN) {
            browser_data_down(b);
            return;
        }

        if (ks == KEY_PGUP) {
            browser_data_pgup(b);
            return;
        }

        if (ks == KEY_PGDN) {
            browser_data_pgdn(b);
            return;
        }

        if (ks == KEY_LEFT) {
            browser_data_left(b);
            return;
        }

        if (ks == KEY_RIGHT) {
            browser_data_right(b);
            return;
        }

        if (ks == KEY_HOME) {
            browser_data_home(b);
            return;
        }

        if (ks == KEY_END) {
            browser_data_end(b);
            return;
        }

        if (ks == 's' || ks == 'S') {
            browser_data_first(b);
            return;
        }

        if (ks == 'e' || ks == 'E') {
            browser_data_last(b);
            return;
        }

        if (ks == 'r' || ks == 'R') {
            browser_data_restore(b);
            return;
        }

        if (ks == 'W' || ks == 'w') {
            browser_data_write(b);
            return;
        }

        if (ks == 'g' || ks == 'G') {
            browser_data_get(b);
            return;
        }

        if (ks == 'f' || ks == 'F') {
            browser_data_find(b);
            return;
        }

        if (ks == 'l' || ks == 'L') {
            browser_data_read(b);
            return;
        }

        if (ks == 'u' || ks == 'U') {
            browser_data_unreplace(b);
            return;
        }

        if (ks == 't' || ks == 'T') {
            browser_data_table(b);
            return;
        }

        if (ks == 'p' || ks == 'P') {
            browser_data_replace(b);
            return;
        }

        if (ks == 'a' || ks == 'A') {
            browser_data_add_col(b);
            return;
        }

        if (ks == 'd' || ks == 'D') {
            browser_data_del_col(b);
            return;
        }

        if (ks == KEY_ESC) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
            return;
        }

    }  // end of cases
    return;
}

void
browser_data_up(Browser *b) {
    if (b->row0 > 0) {
        b->row0--;
        browser_draw_data(*b);
    }
    return;
}

void
browser_data_down(Browser *b) {
    if (b->row0 < (b->maxrow - 1)) {
        b->row0++;
        browser_draw_data(*b);
    }
    return;
}

void
browser_data_pgup(Browser *b) {
    int32 i = b->row0 - b->nrow;
    if (i > 0) {
        b->row0 = i;
    } else {
        b->row0 = 0;
    }
    browser_draw_data(*b);
    return;
}

void
browser_data_pgdn(Browser *b) {
    int32 i = b->row0 + b->nrow;
    if (i < (b->maxrow - 1)) {
        b->row0 = i;
    } else {
        b->row0 = b->maxrow - 1;
    }
    browser_draw_data(*b);
    return;
}

void
browser_data_home(Browser *b) {
    b->row0 = 0;
    b->istart = 0;
    b->iend = b->maxrow;
    browser_draw_data(*b);
    return;
}

void
browser_data_end(Browser *b) {
    b->row0 = b->maxrow - 1;
    browser_draw_data(*b);
    return;
}

void
browser_get_data_xyz(double *x, double *y, double *z, int32 i1, int32 i2,
                     int32 i3, int32 off) {
    int32 in = browser_my.row0 + off;
    *x = browser_my.data[i1][in];
    *y = browser_my.data[i2][in];
    *z = browser_my.data[i3][in];
    return;
}

void
browser_my_get_data(int32 row) {
    browser_my.row0 = row;
    browser_data_get(&browser_my);
    return;
}

void
browser_data_get(Browser *b) {
    int32 in = b->row0;
    set_ivar(0, (double)storage[0][in]);
    for (int32 i = 0; i < NODE; i++) {
        last_ic[i] = (double)storage[i + 1][in];
        set_ivar(i + 1, last_ic[i]);
    }
    for (int32 i = 0; i < NMarkov; i++) {
        last_ic[i + NODE] = (double)storage[i + NODE + 1][in];
        set_ivar(i + 1 + NODE + fix_var, last_ic[i + NODE]);
    }
    for (int32 i = NODE + NMarkov; i < NEQ; i++) {
        set_val(uvar_names[i], storage[i + 1][in]);
    }

    init_conds_redraw_ics();
}

void
browser_data_replace(Browser *b) {
    Window window;
    int32 rev;
    int32 status;
    char var[20];
    char form[80];
    XGetInputFocus(display, &window, &rev);
    strcpy(var, uvar_names[0]);
    strcpy(form, uvar_names[0]);
    status = dialog_box_get("Replace", "Variable:", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status = dialog_box_get("Replace", "Formula:", form, "Replace",
                                "Cancel", 80);
        if (status != 0) {
            browser_replace_column(var, form, b->data, b->maxrow);
        }
        browser_draw_data(*b);
    }

    XSetInputFocus(display, window, rev, CurrentTime);
    return;
}

void
browser_data_unreplace(Browser *b) {
    // browse unreplace column
    int32 n = browser_my.maxrow;
    if (!REPLACE) {
        return;
    }
    for (int32 i = 0; i < n; i++) {
        browser_my.data[R_COL][i] = old_rep[i];
    }
    browser_wipe_rep();

    browser_draw_data(*b);
    return;
}

void
browser_data_table(Browser *b) {
    Window window;
    int32 rev;
    int32 status;

    static char *name[] = {"Variable", "Xlo", "Xhi", "File"};
    char value[LENGTH(name)][MAX_LEN_SBOX];

    double xlo = 0;
    double xhi = 1;
    int32 col;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    snprintf(value[1], sizeof(value[1]), "0.00");
    snprintf(value[2], sizeof(value[2]), "1.00");
    snprintf(value[3], sizeof(value[3]), "%.*s.tab", (int)sizeof(value[0]) - 5,
             value[0]);

    XGetInputFocus(display, &window, &rev);
    status = pop_list_do_string_box(4, 4, 1, "Tabulate", name, value, 40);
    XSetInputFocus(display, window, rev, CurrentTime);
    if (status == 0) {
        return;
    }
    xlo = atof(value[1]);
    xhi = atof(value[2]);
    browser_find_variable(value[0], &col);
    if (col >= 0) {
        browser_make_d_table(xlo, xhi, col, value[3], *b);
    }
    return;
}

void
browser_data_find(Browser *b) {
    Window window;
    int32 rev;
    int32 status;

    static char *name[] = {"*0Variable", "Value"};
    char value[LENGTH(name)][MAX_LEN_SBOX];
    int32 col;
    int32 row = 0;

    double val;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    sprintf(value[1], "0.00");
    XGetInputFocus(display, &window, &rev);
    status = pop_list_do_string_box(2, 2, 1, "Find Data", name, value, 40);

    XSetInputFocus(display, window, rev, CurrentTime);

    if (status == 0) {
        return;
    }
    val = (double)atof(value[1]);
    browser_find_variable(value[0], &col);
    if (col >= 0) {
        browser_find_value(col, val, &row, *b);
    }
    if (row >= 0) {
        b->row0 = row;
        browser_draw_data(*b);
    }
    return;
}

void
browser_open_write_file(FILE **fp, char *fil, int32 *ok) {
    char ans;
    *ok = 0;
    *fp = fopen(fil, "r");
    if (*fp != NULL) {
        fclose(*fp);
        ans = (char)menudrive_two_choice("Yes", "No", "File Exists! Overwrite?",
                                         "yn");
        if (ans != 'y') {
            return;
        }
    }

    *fp = fopen(fil, "w");
    if (*fp == NULL) {
        pop_list_respond_box("Ok", "Cannot open file");
        *ok = 0;
    } else {
        *ok = 1;
    }
    return;
}

void
browser_data_read(Browser *b) {
    int32 status;
    char fil[256];
    char ch;
    FILE *fp;
    int32 len;
    int32 count = 0;
    int32 white = 1;
    double z;

    strcpy(fil, "test.dat");
    status = init_conds_file_selector("Load data", fil, "*.dat");
    if (status == 0) {
        return;
    }
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
        if (isspace((int32)ch) && (1 - white)) {
            white = 1;
        }
    } while (ch != '\n');
    rewind(fp);
    len = 0;
    while (!feof(fp)) {
        for (int32 k = 0; k < count; k++) {
            fscanf(fp, "%lf ", &z);
            if (k < b->maxcol) {
                b->data[k][len] = z;
            }
        }
        ++len;
        if (len >= max_stor) {
            break;
        }
    }
    fclose(fp);
    browser_refresh(len);
    storind = len;
    return;
}

void
browser_data_write(Browser *b) {
    int32 status;
    char fil[256];
    FILE *fp;
    int32 i;
    int32 ok;

    strcpy(fil, "test.dat");

    status = init_conds_file_selector("Write data", fil, "*.dat");
    if (status == 0) {
        return;
    }
    browser_open_write_file(&fp, fil, &ok);
    if (!ok) {
        return;
    }
    for (i = b->istart; i < b->iend; i++) {
        for (int32 j = 0; j < b->maxcol; j++) {
            fprintf(fp, "%.8g ", b->data[j][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return;
}

void
browser_data_left(Browser *b) {
    int32 i = b->col0;
    if (i > 1) {
        b->col0--;
        redraw_browser(*b);
    }
    return;
}

void
browser_data_right(Browser *b) {
    int32 i = b->col0 + b->ncol;
    if (i <= b->maxcol) {
        b->col0++;
        redraw_browser(*b);
    }
    return;
}

void
browser_data_first(Browser *b) {
    b->istart = b->row0;
}

void
browser_data_last(Browser *b) {
    b->iend = b->row0 + 1;
}

void
browser_data_restore(Browser *b) {
    integrate_restore(b->istart, b->iend);
}
