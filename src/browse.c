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

extern char *browse_hint[];
#define xds(a)                                                                 \
    {                                                                          \
        XDrawString(display, w, small_gc, 5, CURY_OFFs, a, strlen(a));         \
        return;                                                                \
    }

#define BMAXCOL 20

extern int32 *my_ode[];

extern int32 *plotlist, N_plist;

extern char *ode_names[MAX_ODE];
extern int32 NEQ, MAXSTOR, NMarkov, FIX_VAR;
extern int32 NEQ_MIN;
extern int32 NODE, NJMP;
extern int32 Xup, TipsFlag;
extern double last_ic[MAX_ODE], DELTA_T;

extern int32 NSYM, NSYM_START, NCON, NCON_START;

extern Display *display;
extern int32 screen, storind;
extern GC gc, small_gc;
extern int32 DCURX, DCURXs, DCURY, DCURYs, CURY_OFFs, CURY_OFF;
extern uint32 MyBackColor, MyForeColor, MyMainWinColor, MyDrawWinColor, GrFore,
    GrBack;

extern Window command_pop;

#define MYMASK                                                                 \
    (ButtonPressMask | ButtonReleaseMask | KeyPressMask | ExposureMask |       \
     StructureNotifyMask | LeaveWindowMask | EnterWindowMask)

#define SIMPMASK                                                               \
    (ButtonPressMask | ButtonReleaseMask | KeyPressMask | ExposureMask |       \
     StructureNotifyMask)

/*  The one and only primitive data browser   */

/*typedef struct {
                Window base,upper;
                Window find,up,down,pgup,pgdn,home,end,left,right;
                Window first,last,restore,write,get,close;
                Window load,repl,unrepl,table,addcol,delcol;
                Window main;
                Window label[BMAXCOL];
                Window time;
                Window hint;
                char hinttxt[256];
                int32 dataflag,xflag;
                int32 col0,row0,ncol,nrow;
                int32 maxrow,maxcol;
                double **data;
                int32 istart,iend;
                } BROWSER;
*/
BROWSER my_browser;

extern int32 noicon;
extern char uvar_names[MAX_ODE][12];
double *old_rep;
int32 REPLACE = 0, R_COL = 0;

extern double **storage;

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
get_data_col(int32 c) {
    return my_browser.data[c];
}

/*Excerpt from the man (Section 2) for  gettimeofday:
"The use of the timezone structure is obsolete; the tz argument should normally
be spec- ified as NULL.  The tz_dsttime field has never been used under Linux;
it has  not  been and will not be supported by libc or glibc.  Each and every
occurrence of this field in the kernel source (other than the declaration) is a
bug."
*/

int32
gettimenow(void) {
    struct timeval now;
    /*struct timezone tz;
    gettimeofday(&now,&tz);
    */
    gettimeofday(&now, NULL);
    return now.tv_usec;
}

void
waitasec(int32 msec) {
    struct timeval tim;
    /*struct timezone tz;*/
    double sec = (double)msec / 1000;
    double t1, t2;
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
get_maxrow_browser(void) {
    return my_browser.maxrow;
}

void
write_mybrowser_data(FILE *fp) {
    write_browser_data(fp, &my_browser);
    return;
}

void
write_browser_data(FILE *fp, BROWSER *b) {
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
check_for_stor(double **data) {
    if (data != storage) {
        err_msg("Only data can be in browser");
        return 0;
    } else
        return 1;
}

void
del_stor_col(char *var, BROWSER *b) {
    int32 nc;
    int32 i, j;

    find_variable(var, &nc);

    if (nc < 0) {
        err_msg("No such column....");
        return;
    }
    if (nc <= NEQ_MIN) { /* NEQ_MIN = NODE+NAUX */
        err_msg("Can't delete that column");
        return;
    }
    if (check_active_plot(nc) == 1) {
        err_msg("This variable is still actively plotted! - Cant delete!");
        return;
    }
    change_plot_vars(nc);
    if (nc < NEQ) {
        for (j = nc; j < NEQ; j++) {
            for (i = 0; i < b->maxrow; i++)
                storage[j][i] = storage[j + 1][i];
            for (i = 0; i < 400; i++)
                my_ode[j - 1 + FIX_VAR][i] = my_ode[j + FIX_VAR][i];
            strcpy(uvar_names[j - 1], uvar_names[j]);
            strcpy(ode_names[j - 1], ode_names[j]);
        }
    }
    free(storage[NEQ + 1]);
    free(ode_names[NEQ]);
    free(my_ode[NEQ + FIX_VAR]);
    NEQ--;
    b->maxcol = NEQ + 1;
    redraw_browser(*b);
}

void
data_del_col(BROWSER *b) {
    /*  this only works with storage  */
    Window w;
    int32 rev;
    if (check_for_stor(b->data) == 0)
        return;
    XGetInputFocus(display, &w, &rev);
    err_msg("Sorry - not working very well yet...");
    return;
}

void
data_add_col(BROWSER *b) {
    Window w;
    int32 rev, status;
    char var[20], form[80];
    if (check_for_stor(b->data) == 0)
        return;
    XGetInputFocus(display, &w, &rev);
    strcpy(var, "");
    strcpy(form, "");
    status = get_dialog("Add Column", "Name", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status =
            get_dialog("Add Column", "Formula:", form, "Add it", "Cancel", 80);
        if (status != 0)
            add_stor_col(var, form, b);
    }
    return;
}

int32
add_stor_col(char *name, char *formula, BROWSER *b) {
    int32 com[4000], i, j;

    if (add_expr(formula, com, &i)) {
        err_msg("Bad Formula .... ");
        return 0;
    }
    if ((my_ode[NEQ + FIX_VAR] = malloc((i + 2)*sizeof(int32))) == NULL) {
        err_msg("Cant allocate formula space");
        return 0;
    }
    if ((storage[NEQ + 1] = malloc(MAXSTOR*sizeof(double))) == NULL) {
        err_msg("Cant allocate space ....");
        free(my_ode[NEQ]);
        return 0;
    }
    if ((ode_names[NEQ] = malloc(80)) == NULL) {
        err_msg("Cannot allocate space ...");
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
chk_seq(char *f, int32 *seq, double *a1, double *a2) {
    int32 i, j = -1;
    char n1[256], n2[256];
    int32 n = strlen(f);
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
replace_column(char *var, char *form, double **dat, int32 n) {
    int32 com[200], i, j;
    int32 intflag = 0;
    int32 dif_var = -1;
    int32 seq = 0;
    double a1, a2, da = 0.0;
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
        find_variable(form, &dif_var);
        if (dif_var < 0) {
            err_msg("No such variable");
            return;
        }
    }

    if (dif_var < 0)
        chk_seq(form, &seq, &a1, &a2);
    if (seq == 1) {
        if (a1 == a2)
            seq = 3;
        else
            da = (a2 - a1) / ((double)(n - 1));
    }
    if (seq == 2)
        da = a2;
    if (seq == 3) {
        err_msg("Illegal sequence");
        return;
    }

    /*  first compile formula ... */

    if (dif_var < 0 && seq == 0) {
        if (add_expr(form, com, &i)) {
            NCON = NCON_START;
            NSYM = NSYM_START;
            err_msg("Illegal formula...");
            return;
        }
    }
    /* next check to see if column is known ... */

    find_variable(var, &i);
    if (i < 0) {
        err_msg("No such column...");
        NCON = NCON_START;
        NSYM = NSYM_START;
        return;
    }
    R_COL = i;

    /* Okay the formula is cool so lets allocate and replace  */

    wipe_rep();
    old_rep = malloc(sizeof(*old_rep)*n);
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
wipe_rep(void) {
    if (!REPLACE)
        return;
    free(old_rep);
    REPLACE = 0;
    return;
}

void
unreplace_column(void) {
    int32 i, n = my_browser.maxrow;
    if (!REPLACE)
        return;
    for (i = 0; i < n; i++)
        my_browser.data[R_COL][i] = old_rep[i];
    wipe_rep();
    return;
}

void
make_d_table(double xlo, double xhi, int32 col, char *filename, BROWSER b) {
    int32 i, npts, ok;
    FILE *fp;
    open_write_file(&fp, filename, &ok);
    if (!ok)
        return;
    npts = b.iend - b.istart;

    fprintf(fp, "%d\n", npts);
    fprintf(fp, "%g\n%g\n", xlo, xhi);
    for (i = 0; i < npts; i++)
        fprintf(fp, "%10.10g\n", b.data[col][i + b.istart]);
    fclose(fp);
    ping();
    return;
}

void
find_value(int32 col, double val, int32 *row, BROWSER b) {
    int32 n = b.maxrow;
    int32 i;
    int32 ihot = 0;
    double err, errm;
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
find_variable(char *s, int32 *col) {
    *col = -1;
    if (strcasecmp("T", s) == 0) {
        *col = 0;
        return;
    }
    *col = find_user_name(2, s);
    if (*col > -1)
        *col = *col + 1;
    return;
}

void
browse_but_on(BROWSER *b, int32 i, Window w, int32 yn) {
    int32 val = 1;
    if (yn)
        val = 2;
    XSetWindowBorderWidth(display, w, val);
    if (yn && TipsFlag && i >= 0) {
        strcpy(b->hinttxt, browse_hint[i]);
        display_browser(b->hint, *b);
    }
    return;
}

void
enter_browser(XEvent ev, BROWSER *b, int32 yn) {
    Window w = ev.xexpose.window;
    if (w == b->find)
        browse_but_on(b, 0, w, yn);
    if (w == b->up)
        browse_but_on(b, 1, w, yn);
    if (w == b->down)
        browse_but_on(b, 2, w, yn);
    if (w == b->pgup)
        browse_but_on(b, 3, w, yn);
    if (w == b->pgdn)
        browse_but_on(b, 4, w, yn);
    if (w == b->left)
        browse_but_on(b, 5, w, yn);
    if (w == b->right)
        browse_but_on(b, 6, w, yn);
    if (w == b->home)
        browse_but_on(b, 7, w, yn);
    if (w == b->end)
        browse_but_on(b, 8, w, yn);
    if (w == b->first)
        browse_but_on(b, 9, w, yn);
    if (w == b->last)
        browse_but_on(b, 10, w, yn);
    if (w == b->restore)
        browse_but_on(b, 11, w, yn);
    if (w == b->write)
        browse_but_on(b, 12, w, yn);
    if (w == b->get)
        browse_but_on(b, 13, w, yn);
    if (w == b->repl)
        browse_but_on(b, 14, w, yn);
    if (w == b->unrepl)
        browse_but_on(b, 15, w, yn);
    if (w == b->table)
        browse_but_on(b, 16, w, yn);
    if (w == b->load)
        browse_but_on(b, 17, w, yn);
    if (w == b->time)
        browse_but_on(b, 18, w, yn);
    if (w == b->addcol)
        browse_but_on(b, 19, w, yn);
    if (w == b->delcol)
        browse_but_on(b, 20, w, yn);
    if (w == b->close)
        browse_but_on(b, -1, w, yn);
    return;
}

void
display_browser(Window w, BROWSER b) {
    int32 i, i0;
    if (w == b.hint) {
        XClearWindow(display, b.hint);
        XDrawString(display, w, small_gc, 8, CURY_OFFs, b.hinttxt,
                    strlen(b.hinttxt));
        return;
    }

    if (w == b.find)
        xds("Find") if (w == b.up) xds("Up") if (w == b.down)
            xds("Down") if (w == b.pgup) xds("PgUp") if (w == b.pgdn)
                xds("PgDn") if (w == b.left) xds("Left") if (w == b.right)
                    xds("Right") if (w == b.home) xds("Home") if (w == b.end)
                        xds("End") if (w == b.first)
                            xds("First") if (w == b.last)
                                xds("Last") if (w == b.restore)
                                    xds("Restore") if (w == b.write)
                                        xds("Write") if (w == b.get)
                                            xds("Get") if (w == b.repl)
                                                xds("Replace");
    if (w == b.unrepl)
        xds("Unrepl");
    if (w == b.table)
        xds("Table");
    if (w == b.load)
        xds("Load");
    if (w == b.time)
        xds("Time") if (w == b.addcol) xds("Add col") if (w == b.close)
            xds("Close") if (w == b.delcol)
                xds("Del col") for (i = 0; i < BMAXCOL; i++) {
            if (w == b.label[i]) {
                i0 = i + b.col0 - 1;
                if (i0 < b.maxcol - 1)
                    XDrawString(display, w, small_gc, 5, CURY_OFFs,
                                uvar_names[i0], strlen(uvar_names[i0]));
            }
        }
    if (w == b.main)
        draw_data(b);
    return;
}

void
redraw_browser(BROWSER b) {
    int32 i, i0;
    Window w;
    draw_data(b);
    for (i = 0; i < BMAXCOL; i++) {
        w = b.label[i];
        i0 = i + b.col0 - 1;
        if (i0 < (b.maxcol - 1)) {
            XClearWindow(display, w);
            XDrawString(display, w, small_gc, 5, CURY_OFFs, uvar_names[i0],
                        strlen(uvar_names[i0]));
        }
    }
    return;
}

void
new_browse_dat(double **new_dat, int32 dat_len) {
    my_browser.data = new_dat;
    refresh_browser(dat_len);
}

void
refresh_browser(int32 length) {
    my_browser.dataflag = 1;
    my_browser.maxrow = length;
    my_browser.iend = length;
    if (Xup && my_browser.xflag == 1)
        draw_data(my_browser);
    return;
}

void
reset_browser(void) {
    my_browser.maxrow = 0;
    my_browser.dataflag = 0;
    return;
}

void
draw_data(BROWSER b) {
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
                        i*drow + DCURYs, string, strlen(string));
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
                            i*drow + DCURYs, string, strlen(string));
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
kill_browser(BROWSER *b) {
    b->xflag = 0;
    waitasec(ClickTime);
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
br_button(Window root, int32 row, int32 col, int32 iflag) {
    Window win;
    int32 dcol = 12*DCURXs;
    int32 drow = (DCURYs + 6);
    int32 width = 8*DCURXs;
    int32 x;
    int32 y;
    if (iflag == 1)
        dcol = 14*DCURXs;
    x = dcol*col + 4;
    y = drow*row + 4;
    win = make_window(root, x, y, width + 5, DCURYs + 1, 1);
    XSelectInput(display, win, MYMASK);
    return win;
}

Window
br_button_data(Window root, int32 row, int32 col, char *name, int32 iflag) {
    Window win;
    int32 dcol = 12*DCURXs;
    int32 drow = (DCURYs + 6);
    int32 width = strlen(name)*DCURXs;

    int32 x;
    int32 y;
    if (iflag == 1)
        dcol = 14*DCURXs;
    x = dcol*col + 4;
    y = drow*row + 4;
    win = make_window(root, x, y, width + 5, DCURYs + 1, 1);
    XSelectInput(display, win, MYMASK);
    return win;
}

void
make_browser(BROWSER *b, char *wname, char *iname, int32 row, int32 col) {
    int32 i;
    int32 ncol = col;
    int32 width, height;
    Window base;
    /* XWMHints wm_hints;
     */
    XTextProperty winname, iconname;
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
    XClassHint class_hints;
    class_hints.res_name = "";
    class_hints.res_class = "";

    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints,
                     NULL, &class_hints);
    make_icon((char *)browse_bits, browse_width, browse_height, base);
    b->upper = make_window(base, 0, 0, width, ystart + drow*6, 1);
    XSetWindowBackground(display, b->upper, MyMainWinColor);
    b->main =
        make_plain_window(base, 0, ystart + drow*6, width, row*drow, 1);
    XSetWindowBackground(display, b->main, MyDrawWinColor);
    b->find = br_button(base, 0, 0, 0);
    b->get = br_button(base, 1, 0, 0);
    b->repl = br_button(base, 2, 0, 0);
    b->restore = br_button(base, 0, 1, 0);
    b->write = br_button(base, 1, 1, 0);
    b->load = br_button(base, 2, 1, 0);
    b->first = br_button(base, 0, 2, 0);
    b->last = br_button(base, 1, 2, 0);
    b->unrepl = br_button(base, 2, 2, 0);
    b->table = br_button(base, 2, 3, 0);
    b->up = br_button(base, 0, 3, 0);
    b->down = br_button(base, 1, 3, 0);
    b->pgup = br_button(base, 0, 4, 0);
    b->pgdn = br_button(base, 1, 4, 0);
    b->left = br_button(base, 0, 5, 0);
    b->right = br_button(base, 1, 5, 0);
    b->home = br_button(base, 0, 6, 0);
    b->end = br_button(base, 1, 6, 0);
    b->addcol = br_button(base, 2, 4, 0);
    b->delcol = br_button(base, 2, 5, 0);
    b->close = br_button(base, 2, 6, 0);
    b->time = br_button(base, 5, 0, 1);
    b->hint = make_window(base, 0, 4*drow, width - 17, drow - 3, 1);
    XSelectInput(display, b->time, SIMPMASK);

    for (i = 0; i < BMAXCOL; i++) {
        b->label[i] = br_button_data(base, 5, i + 1, "1234567890", 1);
        XSelectInput(display, b->label[i], SIMPMASK);
    }
    if (noicon == 0)
        XIconifyWindow(display, base, screen);
    /*  XMapWindow(display,base);  */
    return;
}

/*   These are the global exporters ...   */

void
expose_my_browser(XEvent ev) {
    if (my_browser.xflag == 0)
        return;
    expose_browser(ev, my_browser);
    return;
}

void
enter_my_browser(XEvent ev, int32 yn) {
    if (my_browser.xflag == 0)
        return;
    enter_browser(ev, &my_browser, yn);
    return;
}

void
my_browse_button(XEvent ev) {
    if (my_browser.xflag == 0)
        return;
    browse_button(ev, &my_browser);
    return;
}

void
my_browse_keypress(XEvent ev, int32 *used) {
    if (my_browser.xflag == 0)
        return;
    browse_keypress(ev, used, &my_browser);
    return;
}

void
resize_my_browser(Window win) {
    if (my_browser.xflag == 0)
        return;
    resize_browser(win, &my_browser);
}

void
expose_browser(XEvent ev, BROWSER b) {
    if (my_browser.xflag == 0)
        return;
    if (ev.type != Expose)
        return;
    display_browser(ev.xexpose.window, b);
    return;
}

void
resize_browser(Window win, BROWSER *b) {
    uint32 w, h, hreal;
    int32 dcol = 17*DCURXs, drow = DCURYs + 6;
    int32 i0;
    int32 newrow, newcol;
    if (my_browser.xflag == 0)
        return;
    if (win != b->base)
        return;
    /* w=ev.xconfigure.width;
    h=ev.xconfigure.height; */
    get_new_size(win, &w, &h);
    hreal = h;

    /* first make sure the size is is ok  and an integral
       value of the proper width and height
     */
    i0 = w / dcol;
    if ((w % dcol) > 0)
        i0++;
    if (i0 > b->maxcol)
        i0 = b->maxcol;

    w = i0*dcol;
    if (i0 < 5)
        w = 5*dcol;
    newcol = i0;
    h = hreal - 8 - 5*drow;
    i0 = h / drow;
    if ((h % drow) > 0)
        i0++;
    if (i0 > b->maxrow)
        i0 = b->maxrow;
    h = i0*drow + DCURXs / 2;
    newrow = i0;
    /*  Now resize everything   */
    if (b->ncol == newcol && b->nrow == newrow)
        return;
    b->ncol = newcol;
    b->nrow = newrow;

    XResizeWindow(display, b->base, w - 17, hreal);
    XResizeWindow(display, b->upper, w - 17, 8 + drow*3);
    XResizeWindow(display, b->main, w - 17, h);

    /* Let the browser know how many rows and columns of data  */
    return;
}

/*  if button is pressed in the browser
    then do the following  */

void
browse_button(XEvent ev, BROWSER *b) {
    XEvent zz;
    int32 done = 1;
    Window w = ev.xbutton.window;
    if (my_browser.xflag == 0)
        return;
    if (w == b->up || w == b->down || w == b->pgup || w == b->pgdn ||
        w == b->left || w == b->right) {
        done = 1;
        while (done) {
            if (w == b->up)
                data_up(b);
            if (w == b->down)
                data_down(b);
            if (w == b->pgup)
                data_pgup(b);
            if (w == b->pgdn)
                data_pgdn(b);
            if (w == b->left)
                data_left(b);
            if (w == b->right)
                data_right(b);
            waitasec(100);
            if (XPending(display) > 0) {
                XNextEvent(display, &zz);
                switch (zz.type) {
                case ButtonRelease:
                    done = 0;
                    break;
                }
            }
        }
        return;
    }

    if (w == b->home) {
        data_home(b);
        return;
    }

    if (w == b->end) {
        data_end(b);
        return;
    }

    if (w == b->first) {
        data_first(b);
        return;
    }

    if (w == b->last) {
        data_last(b);
        return;
    }

    if (w == b->restore) {
        data_restore(b);
        return;
    }

    if (w == b->write) {
        data_write(b);
        return;
    }

    if (w == b->get) {
        data_get(b);
        return;
    }

    if (w == b->find) {
        data_find(b);
        return;
    }

    if (w == b->repl) {
        data_replace(b);
        return;
    }

    if (w == b->load) {
        data_read(b);
        return;
    }

    if (w == b->addcol) {
        data_add_col(b);
        return;
    }

    if (w == b->delcol) {
        data_del_col(b);
        return;
    }

    if (w == b->unrepl) {
        data_unreplace(b);
        return;
    }

    if (w == b->table) {
        data_table(b);
        return;
    }

    if (w == b->close) {
        kill_browser(b);
        return;
    }
    return;
}

void
browse_keypress(XEvent ev, int32 *used, BROWSER *b) {
    Window w = ev.xkey.window;

    char ks;
    Window w2;
    int32 rev;

    *used = 0;
    if (my_browser.xflag == 0)
        return;
    XGetInputFocus(display, &w2, &rev);

    if (w == b->main || w == b->base || w == b->upper || w2 == b->base) {
        *used = 1;

        ks = (char)get_key_press(&ev);

        /*
         XLookupString(&ev,buf,maxlen,&ks,&comp);

        */

        if (ks == UP) {
            data_up(b);
            return;
        }

        if (ks == DOWN) {
            data_down(b);
            return;
        }

        if (ks == PGUP) {
            data_pgup(b);
            return;
        }

        if (ks == PGDN) {
            data_pgdn(b);
            return;
        }

        if (ks == LEFT) {
            data_left(b);
            return;
        }

        if (ks == RIGHT) {
            data_right(b);
            return;
        }

        if (ks == HOME) {
            data_home(b);
            return;
        }

        if (ks == END) {
            data_end(b);
            return;
        }

        if (ks == 's' || ks == 'S') {
            data_first(b);
            return;
        }

        if (ks == 'e' || ks == 'E') {
            data_last(b);
            return;
        }

        if (ks == 'r' || ks == 'R') {
            data_restore(b);
            return;
        }

        if (ks == 'W' || ks == 'w') {
            data_write(b);
            return;
        }

        if (ks == 'g' || ks == 'G') {
            data_get(b);
            return;
        }

        if (ks == 'f' || ks == 'F') {
            data_find(b);
            return;
        }

        if (ks == 'l' || ks == 'L') {
            data_read(b);
            return;
        }

        if (ks == 'u' || ks == 'U') {
            data_unreplace(b);
            return;
        }

        if (ks == 't' || ks == 'T') {
            data_table(b);
            return;
        }

        if (ks == 'p' || ks == 'P') {
            data_replace(b);
            return;
        }

        if (ks == 'a' || ks == 'A') {
            data_add_col(b);
            return;
        }

        if (ks == 'd' || ks == 'D') {
            data_del_col(b);
            return;
        }

        if (ks == ESC) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
            return;
        }

    } /* end of cases */
    return;
}

void
data_up(BROWSER *b) {
    if (b->row0 > 0) {
        b->row0--;
        draw_data(*b);
    }
    return;
}

void
data_down(BROWSER *b) {
    if (b->row0 < (b->maxrow - 1)) {
        b->row0++;
        draw_data(*b);
    }
    return;
}

void
data_pgup(BROWSER *b) {
    int32 i = b->row0 - b->nrow;
    if (i > 0)
        b->row0 = i;
    else
        b->row0 = 0;
    draw_data(*b);
    return;
}

void
data_pgdn(BROWSER *b) {
    int32 i = b->row0 + b->nrow;
    if (i < (b->maxrow - 1))
        b->row0 = i;
    else
        b->row0 = b->maxrow - 1;
    draw_data(*b);
    return;
}

void
data_home(BROWSER *b) {
    b->row0 = 0;
    b->istart = 0;
    b->iend = b->maxrow;
    draw_data(*b);
    return;
}

void
data_end(BROWSER *b) {
    b->row0 = b->maxrow - 1;
    draw_data(*b);
    return;
}

void
get_data_xyz(double *x, double *y, double *z, int32 i1, int32 i2, int32 i3,
             int32 off) {
    int32 in = my_browser.row0 + off;
    *x = my_browser.data[i1][in];
    *y = my_browser.data[i2][in];
    *z = my_browser.data[i3][in];
    return;
}

void
data_get_mybrowser(int32 row) {
    my_browser.row0 = row;
    data_get(&my_browser);
    return;
}

void
data_get(BROWSER *b) {
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

    redraw_ics();
}

void
data_replace(BROWSER *b) {
    Window w;
    int32 rev, status;
    char var[20], form[80];
    XGetInputFocus(display, &w, &rev);
    strcpy(var, uvar_names[0]);
    strcpy(form, uvar_names[0]);
    status = get_dialog("Replace", "Variable:", var, "Ok", "Cancel", 20);
    if (status != 0) {
        status =
            get_dialog("Replace", "Formula:", form, "Replace", "Cancel", 80);
        if (status != 0)
            replace_column(var, form, b->data, b->maxrow);
        draw_data(*b);
    }

    XSetInputFocus(display, w, rev, CurrentTime);
    return;
}

void
data_unreplace(BROWSER *b) {
    unreplace_column();
    draw_data(*b);
    return;
}

void
data_table(BROWSER *b) {
    Window w;
    int32 rev, status;

    static char *name[] = {"Variable", "Xlo", "Xhi", "File"};
    char value[4][MAX_LEN_SBOX];

    double xlo = 0, xhi = 1;
    int32 col;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    snprintf(value[1], sizeof(value[0]), "0.00");
    snprintf(value[2], sizeof(value[0]), "1.00");
    snprintf(value[3], sizeof(value[0]), "%s.tab", value[0]);

    XGetInputFocus(display, &w, &rev);
    status = do_string_box(4, 4, 1, "Tabulate", name, value, 40);
    XSetInputFocus(display, w, rev, CurrentTime);
    if (status == 0)
        return;
    xlo = atof(value[1]);
    xhi = atof(value[2]);
    find_variable(value[0], &col);
    if (col >= 0)
        make_d_table(xlo, xhi, col, value[3], *b);
    return;
}

void
data_find(BROWSER *b) {
    Window w;
    int32 rev, status;

    static char *name[] = {"*0Variable", "Value"};
    char value[2][MAX_LEN_SBOX];
    int32 col, row = 0;

    double val;

    strncpy(value[0], uvar_names[0], sizeof(value[0]));
    sprintf(value[1], "0.00");
    XGetInputFocus(display, &w, &rev);
    status = do_string_box(2, 2, 1, "Find Data", name, value, 40);

    XSetInputFocus(display, w, rev, CurrentTime);

    if (status == 0)
        return;
    val = (double)atof(value[1]);
    find_variable(value[0], &col);
    if (col >= 0)
        find_value(col, val, &row, *b);
    if (row >= 0) {
        b->row0 = row;
        draw_data(*b);
    }
    return;
}

void
open_write_file(FILE **fp, char *fil, int32 *ok) {
    char ans;
    *ok = 0;
    *fp = fopen(fil, "r");
    if (*fp != NULL) {
        fclose(*fp);
        ans = (char)TwoChoice("Yes", "No", "File Exists! Overwrite?", "yn");
        if (ans != 'y')
            return;
    }

    *fp = fopen(fil, "w");
    if (*fp == NULL) {
        respond_box("Ok", "Cannot open file");
        *ok = 0;
    } else
        *ok = 1;
    return;
}

void
data_read(BROWSER *b) {
    int32 status;
    char fil[256];
    char ch;
    FILE *fp;
    int32 k;
    int32 len, count = 0, white = 1;
    double z;

    strcpy(fil, "test.dat");
    /*  XGetInputFocus(display,&w,&rev);
    status=get_dialog("Load","Filename:",fil,"Ok","Cancel",40);
    XSetInputFocus(display,w,rev,CurrentTime);
    */
    status = file_selector("Load data", fil, "*.dat");
    if (status == 0)
        return;
    fp = fopen(fil, "r");
    if (fp == NULL) {
        respond_box("Ok", "Cannot open file");
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
    draw_data(*b); */
    return;
}

void
data_write(BROWSER *b) {
    int32 status;
    char fil[256];
    FILE *fp;
    int32 i, j;
    int32 ok;

    strcpy(fil, "test.dat");

    /*
    XGetInputFocus(display,&w,&rev);
     XSetInputFocus(display,command_pop,RevertToParent,CurrentTime);
     strcpy(fil,"test.dat");
     new_string("Write to:",fil);
    */
    /* status=get_dialog("Write","Filename:",fil,"Ok","Cancel",40);

       XSetInputFocus(display,w,rev,CurrentTime); */
    status = file_selector("Write data", fil, "*.dat");
    if (status == 0)
        return;
    open_write_file(&fp, fil, &ok);
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
data_left(BROWSER *b) {
    int32 i = b->col0;
    if (i > 1) {
        b->col0--;
        redraw_browser(*b);
    }
    return;
}

void
data_right(BROWSER *b) {
    int32 i = b->col0 + b->ncol;
    if (i <= b->maxcol) {
        b->col0++;
        redraw_browser(*b);
    }
    return;
}

void
data_first(BROWSER *b) {
    b->istart = b->row0;
}

void
data_last(BROWSER *b) {
    b->iend = b->row0 + 1;
}

void
data_restore(BROWSER *b) {
    restore(b->istart, b->iend);
}

void
get_col_list(char *s, int32 *cl, int32 *n) {
    int32 len, i;

    char sp[256];
    convert(s, sp);
    len = strlen(sp);
    if (len == 0) {
        for (i = 0; i < *n; i++)
            cl[i] = i;
        return;
    }
}
