#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "calc.h"
#include "ggets.h"
#include "init_conds.h"
#include "integers.h"
#include "parserslow.h"
#include "pop_list.h"
#include "xpplim.h"

#define PARAM 1
#define IC 2
#include "browse.h"

extern int32 NCON, NSYM, NCON_START, NSYM_START;

extern Display *display;
extern int32 screen;
extern GC gc, small_gc;
extern int32 DCURX, DCURXs, DCURY, DCURYs, CURY_OFFs, CURY_OFF;

#define MYMASK                                                                 \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask)

#define SIMPMASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

extern double last_ic[MAX_ODE];

struct {
    Window base, quit, answer;
    double last_val;
    int32 use;
} my_calc;

void
draw_calc(Window w) {
    char bob[100];
    if (w == my_calc.answer) {
        XClearWindow(display, w);
        sprintf(bob, "%.16g", my_calc.last_val);
        XDrawString(display, w, small_gc, 0, CURY_OFFs, bob, strlen(bob));
        return;
    }
    if (w == my_calc.quit) {
        XDrawString(display, w, small_gc, 0, CURY_OFFs, "Quit", 4);
        return;
    }
    return;
}

void
make_calc(double z)

{
    int32 width, height;
    static char *name[] = {"Answer"};
    Window base;
    XTextProperty winname;
    XSizeHints size_hints;
    my_calc.last_val = z;

    if (my_calc.use == 0) {
        width = 20 + 24 * DCURXs;
        height = 4 * DCURYs;
        base = make_plain_window(RootWindow(display, screen), 0, 0, width,
                                 height, 4);
        my_calc.base = base;
        XStringListToTextProperty(name, 1, &winname);
        size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
        size_hints.x = 0;
        size_hints.y = 0;
        size_hints.width = width;
        size_hints.height = height;
        size_hints.min_width = width;
        size_hints.min_height = height;
        size_hints.max_width = width;
        size_hints.max_height = height;

        XSetWMProperties(display, base, &winname, &winname, NULL, 0,
                         &size_hints, NULL, NULL);
        my_calc.answer =
            make_window(base, 10, DCURYs / 2, 24 * DCURXs, DCURYs, 0);
        width = (width - 4 * DCURXs) / 2;
        my_calc.quit = make_window(base, width, (int32)(2.5 * DCURYs),
                                   4 * DCURXs, DCURYs, 1);
        XSelectInput(display, my_calc.quit, MYMASK);
        my_calc.use = 1;
    }
    draw_calc(my_calc.answer);
    XFlush(display);
    return;
}

void
quit_calc(void) {
    my_calc.use = 0;
    XSelectInput(display, my_calc.quit, SIMPMASK);
    waitasec(ClickTime);
    XDestroySubwindows(display, my_calc.base);
    XDestroyWindow(display, my_calc.base);
    clr_command();
    return;
}

void
ini_calc_string(char *name, char *value, int32 *pos, int32 *col) {
    strcpy(value, " ");
    strcpy(name, "Formula:");
    *pos = strlen(value);
    *col = (*pos + strlen(name)) * DCURX;
    clr_command();
    display_command(name, value, 2, 0);
    return;
}

void
q_calc(void) {
    char value[80], name[10];
    double z = 0.0;
    XEvent ev;
    int32 done = 0, pos, col, flag;
    my_calc.use = 0;
    make_calc(z);
    ini_calc_string(name, value, &pos, &col);
    while (true) {
        XNextEvent(display, &ev);
        draw_calc(ev.xany.window);
        if (ev.type == ButtonPress)
            if (ev.xbutton.window == my_calc.quit)
                break;
        if (ev.type == EnterNotify && ev.xcrossing.window == my_calc.quit)
            XSetWindowBorderWidth(display, ev.xcrossing.window, 2);

        if (ev.type == LeaveNotify && ev.xcrossing.window == my_calc.quit)
            XSetWindowBorderWidth(display, ev.xcrossing.window, 1);
        edit_command_string(ev, name, value, &done, &pos, &col);
        if (done == 1) {
            flag = do_calc(value, &z);
            if (flag != -1)
                make_calc(z);
            ini_calc_string(name, value, &pos, &col);
            done = 0;
        }
        if (done == -1)
            break;
    }
    quit_calc();
    return;
}

int32
do_calc(char *temp, double *z) {
    char val[15];
    int32 ok;
    int32 i;
    double newz;
    if (strlen(temp) == 0) {
        *z = 0.0;
        return 1;
    }
    if (has_eq(temp, val, &i)) {

        newz = calculate(&temp[i], &ok); /*  calculate quantity  */

        if (ok == 0)
            return -1;
        i = find_user_name(PARAM, val);
        if (i > -1) {
            set_val(val, newz); /* a parameter set to value  */
            *z = newz;
            redraw_params();
        } else {
            i = find_user_name(IC, val);
            if (i < 0) {
                err_msg("No such name!");
                return -1;
            }
            set_val(val, newz);

            last_ic[i] = newz;
            *z = newz;
            redraw_ics();
        }
        return 0;
    }

    newz = calculate(temp, &ok);
    if (ok == 0)
        return -1;
    *z = newz;
    return 1;
}

int32
has_eq(char *z, char *w, int32 *where) {
    int32 i;
    for (i = 0; i < (int32)strlen(z); i++)
        if (z[i] == ':')
            break;
    if (i == (int32)strlen(z))
        return 0;
    strncpy(w, z, i);
    w[i] = 0;
    *where = i + 1;
    return 1;
}

double
calculate(char *expr, int32 *ok) {
    int32 com[400], i;
    double z = 0.0;
    if (add_expr(expr, com, &i)) {
        err_msg("Illegal formula ..");
        *ok = 0;
        goto bye;
    }
    /* fpr_command(com); */
    z = evaluate(com);
    *ok = 1;
bye:
    /* plintf(" old=%d %d  new = %d %d \n",NCON,NSYM,NCON_START,NSYM_START);  */
    NCON = NCON_START;
    NSYM = NSYM_START;
    return z;
}
