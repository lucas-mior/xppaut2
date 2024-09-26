#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "integers.h"
#include "parserslow.h"

#define Param 1
#define IC 2

#define MYMASK                                                                 \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask)

#define SIMPMASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

static struct MyCalc {
    Window base, quit, answer;
    double last_val;
    int32 use;
} my_calc;

static int32 calc_has_eq(char *z, char *w, int32 *where);
static void ini_calc_string(char *name, char *value, int32 *pos, int32 *col);
static void make_calc(double z);
static void draw_calc(Window window);

void
draw_calc(Window window) {
    char bob[100];
    if (window == my_calc.answer) {
        XClearWindow(display, window);
        sprintf(bob, "%.16g", my_calc.last_val);
        XDrawString(display, window, small_gc, 0, CURY_OFFs, bob,
                    (int)strlen(bob));
        return;
    }
    if (window == my_calc.quit) {
        XDrawString(display, window, small_gc, 0, CURY_OFFs, "Quit", 4);
        return;
    }
    return;
}

void
make_calc(double z) {
    int32 width;
    int32 height;
    static char *name[] = {"Answer"};
    Window base;
    XTextProperty winname;
    XSizeHints size_hints;
    my_calc.last_val = z;

    if (my_calc.use == 0) {
        width = 20 + 24*DCURXs;
        height = 4*DCURYs;
        base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0, width,
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
            pop_list_make_window(base, 10, DCURYs / 2, 24*DCURXs, DCURYs, 0);
        width = (width - 4*DCURXs) / 2;
        my_calc.quit = pop_list_make_window(base, width, (int32)(2.5*DCURYs),
                                   4*DCURXs, DCURYs, 1);
        XSelectInput(display, my_calc.quit, MYMASK);
        my_calc.use = 1;
    }
    draw_calc(my_calc.answer);
    XFlush(display);
    return;
}

void
ini_calc_string(char *name, char *value, int32 *pos, int32 *col) {
    strcpy(value, " ");
    strcpy(name, "Formula:");
    *pos = (int32)strlen(value);
    *col = (*pos + (int32)strlen(name))*DCURX;
    ggets_clr_command();
    ggets_display_command(name, value, 2);
    return;
}

void
calc_q_calc(void) {
    char value[80], name[10];
    double z = 0.0;
    XEvent event;
    int32 done = 0, pos, col, flag;
    my_calc.use = 0;
    make_calc(z);
    ini_calc_string(name, value, &pos, &col);
    while (true) {
        XNextEvent(display, &event);
        draw_calc(event.xany.window);
        if (event.type == ButtonPress)
            if (event.xbutton.window == my_calc.quit)
                break;
        if (event.type == EnterNotify && event.xcrossing.window == my_calc.quit)
            XSetWindowBorderWidth(display, event.xcrossing.window, 2);

        if (event.type == LeaveNotify && event.xcrossing.window == my_calc.quit)
            XSetWindowBorderWidth(display, event.xcrossing.window, 1);
        ggets_edit_command_string(event, name, value, &done, &pos, &col);
        if (done == 1) {
            flag = calc_do_calc(value, &z);
            if (flag != -1)
                make_calc(z);
            ini_calc_string(name, value, &pos, &col);
            done = 0;
        }
        if (done == -1)
            break;
    }

    /* quit calc */
    my_calc.use = 0;
    XSelectInput(display, my_calc.quit, SIMPMASK);
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, my_calc.base);
    XDestroyWindow(display, my_calc.base);
    ggets_clr_command();
    return;
}

int32
calc_do_calc(char *temp, double *z) {
    char val[15];
    int32 ok;
    int32 i;
    double newz;
    if (strlen(temp) == 0) {
        *z = 0.0;
        return 1;
    }
    if (calc_has_eq(temp, val, &i)) {
        newz = calc(&temp[i], &ok); /*  calculate quantity  */

        if (ok == 0)
            return -1;
        i = init_conds_find_user_name(Param, val);
        if (i > -1) {
            set_val(val, newz); /* a parameter set to value  */
            *z = newz;
            init_conds_redraw_params();
        } else {
            i = init_conds_find_user_name(IC, val);
            if (i < 0) {
                ggets_err_msg("No such name!");
                return -1;
            }
            set_val(val, newz);

            last_ic[i] = newz;
            *z = newz;
            init_conds_redraw_ics();
        }
        return 0;
    }

    newz = calc(temp, &ok);
    if (ok == 0)
        return -1;
    *z = newz;
    return 1;
}

int32
calc_has_eq(char *z, char *w, int32 *where) {
    int32 i;
    for (i = 0; i < (int32)strlen(z); i++)
        if (z[i] == ':')
            break;
    if (i == (int32)strlen(z))
        return 0;
    strncpy(w, z, (usize)i);
    w[i] = 0;
    *where = i + 1;
    return 1;
}

double
calc(char *expr, int32 *ok) {
    int32 com[400], i;
    double z = 0.0;
    if (parserslow_add_expr(expr, com, &i)) {
        ggets_err_msg("Illegal formula ..");
        *ok = 0;
        goto bye;
    }
    /* fpr_command(com); */
    z = evaluate(com);
    *ok = 1;
bye:
    NCON = NCON_START;
    NSYM = NSYM_START;
    return z;
}
