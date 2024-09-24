#include "integers.h"
#include "functions.h"

#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include "xpplim.h"
#include "info.bitmap"

extern int32 DisplayHeight;
extern int32 DisplayWidth;

#define EV_MASK                                                                \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

#define BUT_MASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     EnterWindowMask | LeaveWindowMask)

extern int32 screen;
extern GC gc;
extern GC small_gc;
extern int32 DCURX;
extern int32 DCURXs;
extern int32 DCURY;
extern int32 DCURYs;
extern int32 CURY_OFFs;
extern int32 CURY_OFF;

extern int32 NUPAR;
extern int32 NODE;
extern int32 NEQ;

extern Window main_win;
extern Window info_pop;
extern int32 TORUS;
extern double TOR_PERIOD;
extern int32 itor[MAX_ODE];

struct {
    Window base, done, cancel;
    Window window[MAX_ODE];
} torbox;

static void do_torus_events(void);
static void make_tor_box(char *title);
static void choose_torus(void);
static void draw_torus_box(Window win);
static void draw_tor_var(int32 i);

void
do_torus_com(int32 c) {
    int32 i;
    TORUS = 0;
    if (c == 0 || c == 2) {
        ggets_new_float("Period :", &TOR_PERIOD);
        if (TOR_PERIOD <= 0.0) {
            ggets_err_msg("Choose positive period");
            return;
        }
        if (c == 0) {
            for (i = 0; i < MAX_ODE; i++)
                itor[i] = 1;
            TORUS = 1;
            return;
        }
        /* Choose them   */
        choose_torus();
        return;
    }
    for (i = 0; i < MAX_ODE; i++)
        itor[i] = 0;
    TORUS = 0;
    return;
}

void
draw_tor_var(int32 i) {
    char strng[sizeof(*uvar_names) + 5];
    XClearWindow(display, torbox.window[i]);
    if (itor[i] == 1)
        snprintf(strng, sizeof(strng), "X  %s", uvar_names[i]);
    else
        snprintf(strng, sizeof(strng), "   %s", uvar_names[i]);
    XDrawString(display, torbox.window[i], small_gc, 0, CURY_OFFs, strng,
                (int)strlen(strng));
    return;
}

void
draw_torus_box(Window win) {
    int32 i;

    if (win == torbox.cancel) {
        XDrawString(display, win, small_gc, 5, CURY_OFFs, "Cancel", 6);
        return;
    }
    if (win == torbox.done) {
        XDrawString(display, win, small_gc, 5, CURY_OFFs, "Done", 4);
        return;
    }

    for (i = 0; i < NEQ; i++) {
        if (win == torbox.window[i])
            draw_tor_var(i);
    }
    return;
}

void
choose_torus(void) {
    int32 i;
    make_tor_box("Fold which");
    do_torus_events();
    for (i = 0; i < NEQ; i++)
        if (itor[i] == 1)
            TORUS = 1;
    return;
}

void
make_tor_box(char *title) {
    int32 ndn, nac, width, height;
    int32 nv;
    /*int32 nh; Not used anywhere*/
    int32 i, i1, j1, xpos, ypos;
    int32 xstart = DCURXs;
    int32 ystart = DCURYs;
    Window base;
    XTextProperty winname;
    XSizeHints size_hints;

    nv = 4*DisplayHeight / (5*(DCURYs + 8));
    /*nh=DisplayWidth/(18*DCURXs);*/

    if (NEQ < nv)
        ndn = NEQ;
    else
        ndn = nv;
    nac = NEQ / ndn;
    if (nac*ndn < NEQ)
        nac++;

    width = 24*DCURXs*nac + 10;
    height = 3*DCURYs + ndn*(DCURYs + 8);

    base =
        make_plain_window(RootWindow(display, screen), 0, 0, width, height, 4);

    torbox.base = base;
    XStringListToTextProperty(&title, 1, &winname);
    size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = width;
    size_hints.min_height = height;
    size_hints.max_width = width;
    size_hints.max_height = height;
    many_pops_make_icon((char *)info_bits, info_width, info_height, base);

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, NULL, NULL, 0, &size_hints,
                         NULL, &class_hints);
    }
    for (i = 0; i < NEQ; i++) {
        i1 = i / nv;
        j1 = i % nv;
        xpos = xstart + 18*DCURXs*i1;
        ypos = ystart + j1*(DCURYs + 8);
        torbox.window[i] =
            make_window(base, xpos, ypos, 15*DCURXs, DCURYs, 1);
    }

    xpos = (width - 16*DCURXs - 10) / 2;
    ypos = height - 3*DCURYs / 2;

    torbox.cancel = make_window(base, xpos, ypos, 8*DCURXs, DCURYs, 1);
    torbox.done =
        make_window(base, xpos + 8*DCURXs + 10, ypos, 8*DCURXs, DCURYs, 1);
    XSelectInput(display, torbox.cancel, BUT_MASK);
    XSelectInput(display, torbox.done, BUT_MASK);
    XRaiseWindow(display, torbox.base);
    return;
}

void
do_torus_events(void) {
    XEvent event;
    int32 status = -1;
    int32 done = 0;
    Window wt;
    int32 i;
    int32 oldit[MAX_ODE];
    for (i = 0; i < NEQ; i++)
        oldit[i] = itor[i];
    while (!done) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:

            many_pops_do_expose(event); /*  menus and graphs etc  */
            draw_torus_box(event.xany.window);
            break;
        case ButtonPress:
            if (event.xbutton.window == torbox.done) {
                status = 1;
                done = 1;
                break;
            }
            if (event.xbutton.window == torbox.cancel) {
                status = -1;
                done = 1;
                break;
            }
            for (i = 0; i < NEQ; i++) {
                if (event.xbutton.window == torbox.window[i]) {
                    itor[i] = 1 - itor[i];
                    draw_tor_var(i);
                    break;
                }
            }
            break;
        case EnterNotify:
            wt = event.xcrossing.window;
            if (wt == torbox.done || wt == torbox.cancel)
                XSetWindowBorderWidth(display, wt, 2);
            break;
        case LeaveNotify:
            wt = event.xcrossing.window;
            if (wt == torbox.done || wt == torbox.cancel)
                XSetWindowBorderWidth(display, wt, 1);
            break;
        default:
            break;
        }
    }

    if (status == -1) {
        for (i = 0; i < NEQ; i++)
            itor[i] = oldit[i];
        TORUS = 0;
    }
    XSelectInput(display, torbox.cancel, EV_MASK);
    XSelectInput(display, torbox.done, EV_MASK);
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, torbox.base);
    XDestroyWindow(display, torbox.base);
    return;
}
