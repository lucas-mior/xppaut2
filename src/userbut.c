#include "integers.h"
#include "functions.h"
#include <stdbool.h>

#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "kbs.h"

#define USERBUTMAX 20

int32 nuserbut = 0;

USERBUT userbut[USERBUTMAX];

extern int32 MyAddedButtonColor;
extern Display *display;
extern Window main_win;
extern int32 DCURYs, DCURXs, CURY_OFFs;
extern GC small_gc;

void
user_button_events(XEvent report) {
    switch (report.type) {
    case Expose:
    case MapNotify:
        user_button_draw(report.xany.window);
        break;
    case EnterNotify:
        user_button_cross(report.xcrossing.window, 2);
        break;
    case LeaveNotify:
        user_button_cross(report.xcrossing.window, 1);
        break;
    case ButtonPress:
        user_button_press(report.xbutton.window);
        break;
    default:
        break;
    }
    return;
}

void
user_button_press(Window w) {
    int32 i;
    for (i = 0; i < nuserbut; i++) {
        if (w == userbut[i].w) {
            run_the_commands(userbut[i].com);
        }
    }
    return;
}

static void
draw_all_user_buttons(void) {
    int32 i = 0;
    for (i = 0; i < nuserbut; i++) {
        user_button_draw(userbut[i].w);
    }
    return;
}

void
user_button_draw(Window w) {
    int32 i;
    for (i = 0; i < nuserbut; i++) {
        if (w == userbut[i].w) {
            XDrawString(display, w, small_gc, 5, CURY_OFFs, userbut[i].bname,
                        (int)strlen(userbut[i].bname));
        }
    }
    return;
}

void
user_button_cross(Window w, int32 b) {
    int32 i;
    for (i = 0; i < nuserbut; i++)
        if (w == userbut[i].w) {
            XSetWindowBorderWidth(display, w, (uint)b);
            return;
        }
    return;
}

int32
find_kbs(char *sc) {
    int32 i = 0;
    while (true) {
        if (strcmp(sc, kbs[i].seq) == 0)
            return kbs[i].com;
        i++;
        if (kbs[i].com == 0)
            return -1;
    }
}

void
create_user_buttons(int32 x0, int32 y0, Window base) {
    int32 i;
    int32 x = x0;
    int32 l;
    if (nuserbut == 0)
        return;
    for (i = 0; i < nuserbut; i++) {
        l = DCURXs*((int32)strlen(userbut[i].bname) + 2);
        userbut[i].w = make_fancy_window(base, x, y0, l, DCURYs, 1);
        x = x + l + DCURXs;
    }
    draw_all_user_buttons();
    return;
}
