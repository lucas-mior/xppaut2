#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"

#define ALL_DONE 2
#define FORGET_ALL 0
#include "struct.h"
#include "integers.h"

static void choice_box_do_checks(ChoiceBox p);
static void choice_box_display(Window window, ChoiceBox p);
static int32 do_choice_box(Window root, char *wname, int32 n, int32 mcc,
                           char **names, int32 *check, int32 type);

void
choice_box_display(Window window, ChoiceBox p) {
    int32 n = p.n;
    XSetFillStyle(display, gc, FillSolid);
    XSetForeground(display, gc, my_fore_color);

    if (window == p.ok) {
        XDrawString(display, window, gc, 0, cury_off, "Ok", 2);
    }
    if (window == p.cancel) {
        XDrawString(display, window, gc, 0, cury_off, "Cancel", 6);
    }
    for (int32 i = 0; i < n; i++) {
        if (window != p.cw[i]) {
            continue;
        }
        XDrawString(display, window, gc, 0, cury_off, p.name[i],
                    (int)strlen(p.name[i]));
        if (p.flag[i] == 1) {
            ggets_set_fore();
        } else {
            ggets_set_back();
        }
        XDrawString(display, window, gc, (p.mc + 1)*dcur_x, cury_off, "X", 1);
    }
    ggets_set_fore();
    return;
}

void
choice_box_do_checks(ChoiceBox p) {
    for (int32 i = 0; i < p.n; i++) {
        if (p.flag[i] == 1) {
            ggets_set_fore();
        } else {
            ggets_set_back();
        }
        XDrawString(display, p.cw[i], gc, (p.mc + 1)*dcur_x, cury_off, "X", 1);
    }
    ggets_set_fore();
    return;
}

void
choice_box_base(char *wname, int32 n, int32 mcc, char **names, int32 *check,
                int32 type) {
    do_choice_box(RootWindow(display, screen), wname, n, mcc, names, check,
                  type);
    return;
}

int32
do_choice_box(Window root, char *wname, int32 n, int32 mcc, char **names,
              int32 *check, int32 type) {
    ChoiceBox p;

    int32 width;
    int32 height;
    int32 maxchar;
    int32 oldcheck[MAXENTRY];
    int32 xpos;
    int32 ypos;
    int32 status;
    int32 xstart;
    int32 ystart;
    XTextProperty winname;
    XSizeHints size_hints;
    Window base;
    maxchar = mcc;
    if (mcc < 10) {
        maxchar = 10;
    }
    width = (maxchar + 5)*dcur_x;
    height = (n + 4)*(dcur_y + 16);
    base = pop_list_make_plain_window(root, 0, 0, width, height, 4);
    XStringListToTextProperty(&wname, 1, &winname);
    size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = width;
    size_hints.min_height = height;
    size_hints.max_width = width;
    size_hints.max_height = height;
    XSetWMProperties(display, base, &winname, NULL, NULL, 0, &size_hints, NULL,
                     NULL);

    ystart = dcur_y;
    xstart = dcur_x;

    p.name = names;
    p.flag = check;
    for (int32 i = 0; i < n; i++) {
        oldcheck[i] = check[i];
        xpos = xstart;
        ypos = ystart + i*(dcur_y + 10);
        p.cw[i] =
            pop_list_make_window(base, xpos, ypos, (mcc + 3)*dcur_x, dcur_y, 1);
    }

    ypos = height - 2*dcur_y;
    xpos = (width - 12*dcur_x) / 2;
    p.ok = pop_list_make_window(base, xpos, ypos, 2*dcur_x, dcur_y, 2);
    p.cancel =
        pop_list_make_window(base, xpos + 4*dcur_x, ypos, 6*dcur_x, dcur_y, 2);
    p.base = base;

    p.n = n;
    p.type = (int16)type;
    p.mc = mcc;
    choice_box_do_checks(p);
    while (true) {
        // choice_box_event_loop
        int32 nn = p.n;
        XEvent event;
        XNextEvent(display, &event);

        status = -1;

        switch (event.type) {
        case ConfigureNotify:
        case Expose:
        case MapNotify:
            choice_box_display(event.xany.window, p);
            break;
        case ButtonPress:
            if (event.xbutton.window == p.ok) {
                ggets_bar(0, 0, 200, 200, p.ok);
                status = ALL_DONE;
            }
            if (event.xbutton.window == p.cancel) {
                ggets_bar(0, 0, 200, 200, p.cancel);
                status = FORGET_ALL;
            }
            for (int32 i = 0; i < nn; i++) {
                if (event.xbutton.window == p.cw[i]) {
                    if (p.type == RADIO) {
                        for (int32 j = 0; j < nn; j++) {
                            p.flag[j] = 0;
                        }
                        p.flag[i] = 1;
                        choice_box_do_checks(p);
                    }
                    if (p.type == CHOICE) {
                        p.flag[i] = 1 - p.flag[i];
                        choice_box_do_checks(p);
                    }
                }
            }

            break;
        case KeyPress:
            break;
        default:
            break;
        }

        if (status != -1) {
            break;
        }
    }
    // choice box destroy
    browser_wait_a_sec(ClickTime);
    XDestroySubwindows(display, p.base);
    XDestroyWindow(display, p.base);

    if (status == FORGET_ALL) {
        for (int32 i = 0; i < n; i++) {
            check[i] = oldcheck[i];
        }
    }
    return status;
}
