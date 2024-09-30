#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/cursorfont.h>

#include "integers.h"
#include "struct.h"
#include "functions.h"

#define DONE_WITH_THIS 1

static void dialog_box_display(Window window, Dialog d, int32 col);
static int32 dialog_box_event_loop(Dialog *d, int32 *pos, int32 *col);

int32
dialog_box_get(char *wname, char *name, char *value, char *ok, char *cancel, int32 max) {
    int32 lm = (int32)strlen(name)*dcur_x;
    int32 lv = max*dcur_x;
    int32 pos;
    int32 colm;
    int32 lo = (int32)strlen(ok)*dcur_x;
    int32 lc = (int32)strlen(cancel)*dcur_x;

    int32 status;
    XTextProperty winname;

    Dialog d;
    strcpy(d.mes_s, name);
    strcpy(d.input_s, value);
    strcpy(d.ok_s, ok);
    strcpy(d.cancel_s, cancel);
    d.base = XCreateSimpleWindow(display, RootWindow(display, screen), 0, 0, (uint)(lm + lv + 20),
                                 (uint)(30 + 2*dcur_y), 2, my_fore_color, my_back_color);
    XStringListToTextProperty(&wname, 1, &winname);

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, d.base, &winname, NULL, NULL, 0, NULL, NULL, &class_hints);
    }

    d.mes = XCreateSimpleWindow(display, d.base, 5, 5, (uint)lm, (uint)dcur_y + 8, 1, my_back_color,
                                my_back_color);
    d.input = XCreateSimpleWindow(display, d.base, 10 + lm, 5, (uint)lv, (uint)dcur_y + 8, 1,
                                  my_back_color, my_back_color);
    d.ok = XCreateSimpleWindow(display, d.base, 5, 10 + dcur_y, (uint)lo + 4, (uint)dcur_y + 8, 1,
                               my_fore_color, my_back_color);
    d.cancel = XCreateSimpleWindow(display, d.base, 5 + lo + 10, 10 + dcur_y, (uint)lc + 4,
                                   (uint)dcur_y + 8, 1, my_fore_color, my_back_color);

    XSelectInput(display, d.base, MASK_EVENT);
    XSelectInput(display, d.input, MASK_EVENT);
    XSelectInput(display, d.mes, MASK_EVENT);
    XSelectInput(display, d.ok, MASK_BUTTON);
    XSelectInput(display, d.cancel, MASK_BUTTON);
    XMapWindow(display, d.base);
    XMapWindow(display, d.mes);
    XMapWindow(display, d.input);
    XMapWindow(display, d.ok);
    XMapWindow(display, d.cancel);
    pos = (int32)strlen(d.input_s);
    colm = dcur_x*pos;
    while (true) {
        status = dialog_box_event_loop(&d, &pos, &colm);
        if (status != -1) {
            break;
        }
    }
    XSelectInput(display, d.cancel, MASK_EVENT);
    XSelectInput(display, d.ok, MASK_EVENT);

    browser_wait_a_sec(CLICK_TIME);
    XDestroySubwindows(display, d.base);
    XDestroyWindow(display, d.base);
    XFlush(display);
    if (status == ALL_DONE || status == DONE_WITH_THIS) {
        strcpy(value, d.input_s);
    }
    return status;
}

int32
dialog_box_event_loop(Dialog *d, int32 *pos, int32 *col) {
    int32 status = -1;
    int32 done = 0;
    int32 ch;
    XEvent event;

    XNextEvent(display, &event);

    switch (event.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        many_pops_do_expose(event);
        dialog_box_display(event.xany.window, *d, *col);
        break;
    case ButtonPress:
        if (event.xbutton.window == d->ok) {
            status = ALL_DONE;
        }
        if (event.xbutton.window == d->cancel) {
            status = ALL_FORGET;
        }
        if (event.xbutton.window == d->input) {
            XSetInputFocus(display, d->input, RevertToParent, CurrentTime);
        }
        break;

    case EnterNotify:
        if (event.xcrossing.window == d->ok || event.xcrossing.window == d->cancel) {
            XSetWindowBorderWidth(display, event.xcrossing.window, 2);
        }
        break;
    case LeaveNotify:
        if (event.xcrossing.window == d->ok || event.xcrossing.window == d->cancel) {
            XSetWindowBorderWidth(display, event.xcrossing.window, 1);
        }
        break;

    case KeyPress:
        ch = ggets_get_key_press(&event);
        ggets_edit_window(d->input, pos, d->input_s, col, &done, ch);
        if (done == -1) {
            status = ALL_FORGET;
        }
        if (done == 1 || done == 2) {
            status = DONE_WITH_THIS;
        }

        break;
    default:
        break;
    }
    return status;
}

void
dialog_box_display(Window window, Dialog d, int32 col) {
    if (window == d.ok) {
        XDrawString(display, window, gc, 0, cury_off + 1, d.ok_s, (int32)strlen(d.ok_s));
    }
    if (window == d.cancel) {
        XDrawString(display, window, gc, 0, cury_off + 1, d.cancel_s, (int32)strlen(d.cancel_s));
    }
    if (window == d.mes) {
        XDrawString(display, window, gc, 0, cury_off + 1, d.mes_s, (int32)strlen(d.mes_s));
    }
    if (window == d.input) {
        XDrawString(display, window, gc, 0, cury_off, d.input_s, (int32)strlen(d.input_s));
        ggets_put_cursor_at(window, col, 0);
    }
    return;
}
/*  Uses Dialog boxes for input of numbers  */
