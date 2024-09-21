#include "functions.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/cursorfont.h>
#include <stdio.h>
#include "integers.h"
#define ALL_DONE 2
#define DONE_WITH_THIS 1
#define FORGET_ALL 0
#include "struct.h"

#define EV_MASK                                                                \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

#define BUT_MASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     EnterWindowMask | LeaveWindowMask)

extern Display *display;
extern Window main_win;
extern uint32 MyBackColor, MyForeColor;
extern int32 screen;
extern GC gc;
extern int32 xor_flag, DCURY, DCURX, CURY_OFF, CURS_X, CURS_Y;

int32
get_dialog(char *wname, char *name, char *value, char *ok, char *cancel,
           int32 max) {
    int32 lm = (int32)strlen(name)*DCURX;
    int32 lv = max*DCURX;
    int32 pos, colm;
    int32 lo = (int32)strlen(ok)*DCURX;
    int32 lc = (int32)strlen(cancel)*DCURX;

    int32 status;
    XTextProperty winname;

    DIALOG d;
    strcpy(d.mes_s, name);
    strcpy(d.input_s, value);
    strcpy(d.ok_s, ok);
    strcpy(d.cancel_s, cancel);
    d.base = XCreateSimpleWindow(display, RootWindow(display, screen), 0, 0,
                                 (uint)(lm + lv + 20), (uint)(30 + 2*DCURY),
                                 2, MyForeColor, MyBackColor);
    XStringListToTextProperty(&wname, 1, &winname);

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, d.base, &winname, NULL, NULL, 0, NULL, NULL,
                         &class_hints);
    }

    d.mes = XCreateSimpleWindow(display, d.base, 5, 5, (uint)lm,
                                (uint)DCURY + 8, 1, MyBackColor, MyBackColor);
    d.input = XCreateSimpleWindow(display, d.base, 10 + lm, 5, (uint)lv,
                                  (uint)DCURY + 8, 1, MyBackColor, MyBackColor);
    d.ok = XCreateSimpleWindow(display, d.base, 5, 10 + DCURY, (uint)lo + 4,
                               (uint)DCURY + 8, 1, MyForeColor, MyBackColor);
    d.cancel = XCreateSimpleWindow(display, d.base, 5 + lo + 10, 10 + DCURY,
                                   (uint)lc + 4, (uint)DCURY + 8, 1,
                                   MyForeColor, MyBackColor);

    XSelectInput(display, d.base, EV_MASK);
    XSelectInput(display, d.input, EV_MASK);
    XSelectInput(display, d.mes, EV_MASK);
    XSelectInput(display, d.ok, BUT_MASK);
    XSelectInput(display, d.cancel, BUT_MASK);
    /* txt=XCreateFontCursor(display,XC_xterm);
    XDefineCursor(display,d.input,txt); */
    XMapWindow(display, d.base);
    XMapWindow(display, d.mes);
    XMapWindow(display, d.input);
    XMapWindow(display, d.ok);
    XMapWindow(display, d.cancel);
    /*  CURS_X=strlen(d.input_s); */
    /* showchar('_', DCURX*CURS_X, 0, d.input); */
    pos = (int32)strlen(d.input_s);
    colm = DCURX*pos;
    while (true) {
        status = dialog_event_loop(&d, &pos, &colm);
        if (status != -1)
            break;
    }
    XSelectInput(display, d.cancel, EV_MASK);
    XSelectInput(display, d.ok, EV_MASK);

    waitasec(ClickTime);
    XDestroySubwindows(display, d.base);
    XDestroyWindow(display, d.base);
    XFlush(display);
    if (status == ALL_DONE || status == DONE_WITH_THIS)
        strcpy(value, d.input_s);
    return status;
}

int32
dialog_event_loop(DIALOG *d, int32 *pos, int32 *col) {
    int32 status = -1;
    int32 done = 0;
    int32 ch;
    XEvent ev;

    XNextEvent(display, &ev);

    switch (ev.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        do_expose(ev);
        display_dialog(ev.xany.window, *d, *col);
        break;
    case ButtonPress:
        if (ev.xbutton.window == d->ok) {
            status = ALL_DONE;
        }
        if (ev.xbutton.window == d->cancel) {
            status = FORGET_ALL;
        }
        if (ev.xbutton.window == d->input)
            XSetInputFocus(display, d->input, RevertToParent, CurrentTime);
        break;

    case EnterNotify:
        if (ev.xcrossing.window == d->ok || ev.xcrossing.window == d->cancel)
            XSetWindowBorderWidth(display, ev.xcrossing.window, 2);
        break;
    case LeaveNotify:
        if (ev.xcrossing.window == d->ok || ev.xcrossing.window == d->cancel)
            XSetWindowBorderWidth(display, ev.xcrossing.window, 1);
        break;

    case KeyPress:
        ch = get_key_press(&ev);
        edit_window(d->input, pos, d->input_s, col, &done, ch);
        if (done == -1)
            status = FORGET_ALL;
        if (done == 1 || done == 2)
            status = DONE_WITH_THIS;

        break;
    default:
        break;
    }
    return status;
}

void
display_dialog(Window window, DIALOG d, int32 col) {
    if (window == d.ok)
        XDrawString(display, window, gc, 0, CURY_OFF + 1, d.ok_s,
                    (int32)strlen(d.ok_s));
    if (window == d.cancel)
        XDrawString(display, window, gc, 0, CURY_OFF + 1, d.cancel_s,
                    (int32)strlen(d.cancel_s));
    if (window == d.mes)
        XDrawString(display, window, gc, 0, CURY_OFF + 1, d.mes_s,
                    (int32)strlen(d.mes_s));
    if (window == d.input) {
        XDrawString(display, window, gc, 0, CURY_OFF, d.input_s,
                    (int32)strlen(d.input_s));
        put_cursor_at(window, col, 0);
        /* showchar('_',DCURX*strlen(d.input_s),0,d.input); */
    }
    return;
}
/*  Uses Dialog boxes for input of numbers  */
