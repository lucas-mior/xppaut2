#include <stdlib.h>
#include "integers.h"
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>

#include "functions.h"

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#include "txtview.bitmap"
#include "mykeydef.h"
#define MAXLINES 5000
#define MAXCOMMENTS 500

#define xds(a)                                                                 \
    do {                                                                       \
        XDrawString(display, window, small_gc, 5, CURY_OFFs, a, strlen(a));    \
        return;                                                                \
    } while (0)

extern char *save_eqn[MAXLINES];
extern ACTION comments[MAXCOMMENTS];
extern int32 n_comments, NLINES;
extern int32 tfBell;
extern Display *display;
extern int32 screen;
extern GC gc, small_gc;
extern int32 DCURX, DCURXs, DCURY, DCURYs, CURY_OFFs, CURY_OFF;

typedef struct {
    Window up, down, pgup, pgdn, kill, home, end, base, text, src, action;
    int32 here, first, hgt, wid, nlines, which;
    int32 dh, dw;
} TXTVIEW;

TXTVIEW txtview;
/*
  [Up]   [Down]  [PgUp]  [PgDn] [Kill]
  [Home] [End]   [Src]   [Actn]
*/

static void redraw_txtview_text(void);
static void txtview_press(Window window, int32 x, int32 y);
static void resize_txtview(int32 w, int32 h);
static void do_txt_action(char *s);
static void enter_txtview(Window window, int32 val);
static void txtview_keypress(XEvent ev);

void
txt_view_events(XEvent ev) {
    int32 x, y;
    if (txtview.here == 0)
        return;

    switch (ev.type) {
    case Expose:
    case MapNotify:
        redraw_txtview(ev.xany.window);
        break;
    case ConfigureNotify:
        if (ev.xconfigure.window != txtview.base)
            return;
        x = ev.xconfigure.width;
        y = ev.xconfigure.height;
        resize_txtview(x, y);
        break;
    case EnterNotify:
        enter_txtview(ev.xexpose.window, 2);
        break;
    case LeaveNotify:
        enter_txtview(ev.xexpose.window, 1);
        break;
    case ButtonPress:
        txtview_press(ev.xbutton.window, ev.xbutton.x, ev.xbutton.y);
        break;
    case KeyPress:
        txtview_keypress(ev);
        break;
    default:
        break;
    }
    return;
}

void
txtview_keypress(XEvent ev) {
    Window window = ev.xkey.window;
    char ks;
    if (window == txtview.base || window == txtview.text) {
        ks = (char)get_key_press(&ev);
        if (ks == KEY_UP) {
            txtview_press(txtview.up, 0, 0);
            return;
        }
        if (ks == KEY_DOWN) {
            txtview_press(txtview.down, 0, 0);
            return;
        }
        if (ks == KEY_PGUP) {
            txtview_press(txtview.pgup, 0, 0);
            return;
        }
        if (ks == KEY_PGDN) {
            txtview_press(txtview.pgdn, 0, 0);
            return;
        }
        if (ks == KEY_HOME) {
            txtview_press(txtview.home, 0, 0);
            return;
        }
        if (ks == KEY_END) {
            txtview_press(txtview.end, 0, 0);
            return;
        }
    }
    return;
}

void
enter_txtview(Window window, int32 val) {
    Window w = window;
    if (w == txtview.up || w == txtview.down || w == txtview.pgup ||
        w == txtview.pgdn || w == txtview.home || w == txtview.end ||
        w == txtview.src || w == txtview.action || w == txtview.kill)
        XSetWindowBorderWidth(display, w, (uint)val);
    return;
}

void
do_txt_action(char *s) {
    int32 tb = tfBell;
    tfBell = 1;
    get_graph();
    extract_action(s);
    ping();
    tfBell = tb;
    chk_delay();
    redraw_params();
    redraw_ics();
    reset_graph();
    return;
}

void
resize_txtview(int32 w, int32 h) {
    int32 hgt = h - 8 - 3*DCURYs;
    XMoveResizeWindow(display, txtview.text, 2, 3*DCURYs + 5, (uint)w - 4,
                      (uint)hgt);
    txtview.nlines = (int32)(hgt / DCURY);
    return;
}

void
txtview_press(Window window, int32 x, int32 y) {
    int32 j;
    int32 nt;
    if (txtview.which == 1)
        nt = n_comments;
    else
        nt = NLINES;

    if (window == txtview.text) {
        if (txtview.which == 0)
            return;
        if (x > (2*txtview.dw))
            return;
        j = txtview.first + y / txtview.dh;
        if ((j < n_comments) && (comments[j].aflag > 0))
            do_txt_action(comments[j].action);
        return;
    }

    if (window == txtview.up) {
        if (txtview.first > 0) {
            txtview.first -= 1;
            redraw_txtview_text();
        }
    }
    if (window == txtview.down) {
        j = txtview.first + 1 + txtview.nlines;
        if (j <= nt) {
            txtview.first += 1;
            redraw_txtview_text();
        }
    }
    if (window == txtview.home) {
        txtview.first = 0;
        redraw_txtview_text();
    }
    if (window == txtview.end) {
        j = nt - txtview.nlines;
        if (j >= 0) {
            txtview.first = j;
            redraw_txtview_text();
        }
    }
    if (window == txtview.kill) {
        txtview.here = 0;
        waitasec(ClickTime);
        XDestroySubwindows(display, txtview.base);
        XDestroyWindow(display, txtview.base);
    }
    if (window == txtview.pgup) {
        j = txtview.first - txtview.nlines;
        if (j < 0)
            j = 0;
        txtview.first = j;
        redraw_txtview_text();
    }

    if (window == txtview.pgdn) {
        j = txtview.first + txtview.nlines;
        if (j < nt) {
            txtview.first = j;
            redraw_txtview_text();
        }
    }

    if (window == txtview.src) {
        txtview.which = 0;
        redraw_txtview_text();
    }

    if (window == txtview.action) {
        if (n_comments > 0) {
            txtview.which = 1;
            redraw_txtview_text();
        }
    }
    return;
}

void
redraw_txtview(Window window) {
    if (window == txtview.text)
        redraw_txtview_text();
    if (window == txtview.up)
        xds("Up");
    if (window == txtview.down)
        xds("Down");
    if (window == txtview.pgup)
        xds("PgUp");
    if (window == txtview.pgdn)
        xds("PgDn");
    if (window == txtview.kill)
        xds("Kill");
    if (window == txtview.home)
        xds("Home");
    if (window == txtview.end)
        xds("End");
    if (window == txtview.src)
        xds("Source");
    if (window == txtview.action)
        xds("Action");
    return;
}

void
redraw_txtview_text(void) {
    int32 i, j;
    XClearWindow(display, txtview.text);
    for (i = 0; i < txtview.nlines; i++) {
        j = i + txtview.first;
        switch (txtview.which) {
        case 0:
            if (j < NLINES) {
                XDrawString(display, txtview.text, gc, txtview.dw,
                            i*txtview.dh + CURY_OFFs, save_eqn[j],
                            (int)strlen(save_eqn[j]));
            }
            break;
        case 1:
            if (j < n_comments)
                XDrawString(display, txtview.text, gc, txtview.dw,
                            i*DCURY + CURY_OFFs, comments[j].text,
                            (int)strlen(comments[j].text));
            break;
        default:
            fprintf(stderr, "Unexpected case in %s.\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
    return;
}

void
init_txtview(void) {
    txtview.here = 0;
    txtview.dh = DCURY;
    txtview.dw = DCURX;
    txtview.which = 0;
    txtview.first = 0;
    return;
}

void
make_txtview(void) {
    int32 minwid = DCURXs*60, minlen = 3*DCURYs + 8 + 10*DCURY;
    Window base;
    int32 ww = 9*DCURXs, hh = DCURYs + 4;
    static char *wname[] = {"Text Viewer"}, *iname[] = {"Txtview"};

    /*XWMHints wm_hints;
     */
    XTextProperty winname, iconname;
    XSizeHints size_hints;
    if (txtview.here == 1)
        return;
    base =
        make_plain_window(RootWindow(display, screen), 0, 0, minwid, minlen, 4);
    txtview.base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask);

    XStringListToTextProperty(wname, 1, &winname);
    XStringListToTextProperty(iname, 1, &iconname);
    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.min_width = minwid;
    size_hints.min_height = minlen;
    size_hints.x = 0;
    size_hints.y = 0;
    /*wm_hints.initial_state=IconicState;
    wm_hints.flags=StateHint;
    */
    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints,
                     NULL, NULL);
    make_icon((char *)txtview_bits, txtview_width, txtview_height, base);
    txtview.up = make_window(base, DCURXs, 2, 8*DCURXs, DCURYs, 1);
    txtview.down = make_window(base, DCURXs + ww, 2, 8*DCURXs, DCURYs, 1);
    txtview.pgup = make_window(base, DCURXs + 2*ww, 2, 8*DCURXs, DCURYs, 1);
    txtview.pgdn = make_window(base, DCURXs + 3*ww, 2, 8*DCURXs, DCURYs, 1);
    txtview.kill = make_window(base, DCURXs + 4*ww, 2, 8*DCURXs, DCURYs, 1);
    txtview.home = make_window(base, DCURXs, 2 + hh, 8*DCURXs, DCURYs, 1);
    txtview.end = make_window(base, DCURXs + ww, 2 + hh, 8*DCURXs, DCURYs, 1);
    txtview.src =
        make_window(base, DCURXs + 2*ww, 2 + hh, 8*DCURXs, DCURYs, 1);
    txtview.action =
        make_window(base, DCURXs + 3*ww, 2 + hh, 8*DCURXs, DCURYs, 1);
    txtview.text =
        make_plain_window(base, 2, 3*DCURYs + 5, minwid - 4, 10*DCURY, 1);
    txtview.here = 1;
    txtview.nlines = 10;
    txtview.which = 0;
    txtview.first = 0;
    txtview.dh = DCURY;
    txtview.dw = DCURX;
    return;
}
