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

#define XDS(a)                                                                                     \
    do {                                                                                           \
        XDrawString(display, window, small_gc, 5, cury_offs, a, strlen(a));                        \
        return;                                                                                    \
    } while (0)

static struct TxtView {
    Window up;
    Window down;
    Window pgup;
    Window pgdn;
    Window kill;
    Window home;
    Window end;
    Window base;
    Window text;
    Window src;
    Window action;
    int32 here;
    int32 first;
    int32 hgt;
    int32 wid;
    int32 nlines;
    int32 which;
    int32 dh;
    int32 dw;
} txtview;

/*
  [Up]   [Down]  [PgUp]  [PgDn] [Kill]
  [Home] [End]   [Src]   [Actn]
*/

static void redraw_txtview_text(void);
static void txtview_press(Window window, int32 x, int32 y);
static void resize_txtview(int32 w, int32 h);
static void do_txt_action(char *s);
static void enter_txtview(Window window, int32 val);
static void txtview_keypress(XEvent event);

void
txt_view_events(XEvent event) {
    int32 x;
    int32 y;

    if (txtview.here == 0) {
        return;
    }

    switch (event.type) {
    case Expose:
    case MapNotify:
        txt_redraw_view(event.xany.window);
        break;
    case ConfigureNotify:
        if (event.xconfigure.window != txtview.base) {
            return;
        }
        x = event.xconfigure.width;
        y = event.xconfigure.height;
        resize_txtview(x, y);
        break;
    case EnterNotify:
        enter_txtview(event.xexpose.window, 2);
        break;
    case LeaveNotify:
        enter_txtview(event.xexpose.window, 1);
        break;
    case ButtonPress:
        txtview_press(event.xbutton.window, event.xbutton.x, event.xbutton.y);
        break;
    case KeyPress:
        txtview_keypress(event);
        break;
    default:
        break;
    }
    return;
}

void
txtview_keypress(XEvent event) {
    Window window = event.xkey.window;
    char ks;
    if (window == txtview.base || window == txtview.text) {
        ks = (char)ggets_get_key_press(&event);
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
    if (w == txtview.up || w == txtview.down || w == txtview.pgup || w == txtview.pgdn ||
        w == txtview.home || w == txtview.end || w == txtview.src || w == txtview.action ||
        w == txtview.kill) {
        XSetWindowBorderWidth(display, w, (uint)val);
    }
    return;
}

void
do_txt_action(char *s) {
    int32 tb = tfBell;
    tfBell = 1;
    graphics_get_graph();
    load_eqn_extract_action(s);
    ggets_ping();
    tfBell = tb;
    numerics_chk_delay();
    init_conds_redraw_params();
    init_conds_redraw_ics();
    graphics_reset_graph();
    return;
}

void
resize_txtview(int32 w, int32 h) {
    int32 hgt = h - 8 - 3*dcur_ys;
    XMoveResizeWindow(display, txtview.text, 2, 3*dcur_ys + 5, (uint)w - 4, (uint)hgt);
    txtview.nlines = (int32)(hgt / dcur_y);
    return;
}

void
txtview_press(Window window, int32 x, int32 y) {
    int32 j;
    int32 nt;
    if (txtview.which == 1) {
        nt = n_comments;
    } else {
        nt = NLINES;
    }

    if (window == txtview.text) {
        if (txtview.which == 0) {
            return;
        }
        if (x > (2*txtview.dw)) {
            return;
        }
        j = txtview.first + y / txtview.dh;
        if ((j < n_comments) && (comments[j].aflag > 0)) {
            do_txt_action(comments[j].action);
        }
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
        browser_wait_a_sec(CLICK_TIME);
        XDestroySubwindows(display, txtview.base);
        XDestroyWindow(display, txtview.base);
    }
    if (window == txtview.pgup) {
        j = txtview.first - txtview.nlines;
        if (j < 0) {
            j = 0;
        }
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
txt_redraw_view(Window window) {
    if (window == txtview.text) {
        redraw_txtview_text();
    }
    if (window == txtview.up) {
        XDS("Up");
    }
    if (window == txtview.down) {
        XDS("Down");
    }
    if (window == txtview.pgup) {
        XDS("PgUp");
    }
    if (window == txtview.pgdn) {
        XDS("PgDn");
    }
    if (window == txtview.kill) {
        XDS("Kill");
    }
    if (window == txtview.home) {
        XDS("Home");
    }
    if (window == txtview.end) {
        XDS("End");
    }
    if (window == txtview.src) {
        XDS("Source");
    }
    if (window == txtview.action) {
        XDS("Action");
    }
    return;
}

void
redraw_txtview_text(void) {
    int32 j;
    XClearWindow(display, txtview.text);
    for (int32 i = 0; i < txtview.nlines; i++) {
        j = i + txtview.first;
        switch (txtview.which) {
        case 0:
            if (j < NLINES) {
                XDrawString(display, txtview.text, gc, txtview.dw, i*txtview.dh + cury_offs,
                            save_eqn[j], (int)strlen(save_eqn[j]));
            }
            break;
        case 1:
            if (j < n_comments) {
                XDrawString(display, txtview.text, gc, txtview.dw, i*dcur_y + cury_offs,
                            comments[j].text, (int)strlen(comments[j].text));
            }
            break;
        default:
            fprintf(stderr, "Unexpected case in %s.\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
    return;
}

void
txt_init_view(void) {
    txtview.here = 0;
    txtview.dh = dcur_y;
    txtview.dw = dcur_x;
    txtview.which = 0;
    txtview.first = 0;
    return;
}

void
txt_make_view(void) {
    int32 minwid = dcur_xs*60, minlen = 3*dcur_ys + 8 + 10*dcur_y;
    Window base;
    int32 ww = 9*dcur_xs, hh = dcur_ys + 4;
    static char *wname[] = {"Text Viewer"}, *iname[] = {"Txtview"};

    XTextProperty winname;
    XTextProperty iconname;
    XSizeHints size_hints;
    if (txtview.here == 1) {
        return;
    }
    base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0, minwid, minlen, 4);
    txtview.base = base;
    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask | StructureNotifyMask);

    XStringListToTextProperty(wname, 1, &winname);
    XStringListToTextProperty(iname, 1, &iconname);
    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.min_width = minwid;
    size_hints.min_height = minlen;
    size_hints.x = 0;
    size_hints.y = 0;
    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints, NULL, NULL);
    many_pops_make_icon((char *)txtview_bits, txtview_width, txtview_height, base);
    txtview.up = pop_list_make_window(base, dcur_xs, 2, 8*dcur_xs, dcur_ys, 1);
    txtview.down = pop_list_make_window(base, dcur_xs + ww, 2, 8*dcur_xs, dcur_ys, 1);
    txtview.pgup = pop_list_make_window(base, dcur_xs + 2*ww, 2, 8*dcur_xs, dcur_ys, 1);
    txtview.pgdn = pop_list_make_window(base, dcur_xs + 3*ww, 2, 8*dcur_xs, dcur_ys, 1);
    txtview.kill = pop_list_make_window(base, dcur_xs + 4*ww, 2, 8*dcur_xs, dcur_ys, 1);
    txtview.home = pop_list_make_window(base, dcur_xs, 2 + hh, 8*dcur_xs, dcur_ys, 1);
    txtview.end = pop_list_make_window(base, dcur_xs + ww, 2 + hh, 8*dcur_xs, dcur_ys, 1);
    txtview.src = pop_list_make_window(base, dcur_xs + 2*ww, 2 + hh, 8*dcur_xs, dcur_ys, 1);
    txtview.action = pop_list_make_window(base, dcur_xs + 3*ww, 2 + hh, 8*dcur_xs, dcur_ys, 1);
    txtview.text = pop_list_make_plain_window(base, 2, 3*dcur_ys + 5, minwid - 4, 10*dcur_y, 1);
    txtview.here = 1;
    txtview.nlines = 10;
    txtview.which = 0;
    txtview.first = 0;
    txtview.dh = dcur_y;
    txtview.dw = dcur_x;
    return;
}
