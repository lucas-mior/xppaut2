#include "functions.h"
#include "integers.h"
#include "auto_nox.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include "xpplim.h"
#include "eqns.bitmap"
#include "equilib.bitmap"

#include "mykeydef.h"
#define XDS(a)                                                                 \
    do {                                                                       \
        XDrawString(display, window, small_gc, 5, cury_offs, a, strlen(a));    \
        return;                                                                \
    } while (0)

#define MYMASK                                                                 \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask)

static struct EqBox {
    Window base;
    Window stab;
    Window rest;
    Window top;
    Window close;
    Window import;
    double y[MAX_ODE];
    double ev[MAX_ODE + MAX_ODE];
    int32 n;
    int32 flag;
    int32 info[5];
    char type[15];
} eq_box;

static struct EqList {
    Window base;
    Window up;
    Window down;
    Window list;
    Window main;
    Window close;
    int32 istart;
    int32 nlines;
    int32 flag;
} eq_list;

static int32 sparity;

static void eig_list_eq_list_down(void);
static void eig_list_eq_list_up(void);

void
eig_list_draw_eq_list(Window window) {
    char bob[300];
    char fstr[15];
    if (eq_list.flag == 0) {
        return;
    }
    if (window == eq_list.up) {
        XDS("Up");
    }
    if (window == eq_list.down) {
        XDS("Down");
    }
    if (window == eq_list.close) {
        XDS("Close");
    }
    if (window == eq_list.list) {
        for (int32 i = eq_list.istart; i < eq_list.istart + eq_list.nlines;
             i++) {
            if (i >= NEQ) {
                break;
            }
            if (i < NODE && METHOD > 0) {
                strcpy(fstr, "d%s/dT=%s");
            }
            if (i < NODE && METHOD == 0) {
                strcpy(fstr, "%s(n+1)=%s");
            }
            if (i < NODE && eq_type[i] == 1) {
                strcpy(fstr, "%s(t)=%s");
            }
            if (i >= NODE) {
                strcpy(fstr, "%s=%s");
            }
            sprintf(bob, fstr, uvar_names[i], ode_names[i]);

            bob[299] = 0;
            XDrawString(display, window, small_gc, 0,
                        cury_offs + (i - eq_list.istart)*(dcur_ys + 2), bob,
                        (int32)strlen(bob));
        }
    }
    return;
}

void
eig_list_create_eq_list(void) {
    int32 width;
    int32 height;
    int32 hlist;
    int32 hmain;
    Window base;
    static char *wname[] = {"Equations"};
    static char *iname[] = {"Eqns"};
    XTextProperty winname;
    XTextProperty iconame;
    XSizeHints size_hints;

    if (eq_list.flag == 1) {
        XRaiseWindow(display, eq_list.base);
        return;
    }

    eq_list.flag = 0;  //  this is to tell that no eq_box is here

    hmain = 3*dcur_ys;
    hlist = NEQ*(dcur_ys + 2);
    height = hlist + hmain;
    if (height > 300) {
        height = 300;
    }
    eq_list.istart = 0;
    eq_list.nlines = (height - hmain) / (dcur_ys + 2);

    width = 300;
    base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0, width,
                                      height, 4);
    eq_list.base = base;

    XStringListToTextProperty(wname, 1, &winname);
    XStringListToTextProperty(iname, 1, &iconame);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = width;
    size_hints.min_height = height;

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, &iconame, NULL, 0,
                         &size_hints, NULL, &class_hints);
    }
    many_pops_make_icon((char *)eqns_bits, eqns_width, eqns_height, base);
    eq_list.main = pop_list_make_plain_window(base, 0, 0, width, hmain, 1);
    eq_list.list = pop_list_make_plain_window(base, 0, hmain, width, hlist, 1);
    eq_list.close =
        pop_list_make_window(eq_list.main, 10, 5, 7*dcur_xs, dcur_ys + 2, 1);
    eq_list.up = pop_list_make_window(eq_list.main, 10 + 7*dcur_xs + 14, 5,
                                      7*dcur_xs, dcur_ys + 2, 1);
    eq_list.down = pop_list_make_window(eq_list.main, 10 + 14*dcur_xs + 28, 5,
                                        7*dcur_xs, dcur_ys + 2, 1);

    XSelectInput(display, eq_list.up, MYMASK);
    XSelectInput(display, eq_list.down, MYMASK);
    XSelectInput(display, eq_list.close, MYMASK);
    eq_list.flag = 1;
    return;
}

void
eig_list_eq_list_keypress(XEvent event, int32 *used) {
    Window window = event.xkey.window;

    char ks;

    *used = 0;

    if (eq_list.flag == 0) {
        return;
    }
    if (window == eq_list.main || window == eq_list.base ||
        window == eq_list.list) {
        *used = 1;
        ks = (char)ggets_get_key_press(&event);

        if (ks == KEY_UP) {
            eig_list_eq_list_up();
            return;
        }

        if (ks == KEY_DOWN) {
            eig_list_eq_list_down();
            return;
        }
    }
    return;
}

void
eig_list_enter_eq_stuff(Window window, int32 b) {
    if (eq_list.flag == 1) {
        if (window == eq_list.close || window == eq_list.up ||
            window == eq_list.down) {
            XSetWindowBorderWidth(display, window, (uint)b);
        }
    }
    if (eq_box.flag == 1 &&
        (window == eq_box.close || window == eq_box.import)) {
        XSetWindowBorderWidth(display, window, (uint)b);
    }
    return;
}

void
eig_list_eq_list_button(XEvent event) {
    Window window = event.xbutton.window;
    // pure laziness here - use this to go to eq_box
    do {
        // eig_list_eq_box_button
        if (eq_box.flag == 0) {
            break;
        }
        if (window == eq_box.import) {
            // eig list eq box import
            int32 n = eq_box.n;
            for (int32 i = 0; i < n; i++) {
                last_ic[i] = eq_box.y[i];
            }

            if (n < 20) {
                if (sparity == 0) {
                    for (int32 i = 0; i < n; i++) {
                        homo_l[i] = eq_box.y[i];
                    }
                    printf("Saved to left equilibrium\n");
                }
                if (sparity == 1) {
                    for (int32 i = 0; i < n; i++) {
                        homo_r[i] = eq_box.y[i];
                    }
                    printf("Saved to right equilibrium\n");
                }
                sparity = 1 - sparity;
            }
            init_conds_redraw_ics();
            break;
        }
        if (eq_box.close == window) {
            eq_box.flag = 0;
            browser_wait_a_sec(ClickTime);
            XDestroySubwindows(display, eq_box.base);
            XDestroyWindow(display, eq_box.base);
        }
    } while (0);
    if (eq_list.flag == 0) {
        return;
    }

    if (window == eq_list.up) {
        eig_list_eq_list_up();
        return;
    }
    if (window == eq_list.down) {
        eig_list_eq_list_down();
        return;
    }
    if (window == eq_list.close) {
        eq_list.flag = 0;
        browser_wait_a_sec(2*ClickTime);
        XDestroySubwindows(display, eq_list.base);
        XDestroyWindow(display, eq_list.base);
    }
    return;
}

void
eig_list_eq_list_up(void) {
    if (eq_list.istart > 0) {
        eq_list.istart--;
        XClearWindow(display, eq_list.list);
        eig_list_draw_eq_list(eq_list.list);
    }
    return;
}

void
eig_list_eq_list_down(void) {
    if (eq_list.istart < (NEQ - 1)) {
        eq_list.istart++;
        XClearWindow(display, eq_list.list);
        eig_list_draw_eq_list(eq_list.list);
    }
    return;
}

void
eig_list_get_new_size(Window win, uint32 *wid, uint32 *hgt) {
    int32 x;
    int32 y;
    uint32 bw;
    uint32 de;
    Window root;
    XGetGeometry(display, win, &root, &x, &y, wid, hgt, &bw, &de);
    return;
}

void
eig_list_resize_eq_list(Window win) {
    int32 nlines;
    uint32 w;
    uint32 h;
    if (eq_list.flag == 0) {
        return;
    }
    if (win != eq_list.base) {
        return;
    }
    eig_list_get_new_size(win, &w, &h);
    nlines = ((int32)h - cury_offs - 2*dcur_ys) / (dcur_ys + 2);
    eq_list.nlines = nlines;
    XResizeWindow(display, eq_list.base, w, h);
    XResizeWindow(display, eq_list.list, w, h - (uint)(2*dcur_ys));
    XResizeWindow(display, eq_list.main, w, (uint)(2*dcur_ys));
    return;
}

void
eig_list_create_eq_box(int32 cp, int32 cm, int32 rp, int32 rm, int32 im,
                       double *y, int32 n) {
    int32 width;
    int32 hstab;
    int32 hequil;
    int32 height;
    static char *name[] = {"Equilibria"};
    static char *iname[] = {"Equil"};
    int32 tpos;
    int32 tpos2;
    Window base;
    XTextProperty winname;
    XTextProperty iconame;
    XSizeHints size_hints;
    //    Do this every time
    init_conds_redraw_ics();
    for (int32 i = 0; i < n; i++) {
        eq_box.y[i] = y[i];
    }
    eq_box.n = n;
    eq_box.info[0] = cp;
    eq_box.info[1] = cm;
    eq_box.info[2] = im;
    eq_box.info[3] = rp;
    eq_box.info[4] = rm;
    if (cp > 0 || rp > 0) {
        sprintf(eq_box.type, "UNSTABLE");
    } else if (im > 0) {
        sprintf(eq_box.type, "NEUTRAL");
    } else {
        sprintf(eq_box.type, "STABLE");
    }

    if (eq_box.flag == 0) {  //   the box is not made yet
        width = (30 + 30*(int32)(n / 20))*dcur_xs;
        if (n >= 20) {
            hequil = 20*(dcur_ys + 4);
        } else {
            hequil = n*(dcur_ys + 4) + 10;
        }
        hstab = 2*dcur_y + 4*dcur_ys;
        height = hequil + hstab;
        tpos = (width - 8*dcur_x) / 2;
        tpos2 = tpos + 9*dcur_x;
        base = pop_list_make_plain_window(RootWindow(display, screen), 0, 0,
                                          width, height, 4);

        eq_box.base = base;

        XStringListToTextProperty(name, 1, &winname);
        XStringListToTextProperty(iname, 1, &iconame);
        size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
        size_hints.x = 0;
        size_hints.y = 0;
        size_hints.width = width;
        size_hints.height = height;
        size_hints.min_width = width;
        size_hints.min_height = height;
        size_hints.max_width = width;
        size_hints.max_height = height;

        XSetWMProperties(display, eq_box.base, &winname, &iconame, NULL, 0,
                         &size_hints, NULL, NULL);
        many_pops_make_icon((char *)equilib_bits, equilib_width, equilib_height,
                            base);
        eq_box.stab =
            pop_list_make_plain_window(eq_box.base, 0, 0, width, hstab, 1);
        eq_box.rest =
            pop_list_make_plain_window(eq_box.base, 0, hstab, width, hequil, 1);
        eq_box.top =
            pop_list_make_window(eq_box.stab, tpos, 2, 8*dcur_x, dcur_y + 5, 1);
        eq_box.close =
            pop_list_make_window(eq_box.base, 2, 2, 8*dcur_xs, dcur_ys + 4, 1);
        eq_box.import = pop_list_make_window(eq_box.base, tpos2, 2, 8*dcur_xs,
                                             dcur_ys + 4, 1);
        eq_box.flag = 1;
    } else {  //   Already it has been created so we are updating it
        XClearWindow(display, eq_box.top);
        eig_list_draw_eq_box(eq_box.top);
        XClearWindow(display, eq_box.stab);
        eig_list_draw_eq_box(eq_box.stab);
        XClearWindow(display, eq_box.rest);
        eig_list_draw_eq_box(eq_box.rest);
    }
    return;
}

void
eig_list_draw_eq_box(Window window) {
    int32 ncol;
    int32 n = eq_box.n;
    int32 nrow;
    int32 in;
    char temp[50];
    if (eq_box.flag == 0) {
        return;
    }
    if (window == eq_box.close) {
        XDS("Close");
    }
    if (window == eq_box.import) {
        XDS("Import");
    }
    if (window == eq_box.top) {
        XDrawString(display, eq_box.top, gc, 5, cury_off, eq_box.type,
                    (int)strlen(eq_box.type));
        return;
    }
    if (window == eq_box.stab) {
        sprintf(temp, "c+ = %d", eq_box.info[0]);
        XDrawString(display, eq_box.stab, small_gc, 2, 2*dcur_y + 6, temp,
                    (int)strlen(temp));
        sprintf(temp, "c- = %d", eq_box.info[1]);
        XDrawString(display, eq_box.stab, small_gc, 2 + 9*dcur_xs,
                    2*dcur_y + 6, temp, (int)strlen(temp));
        sprintf(temp, "im = %d", eq_box.info[2]);
        XDrawString(display, eq_box.stab, small_gc, 2 + 18*dcur_xs,
                    2*dcur_y + 6, temp, (int)strlen(temp));
        sprintf(temp, "r+ = %d", eq_box.info[3]);
        XDrawString(display, eq_box.stab, small_gc, 2,
                    2*dcur_y + 2*dcur_ys + 6, temp, (int)strlen(temp));
        sprintf(temp, "r- = %d", eq_box.info[4]);
        XDrawString(display, eq_box.stab, small_gc, 2 + 9*dcur_xs,
                    2*dcur_y + 2*dcur_ys + 6, temp, (int)strlen(temp));
        return;
    }
    if (window == eq_box.rest) {
        if (n >= 20) {
            nrow = 20;
        } else {
            nrow = n;
        }

        ncol = 1 + n / 3;

        for (int32 j = 0; j < ncol; j++) {
            for (int32 i = 0; i < nrow; i++) {
                in = j*20 + i;
                if (in >= n) {
                    continue;
                }
                sprintf(temp, "%s=%.5g", uvar_names[in], eq_box.y[in]);
                XDrawString(display, eq_box.rest, small_gc, j*28*dcur_xs + 8,
                            i*(dcur_ys + 3) + 13, temp, (int)strlen(temp));
            }
        }
        return;
    }
}
