#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "functions.h"
#include "autlim.h"
#include "auto.bitmap"
#include "auto_nox.h"
#include "integers.h"
#include "math.h"
#include "mykeydef.h"
#include "parserslow.h"

#define xds(a)                                                                 \
    do {                                                                       \
        XDrawString(display, window, gc, 5, CURY_OFFb, a, strlen(a));          \
        return;                                                                \
    } while (0)

#define SBW XSetWindowBorderWidth(display, window, 1)

#define MYMASK                                                                 \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask | ButtonMotionMask)

extern int32 TrueColorFlag;
extern uint32 MyBackColor, MyForeColor, MyMainWinColor, MyDrawWinColor;
int32 AutoRedrawFlag = 1;

extern int32 screen, storind, NODE;
extern GC gc, small_gc;
extern int32 DCURX, DCURXs, DCURY, DCURYs, CURY_OFFs, CURY_OFFb, CURY_OFF;
static int32 STD_HGT_var = 0;
static int32 STD_WID_var = 0;
static int32 Auto_extra_wid, Auto_extra_hgt;
static int32 Auto_x0, Auto_y0;

/* stuff for marking a branch  */
int32 mark_flag = 0;
int32 mark_ibrs, mark_ibre;
int32 mark_ipts, mark_ipte;
int32 mark_ixs, mark_ixe, mark_iys, mark_iye;
extern Window command_pop;

extern double TEND;

extern int32 xorfix;

extern int32 TipsFlag;
extern uint32 MyBackColor, MyForeColor, MyMainWinColor, MyDrawWinColor, GrFore,
    GrBack;

extern char *auto_hint[], *aaxes_hint[], *afile_hint[], *arun_hint[];
extern char *no_hint[], *aspecial_hint[];

extern double constants[];

static uint32 DONT_XORCross = 0;

static struct {
    Window canvas, axes, numerics, grab, next, run, clear, redraw, base, per;
    Window info, param, file, abort, stab, hint, kill;
} auto_win;

static Diagram *CUR_DIAGRAM;

static void MarkAuto(int32 x, int32 y);
static int32 query_special(char *title, char *nsymb);
static void clear_msg(void);
static void find_point(int32 ibr, int32 pt);
static void RedrawMark(void);
static void auto_kill(void);
static void a_msg(int32 i, int32 v);
static Window lil_button(Window root, int32 x, int32 y);

void
auto_x11_line(int32 a, int32 b, int32 c, int32 d) {
    XDrawLine(display, auto_win.canvas, small_gc, (a), (b), (c), (d));
    return;
}

void
auto_x11_line_trans(double a, double b, double c, double d) {
    auto_x11_line(IXVal(a), IYVal(b), IXVal(c), IYVal(d));
    return;
}

void
auto_x11_text(int32 a, int32 b, char *c) {
    XDrawString(display, auto_win.canvas, small_gc, (a), (b), (c),
                (int)strlen(c));
    return;
}

void
clr_stab(void) {
    int32 r = Auto.st_wid / 4;
    XClearWindow(display, auto_win.stab);
    XDrawArc(display, auto_win.stab, small_gc, r, r, (uint)(2*r),
             (uint)(2*r), 0, 360*64);
    return;
}

void
auto_stab_line(int32 x, int32 y, int32 xp, int32 yp) {
    XDrawLine(display, auto_win.stab, small_gc, x, y, xp, yp);
    return;
}

void
clear_auto_plot(void) {
    XClearWindow(display, auto_win.canvas);
    redraw_auto_menus();
}

void
redraw_auto_menus(void) {
    display_auto(auto_win.axes);
    display_auto(auto_win.numerics);
    display_auto(auto_win.grab);
    display_auto(auto_win.run);
    display_auto(auto_win.redraw);
    display_auto(auto_win.clear);
    display_auto(auto_win.per);
    display_auto(auto_win.param);
    display_auto(auto_win.kill);
    display_auto(auto_win.file);
    display_auto(auto_win.abort);
    return;
}

int32
query_special(char *title, char *nsymb) {
    int32 status = 1;
    static char *m[] = {"BP", "EP", "HB", "LP", "MX", "PD", "TR", "UZ"};
    static char key[] = "behlmptu";
    int32 ch = (char)auto_pop_up_list(title, m, key, 8, 11, 1, 10, 10,
                                      aspecial_hint, Auto.hinttxt);
    if (ch == 'b') {
        sprintf(nsymb, "BP");
    } else if (ch == 'e') {
        sprintf(nsymb, "EP");
    } else if (ch == 'h') {
        sprintf(nsymb, "HB");
    } else if (ch == 'l') {
        sprintf(nsymb, "LP");
    } else if (ch == 'm') {
        sprintf(nsymb, "MX");
    } else if (ch == 'p') {
        sprintf(nsymb, "PD");
    } else if (ch == 't') {
        sprintf(nsymb, "TR");
    } else if (ch == 'u') {
        sprintf(nsymb, "UZ");
    } else {
        status = 0;
        sprintf(nsymb, "  ");
    }
    redraw_auto_menus();
    return status;
}

void
do_auto_range(void) {
    double t = TEND;

    if (mark_flag == 2)
        do_auto_range_go();
    TEND = t;
    return;
}

void
auto_get_info(int32 *n, char *pname) {
    int32 i1, i2, ibr;
    Diagram *d, *dnew;

    if (mark_flag == 2) {
        i1 = abs(mark_ipts);
        ibr = mark_ibrs;
        i2 = abs(mark_ipte);
        *n = abs(i2 - i1);
        d = bifd;
        while (true) {
            if (d->ibr == ibr && ((d->ntot == i1) || (d->ntot == (-i1)))) {
                strcpy(pname, upar_names[AutoPar[d->icp1]]);
                break;
            }
            dnew = d->next;
            if (dnew == NULL) {
                break;
            }
            d = dnew;
        }
    }
    return;
}

void
auto_set_mark(int32 i) {
    int32 pt, ibr;
    if (mark_flag == 2) {
        ibr = mark_ibrs;
        if (abs(mark_ipts) < abs(mark_ipte))
            pt = abs(mark_ipts) + i;
        else
            pt = abs(mark_ipte) + i;
        find_point(ibr, pt);
    }
    return;
}

void
find_point(int32 ibr, int32 pt) {
    int32 i;
    Diagram *d, *dnew;
    if (NBifs < 2)
        return;
    d = bifd;
    while (true) {
        if (d->ibr == ibr &&
            ((d->ntot == pt) ||
             (d->ntot ==
              (-pt)))) { /* need to look at both signs to ignore stability */
            /* now we use this info to set parameters and init data */
            for (i = 0; i < NODE; i++)
                set_ivar(i + 1, d->u0[i]);
            get_ic(0, d->u0);
            for (i = 0; i < NAutoPar; i++)
                constants[Auto_index_to_array[i]] = d->par[i];
            evaluate_derived();
            redo_all_fun_tables();
            redraw_params();
            redraw_ics();
            if ((d->per) > 0)
                set_total(d->per);
            break;
        }
        dnew = d->next;
        if (dnew == NULL) {
            break;
        }
        d = dnew;
    }
    return;
}

void
traverse_diagram(void) {
    Diagram *d, *dnew, *dold;
    int32 done = 0;
    int32 ix, iy, i;
    int32 lalo;
    XEvent ev;
    int32 kp;
    mark_flag = 0;
    if (NBifs < 2)
        return;

    d = bifd;
    DONT_XORCross = 0;
    traverse_out(d, &ix, &iy, 1);

    while (done == 0) {
        XNextEvent(display, &ev);
        if (ev.type == ButtonPress) {
            int32 mindex;
            double dist;
            double ndist;
            int32 xm = ev.xmotion.x;
            int32 ym = ev.xmotion.y;

            Window window = ev.xmotion.window;

            if (window == auto_win.canvas) {
                clear_msg();
                XORCross(ix, iy);
                DONT_XORCross = 1;
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 0);

                mindex = 0;
                ndist = Auto.wid*Auto.hgt;
                XORCross(ix, iy);
                lalo = load_all_labeled_orbits;
                load_all_labeled_orbits = 0;
                while (true) {
                    dist = sqrt(((double)(xm - ix))*((double)(xm - ix)) +
                                ((double)(ym - iy))*((double)(ym - iy)));
                    if (dist < ndist) {
                        ndist = dist;
                        mindex = d->index;
                    }
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                    traverse_out(
                        d, &ix, &iy,
                        0); /*Need this each time to update the distance calc*/
                }
                d = dnew;
                CUR_DIAGRAM = d;
                load_all_labeled_orbits = lalo;
                traverse_out(d, &ix, &iy, 0);

                XORCross(ix, iy);
                while (true) {
                    if (d->index == mindex) {
                        dnew = d;
                        break;
                    }
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                DONT_XORCross = 0;
                traverse_out(d, &ix, &iy, 1);
            }
        } else if (ev.type == KeyPress) {
            int32 found = 0;
            char symb[3], nsymb[3];
            clear_msg();
            kp = get_key_press(&ev);

            switch (kp) {
            case KEY_RIGHT:
                dnew = d->next;
                if (dnew == NULL)
                    dnew = bifd;
                XORCross(ix, iy);
                d = dnew;
                CUR_DIAGRAM = dnew;
                traverse_out(d, &ix, &iy, 1);
                break;

            case KEY_LEFT:
                dnew = d->prev;
                if (dnew == NULL)
                    dnew = bifd;
                XORCross(ix, iy);
                d = dnew;
                CUR_DIAGRAM = dnew;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_UP:
                if (!query_special("Next...", nsymb)) {
                    break;
                }
                XORCross(ix, iy);
                found = 0;
                dold = d;
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    get_bif_sym(symb, dnew->itp);
                    if (strcmp(symb, nsymb) == 0) {
                        d = dnew;
                        found = 1;
                        break;
                    }
                    d = dnew;
                }
                if (found) {
                    d = dnew;
                } else {
                    snprintf(Auto.hinttxt, 255, "  Higher %s not found", nsymb);
                    display_auto(auto_win.hint);
                    d = dold;
                }
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_DOWN:
                if (!query_special("Previous...", nsymb)) {
                    break;
                }
                XORCross(ix, iy);
                found = 0;
                dold = d;
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    get_bif_sym(symb, dnew->itp);
                    if (strcmp(symb, nsymb) == 0) {
                        d = dnew;
                        found = 1;
                        break;
                    }
                    d = dnew;
                }
                if (found) {
                    d = dnew;
                } else {
                    snprintf(Auto.hinttxt, 255, "  Lower %s not found", nsymb);
                    display_auto(auto_win.hint);
                    d = dold;
                }
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_TAB:
                XORCross(ix, iy);
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = bifd;
                        break;
                    }
                    d = dnew;
                    if (d->lab != 0)
                        break;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case 's': /* mark the start of a branch */
                if (mark_flag == 0) {
                    MarkAuto(ix, iy);
                    mark_ibrs = d->ibr;
                    mark_ipts = d->ntot;
                    mark_flag = 1;
                    mark_ixs = ix;
                    mark_iys = iy;
                }
                break;
            case 'e': /* mark end of branch */
                if (mark_flag == 1) {
                    MarkAuto(ix, iy);
                    mark_ibre = d->ibr;
                    mark_ipte = d->ntot;
                    mark_flag = 2;
                    mark_ixe = ix;
                    mark_iye = iy;
                }
                break;
            case KEY_END: /*All the way to end*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_HOME: /*All the way to beginning*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_PGUP: /*Same as KEY_TAB except we don't wrap*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                    if (d->lab != 0)
                        break;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case KEY_PGDN: /*REVERSE KEY_TAB*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    d = dnew;
                    if (d->lab != 0)
                        break;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;

            case KEY_FINE:
                done = 1;
                XORCross(ix, iy);
                /*Cross should be erased now that we have made our selection.*/
                /*Seems XORing it with new draw can tend to bring it back
                randomly depending on the order of window expose events.  Best
                not to do the XORCross function at all.*/
                DONT_XORCross = 1;
                redraw_diagram();
                RedrawMark();
                break;
            case KEY_ESC:
                done = -1;
                break;
            default:
                fprintf(stderr, "Unexpected case in %s.\n", __func__);
                exit(EXIT_FAILURE);
            }
        }
    }
    /* check mark_flag branch similarity */
    if (mark_flag == 2) {
        if (mark_ibrs != mark_ibre)
            mark_flag = 0;
    }
    if (done == 1) {
        grabpt.ibr = d->ibr;
        grabpt.lab = d->lab;
        for (i = 0; i < 8; i++)
            grabpt.par[i] = d->par[i];
        grabpt.icp1 = d->icp1;
        grabpt.icp2 = d->icp2;
        grabpt.per = d->per;
        grabpt.torper = d->torper;
        for (i = 0; i < NODE; i++) {
            grabpt.uhi[i] = d->uhi[i];
            grabpt.ulo[i] = d->ulo[i];
            grabpt.u0[i] = d->u0[i];
            grabpt.ubar[i] = d->ubar[i];
            set_ivar(i + 1, grabpt.u0[i]);
        }
        get_ic(0, grabpt.u0);
        grabpt.flag = 1;
        grabpt.itp = d->itp;
        grabpt.ntot = d->ntot;
        grabpt.nfpar = d->nfpar;
        grabpt.index = d->index;
        for (i = 0; i < NAutoPar; i++)
            constants[Auto_index_to_array[i]] = grabpt.par[i];
    }
    evaluate_derived();
    redo_all_fun_tables();
    redraw_params();
    redraw_ics();
}

void
clear_auto_info(void) {
    XClearWindow(display, auto_win.info);
    return;
}

void
draw_auto_info(char *bob, int32 x, int32 y) {
    XDrawString(display, auto_win.info, small_gc, x, y, bob, (int)strlen(bob));
    return;
}

void
refreshdisplay(void) {
    XFlush(display);
    return;
}

int32
byeauto_(int32 *iflag) {
    XEvent event;
    Window window;
    char ch;
    if (Auto.exist == 0)
        return 1;
    *iflag = 0;
    while (XPending(display) > 0) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:
            do_expose(event);
            break;
        case ButtonPress:
            window = event.xbutton.window;
            if (window == auto_win.abort) {
                SBW;
                *iflag = 1;
                return 1;
            }
            break;
        case KeyPress:
            ch = (char)(get_key_press(&event));
            if (ch == KEY_ESC) {
                *iflag = 1;
                return 0;
            }
            break;
        default:
            break;
        }
    }

    return 0;
}

void
Circle(int32 x, int32 y, int32 r) {
    XDrawArc(display, auto_win.canvas, small_gc, x - r, y - r, (uint)r << 1,
             (uint)r << 1, 0, 360*64);
    return;
}

void
autocol(int32 col) {
    set_scolor(col);
    return;
}

void
autobw(void) {
    XSetBackground(display, small_gc, MyBackColor);
    XSetForeground(display, small_gc, MyForeColor);
    return;
}

int32
auto_rubber(int32 *i1, int32 *j1, int32 *i2, int32 *j2, int32 flag) {
    return rubber(i1, j1, i2, j2, auto_win.canvas, flag);
}

int32
auto_pop_up_list(char *title, char **list, char *key, int32 n, int32 max,
                 int32 def, int32 x, int32 y, char **hints, char *httxt) {
    Window temp = auto_win.base;
    int32 value = pop_up_list(&temp, title, list, key, n, max, def, x, y, hints,
                              auto_win.hint, httxt);
    return value;
}

void
RedrawMark(void) {
    if (mark_flag == 2) {
        MarkAuto(mark_ixs, mark_iys);
        MarkAuto(mark_ixe, mark_iye);
    }
    return;
}

void
MarkAuto(int32 x, int32 y) {
    LineWidth(2);
    auto_x11_line(x - 8, y - 8, x + 8, y + 8);
    auto_x11_line(x + 8, y - 8, x - 8, y + 8);
    LineWidth(1);
    return;
}

void
XORCross(int32 x, int32 y) {
    if (DONT_XORCross)
        return;

    if (xorfix) {
        XSetForeground(display, small_gc, MyDrawWinColor);
        XSetBackground(display, small_gc, MyForeColor);
    }

    XSetFunction(display, small_gc, GXxor);
    LineWidth(2);
    auto_x11_line(x - 8, y, x + 8, y);
    auto_x11_line(x, y + 8, x, y - 8);
    XSetFunction(display, small_gc, GXcopy);
    LineWidth(1);
    if (xorfix) {
        XSetForeground(display, small_gc, MyForeColor);
        XSetBackground(display, small_gc, MyDrawWinColor);
    }

    XFlush(display);
    return;
}

void
FillCircle(int32 x, int32 y, int32 r) {
    int32 r2 = (int32)(r / 1.41421356 + 0.5);
    uint32 wh = (uint32)(2*r2);

    XFillArc(display, auto_win.canvas, small_gc, x - r2, y - r2, wh, wh, 0,
             360*64);
    return;
}

static void
auto_update_view(double xlo, double xhi, double ylo, double yhi) {
    Auto.xmin = xlo;
    Auto.ymin = ylo;
    Auto.xmax = xhi;
    Auto.ymax = yhi;
    redraw_diagram();
}

void
auto_scroll_window(void) {
    XEvent ev;
    int32 i = 0, j = 0;
    int32 i0 = 0, j0 = 0;
    int32 state = 0;
    double xlo = Auto.xmin;
    double ylo = Auto.ymin;
    double xhi = Auto.xmax;
    double yhi = Auto.ymax;
    double dx = 0, dy = 0;
    int32 alldone = 0;
    XSelectInput(display, auto_win.canvas,
                 KeyPressMask | ButtonPressMask | ButtonReleaseMask |
                     PointerMotionMask | ButtonMotionMask | ExposureMask);
    while (!alldone) {
        XNextEvent(display, &ev);
        switch (ev.type) {
        case KeyPress:
            alldone = 1;
            break;
        case Expose:
            do_expose(ev);
            break;
        case ButtonPress:
            if (state == 0) {
                i0 = ev.xkey.x;
                j0 = ev.xkey.y;

                state = 1;
            }
            break;
        case MotionNotify:
            if (state == 1) {
                i0 = ev.xmotion.x;
                j0 = ev.xmotion.y;
                dx = 0.0;
                dy = 0.0;
                auto_update_view(xlo + dx, xhi + dx, ylo + dy, yhi + dy);

                state = 2;
                break;
            }
            if (state == 2) {
                i = ev.xmotion.x;
                j = ev.xmotion.y;
                dx = (double)(i0 - i)*(Auto.xmax - Auto.xmin) /
                     (double)Auto.wid;
                dy = (double)(j - j0)*(Auto.ymax - Auto.ymin) /
                     (double)Auto.hgt;
                auto_update_view(xlo + dx, xhi + dx, ylo + dy, yhi + dy);
            }
            break;
        case ButtonRelease:
            state = 0;
            xlo = xlo + dx;
            xhi = xhi + dx;
            ylo = ylo + dy;
            yhi = yhi + dy;
            break;
        default:
            break;
        }
    }
    return;
}

void
LineWidth(int32 wid) {
    int32 ls = LineSolid;
    int32 cs = CapButt;
    int32 js = JoinRound;
    XSetLineAttributes(display, small_gc, (uint)wid, ls, cs, js);
    return;
}

void
auto_motion(XEvent ev) {
    int32 i = ev.xmotion.x;
    int32 j = ev.xmotion.y;
    double x, y;
    Window window = ev.xmotion.window;
    if (Auto.exist == 0)
        return;
    if (window == auto_win.canvas) {
        x = Auto.xmin +
            (double)(i - Auto.x0)*(Auto.xmax - Auto.xmin) / (double)Auto.wid;
        y = Auto.ymin + (double)(Auto.y0 - j + Auto.hgt) *
                            (Auto.ymax - Auto.ymin) / (double)Auto.hgt;
        sprintf(Auto.hinttxt, "x=%g,y=%g", x, y);
        storeautopoint(x, y);
        display_auto(auto_win.hint);
    }
    return;
}

void
display_auto(Window window) {
    int32 ix, iy;
    if (Auto.exist == 0)
        return;
    if (window == auto_win.canvas) {
        if (AutoRedrawFlag == 1)
            redraw_diagram();
    }
    if (window == auto_win.stab) {
        int32 r = Auto.st_wid / 4;
        XFlush(display);
        XDrawArc(display, auto_win.stab, small_gc, r, r, (uint)(2*r),
                 (uint)(2*r), 0, 360*64);
        if (CUR_DIAGRAM != NULL) {
            traverse_out(CUR_DIAGRAM, &ix, &iy, 1);
        }
        XFlush(display);
    }
    if (window == auto_win.axes)
        xds("Axes");
    if (window == auto_win.numerics)
        xds("Numerics");
    if (window == auto_win.grab)
        xds("Grab");
    if (window == auto_win.run)
        xds("Run");
    if (window == auto_win.redraw)
        xds("reDraw");
    if (window == auto_win.clear)
        xds("Clear");
    if (window == auto_win.per)
        xds("Usr period");
    if (window == auto_win.kill)
        xds("Close");
    if (window == auto_win.param)
        xds("Parameter");
    if (window == auto_win.file)
        xds("File");
    if (window == auto_win.abort)
        xds("ABORT");
    if (window == auto_win.hint) {
        XClearWindow(display, window);
        XDrawString(display, window, gc, 8, CURY_OFF, Auto.hinttxt,
                    (int)strlen(Auto.hinttxt));
        return;
    }
    return;
}

Window
lil_button(Window root, int32 x, int32 y) {
    Window win;
    int32 width = 12*DCURX;
    win = make_window(root, x, y, width, DCURY + 1, 1);
    XSelectInput(display, win, MYMASK);
    return win;
}

void
make_auto(char *wname, char *iname) {
    int32 x, y, wid, hgt;
    int32 addwid = 16*DCURX, addhgt = 3*DCURY, hinthgt = DCURY + 6;
    Window base = 0;
    int32 dely = DCURY + 5;
    int32 ymargin = 4*DCURYs, xmargin = 12*DCURXs;
    XTextProperty winname, iconname;
    XSizeHints size_hints;

    STD_HGT_var = 20*DCURY;
    STD_WID_var = 67*DCURX;
    Auto_extra_wid = 10 + addwid;
    Auto_extra_hgt = addhgt + 2*DCURY + hinthgt;
    wid = 10 + addwid + STD_WID_var + xmargin;
    hgt = addhgt + 2*DCURY + STD_HGT_var + ymargin + hinthgt;
    x = addwid + 5;
    y = DCURY;
    Auto_x0 = x;
    Auto_y0 = y;
    base = make_plain_window(RootWindow(display, screen), 0, 0, wid, hgt, 4);
    XSetWindowBackground(display, base, MyMainWinColor);
    auto_win.base = base;

    strcpy(Auto.hinttxt, "hint");

    XSelectInput(display, base,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask);

    XStringListToTextProperty(&wname, 1, &winname);
    XStringListToTextProperty(&iname, 1, &iconname);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.min_width = wid;
    size_hints.min_height = hgt;

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, &iconname, NULL, 0,
                         &size_hints, NULL, &class_hints);
    }

    make_icon((char *)auto_bits, auto_width, auto_height, base);

    auto_win.canvas = make_plain_window(base, x, y, STD_WID_var + xmargin,
                                        STD_HGT_var + ymargin, 1);
    XSetWindowBackground(display, auto_win.canvas, MyDrawWinColor);
    XSelectInput(display, auto_win.canvas, MYMASK);

    x = DCURX;
    y = DCURY + STD_HGT_var + ymargin - 8*DCURX;
    auto_win.stab = make_plain_window(base, x, y, 12*DCURX, 12*DCURX, 2);
    Auto.st_wid = 12*DCURX;
    x = DCURX + 2;
    y = 2*DCURY;
    Auto.hgt = STD_HGT_var;
    Auto.wid = STD_WID_var;
    Auto.x0 = 10*DCURXs;
    Auto.y0 = 2*DCURYs;
    auto_win.kill = lil_button(base, 2, 2);
    auto_win.param = lil_button(base, x, y);
    y += dely;
    auto_win.axes = lil_button(base, x, y);
    y += dely;
    auto_win.numerics = lil_button(base, x, y);
    y += dely;
    auto_win.run = lil_button(base, x, y);
    y += dely;
    auto_win.grab = lil_button(base, x, y);
    y += dely;
    auto_win.per = lil_button(base, x, y);
    y += dely;
    auto_win.clear = lil_button(base, x, y);
    y += dely;
    auto_win.redraw = lil_button(base, x, y);
    y += dely;
    auto_win.file = lil_button(base, x, y);

    y += 3*dely;
    auto_win.abort = lil_button(base, x, y);

    y = DCURY + STD_HGT_var + ymargin + 5;
    x = addwid + 5;
    auto_win.info =
        make_plain_window(base, x, y, STD_WID_var + xmargin, addhgt, 2);
    auto_win.hint = make_plain_window(base, x, y + addhgt + 6,
                                      STD_WID_var + xmargin, DCURY + 2, 2);

    draw_bif_axes();
    return;
}

void
resize_auto_window(XEvent ev) {
    int32 wid, hgt;
    int32 addhgt = (int32)(3.5*DCURY);
    int32 ymargin = 4*DCURYs, xmargin = 12*DCURXs;
    STD_HGT_var = 20*DCURY;
    STD_WID_var = 50*DCURX;

    if (ev.xconfigure.window == auto_win.base) {
        Window root;
        int32 xloc;
        int32 yloc;
        uint32 cwid;
        uint32 chgt;
        uint32 cbwid;
        uint32 cdepth;
        int32 ix, iy;

        wid = ev.xconfigure.width - Auto_extra_wid;
        hgt = ev.xconfigure.height - Auto_extra_hgt;

        addhgt = 3*DCURY;

        XResizeWindow(display, auto_win.canvas, (uint)wid, (uint)hgt);

        XGetGeometry(display, auto_win.canvas, &root, &xloc, &yloc, &cwid,
                     &chgt, &cbwid, &cdepth);

        Auto.hgt = (int32)chgt - ymargin;
        Auto.wid = (int32)cwid - xmargin;
        if (TrueColorFlag > 0) {
            XMoveResizeWindow(display, auto_win.info, xloc,
                              yloc + (int32)chgt + 4, (uint)wid, (uint)addhgt);

            XMoveResizeWindow(display, auto_win.hint, xloc,
                              yloc + (int32)chgt + addhgt + 10, (uint)wid,
                              (uint)DCURY + 2);
        }

        if (NBifs < 2)
            return;
        traverse_out(CUR_DIAGRAM, &ix, &iy, 1);
    }
    return;
}

void
a_msg(int32 i, int32 v) {
    if (v == 0 || TipsFlag == 0)
        return;
    strncpy(Auto.hinttxt, auto_hint[i], sizeof(Auto.hinttxt));
    display_auto(auto_win.hint);
    return;
}

void
clear_msg(void) {
    Auto.hinttxt[0] = '\0';
    display_auto(auto_win.hint);
    return;
}

void
auto_enter(Window window, int32 v) {
    if (Auto.exist == 0)
        return;
    if (window == auto_win.axes) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(1, v);
        return;
    }
    if (window == auto_win.numerics) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(2, v);
        return;
    }
    if (window == auto_win.grab) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(4, v);
        return;
    }
    if (window == auto_win.run) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(3, v);
        return;
    }
    if (window == auto_win.redraw) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(7, v);
        return;
    }
    if (window == auto_win.clear) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(6, v);
        return;
    }
    if (window == auto_win.per) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(5, v);
        return;
    }
    if (window == auto_win.param) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(0, v);
        return;
    }
    if (window == auto_win.kill) {
        XSetWindowBorderWidth(display, window, (uint)v);
        return;
    }
    if (window == auto_win.file) {
        XSetWindowBorderWidth(display, window, (uint)v);
        a_msg(8, v);
        return;
    }
    return;
}

void
auto_button(XEvent ev) {
    Window window = ev.xbutton.window;
    if (Auto.exist == 0)
        return;
    if (window == auto_win.axes) {
        SBW;
        auto_plot_par();
        return;
    }
    if (window == auto_win.numerics) {
        SBW;
        auto_num_par();
        return;
    }
    if (window == auto_win.grab) {
        SBW;
        auto_grab();
        return;
    }
    if (window == auto_win.run) {
        SBW;
        auto_run();
        return;
    }
    if (window == auto_win.redraw) {
        SBW;
        redraw_diagram();
        return;
    }
    if (window == auto_win.clear) {
        SBW;
        draw_bif_axes();
        return;
    }
    if (window == auto_win.per) {
        SBW;
        auto_per_par();
        return;
    }
    if (window == auto_win.param) {
        SBW;
        auto_params();
        return;
    }
    if (window == auto_win.kill) {
        SBW;
        auto_kill();
        return;
    }
    if (window == auto_win.file) {
        SBW;
        auto_file();
        return;
    }
    return;
}

void
auto_kill(void) {
    Auto.exist = 0;
    waitasec(ClickTime);
    XDestroySubwindows(display, auto_win.base);
    XDestroyWindow(display, auto_win.base);
    return;
}

void
auto_keypress(XEvent ev, int32 *used) {
    Window window = ev.xkey.window;
    char ks;
    Window w2;
    int32 rev;

    *used = 0;
    if (Auto.exist == 0)
        return;
    XGetInputFocus(display, &w2, &rev);

    if (window == auto_win.base || window == auto_win.canvas ||
        w2 == auto_win.base) {
        *used = 1;
        ks = (char)get_key_press(&ev);

        if (ks == 'a' || ks == 'A') {
            auto_plot_par();
            return;
        }
        if (ks == 'n' || ks == 'N') {
            auto_num_par();
            return;
        }
        if (ks == 'G' || ks == 'g') {
            auto_grab();
            return;
        }
        if (ks == 'R' || ks == 'r') {
            auto_run();
            return;
        }
        if (ks == 'D' || ks == 'd') {
            redraw_diagram();
            return;
        }
        if (ks == 'C' || ks == 'c') {
            draw_bif_axes();
            return;
        }
        if (ks == 'U' || ks == 'u') {
            auto_per_par();
            return;
        }
        if (ks == 'P' || ks == 'p') {
            auto_params();
            return;
        }
        if (ks == 'F' || ks == 'f') {
            auto_file();
            return;
        }

        if (ks == KEY_ESC) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
            return;
        }
    }
    return;
}
