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

#define RUBBOX 0
#define RUBLINE 1

#define RIGHT 6
#define LEFT 2
#define ESC 27
#define TAB 10
#define BAD 0
#define FINE 13

#define STD_WID 460 /* golden mean  */
#define STD_HGT 284

#define xds(a)                                                                 \
    do {                                                                          \
        XDrawString(display, w, gc, 5, CURY_OFFb, a, strlen(a));               \
        return;                                                                \
    } while (0)

#define SBW XSetWindowBorderWidth(display, w, 1)

#define MAX_AUT_PER 10

#define MYMASK                                                                 \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask | ButtonMotionMask)

#define SIMPMASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

extern Display *display;
extern int32 TrueColorFlag;
extern uint32 MyBackColor, MyForeColor, MyMainWinColor, MyDrawWinColor;
int32 AutoRedrawFlag = 1;

extern int32 screen, storind, NODE;
extern GC gc, small_gc;
extern int32 DCURX, DCURXs, DCURY, DCURYs, CURY_OFFs, CURY_OFFb, CURY_OFF;
int32 STD_HGT_var = 0;
int32 STD_WID_var = 0;
int32 Auto_extra_wid, Auto_extra_hgt;
int32 Auto_x0, Auto_y0;
/* stuff for marking a branch  */
int32 mark_flag = 0;
int32 mark_ibrs, mark_ibre;
int32 mark_ipts, mark_ipte;
int32 mark_ixs, mark_ixe, mark_iys, mark_iye;
extern Window command_pop;

extern double TEND;
extern int32 Auto_index_to_array[8];
extern int32 AutoPar[8];

extern int32 xorfix;

extern int32 TipsFlag;
extern uint32 MyBackColor, MyForeColor, MyMainWinColor, MyDrawWinColor, GrFore,
    GrBack;

extern char *auto_hint[], *aaxes_hint[], *afile_hint[], *arun_hint[],
    *no_hint[], *aspecial_hint[];

extern double constants[];

static uint32 DONT_XORCross = 0;

AUTOWIN AutoW;

extern BIFUR Auto;

extern GRABPT grabpt;

extern DIAGRAM *bifd;
DIAGRAM *CUR_DIAGRAM;

extern int32 NBifs;

/* ****************************************************
   Code here
*****************************************************/

void
ALINE(int32 a, int32 b, int32 c, int32 d) {
    XDrawLine(display, AutoW.canvas, small_gc, (a), (b), (c), (d));
    return;
}

void
DLINE(double a, double b, double c, double d) {
    ALINE(IXVal(a), IYVal(b), IXVal(c), IYVal(d));
    return;
}

void
ATEXT(int32 a, int32 b, char *c) {
    XDrawString(display, AutoW.canvas, small_gc, (a), (b), (c), strlen(c));
    return;
}

void
clr_stab(void) {
    int32 r = Auto.st_wid / 4;
    XClearWindow(display, AutoW.stab);
    XDrawArc(display, AutoW.stab, small_gc, r, r, 2*r, 2*r, 0, 360*64);
    return;
}

void
auto_stab_line(int32 x, int32 y, int32 xp, int32 yp) {
    XDrawLine(display, AutoW.stab, small_gc, x, y, xp, yp);
    return;
}

void
clear_auto_plot(void) {
    XClearWindow(display, AutoW.canvas);
    redraw_auto_menus();
}

void
redraw_auto_menus(void) {
    display_auto(AutoW.axes);
    display_auto(AutoW.numerics);
    display_auto(AutoW.grab);
    display_auto(AutoW.run);
    display_auto(AutoW.redraw);
    display_auto(AutoW.clear);
    display_auto(AutoW.per);
    display_auto(AutoW.param);
    display_auto(AutoW.kill);
    display_auto(AutoW.file);
    display_auto(AutoW.abort);
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
    DIAGRAM *d, *dnew;

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
    DIAGRAM *d, *dnew;
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
    DIAGRAM *d, *dnew, *dold;
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

            Window w = ev.xmotion.window;

            if (w == AutoW.canvas) {
                clear_msg();
                /*
                GO HOME
                */
                XORCross(ix, iy);
                DONT_XORCross = 1;
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    /*bifd = dnew;*/
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 0);
                /*
                END GO HOME
                */

                /*
                GO END
                */
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
                /*
                END GO END
                */

                /*
                GO HOME
                */
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
                    /*bifd = dnew;*/
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                DONT_XORCross = 0;
                traverse_out(d, &ix, &iy, 1);
                /*
                END GO HOME
                */
            }
        } else if (ev.type == KeyPress) {
            int32 found = 0;
            char symb[3], nsymb[3];
            clear_msg();
            kp = get_key_press(&ev);

            switch (kp) {
            case RIGHT:
                dnew = d->next;
                if (dnew == NULL)
                    dnew = bifd;
                XORCross(ix, iy);
                d = dnew;
                CUR_DIAGRAM = dnew;
                traverse_out(d, &ix, &iy, 1);
                break;

            case LEFT:
                dnew = d->prev;
                if (dnew == NULL)
                    dnew = bifd;
                XORCross(ix, iy);
                d = dnew;
                CUR_DIAGRAM = dnew;
                traverse_out(d, &ix, &iy, 1);
                break;
            case UP:
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
                    /*if(d->lab==0)break;*/
                }
                if (found) {
                    d = dnew;
                } else {
                    snprintf(Auto.hinttxt, 255, "  Higher %s not found", nsymb);
                    display_auto(AutoW.hint);
                    d = dold;
                }
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case DOWN:
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
                    display_auto(AutoW.hint);
                    d = dold;
                }
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case TAB:
                XORCross(ix, iy);
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = bifd;
                        break;
                    } /*TAB wraps*/
                    d = dnew;
                    if (d->lab != 0)
                        break;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
                /* New code */
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
            case END: /*All the way to end*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->next;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    /*bifd = dnew;*/
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case HOME: /*All the way to beginning*/
                XORCross(ix, iy);
                while (true) {
                    dnew = d->prev;
                    if (dnew == NULL) {
                        dnew = d;
                        break;
                    }
                    /*bifd = dnew;*/
                    d = dnew;
                }
                d = dnew;
                CUR_DIAGRAM = d;
                traverse_out(d, &ix, &iy, 1);
                break;
            case PGUP: /*Same as TAB except we don't wrap*/
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
            case PGDN: /*REVERSE TAB*/
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

            case FINE:
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
            case ESC:
                done = -1;
                break;
            default:
                fprintf(stderr, "Unexpected case in %s.\n", __func__);
                exit(EXIT_FAILURE);
            }
        }
    }
    /*XORCross(ix,iy);
     */
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
    XClearWindow(display, AutoW.info);
    return;
}

void
draw_auto_info(char *bob, int32 x, int32 y) {
    XDrawString(display, AutoW.info, small_gc, x, y, bob, strlen(bob));
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
    Window w;
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
            w = event.xbutton.window;
            if (w == AutoW.abort) {
                SBW;
                *iflag = 1;
                return 1;
            }
            break;
        case KeyPress:
            ch = (char)(get_key_press(&event));
            if (ch == ESC) {
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
    XDrawArc(display, AutoW.canvas, small_gc, x - r, y - r, r << 1, r << 1, 0,
             360*64);
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
    return rubber(i1, j1, i2, j2, AutoW.canvas, flag);
}

int32
auto_pop_up_list(char *title, char **list, char *key, int32 n, int32 max,
                 int32 def, int32 x, int32 y, char **hints, char *httxt) {
    Window temp = AutoW.base;
    int32 value = pop_up_list(&temp, title, list, key, n, max, def, x, y, hints,
                              AutoW.hint, httxt);
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
    ALINE(x - 8, y - 8, x + 8, y + 8);
    ALINE(x + 8, y - 8, x - 8, y + 8);
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
    ALINE(x - 8, y, x + 8, y);
    ALINE(x, y + 8, x, y - 8);
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
    int32 wh = 2*r2;

    XFillArc(display, AutoW.canvas, small_gc, x - r2, y - r2, wh, wh, 0,
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
    /*    printf("xin: %g %g %g %g\n",xlo,xhi,ylo,yhi); */
    XSelectInput(display, AutoW.canvas,
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

                /*  x0=Auto.xmin+(double)(i-Auto.x0)*(Auto.xmax-Auto.xmin)/(double)Auto.wid;
                    y0=Auto.ymin+(double)(Auto.y0-j+Auto.hgt)*(Auto.ymax-Auto.ymin)/(double)Auto.hgt;
                    printf("%d %d %g %g \n",i,j,x0,y0); */
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
                /* x=Auto.xmin+(double)(i-Auto.x0)*(Auto.xmax-Auto.xmin)/(double)Auto.wid;
               y=Auto.ymin+(double)(Auto.y0-j+Auto.hgt)*(Auto.ymax-Auto.ymin)/(double)Auto.hgt;
               printf("%d %d %g %g \n",i,j,x,y,x0,y0);
                dx=-(x-x0)/2;
                dy=-(y-y0)/2; */
                dx = (double)(i0 - i)*(Auto.xmax - Auto.xmin) /
                     (double)Auto.wid;
                dy = (double)(j - j0)*(Auto.ymax - Auto.ymin) /
                     (double)Auto.hgt;
                /*    printf("%d %d %d %d %g %g\n",i,j,i0,j0,dx,dy); */
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
    XSetLineAttributes(display, small_gc, wid, ls, cs, js);
    return;
}

void
auto_motion(XEvent ev) {
    int32 i = ev.xmotion.x;
    int32 j = ev.xmotion.y;
    double x, y;
    Window w = ev.xmotion.window;
    if (Auto.exist == 0)
        return;
    if (w == AutoW.canvas) {
        x = Auto.xmin +
            (double)(i - Auto.x0)*(Auto.xmax - Auto.xmin) / (double)Auto.wid;
        y = Auto.ymin + (double)(Auto.y0 - j + Auto.hgt) *
                            (Auto.ymax - Auto.ymin) / (double)Auto.hgt;
        sprintf(Auto.hinttxt, "x=%g,y=%g", x, y);
        storeautopoint(x, y);
        display_auto(AutoW.hint);
    }
    return;
}

void
display_auto(Window w) {
    int32 ix, iy;
    if (Auto.exist == 0)
        return;
    if (w == AutoW.canvas) {
        if (AutoRedrawFlag == 1)
            redraw_diagram();
    }
    if (w == AutoW.stab) {
        int32 r = Auto.st_wid / 4;
        XFlush(display);
        XDrawArc(display, AutoW.stab, small_gc, r, r, 2*r, 2*r, 0,
                 360*64);
        if (CUR_DIAGRAM != NULL) {
            traverse_out(CUR_DIAGRAM, &ix, &iy, 1); /*clr_stab();*/
        }
        XFlush(display);
    }
    if (w == AutoW.axes)
        xds("Axes");
    if (w == AutoW.numerics)
        xds("Numerics");
    if (w == AutoW.grab)
        xds("Grab");
    if (w == AutoW.run)
        xds("Run");
    if (w == AutoW.redraw)
        xds("reDraw");
    if (w == AutoW.clear)
        xds("Clear");
    if (w == AutoW.per)
        xds("Usr period");
    if (w == AutoW.kill)
        xds("Close");
    if (w == AutoW.param)
        xds("Parameter");
    if (w == AutoW.file)
        xds("File");
    if (w == AutoW.abort)
        xds("ABORT");
    if (w == AutoW.hint) {
        XClearWindow(display, w);
        XDrawString(display, w, gc, 8, CURY_OFF, Auto.hinttxt,
                    strlen(Auto.hinttxt));
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

static void
aw(void) {
    XFlush(display);
    sleep(5);
    return;
}

void
make_auto(char *wname, char *iname) {
    int32 x, y, wid, hgt, addwid = 16*DCURX, addhgt = 3.0*DCURY,
                          hinthgt = DCURY + 6;
    Window base = 0;
    int32 dely = DCURY + 5;
    STD_HGT_var = 20*DCURY;
    /*STD_WID_var =1.62*STD_HGT_var;*/
    STD_WID_var = 67*DCURX;
    int32 ymargin = 4*DCURYs, xmargin = 12*DCURXs;
    XTextProperty winname, iconname;
    XSizeHints size_hints;
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
    AutoW.base = base;

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

    XClassHint class_hints;
    class_hints.res_name = "";
    class_hints.res_class = "";

    XSetWMProperties(display, base, &winname, &iconname, NULL, 0, &size_hints,
                     NULL, &class_hints);

    make_icon((char *)auto_bits, auto_width, auto_height, base);

    AutoW.canvas = make_plain_window(base, x, y, STD_WID_var + xmargin,
                                     STD_HGT_var + ymargin, 1);
    XSetWindowBackground(display, AutoW.canvas, MyDrawWinColor);
    XSelectInput(display, AutoW.canvas, MYMASK);

    x = DCURX;
    y = DCURY + STD_HGT_var + ymargin - 8*DCURX;
    AutoW.stab = make_plain_window(base, x, y, 12*DCURX, 12*DCURX, 2);
    Auto.st_wid = 12*DCURX;
    x = DCURX + 2;
    y = 2*DCURY;
    Auto.hgt = STD_HGT_var;
    Auto.wid = STD_WID_var;
    Auto.x0 = 10*DCURXs;
    Auto.y0 = 2*DCURYs;
    AutoW.kill = lil_button(base, 2, 2);
    AutoW.param = lil_button(base, x, y);
    y += dely;
    AutoW.axes = lil_button(base, x, y);
    y += dely;
    AutoW.numerics = lil_button(base, x, y);
    y += dely;
    AutoW.run = lil_button(base, x, y);
    y += dely;
    AutoW.grab = lil_button(base, x, y);
    y += dely;
    AutoW.per = lil_button(base, x, y);
    y += dely;
    AutoW.clear = lil_button(base, x, y);
    y += dely;
    AutoW.redraw = lil_button(base, x, y);
    y += dely;
    AutoW.file = lil_button(base, x, y);

    y += 3*dely;
    AutoW.abort = lil_button(base, x, y);

    y = DCURY + STD_HGT_var + ymargin + 5;
    x = addwid + 5;
    AutoW.info =
        make_plain_window(base, x, y, STD_WID_var + xmargin, addhgt, 2);
    AutoW.hint = make_plain_window(base, x, y + addhgt + 6,
                                   STD_WID_var + xmargin, DCURY + 2, 2);

    draw_bif_axes();
    return;
}

void
resize_auto_window(XEvent ev) {
    int32 wid, hgt, addhgt = 3.5*DCURY;
    STD_HGT_var = 20*DCURY;
    /*STD_WID_var =1.62*STD_HGT_var;*/
    STD_WID_var = 50*DCURX;
    int32 ymargin = 4*DCURYs, xmargin = 12*DCURXs;
    if (ev.xconfigure.window == AutoW.base) {
        wid = ev.xconfigure.width - Auto_extra_wid;
        hgt = ev.xconfigure.height - Auto_extra_hgt;

        addhgt = 3.0*DCURY;

        XResizeWindow(display, AutoW.canvas, wid, hgt);
        Window root;
        int32 xloc;
        int32 yloc;
        uint32 cwid;
        uint32 chgt;
        uint32 cbwid;
        uint32 cdepth;

        XGetGeometry(display, AutoW.canvas, &root, &xloc, &yloc, &cwid, &chgt,
                     &cbwid, &cdepth);

        Auto.hgt = chgt - ymargin;
        Auto.wid = cwid - xmargin;
        /*  printf("%l %d %d %d %d\n", AutoW.info,xloc,yloc+chgt+4,wid,addhgt);
            printf("%l %d %d %d %d\n",
           AutoW.hint,xloc,yloc+chgt+addhgt+10,wid,DCURY+2);  */
        /* XMoveWindow(display,AutoW.info,xloc,yloc+chgt+4); */
        if (TrueColorFlag > 0) {
            XMoveResizeWindow(display, AutoW.info, xloc, yloc + chgt + 4, wid,
                              addhgt);

            XMoveResizeWindow(display, AutoW.hint, xloc,
                              yloc + chgt + addhgt + 10, wid, DCURY + 2);
        }

        int32 ix, iy;

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
    display_auto(AutoW.hint);
    return;
}

void
clear_msg(void) {
    Auto.hinttxt[0] = '\0';
    display_auto(AutoW.hint);
    return;
}

/*  Auto event handlers   */

void
auto_enter(Window w, int32 v) {
    if (Auto.exist == 0)
        return;
    if (w == AutoW.axes) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(1, v);
        return;
    }
    if (w == AutoW.numerics) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(2, v);
        return;
    }
    if (w == AutoW.grab) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(4, v);
        return;
    }
    if (w == AutoW.run) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(3, v);
        return;
    }
    if (w == AutoW.redraw) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(7, v);
        return;
    }
    if (w == AutoW.clear) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(6, v);
        return;
    }
    if (w == AutoW.per) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(5, v);
        return;
    }
    if (w == AutoW.param) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(0, v);
        return;
    }
    if (w == AutoW.kill) {
        XSetWindowBorderWidth(display, w, v);
        return;
    }
    if (w == AutoW.file) {
        XSetWindowBorderWidth(display, w, v);
        a_msg(8, v);
        return;
    }
    return;
}

void
auto_button(XEvent ev) {
    Window w = ev.xbutton.window;
    if (Auto.exist == 0)
        return;
    if (w == AutoW.axes) {
        SBW;
        auto_plot_par();
        return;
    }
    if (w == AutoW.numerics) {
        SBW;
        auto_num_par();
        return;
    }
    if (w == AutoW.grab) {
        SBW;
        auto_grab();
        return;
    }
    if (w == AutoW.run) {
        SBW;
        auto_run();
        return;
    }
    if (w == AutoW.redraw) {
        SBW;
        redraw_diagram();
        return;
    }
    if (w == AutoW.clear) {
        SBW;
        draw_bif_axes();
        return;
    }
    if (w == AutoW.per) {
        SBW;
        auto_per_par();
        return;
    }
    if (w == AutoW.param) {
        SBW;
        auto_params();
        return;
    }
    if (w == AutoW.kill) {
        SBW;
        auto_kill();
        return;
    }
    if (w == AutoW.file) {
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
    XDestroySubwindows(display, AutoW.base);
    XDestroyWindow(display, AutoW.base);
    return;
}

void
auto_keypress(XEvent ev, int32 *used) {
    Window w = ev.xkey.window;
    /*
     int32 maxlen=64;
     char buf[65];
     XComposeStatus comp;
     KeySym ks;  */
    char ks;
    Window w2;
    int32 rev;

    *used = 0;
    if (Auto.exist == 0)
        return;
    XGetInputFocus(display, &w2, &rev);

    if (w == AutoW.base || w == AutoW.canvas || w2 == AutoW.base) {
        *used = 1;
        ks = (char)get_key_press(&ev);
        /* XLookupString(&ev,buf,maxlen,&ks,&comp); */

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

        if (ks == ESC) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
            return;
        }
    }
}
