#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

#include "integers.h"

#include "functions.h"

extern Window draw_win;
extern int32 screen;
extern int32 xor_flag;
extern int32 xorfix;
extern GC gc, gc_graph;
extern uint32 MyBackColor;
extern uint32 MyForeColor;
extern uint32 MyMainWinColor;
extern uint32 MyDrawWinColor;
extern uint32 GrFore, GrBack;

static void rbox(int32 i1, int32 j1, int32 i2, int32 j2, Window window,
                 int32 f);

int32
rubber(int32 *x1, int32 *y1, int32 *x2, int32 *y2, Window window, int32 f) {
    XEvent ev;
    int32 there = 0;
    int32 error = 0;
    int32 dragx = 0, dragy = 0;
    int32 oldx = 0, oldy = 0;
    int32 state = 0;
    xor_flag = 1;
    XFlush(display);
    chk_xor();
    if (xorfix) {
        XSetForeground(display, gc, MyDrawWinColor);
        XSetBackground(display, gc, MyForeColor);
        /*XSetForeground(display,gc,GrFore);*/
    }

    XSelectInput(display, window,
                 KeyPressMask | ButtonPressMask | ButtonReleaseMask |
                     PointerMotionMask | ButtonMotionMask | ExposureMask);
    while (!there) {
        XNextEvent(display, &ev);
        switch (ev.type) {
        case Expose:
            do_expose(ev);
            xor_flag = 1;
            chk_xor();
            if (xorfix) {
                XSetForeground(display, gc, MyDrawWinColor);
                XSetBackground(display, gc, MyForeColor);
                /*XSetForeground(display,gc,GrFore);*/
            }
            break;

        case KeyPress:
            if (state > 0)
                break; /* too late Bozo   */
            there = 1;
            error = 1;
            break;
        case ButtonPress:
            if (state > 0)
                break;
            state = 1;
            dragx = ev.xkey.x;
            dragy = ev.xkey.y;
            oldx = dragx;
            oldy = dragy;
            rbox(dragx, dragy, oldx, oldy, window, f);
            break;
        case MotionNotify:
            if (state == 0)
                break;
            rbox(dragx, dragy, oldx, oldy, window, f);
            oldx = ev.xmotion.x;
            oldy = ev.xmotion.y;
            rbox(dragx, dragy, oldx, oldy, window, f);
            break;
        case ButtonRelease:
            if (state == 0)
                break;
            there = 1;
            rbox(dragx, dragy, oldx, oldy, window, f);
            break;
        default:
            break;
        }
    }
    xor_flag = 0;
    chk_xor();

    if (xorfix) {
        /*XSetForeground(display,gc,GrBack); */
        XSetForeground(display, gc, MyForeColor);
        XSetBackground(display, gc, MyDrawWinColor);
    }

    if (!error) {
        *x1 = dragx;
        *y1 = dragy;
        *x2 = oldx;
        *y2 = oldy;
    }

    XSelectInput(display, window,
                 KeyPressMask | ButtonPressMask | ExposureMask |
                     ButtonReleaseMask | ButtonMotionMask);
    if (error)
        return 0;
    return 1;
}

void
rbox(int32 i1, int32 j1, int32 i2, int32 j2, Window window, int32 f) {
    int32 x1 = i1, x2 = i2, y1 = j1, y2 = j2;
    if (f == RUBLINE) {
        XDrawLine(display, window, gc, i1, j1, i2, j2);
        return;
    }
    if (x1 > x2) {
        x2 = i1;
        x1 = i2;
    }
    if (y1 > y2) {
        y1 = j2;
        y2 = j1;
    }
    rectangle(x1, y1, x2, y2, window);
    return;
}
