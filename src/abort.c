#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "functions.h"
#include "integers.h"

extern Window command_pop;
extern GC gc, gc_graph, small_gc;
extern Display *display;

extern int32 DCURYb, DCURXb, CURY_OFFb;
extern int32 DCURYs, DCURXs, CURY_OFFs;
extern int32 DCURY, DCURX, CURY_OFF;

void
plot_command(int32 nit, int32 icount, int32 cwidth) {
    double dx;

    if (nit == 0)
        return;

    dx = (double)icount*(double)cwidth / (double)nit;
    XDrawPoint(display, command_pop, gc, (int32)dx, 5);
    return;
}

int32
my_abort(void) {
    int32 ch;

    while (XPending(display) > 0) {
        XEvent event;
        XNextEvent(display, &event);

        if (check_ani_pause(event) == 27)
            return 27;

        switch (event.type) {
        case Expose:
            do_expose(event);
            break;
        case ButtonPress:
            break;
        case KeyPress:
            ch = get_key_press(&event);
            return ch;
        default:
            break;
        }

        return 0;
    }

    return 64;
}
