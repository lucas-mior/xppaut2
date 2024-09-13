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
    int32 i;
    float dx;
    if (nit == 0)
        return;
    dx = (float)icount * (float)cwidth / (float)nit;
    i = (int32)dx;

    XDrawPoint(display, command_pop, gc, i, 5);
    return;
}

int32
my_abort(void) {
    int32 ch;
    XEvent event;
    while (XPending(display) > 0) {
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
        }
        return 0;
    }

    return 64;
}
