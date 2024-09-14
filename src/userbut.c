#include "integers.h"
#include "functions.h"
#include <stdbool.h>

#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "kbs.h"

#define USERBUTCOLOR 24
#define USERBUTMAX 20

int32 nuserbut = 0;

USERBUT userbut[USERBUTMAX];

extern int32 MyAddedButtonColor;
extern Display *display;
extern Window main_win;
extern int32 DCURYs, DCURXs, CURY_OFFs;
extern GC small_gc;

void
user_button_events(XEvent report) {
    switch (report.type) {
    case Expose:
    case MapNotify:
        user_button_draw(report.xany.window);
        break;
    case EnterNotify:
        user_button_cross(report.xcrossing.window, 2);
        break;
    case LeaveNotify:
        user_button_cross(report.xcrossing.window, 1);
        break;
    case ButtonPress:
        user_button_press(report.xbutton.window);
        break;
    }
    return;
}

void
user_button_press(Window w) {
    int32 i;
    for (i = 0; i < nuserbut; i++) {
        if (w == userbut[i].w) {
            run_the_commands(userbut[i].com);
        }
    }
    return;
}

void
draw_all_user_buttons(void) {
    int32 i = 0;
    for (i = 0; i < nuserbut; i++) {
        user_button_draw(userbut[i].w);
    }
    return;
}

void
user_button_draw(Window w) {
    int32 i;
    for (i = 0; i < nuserbut; i++) {
        if (w == userbut[i].w) {
            XDrawString(display, w, small_gc, 5, CURY_OFFs, userbut[i].bname,
                        strlen(userbut[i].bname));
        }
    }
    return;
}

void
user_button_cross(Window w, int32 b) {
    int32 i;
    for (i = 0; i < nuserbut; i++)
        if (w == userbut[i].w) {
            XSetWindowBorderWidth(display, w, b);
            return;
        }
    return;
}

int32
get_button_info(char *s, char *bname, char *sc) {
    int32 i = 0, j = 0, f = 0, n = strlen(s);
    char c;
    if (n == 0)
        return -1;
    bname[0] = 0;
    sc[0] = 0;
    while (true) {
        if (i == n)
            break;
        c = s[i];
        if (c == ':') {
            f = 1;
            bname[j] = 0;
            j = 0;
            i++;
        } else {
            if (f == 0) {
                bname[j] = c;
                j++;
            } else {
                sc[j] = c;
                j++;
            }
            i++;
        }
    }
    sc[j] = 0;

    return 1;
}

int32
find_kbs(char *sc) {
    int32 i = 0;
    while (true) {
        if (strcmp(sc, kbs[i].seq) == 0)
            return kbs[i].com;
        i++;
        if (kbs[i].com == 0)
            return -1;
    }
    exit(EXIT_FAILURE);
}

void
add_user_button(char *s) {
    char bname[10], sc[10];
    int32 z;
    if (nuserbut >= USERBUTMAX)
        return;
    if (strlen(s) == 0)
        return;
    get_button_info(s, bname, sc);
    if (strlen(bname) == 0 || strlen(sc) == 0)
        return;
    z = find_kbs(sc);
    if (z == -1) {
        plintf("%s - not implemented\n", sc);
        return;
    }
    /*Don't add buttons with same functionality twice*/
    int32 i;
    for (i = 0; i < nuserbut; i++) {
        if (userbut[i].com == z) {
            /*		plintf("But=%s:%s already implemented as button
             * '%s'\n",bname,sc,userbut[i].bname); */
            return;
        }
    }
    userbut[nuserbut].com = z;
    strcpy(userbut[nuserbut].bname, bname);
    plintf(" added button(%d)  -- %s %d\n", nuserbut, userbut[nuserbut].bname,
           userbut[nuserbut].com);
    nuserbut++;
    return;
}

void
create_user_buttons(int32 x0, int32 y0, Window base) {
    int32 i;
    int32 x = x0;
    int32 l;
    if (nuserbut == 0)
        return;
    for (i = 0; i < nuserbut; i++) {
        l = DCURXs*(strlen(userbut[i].bname) + 2);
        userbut[i].w = make_fancy_window(base, x, y0, l, DCURYs, 1);
        x = x + l + DCURXs;
    }
    draw_all_user_buttons();
    return;
}
