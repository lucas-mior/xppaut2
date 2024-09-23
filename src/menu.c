#include "menus.h"
#include "functions.h"
#include "integers.h"

#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>

int32 help_menu;
MENUDEF my_menus[3];
extern int32 tfBell;
extern int32 TipsFlag;
extern int32 DCURY;
extern int32 DCURX;
extern int32 CURY_OFF;
extern int32 DCURYs;
extern int32 DCURYb;
extern GC gc;

static void unshow_menu(int32 j);
static void show_menu(int32 j);

void
add_menu(Window base, int32 j, int32 n, char **names, char *key, char **hint) {
    Window window;
    int32 i;
    Cursor cursor;
    cursor = XCreateFontCursor(display, XC_hand2);
    window = make_plain_unmapped_window(base, 0, DCURYs + DCURYb + 10,
                                        16*DCURX, 21*(DCURY + 2) - 3, 1);
    my_menus[j].base = window;
    XDefineCursor(display, window, cursor);
    my_menus[j].names = names;
    my_menus[j].n = n;
    my_menus[j].hints = hint;
    strcpy(my_menus[j].key, key);
    my_menus[j].title =
        make_unmapped_window(window, 0, 0, 16*DCURX, DCURY, 1);
    for (i = 0; i < n; i++) {
        my_menus[j].window[i] = make_unmapped_window(
            window, 0, (i + 1)*(DCURY + 2), 16*DCURX, DCURY, 0);
    }
    my_menus[j].visible = 0;
    XMapRaised(display, my_menus[j].base);
    XMapSubwindows(display, my_menus[j].base);
    return;
}

void
create_the_menus(Window base) {
    char key[30];
    strcpy(key, "icndwakgufpemtsvxr3b");
    add_menu(base, MAIN_MENU, MAIN_ENTRIES, main_menu, key, main_hint);
    strcpy(key, "tsrdniobmechpukva");
    key[17] = 27;
    key[18] = 0;
    add_menu(base, NUM_MENU, NUM_ENTRIES, num_menu, key, num_hint);
    /* CLONE */
    strcpy(key, "pwracesbhqtiglxu");
    add_menu(base, FILE_MENU, FILE_ENTRIES, fileon_menu, key, file_hint);
    help_menu = -1;
    return;
}

void
show_menu(int32 j) {
    /*  XMapRaised(display,my_menus[j].base);
    XMapSubwindows(display,my_menus[j].base);
    */
    XRaiseWindow(display, my_menus[j].base);
    my_menus[j].visible = 1;
    help_menu = j;
    return;
}

void
unshow_menu(int32 j) {
    if (j < 0)
        return;
    my_menus[j].visible = 0;
    /* XUnmapSubwindows(display,my_menus[j].base);
     XUnmapWindow(display,my_menus[j].base); */
    return;
}

void
help(void) {
    unshow_menu(help_menu);
    show_menu(MAIN_MENU);
    return;
}

void
help_num(void) {
    unshow_menu(help_menu);
    show_menu(NUM_MENU);
    return;
}

void
help_file(void) {
    if (tfBell)
        my_menus[FILE_MENU].names = fileon_menu;
    else
        my_menus[FILE_MENU].names = fileoff_menu;
    unshow_menu(help_menu);
    show_menu(FILE_MENU);
    return;
}

void
menu_crossing(Window win, int32 yn) {
    int32 i, n, j = help_menu;
    char **z;
    if (j < 0)
        return;
    if (my_menus[j].visible == 0)
        return;
    n = my_menus[j].n;
    z = my_menus[j].hints;
    for (i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            XSetWindowBorderWidth(display, win, (uint)yn);
            if (yn && TipsFlag)
                bottom_msg(z[i]);
            return;
        }
    }
    return;
}

void
menu_expose(Window win) {
    int32 i, n, j = help_menu;
    char **z;
    if (j < 0)
        return;
    if (my_menus[j].visible == 0)
        return;
    n = my_menus[j].n;
    z = my_menus[j].names;
    if (win == my_menus[j].title) {
        ggets_set_fore();
        bar(0, 0, 16*DCURX, DCURY, win);
        set_back();
        XDrawString(display, win, gc, DCURX / 2 + 5, CURY_OFF, z[0],
                    (int)strlen(z[0]));
        ggets_set_fore();
        /* BaseCol();
        XDrawString(display,win,gc,0,CURY_OFF,z[0],strlen(z[0]));
        */
        return;
    }
    for (i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            BaseCol();
            XDrawString(display, win, gc, 5, CURY_OFF, z[i + 1],
                        (int)strlen(z[i + 1]));
            return;
        }
    }
    return;
}

void
menu_button(Window win) {
    int32 i, n, j = help_menu;
    if (j < 0)
        return;
    if (my_menus[j].visible == 0)
        return;
    n = my_menus[j].n;
    for (i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            XSetWindowBorderWidth(display, win, 0);
            commander(my_menus[j].key[i]);
            return;
        }
    }
    return;
}

void
draw_help(void) {
    int32 i, j = help_menu, n;
    /*char **z;
     */
    if (j < 0)
        return;
    if (my_menus[j].visible == 0)
        return;
    n = my_menus[j].n;
    /*z=my_menus[j].names;
     */
    menu_expose(my_menus[j].title);
    for (i = 0; i < n; i++)
        menu_expose(my_menus[j].window[i]);
}
/*
menu_events(ev)
     XEvent event;
{
  switch(ev.type){
  case Expose:
  case MapNotify:
    menu_expose(ev.xexpose.window);
    break;
  case ButtonPress:
    menu_button(ev.xbutton.window);
    break;
  case EnterNotify:
    menu_crossing(ev.xcrossing.window,1);
    break;
  case LeaveNotify:
    menu_crossing(ev.xcrossing.window,0);
    break;
  }
}
*/
