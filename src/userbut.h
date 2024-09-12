#ifndef _userbut_h_
#define _userbut_h_
#include "integers.h"

#include <X11/Xlib.h>

typedef struct {
    Window w;
    char bname[10];
    int32 com;
} USERBUT;

void user_button_events(XEvent report);
void user_button_press(Window w);
void user_button_draw(Window w);
void user_button_cross(Window w, int32 b);
int32 get_button_info(char *s, char *bname, char *sc);
int32 find_kbs(char *sc);
void add_user_button(char *s);
void create_user_buttons(int32 x0, int32 y0, Window base);

#endif
