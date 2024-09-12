#ifndef _torus_h_
#define _torus_h_
#include "integers.h"

#include <X11/Xlib.h>

void do_torus_com(int32 c);
void draw_tor_var(int32 i);
void draw_torus_box(Window win);
void choose_torus(void);
void make_tor_box(char *title);
void do_torus_events(void);

#endif
