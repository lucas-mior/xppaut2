#ifndef _eig_list_h_
#define _eig_list_h_

#include <X11/Xlib.h>
#include "integers.h"

void draw_eq_list(Window w);
void create_eq_list(void);
void eq_list_keypress(XEvent ev, int32 *used);
void enter_eq_stuff(Window w, int32 b);
void eq_list_button(XEvent ev);
void eq_list_up(void);
void eq_list_down(void);
void eq_box_import(void);
void get_new_size(Window win, uint32 *wid, uint32 *hgt);
void resize_eq_list(Window win);
void eq_box_button(Window w);
void create_eq_box(int32 cp, int32 cm, int32 rp, int32 rm, int32 im, double *y,
                   double *ev, int32 n);
void draw_eq_box(Window w);

#endif
