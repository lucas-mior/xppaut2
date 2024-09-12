#ifndef _txtread_h
#define _txtread_h
#include "integers.h"

void txt_view_events(XEvent ev);
void txtview_keypress(XEvent ev);
void enter_txtview(Window w, int32 val);
void do_txt_action(char *s);
void resize_txtview(int32 w, int32 h);
void txtview_press(Window w, int32 x, int32 y);
void redraw_txtview(Window w);
void redraw_txtview_text(void);
void init_txtview(void);
void make_txtview(void);

#endif
