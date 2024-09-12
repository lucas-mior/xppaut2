#ifndef _auto_x11_h_
#define _auto_x11_h_
#include "integers.h"

#include <X11/Xlib.h>

typedef struct {
    Window canvas, axes, numerics, grab, next, run, clear, redraw, base, per;
    Window info, param, file, abort, stab, hint, kill;
} AUTOWIN;

void ALINE(int32 a, int32 b, int32 c, int32 d);
void DLINE(double a, double b, double c, double d);
void ATEXT(int32 a, int32 b, char *c);
void clr_stab(void);
void auto_stab_line(int32 x, int32 y, int32 xp, int32 yp);
void clear_auto_plot(void);
void redraw_auto_menus(void);
void traverse_diagram(void);
void clear_auto_info(void);
void draw_auto_info(char *bob, int32 x, int32 y);
void refreshdisplay(void);
int32 byeauto_(int32 *iflag);
void Circle(int32 x, int32 y, int32 r);
void autocol(int32 col);
void autobw(void);
int32 auto_rubber(int32 *i1, int32 *j1, int32 *i2, int32 *j2, int32 flag);
int32 auto_pop_up_list(char *title, char **list, char *key, int32 n, int32 max,
                       int32 def, int32 x, int32 y, char **hints, char *httxt);
void MarkAuto(int32 x, int32 y);
void XORCross(int32 x, int32 y);
void FillCircle(int32 x, int32 y, int32 r);
void LineWidth(int32 wid);
void auto_motion(XEvent ev);
void display_auto(Window w);
Window lil_button(Window root, int32 x, int32 y, char *name);
void make_auto(char *wname, char *iname);
void resize_auto_window(XEvent ev);
void a_msg(int32 i, int32 v);
void auto_enter(Window w, int32 v);
void auto_button(XEvent ev);
void auto_kill(void);
void auto_keypress(XEvent ev, int32 *used);
int32 query_special(char *title, char *nsymb);
void clear_msg();
void find_point(int32 ibr, int32 pt);
void auto_get_info(int32 *n, char *pname);
void auto_set_mark(int32 i);
void do_auto_range();
void RedrawMark();

#endif
