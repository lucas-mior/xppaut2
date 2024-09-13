#ifndef _ggets_h
#define _ggets_h
#include "integers.h"

#include <X11/Xlib.h>

#define MaxIncludeFiles 10
#define ClickTime 200

void ping(void);
void reset_graphics(void);
void blank_screen(Window w);
void set_fore(void);
void set_back(void);
void showchar(int32 ch, int32 col, int32 row, Window or);
void chk_xor(void);
void set_gcurs(int32 y, int32 x);
void clr_command(void);
void draw_info_pop(Window win);
void bottom_msg(int32 line, char *msg);
void gputs(char *string, Window win);
void err_msg(char *string);
int32 plintf(char *fmt, ...);
int32 show_position(XEvent ev, int32 *com);
void gpos_prn(char *string, int32 row, int32 col);
void put_command(char *string);
int32 get_key_press(XEvent *ev);
void cput_text(void);
int32 get_mouse_xy(int32 *x, int32 *y, Window w);
void Ftext(int32 x, int32 y, char *string, Window o);
void bar(int32 x, int32 y, int32 x2, int32 y2, Window w);
void rectangle(int32 x, int32 y, int32 x2, int32 y2, Window w);
void setfillstyle(int32 type, int32 color);
void circle(int32 x, int32 y, int32 radius, Window w);
void xline(int32 x0, int32 y0, int32 x1, int32 y1, Window w);
int32 new_float(char *name, double *value);
int32 new_int(char *name, int32 *value);
void display_command(char *name, char *value, int32 pos, int32 col);
void clr_line_at(Window w, int32 col0, int32 pos, int32 n);
void put_cursor_at(Window w, int32 col0, int32 pos);
void put_string_at(Window w, int32 col, char *s, int32 off);
void movmem(char *s1, char *s2, int32 len);
void memmov(char *s1, char *s2, int32 len);
void edit_window(Window w, int32 *pos, char *value, int32 *col, int32 *done,
                 int32 ch);
void do_backspace(int32 *pos, char *value, int32 *col, Window w);
void edit_command_string(XEvent ev, char *name, char *value, int32 *done,
                         int32 *pos, int32 *col);
int32 new_string(char *name, char *value);

#endif
