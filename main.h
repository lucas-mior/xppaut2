#ifndef _main_h__

#define _main_h__
#include "integers.h"

#include <X11/Xlib.h>

void do_main(int32 argc, char **argv);
void check_for_quiet(int32 argc, char **argv);
void do_vis_env(void);
void init_X(void);
void set_big_font(void);
void set_small_font(void);
void xpp_events(XEvent report, int32 min_wid, int32 min_hgt);
void do_events(uint32 min_wid, uint32 min_hgt);
void bye_bye(void);
void clr_scrn(void);
void redraw_all(void);
void commander(int32 ch);
Window init_win(uint32 bw, char *icon_name, char *win_name, int32 x, int32 y,
                uint32 min_wid, uint32 min_hgt, int32 argc,
                char **argv);
void top_button_draw(Window w);
void top_button_cross(Window w, int32 b);
void top_button_press(Window w);
void top_button_events(XEvent report);
void make_top_buttons(void);
void getGC(GC *gc);
void load_fonts(void);
void make_pops(void);
void FixWindowSize(Window w, int32 width, int32 height, int32 flag);
int32 getxcolors(XWindowAttributes *win_info, XColor **colors);
void test_color_info(void);

#endif
