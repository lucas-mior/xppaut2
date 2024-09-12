#ifndef _many_pops_h
#define _many_pops_h
#include "integers.h"

#include <X11/Xlib.h>

int32 select_table(void);
void get_intern_set(void);
void make_icon(char *icon, int32 wid, int32 hgt, Window w);
void title_text(char *string);
void gtitle_text(char *string, Window win);
void restore_off(void);
void restore_on(void);
void add_label(char *s, int32 x, int32 y, int32 size, int32 font);
void draw_marker(double x, double y, double size, int32 type);
void draw_grob(int32 i);
void arrow_head(double xs, double ys, double xe, double ye, double size);
void destroy_grob(Window w);
void destroy_label(Window w);
void draw_label(Window w);
void add_grob(double xs, double ys, double xe, double ye, double size,
              int32 type, int32 color);
int32 select_marker_type(int32 *type);
int32 man_xy(float *xe, float *ye);
int32 get_marker_info(void);
int32 get_markers_info(void);
void add_marker(void);
void add_marker_old(void);
void add_markers(void);
void add_markers_old(void);
void add_pntarr(int32 type);
void edit_object_com(int32 com);
void do_gr_objs_com(int32 com);
void do_windows_com(int32 c);
void set_restore(int32 flag);
int32 is_col_plotted(int32 nc);
void destroy_a_pop(void);
void init_grafs(int32 x, int32 y, int32 w, int32 h);
void ps_restore(void);
void svg_restore(void);
int32 rotate3dcheck(XEvent ev);
void do_motion_events(XEvent ev);
void do_expose(XEvent ev);
void resize_all_pops(int32 wid, int32 hgt);
void kill_all_pops(void);
void create_a_pop(void);
void GrCol(void);
void BaseCol(void);
void SmallGr(void);
void SmallBase(void);
void change_plot_vars(int32 k);
int32 check_active_plot(int32 k);
int32 graph_used(int32 i);
void make_active(int32 i, int32 flag);
void select_window(Window w);
void set_gr_fore(void);
void set_gr_back(void);
void hi_lite(Window wi);
void lo_lite(Window wi);
void select_sym(Window w);
void canvas_xy(char *buf);
void check_draw_button(XEvent ev);
void set_active_windows(void);

#endif
