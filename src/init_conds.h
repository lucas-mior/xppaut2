#ifndef _init_conds_h_
#define _init_conds_h_
#include "integers.h"

#include <X11/Xlib.h>
#include "read_dir.h"

#define FILESELNWIN 10
typedef struct {
    int32 n, n0, here;
    Window base, cancel, ok, up, dn, pgup, pgdn, file, wild, w[FILESELNWIN],
        dir, home, start;
    Window fw, ww;
    char wildtxt[256], filetxt[256];
    int32 nwin, minwid, minhgt;
    int32 off, pos, hot;
    char title[256];
} FILESEL;

typedef struct {
    int32 pos, n, n0, npos;
    int32 ihot, twid;
    int32 max;
    char **v;
    Window side, up, down, text;
} SCROLL_LIST;

typedef struct {
    int32 use, pos, l;
    char parname[20];
    double lo, hi, val;
    int32 hgt;
    int32 type, index;
    Window left, right, top, main, slide, go;
} PAR_SLIDER;

typedef struct {
    int32 use, type;
    int32 n;
    Window base;
    Window cancel, ok, def, go;
    Window *w;
    Window *we;
    char **value;
    int32 mc, *off, *pos;
} BoxListold;

typedef struct {
    int32 use, type, xuse;
    int32 n, n0;
    int32 nwin, minwid, minhgt;
    Window up, dn;
    Window pgup, pgdn;
    Window base;
    Window cancel, ok, def, go, close;
    Window xvt, pp, arr;
    Window *w;
    Window *we;
    Window *ck;
    char **value, *iname, *wname;
    int32 *isck;
    int32 mc, *off, *pos;
} BoxList;

void create_scroll_list(Window base, int32 x, int32 y, int32 width,
                        int32 height, SCROLL_LIST *sl);
void free_scroll_list(SCROLL_LIST *sl);
void add_scroll_item(char *v, SCROLL_LIST *sl);
int32 expose_scroll_list(Window w, SCROLL_LIST sl);
void redraw_scroll_list(SCROLL_LIST sl);
void c_hints(void);
void clone_ode(void);
int32 find_user_name(int32 type, char *oname);
void create_par_sliders(Window base, int32 x0, int32 h0);
void resize_par_slides(int32 h);
void slide_button_press(Window w);
void do_slide_button(int32 w, PAR_SLIDER *p);
void expose_selector(Window w);
void redraw_directory(void);
void redraw_file_list(void);
void redraw_fs_text(char *string, Window w, int32 flag);
void display_file_sel(FILESEL f, Window w);
void new_wild(void);
void fs_scroll(int32 i);
int32 button_selector(Window w);
void crossing_selector(Window w, int32 c);
int32 do_file_select_events(void);
void create_file_selector(char *title, char *file, char *wild);
void stringintersect(char *target, char *sother);
int32 edit_fitem(int32 ch, char *string, Window w, int32 *off1, int32 *pos1,
                 int32 mc);
int32 selector_key(XEvent ev);
void destroy_selector(void);
int32 file_selector(char *title, char *file, char *wild);
void reset_sliders(void);
void redraw_slide(PAR_SLIDER *p);
void set_slide_pos(PAR_SLIDER *p);
void slide_release(Window w);
void do_slide_release(int32 w, PAR_SLIDER *p);
void slider_motion(XEvent ev);
void do_slide_motion(Window w, int32 x, PAR_SLIDER *p, int32 state);
void enter_slides(Window w, int32 val);
void enter_slider(Window w, PAR_SLIDER *p, int32 val);
void expose_slides(Window w);
void expose_slider(Window w, PAR_SLIDER *p);
void draw_slider(Window w, int32 x, int32 hgt, int32 l);
void make_par_slider(Window base, int32 x, int32 y, int32 width, int32 index);
void make_new_ic_box(void);
void make_new_bc_box(void);
void make_new_delay_box(void);
void make_new_param_box(void);
void initialize_box(void);
void resize_par_box(Window win);
void get_nrow_from_hgt(int32 h, int32 *n, int32 *w);
void destroy_box(BoxList *b);
void make_box_list_window(BoxList *b, int32 type);
void make_box_list(BoxList *b, char *wname, char *iname, int32 n, int32 type,
                   int32 use);
void do_box_expose(Window w);
void justify_string(Window w1, char *s1);
void draw_one_box(BoxList b, int32 index);
void redraw_params(void);
void redraw_delays(void);
void redraw_ics(void);
void redraw_bcs(void);
void display_box(BoxList b, Window w);
void box_enter_events(Window w, int32 yn);
void box_enter(BoxList b, Window w, int32 val);
int32 find_the_box(BoxList b, Window w, int32 *index);
void set_up_xvt(void);
void set_up_pp(void);
void set_up_arry(void);
void redraw_entire_box(BoxList *b);
void do_box_button(BoxList *b, Window w);
void box_list_scroll(BoxList *b, int32 i);
void box_buttons(Window w);
void box_keypress(XEvent ev, int32 *used);
void do_box_key(BoxList *b, XEvent ev, int32 *used);
void man_ic(void);
void new_parameter(void);
void redo_stuff(void);
void set_default_ics(void);
void set_default_params(void);
void draw_editable(Window win, char *string, int32 off, int32 cursor, int32 mc);
void put_edit_cursor(Window w, int32 pos);
int32 edit_bitem(BoxList *b, int32 i, int32 ch);
void add_edit_float(BoxList *b, int32 i, double z);
void set_edit_params(BoxList *b, int32 i, char *string);
void add_editval(BoxList *b, int32 i, char *string);
void check_box_cursor(void);
void prt_focus(void);
int32 to_float(char *s, double *z);
void set_value_from_box(BoxList *b, int32 i);
void load_entire_box(BoxList *b);

#endif
