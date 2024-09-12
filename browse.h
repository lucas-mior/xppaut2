#ifndef _browse_h_
#define _browse_h_
#include "integers.h"
#include "max_len_sbox.h"

#define BMAXCOL 20

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct {
    Window base, upper;
    Window find, up, down, pgup, pgdn, home, end, left, right;
    Window first, last, restore, write, get, close;
    Window load, repl, unrepl, table, addcol, delcol;
    Window main;
    Window label[BMAXCOL];
    Window time;
    Window hint;
    char hinttxt[256];
    int32 dataflag, xflag;
    int32 col0, row0, ncol, nrow;
    int32 maxrow, maxcol;
    float **data;
    int32 istart, iend;
} BROWSER;

/*extern BROWSER my_browser;
 */

float **get_browser_data(void);
void set_browser_data(float **data, int32 col0);
float *get_data_col(int32 c);
int32 gettimenow(void);
void waitasec(int32 msec);
int32 get_maxrow_browser(void);
void write_mybrowser_data(FILE *fp);
void write_browser_data(FILE *fp, BROWSER *b);
int32 check_for_stor(float **data);
void del_stor_col(char *var, BROWSER *b);
void data_del_col(BROWSER *b);
void data_add_col(BROWSER *b);
int32 add_stor_col(char *name, char *formula, BROWSER *b);
void chk_seq(char *f, int32 *seq, double *a1, double *a2);
void replace_column(char *var, char *form, float **dat, int32 n);
void wipe_rep(void);
void unreplace_column(void);
void make_d_table(double xlo, double xhi, int32 col, char *filename, BROWSER b);
void find_value(int32 col, double val, int32 *row, BROWSER b);
void find_variable(char *s, int32 *col);
void browse_but_on(BROWSER *b, int32 i, Window w, int32 yn);
void enter_browser(XEvent ev, BROWSER *b, int32 yn);
void display_browser(Window w, BROWSER b);
void redraw_browser(BROWSER b);
void new_browse_dat(float **new_dat, int32 dat_len);
void refresh_browser(int32 length);
void reset_browser(void);
void draw_data(BROWSER b);
void init_browser(void);
void kill_browser(BROWSER *b);
void make_new_browser(void);
Window br_button(Window root, int32 row, int32 col, char *name, int32 iflag);
Window br_button_data(Window root, int32 row, int32 col, char *name,
                      int32 iflag);
void make_browser(BROWSER *b, char *wname, char *iname, int32 row, int32 col);
void expose_my_browser(XEvent ev);
void enter_my_browser(XEvent ev, int32 yn);
void my_browse_button(XEvent ev);
void my_browse_keypress(XEvent ev, int32 *used);
void resize_my_browser(Window win);
void expose_browser(XEvent ev, BROWSER b);
void resize_browser(Window win, BROWSER *b);
void browse_button(XEvent ev, BROWSER *b);
void browse_keypress(XEvent ev, int32 *used, BROWSER *b);
void data_up(BROWSER *b);
void data_down(BROWSER *b);
void data_pgup(BROWSER *b);
void data_pgdn(BROWSER *b);
void data_home(BROWSER *b);
void data_end(BROWSER *b);
void get_data_xyz(float *x, float *y, float *z, int32 i1, int32 i2, int32 i3,
                  int32 off);
void data_get(BROWSER *b);
void data_replace(BROWSER *b);
void data_unreplace(BROWSER *b);
void data_table(BROWSER *b);
void data_find(BROWSER *b);
void open_write_file(FILE **fp, char *fil, int32 *ok);
void data_read(BROWSER *b);
void data_write(BROWSER *b);
void data_left(BROWSER *b);
void data_right(BROWSER *b);
void data_first(BROWSER *b);
void data_last(BROWSER *b);
void data_restore(BROWSER *b);
void get_col_list(char *s, int32 *cl, int32 *n);

#endif
