#ifndef _aniparse_h_
#define _aniparse_h_
#include "integers.h"

#include <X11/Xlib.h>
#include <stdio.h>

/**************  New stuff for the Grabber ***************************/
#define MAX_GEVENTS 20  /* maximum variables you can change per grabbable */
#define MAX_ANI_GRAB 50 /* max grabbable objects  */

typedef struct { /* tasks have the form {name1=formula1;name2=formula2;...} */

    double vrhs[MAX_GEVENTS];
    char lhsname[MAX_GEVENTS][11];
    int32 lhsivar[MAX_GEVENTS];
    int32 *comrhs[MAX_GEVENTS];
    int32 runnow;
    int32 n; /* number of tasks <= MAX_GEVENTS */
} GRAB_TASK;

typedef struct {
    int32 ok;
    double zx, zy, tol;
    int32 *x, *y;
    GRAB_TASK start, end;
} ANI_GRAB;

/***************  End of grabber stuff  in header **************/

typedef struct {
    int32 flag;
    int32 skip;
    char root[100];
    char filter[256];
    int32 aviflag, filflag;
} MPEG_SAVE;

typedef struct {
    int32 n;
    int32 *x, *y, *col;
    int32 i;
} Comet;

typedef struct {
    Comet c;
    int32 type, flag;
    int32 *col, *x1, *y1, *x2, *y2, *who;
    double zcol, zx1, zy1, zx2, zy2, zrad, zval;
    int32 zthick, tfont, tsize, tcolor;
} ANI_COM;

void new_vcr(void);
void create_vcr(char *name);
void ani_border(Window w, int32 i);
void destroy_vcr(void);
void do_ani_events(XEvent ev);
void ani_motion_stuff(Window w, int32 x, int32 y);
double get_current_time(void);
void update_ani_motion_stuff(int32 x, int32 y);
void ani_buttonx(XEvent ev, int32 flag);
void ani_button(Window w);
void ani_create_mpeg(void);
void ani_expose(Window w);
void ani_resize(int32 x, int32 y);
void ani_newskip(void);
void check_on_the_fly(void);
void on_the_fly(int32 task);
void ani_frame(int32 task);
void set_to_init_data(void);
void set_from_init_data(void);
void ani_flip1(int32 n);
void ani_flip(void);
void ani_disk_warn(void);
int32 getppmbits(Window window, int32 *wid, int32 *hgt, unsigned char *out);
int32 writeframe(char *filename, Window window, int32 wid, int32 hgt);
void ani_zero(void);
void get_ani_file(char *fname);
int32 ani_new_file(char *filename);
int32 load_ani_file(FILE *fp);
int32 parse_ani_string(char *s, FILE *fp);
void set_ani_dimension(char *x1, char *y1, char *x2, char *y2);
int32 add_ani_com(int32 type, char *x1, char *y1, char *x2, char *y2, char *col,
                  char *thick);
void init_ani_stuff(void);
void free_ani(void);
int32 chk_ani_color(char *s, int32 *index);
int32 add_ani_expr(char *x, int32 *c);
int32 add_ani_rline(ANI_COM *a, char *x1, char *y1, char *col, char *thick);
void reset_comets(void);
void roll_comet(ANI_COM *a, int32 xn, int32 yn, int32 col);
int32 add_ani_comet(ANI_COM *a, char *x1, char *y1, char *x2, char *col,
                    char *thick);
int32 add_ani_line(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                   char *col, char *thick);
int32 add_ani_null(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                   char *col, char *who);
int32 add_ani_rect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                   char *col, char *thick);
int32 add_ani_frect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                    char *col, char *thick);
int32 add_ani_ellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                    char *col, char *thick);
int32 add_ani_fellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
                     char *col, char *thick);
int32 add_ani_circle(ANI_COM *a, char *x1, char *y1, char *x2, char *col,
                     char *thick);
int32 add_ani_fcircle(ANI_COM *a, char *x1, char *y1, char *x2, char *col,
                      char *thick);
int32 add_ani_text(ANI_COM *a, char *x1, char *y1, char *y2);
int32 add_ani_vtext(ANI_COM *a, char *x1, char *y1, char *x2, char *y2);
int32 add_ani_settext(ANI_COM *a, char *x1, char *y1, char *col);
void render_ani(void);
void set_ani_perm(void);
void eval_ani_color(int32 j);
void eval_ani_com(int32 j);
void set_ani_thick(int32 t);
void set_ani_font_stuff(int32 size, int32 font, int32 color);
void set_ani_col(int32 j);
void xset_ani_col(int32 icol);
void ani_rad2scale(double rx, double ry, int32 *ix, int32 *iy);
void ani_radscale(double rad, int32 *ix, int32 *iy);
void ani_ij_to_xy(int32 ix, int32 iy, double *x, double *y);
void ani_xyscale(double x, double y, int32 *ix, int32 *iy);
void draw_ani_comet(int32 j);
void draw_ani_null(int32 j, int32 id);
void draw_ani_line(int32 j);
void draw_ani_rline(int32 j);
void draw_ani_circ(int32 j);
void draw_ani_fcirc(int32 j);
void draw_ani_rect(int32 j);
void draw_ani_frect(int32 j);
void draw_ani_ellip(int32 j);
void draw_ani_fellip(int32 j);
void draw_ani_text(int32 j);
void draw_ani_vtext(int32 j);
void tst_pix_draw(void);
void read_ani_line(FILE *fp, char *s);
void de_space(char *s);
int32 add_grab_command(char *xs, char *ys, char *ts, FILE *fp);
void info_grab_stuff(void);
int32 ani_grab_tasks(char *line, int32 igrab, int32 which);
int32 run_now_grab(void);
int32 search_for_grab(double x, double y);
void do_grab_tasks(int32 which);
int32 add_grab_task(char *lhs, char *rhs, int32 igrab, int32 which);
void draw_grab_points(void);
void free_grabber(void);
int32 check_ani_pause(XEvent ev);
void do_ani_slider_motion(Window w, int32 x);
void draw_ani_slider(Window w, int32 x);
void redraw_ani_slider(void);
#endif
