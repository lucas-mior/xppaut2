#ifndef aniparse_avi_h
#define aniparse_avi_h
#include "integers.h"
#include "X11/Xlib.h"
#include <stdio.h>

typedef struct {
    int32 nframe, wid, hgt, fps;
    uchar *image;
    int32 cur_frame;
    int32 task;
} AVI_INFO;

AVI_INFO avi_info;
typedef struct {
    int32 flag;
    int32 skip;
    char root[100];
    char filter[256];
    int32 aviflag, filflag;
} MpegSave;

typedef struct {
    int32 n;
    int32 *x, *y, *col;
    int32 i;
} Comet;

typedef struct {
    Comet c;
    int32 type, flag;
    int32 *col, *x1, *y1, *x2, *y2;
    double zcol, zx1, zy1, zx2, zy2, zrad, zval;
    int32 zthick, tfont, tsize, tcolor;
} AniCom;

int32 new_vcr(void);
int32 create_vcr(char *name);
int32 ani_border(Window w, int32 i);
int32 do_ani_events(XEvent ev);
int32 ani_button(Window w);
int32 ani_create_mpeg(void);
int32 ani_expose(Window w);
int32 ani_resize(int32 x, int32 y);
int32 ani_newskip(void);
int32 ani_flip1(int32 n);
int32 ani_flip(void);
int32 ani_disk_warn(void);
int32 getppmbits(Window window, int32 *wid, int32 *hgt, uchar *out);
int32 writeframe(char *filename, Window window, int32 wid, int32 hgt);
int32 ani_zero(void);
int32 get_ani_file(void);
int32 ani_new_file(char *filename);
int32 load_ani_file(FILE *fp);
int32 parse_ani_string(char *s);
int32 set_ani_dimension(char *x1, char *y1, char *x2, char *y2);
int32 add_ani_com(int32 type, char *x1, char *y1, char *x2, char *y2, char *col,
                  char *thick);
int32 init_ani_stuff(void);
int32 free_ani(void);
int32 chk_ani_color(char *s, int32 *index);
int32 add_ani_expr(char *x, int32 *c);
int32 add_ani_rline(AniCom *a, char *x1, char *y1, char *col, char *thick);
int32 reset_comets(void);
int32 roll_comet(AniCom *a, int32 xn, int32 yn, int32 col);
int32 add_ani_comet(AniCom *a, char *x1, char *y1, char *x2, char *y2,
                    char *col, char *thick);
int32 add_ani_line(AniCom *a, char *x1, char *y1, char *x2, char *y2, char *col,
                   char *thick);
int32 add_ani_rect(AniCom *a, char *x1, char *y1, char *x2, char *y2, char *col,
                   char *thick);
int32 add_ani_frect(AniCom *a, char *x1, char *y1, char *x2, char *y2,
                    char *col, char *thick);
int32 add_ani_ellip(AniCom *a, char *x1, char *y1, char *x2, char *y2,
                    char *col, char *thick);
int32 add_ani_fellip(AniCom *a, char *x1, char *y1, char *x2, char *y2,
                     char *col, char *thick);
int32 add_ani_circle(AniCom *a, char *x1, char *y1, char *x2, char *col,
                     char *thick);
int32 add_ani_fcircle(AniCom *a, char *x1, char *y1, char *x2, char *col,
                      char *thick);
int32 add_ani_text(AniCom *a, char *x1, char *y1, char *y2);
int32 add_ani_vtext(AniCom *a, char *x1, char *y1, char *x2, char *y2);
int32 add_ani_settext(AniCom *a, char *x1, char *y1, char *col);
int32 render_ani(void);
int32 set_ani_perm(void);
int32 eval_ani_color(int32 j);
int32 eval_ani_com(int32 j);
int32 set_ani_thick(int32 t);
int32 set_ani_font_stuff(int32 size, int32 font, int32 color);
int32 set_ani_col(int32 j);
int32 xset_ani_col(int32 icol);
int32 ani_rad2scale(double rx, double ry, int32 *ix, int32 *iy);
int32 ani_radscale(double rad, int32 *ix, int32 *iy);
int32 ani_xyscale(double x, double y, int32 *ix, int32 *iy);
int32 draw_ani_comet(int32 j);
int32 draw_ani_line(int32 j);
int32 draw_ani_rline(int32 j);
int32 draw_ani_circ(int32 j);
int32 draw_ani_fcirc(int32 j);
int32 draw_ani_rect(int32 j);
int32 draw_ani_frect(int32 j);
int32 draw_ani_ellip(int32 j);
int32 draw_ani_fellip(int32 j);
int32 draw_ani_text(int32 j);
int32 draw_ani_vtext(int32 j);
int32 tst_pix_draw(void);
int32 read_ani_line(FILE *fp, char *s);
int32 de_space(char *s);

#endif
