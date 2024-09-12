#ifndef _graf_par_h_
#define _graf_par_h_
#include "integers.h"

#define RUBBOX 0
#define RUBLINE 1

#define SCRNFMT 0
#define PSFMT 1
#define SVGFMT 2

#define REAL_SMALL 1.e-6
#define MAX_LEN_SBOX 25
#define MAXBIFCRV 100
#define lmax(a, b) ((a > b) ? a : b)

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct {
    char angle[20];
    char yes[3];
    double start;
    double incr;
    int32 nclip;
} MOV3D;

typedef struct {
    float *x[MAXBIFCRV], *y[MAXBIFCRV];
    int32 color[MAXBIFCRV], npts[MAXBIFCRV], nbifcrv;
    Window w;
} BD;

void change_view_com(int32 com);
void ind_to_sym(int32 ind, char *str);
void check_flags(void);
void get_2d_view(int32 ind);
void axes_opts(void);
void get_3d_view(int32 ind);
void check_val(double *x1, double *x2, double *xb, double *xd);
void get_max(int32 index, double *vmin, double *vmax);
void pretty(double *x1, double *x2);
void corner_cube(double *xlo, double *xhi, double *ylo, double *yhi);
void fit_window(void);
void check_windows(void);
void user_window(void);
void xi_vs_t(void);
void redraw_the_graph(void);
void movie_rot(double start, double increment, int32 nclip, int32 angle);
void test_rot(void);
void get_3d_par_com(void);
void get_3d_par_noper(void);
void window_zoom_com(int32 c);
void zoom_in(int32 i1, int32 j1, int32 i2, int32 j2);
void zoom_out(int32 i1, int32 j1, int32 i2, int32 j2);
void graph_all(int32 *list, int32 n, int32 type);
int32 find_color(int32 in);
int32 alter_curve(char *title, int32 in_it, int32 n);
void edit_curve(void);
void new_curve(void);
void create_ps(void);
void change_cmap_com(int32 i);
void freeze_com(int32 c);
void set_key(int32 x, int32 y);
void draw_freeze_key(void);
void key_frz_com(int32 c);
void edit_frz(void);
void delete_frz_crv(int32 i);
void delete_frz(void);
void kill_frz(void);
int32 freeze_crv(int32 ind);
void auto_freeze_it(void);
int32 create_crv(int32 ind);
void edit_frz_crv(int32 i);
void draw_frozen_cline(int32 index, Window w);
void draw_freeze(Window w);
void init_bd(void);
void draw_bd(Window w);
void free_bd(void);
void add_bd_crv(float *x, float *y, int32 len, int32 type, int32 ncrv);
void frz_bd(void);
void read_bd(FILE *fp);
int32 get_frz_index(Window w);
void export_graf_data(void);
void add_a_curve_com(int32 c);
void default_window();
void dump_ps(int32 i);

#endif
