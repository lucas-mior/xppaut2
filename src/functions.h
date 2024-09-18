#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "integers.h"

#define MAIN_MENU 0
#define FILE_MENU 1
#define NUM_MENU 2

#define SETVAR(i, x)                                                           \
    do {                                                                       \
        if ((i) < NVAR)                                                        \
            variables[(i)] = (x);                                              \
    } while (0)
#define GETVAR(i) (i) < NVAR ? variables[(i)] : 0.0

#define MAX_LEN_SBOX 30

#ifndef ABORT_H
#define ABORT_H

void plot_command(int32 nit, int32 icount, int32 cwidth);
int32 my_abort(void);

#endif

#ifndef ADJ2_H
#define ADJ2_H

#include <stdio.h>

extern int32 AdjRange;

void init_trans(void);
void dump_transpose_info(FILE *fp, int32 f);
int32 do_transpose(void);
int32 create_transpose(void);
void alloc_h_stuff(void);
void data_back(void);
void adj_back(void);
void h_back(void);
void make_adj_com(int32 com);
void adjoint_parameters(void);
void new_h_fun(int32 silent);
void dump_h_stuff(FILE *fp, int32 f);
int32 make_h(double **orb, double **adj, int32 nt, int32 node, int32 silent);
void new_adjoint(void);
void test_test(void);
void compute_one_orbit(double *ic, double per);
int32 adjoint(double **orbit, double **adjnt, int32 nt, double dt, double eps,
              double minerr, int32 maxit, int32 node);
void eval_rhs(double **jac, int32 k1, int32 k2, double t, double *y, double *yp,
              int32 node);
int32 rk_interp(double **jac, int32 k1, int32 k2, double *y, double *work,
                int32 neq, double del, int32 nstep);
int32 step_eul(double **jac, int32 k, int32 k2, double *yold, double *work,
               int32 node, double dt);
void do_liapunov(void);
void alloc_liap(int32 n);
void do_this_liaprun(int32 i, double p);
void norm_vec(double *v, double *mu, int32 n);
int32 hrw_liapunov(double *liap, int32 batch, double eps);
void new_adjoint(void);

#endif

#ifndef ANIPARSE_H
#define ANIPARSE_H

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
} GrabTask;

typedef struct {
    int32 ok;
    double zx, zy, tol;
    int32 *x, *y;
    GrabTask start, end;
} AniGrab;

/***************  End of grabber stuff  in header **************/

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
    int32 *col, *x1, *y1, *x2, *y2, *who;
    double zcol, zx1, zy1, zx2, zy2, zrad, zval;
    int32 zthick, tfont, tsize, tcolor;
} AniCom;

extern int32 animation_on_the_fly;

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
int32 getppmbits(Window window, int32 *wid, int32 *hgt, uchar *out);
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
int32 add_ani_rline(AniCom *a, char *x1, char *y1, char *col, char *thick);
void reset_comets(void);
void roll_comet(AniCom *a, int32 xn, int32 yn, int32 col);
int32 add_ani_comet(AniCom *a, char *x1, char *y1, char *x2, char *col,
                    char *thick);
int32 add_ani_line(AniCom *a, char *x1, char *y1, char *x2, char *y2, char *col,
                   char *thick);
int32 add_ani_null(AniCom *a, char *x1, char *y1, char *x2, char *y2, char *col,
                   char *who);
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

#ifndef ARRAY_PLOT_H
#define ARRAY_PLOT_H

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct {
    Window base, wclose, wedit, wprint, wstyle, wscale, wmax, wmin, wplot,
        wredraw, wtime, wgif, wrange, wfit;
    int32 index0, indexn, alive, nacross, ndown, plotdef;
    int32 height, width, ploth, plotw;
    int32 nstart, nskip, ncskip;
    char name[20];
    double tstart, tend, zmin, zmax, dt;
    char xtitle[256], ytitle[256], filename[256], bottom[256];
    int32 type;
} ArrayPlot;

extern int32 aplot_range;

void set_acolor(int32);
void tag_aplot(char *);
void close_aplot_files(void);
void draw_one_array_plot(char *);
void gif_aplot_all(char *, int32);
void optimize_aplot(int32 *plist);
void make_my_aplot(char *name);
void scale_aplot(ArrayPlot *ap, double *zmax, double *zmin);
void init_array_plot(ArrayPlot *ap);
void expose_aplot(Window w);
void do_array_plot_events(XEvent ev);
void wborder(Window w, int32 i, ArrayPlot ap);
void destroy_aplot(void);
void init_my_aplot(void);
void create_array_plot(ArrayPlot *ap, char *wname, char *iname);
void print_aplot(ArrayPlot *ap);
void apbutton(Window w);
void draw_scale(ArrayPlot ap);
void draw_aplot(ArrayPlot ap);
void edit_aplot(void);
void get_root(char *s, char *sroot, int32 *num);
void reset_aplot_axes(ArrayPlot ap);
void dump_aplot(FILE *fp, int32 f);
int32 editaplot(ArrayPlot *ap);
void gif_aplot(void);
void grab_aplot_screen(ArrayPlot ap);
void redraw_aplot(ArrayPlot ap);
void display_aplot(Window w, ArrayPlot ap);

#endif

#ifndef ARRAY_PRINT_H
#define ARRAY_PRINT_H

int32 array_print(char *filename, char *xtitle, char *ytitle, char *bottom,
                  int32 nacross, int32 ndown, int32 col0, int32 row0,
                  int32 nskip, int32 ncskip, int32 maxrow, int32 maxcol,
                  double **data, double zmin, double zmax, double tlo,
                  double thi, int32 type);
void ps_replot(double **z, int32 col0, int32 row0, int32 nskip, int32 ncskip,
               int32 maxrow, int32 maxcol, int32 nacross, int32 ndown,
               double zmin, double zmax, int32 type);
void ps_begin(double xlo, double ylo, double xhi, double yhi, double sx,
              double sy);
void ps_convert(double x, double y, double *xs, double *ys);
void ps_col_scale(double y0, double x0, double dy, double dx, int32 n,
                  double zlo, double zhi, int32 type);
void ps_boxit(double tlo, double thi, double jlo, double jhi, double zlo,
              double zhi, char *sx, char *sy, char *sb, int32 type);
void ps_close(void);
void ps_setline(double fill, int32 thick);
void ps_put_char(int32 ch, double *x, double *y);
void ps_text2(char *str, double xr, double yr, int32 icent);
void ps_line2(double x1r, double y1r, double x2r, double y2r);
void ps_set_text(double angle, double slant, double x_size, double y_size);
void ps_rect(double x, double y, double wid, double len);
void ps_bar(double x, double y, double wid, double len, double fill,
            int32 flag);
void ps_rgb_bar(double x, double y, double wid, double len, double fill,
                int32 flag, int32 rgb);
void ps_hsb_bar(double x, double y, double wid, double len, double fill,
                int32 flag);

#endif

#ifndef auto_x11_h
#define auto_x11_h

#include <X11/Xlib.h>

typedef struct {
    Window canvas, axes, numerics, grab, next, run, clear, redraw, base, per;
    Window info, param, file, abort, stab, hint, kill;
} AUTOWIN;

extern int32 AutoRedrawFlag;

extern int32 mark_flag;
extern int32 mark_ibrs, mark_ibre;
extern int32 mark_ipts, mark_ipte;
extern int32 mark_ixs, mark_ixe, mark_iys, mark_iye;

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
void auto_scroll_window(void);
int32 auto_rubber(int32 *i1, int32 *j1, int32 *i2, int32 *j2, int32 flag);
int32 auto_pop_up_list(char *title, char **list, char *key, int32 n, int32 max,
                       int32 def, int32 x, int32 y, char **hints, char *httxt);
void MarkAuto(int32 x, int32 y);
void XORCross(int32 x, int32 y);
void FillCircle(int32 x, int32 y, int32 r);
void LineWidth(int32 wid);
void auto_motion(XEvent ev);
void display_auto(Window w);
Window lil_button(Window root, int32 x, int32 y);
void make_auto(char *wname, char *iname);
void resize_auto_window(XEvent ev);
void a_msg(int32 i, int32 v);
void auto_enter(Window w, int32 v);
void auto_button(XEvent ev);
void auto_kill(void);
void auto_keypress(XEvent ev, int32 *used);
int32 query_special(char *title, char *nsymb);
void clear_msg(void);
void find_point(int32 ibr, int32 pt);
void auto_get_info(int32 *n, char *pname);
void auto_set_mark(int32 i);
void do_auto_range(void);
void RedrawMark(void);

#endif

#ifndef axes2_h
#define axes2_h

#include <X11/Xlib.h>

extern int32 DOING_AXES;
extern int32 DOING_BOX_AXES;

void re_title(void);
void get_title_str(char *s1, char *s2, char *s3);
void make_title(char *str);
double dbl_raise(double x, int32 y);
double make_tics(double tmin, double tmax);
void find_max_min_tic(double *tmin, double *tmax, double tic);
void redraw_cube_pt(double theta, double phi);
void do_axes(void);
void redraw_cube(double theta, double phi);
void draw_unit_cube(void);
void Frame_3d(void);
void Box_axis(double x_min, double x_max, double y_min, double y_max, char *sx,
              char *sy, int32 flag);
void draw_ytics(char *s1, double start, double incr, double end);
void draw_xtics(char *s2, double start, double incr, double end);

#endif

/******************************************************************
 *                                                                *
 * File          : band.h                                         *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a generic BAND linear solver       *
 * package. There are two sets of band solver routines listed in  *
 * this file: one set uses type BandMat defined below and the     *
 * other set uses the type double ** for band matrix arguments.     *
 * The two sets of band solver routines make it easy to work      *
 * with two types of band matrices:                               *
 *                                                                *
 * (1) The BandMat type is intended for use with large            *
 *     band matrices whose elements/columns may be stored in      *
 *     non-contiguous memory locations or even distributed into   *
 *     different processor memories. This type may be modified to *
 *     include such distribution information. If this is done,    *
 *     then all the routines that use BandMat must be modified to *
 *     reflect the new data structure.                            *
 *                                                                *
 * (2) The set of routines that use double ** (and NOT the BandMat  *
 *     type) is intended for use with small matrices which can    *
 *     easily be allocated within a contiguous block of memory    *
 *     on a single processor.                                     *
 *                                                                *
 * Routines that work with the type BandMat begin with "Band".    *
 * The BandAllocMat function allocates a band matrix for use in   *
 * the other matrix routines listed in this file. Matrix storage  *
 * details are given in the documentation for the type BandMat.   *
 * The BandAllocPiv function allocates memory for pivot           *
 * information. The storage allocated by BandAllocMat and         *
 * BandAllocPiv is deallocated by the routines BandFreeMat and    *
 * BandFreePiv, respectively. The BandFactor and BandBacksolve    *
 * routines perform the actual solution of a band linear system.  *
 * Note that the BandBacksolve routine has a parameter b of type  *
 * N_Vector. The current implementation makes use of a machine    *
 * environment specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type N_Vector. Thus, the          *
 * implementation of BandBacksolve may need to change if the      *
 * type N_Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "band" (except for  *
 * the factor and solve routines which are called gbfa and gbsl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for bandalloc.                               *
 *                                                                *
 ******************************************************************/

#ifndef BAND_H
#define BAND_H

#include "llnltyps.h"
#include "vector.h"

/******************************************************************
 *                                                                *
 * Type: BandMat                                                  *
 *----------------------------------------------------------------*
 * The type BandMat is the type of a large (possibly distributed) *
 * band matrix. It is defined to be a pointer to a structure      *
 * with the following fields:                                     *
 *                                                                *
 * size is the number of columns (== number of rows)              *
 *                                                                *
 * mu   is the upper bandwidth, 0 <= mu <= size-1                 *
 *                                                                *
 * ml   is the lower bandwidth, 0 <= ml <= size-1                 *
 *                                                                *
 * smu  is the storage upper bandwidth, mu <= smu <= size-1.      *
 *         The BandFactor routine writes the LU factors           *
 *         into the storage for A. The upper triangular factor U, *
 *         however, may have an upper bandwidth as big as         *
 *         MIN(size-1,mu+ml) because of partial pivoting. The smu *
 *         field holds the upper bandwidth allocated for A.       *
 *                                                                *
 * data is a two dimensional array used for component storage.    *
 *         The elements of a band matrix of type BandMat are      *
 *         stored columnwise (i.e. columns are stored one on top  *
 *         of the other in memory). Only elements within the      *
 *         specified bandwidths are stored.                       *
 *                                                                *
 * If we number rows and columns in the band matrix starting      *
 * from 0, then                                                   *
 *                                                                *
 * data[0] is a pointer to (smu+ml+1)*size contiguous locations   *
 *            which hold the elements within the band of A        *
 *                                                                *
 * data[j] is a pointer to the uppermost element within the band  *
 *            in the jth column. This pointer may be treated as   *
 *            an array indexed from smu-mu (to access the         *
 *            uppermost element within the band in the jth        *
 *            column) to smu+ml (to access the lowest element     *
 *            within the band in the jth column). (Indices from 0 *
 *            to smu-mu-1 give access to extra storage elements   *
 *            required by BandFactor.)                            *
 *                                                                *
 * data[j][i-j+smu] is the (i,j)th element, j-mu <= i <= j+ml.    *
 *                                                                *
 * The macros below allow a user to access individual matrix      *
 * elements without writing out explicit data structure           *
 * references and without knowing too much about the underlying   *
 * element storage. The only storage assumption needed is that    *
 * elements are stored columnwise and that a pointer into the jth *
 * column of elements can be obtained via the BAND_COL macro. The *
 * BAND_COL_ELEM macro selects an element from a column which has *
 * already been isolated via BAND_COL. BAND_COL_ELEM allows the   *
 * user to avoid the translation from the matrix location (i,j)   *
 * to the index in the array returned by BAND_COL at which the    *
 * (i,j)th element is stored. See the documentation for BAND_COL  *
 * and BAND_COL_ELEM for usage details. Users should use these    *
 * macros whenever possible.                                      *
 *                                                                *
 ******************************************************************/

typedef struct {
    int64 size;
    int64 mu, ml, smu;
    double **data;
} *BandMat;

/* BandMat accessor macros */

/******************************************************************
 *                                                                *
 * Macro : BAND_ELEM                                              *
 * Usage : BAND_ELEM(A,i,j) = a_ij;  OR                           *
 *         a_ij = BAND_ELEM(A,i,j);                               *
 *----------------------------------------------------------------*
 * BAND_ELEM(A,i,j) references the (i,j)th element of the         *
 * N by N band matrix A, where 0 <= i,j <= N-1. The location      *
 * (i,j) should further satisfy j-(A->mu) <= i <= j+(A->ml).      *
 *                                                                *
 ******************************************************************/

#define BAND_ELEM(A, i, j) ((A->data)[j][i - j + (A->smu)])

/******************************************************************
 *                                                                *
 * Macro : BAND_COL                                               *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *----------------------------------------------------------------*
 * BAND_COL(A,j) references the diagonal element of the jth       *
 * column of the N by N band matrix A, 0 <= j <= N-1. The type of *
 * the expression BAND_COL(A,j) is double *. The pointer returned   *
 * by the call BAND_COL(A,j) can be treated as an array which is  *
 * indexed from -(A->mu) to (A->ml).                              *
 *                                                                *
 ******************************************************************/

#define BAND_COL(A, j) (((A->data)[j]) + (A->smu))

/******************************************************************
 *                                                                *
 * Macro : BAND_COL_ELEM                                          *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *         BAND_COL_ELEM(col_j,i,j) = a_ij;  OR                   *
 *         a_ij = BAND_COL_ELEM(col_j,i,j);                       *
 *----------------------------------------------------------------*
 * This macro references the (i,j)th entry of the band matrix A   *
 * when used in conjunction with BAND_COL as shown above. The     *
 * index (i,j) should satisfy j-(A->mu) <= i <= j+(A->ml).        *
 *                                                                *
 ******************************************************************/

#define BAND_COL_ELEM(col_j, i, j) (col_j[i - j])

/* Functions that use the BandMat representation for a band matrix */

/******************************************************************
 *                                                                *
 * Function : BandAllocMat                                        *
 * Usage    : A = BandAllocMat(N, mu, ml, smu);                   *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * BandAllocMat allocates memory for an N by N band matrix with   *
 * upper bandwidth mu, lower bandwidth ml, and storage upper      *
 * bandwidth smu. Pass smu as follows depending on whether A will *
 * be factored by BandFactor:                                     *
 *                                                                *
 * (1) Pass smu = mu if A will not be factored.                   *
 *                                                                *
 * (2) Pass smu = MIN(N-1,mu+ml) if A will be factored.           *
 *                                                                *
 * BandAllocMat returns the storage allocated (type BandMat) or   *
 * NULL if the request for matrix storage cannot be satisfied.    *
 * See the documentation for the type BandMat for matrix storage  *
 * details.                                                       *
 *                                                                *
 ******************************************************************/

BandMat BandAllocMat(int64 N, int64 mu, int64 ml, int64 smu);

/******************************************************************
 *                                                                *
 * Function : BandAllocPiv                                        *
 * Usage    : p = BandAllocPiv(N);                                *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * BandAllocPiv allocates memory for pivot information to be      *
 * filled in by the BandFactor routine during the factorization   *
 * of an N by N band matrix. The underlying type for pivot        *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, BandAllocPiv returns NULL.  *
 *                                                                *
 ******************************************************************/

int64 *BandAllocPiv(int64 N);

/******************************************************************
 *                                                                *
 * Function : BandFactor                                          *
 * Usage    : ier = BandFactor(A, p);                             *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * BandFactor performs the LU factorization of the N by N band    *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the "matrix" A and the    *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 *                                                                *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * BandFactor returns 0 if successful. Otherwise it encountered   *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 * Important Note: A must allocated to accomodate the increase    *
 * in upper bandwidth that occurs during factorization. If,       *
 * mathematically, A is a band matrix with upper bandwidth mu and *
 * lower bandwidth ml, then the upper triangular factor U can     *
 * have upper bandwidth as big as smu=MIN(n-1,mu+ml). The lower   *
 * triangular factor L has lower bandwidth ml. Allocate A with    *
 * call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are  *
 * as defined above. The user does not have to zero the "extra"   *
 * storage allocated for the purpose of factorization. This will  *
 * handled by the BandFactor routine.                             *
 *                                                                *
 ******************************************************************/

int64 BandFactor(BandMat A, int64 *p);

/******************************************************************
 *                                                                *
 * Function : BandBacksolve                                       *
 * Usage    : BandBacksolve(A, p, b);                             *
 *----------------------------------------------------------------*
 * BandBacksolve solves the N-dimensional system A x = b using    *
 * the LU factorization in A and the pivot information in p       *
 * computed in BandFactor. The solution x is returned in b. This  *
 * routine cannot fail if the corresponding call to BandFactor    *
 * did not fail.                                                  *
 *                                                                *
 ******************************************************************/

void BandBacksolve(BandMat A, int64 *p, N_Vector b);

/******************************************************************
 *                                                                *
 * Function : BandZero                                            *
 * Usage    : BandZero(A);                                        *
 *----------------------------------------------------------------*
 * A(i,j) <- 0.0,    j-(A->mu) <= i <= j+(A->ml).                 *
 *                                                                *
 ******************************************************************/

void BandZero(BandMat A);

/******************************************************************
 *                                                                *
 * Function : BandCopy                                            *
 * Usage    : BandCopy(A, B, copymu, copyml);                     *
 *----------------------------------------------------------------*
 * BandCopy copies the submatrix with upper and lower bandwidths  *
 * copymu, copyml of the N by N band matrix A into the N by N     *
 * band matrix B.                                                 *
 *                                                                *
 ******************************************************************/

void BandCopy(BandMat A, BandMat B, int64 copymu, int64 copyml);

/******************************************************************
 *                                                                *
 * Function: BandScale                                            *
 * Usage   : BandScale(c, A);                                     *
 *----------------------------------------------------------------*
 * A(i,j) <- c*A(i,j),   j-(A->mu) <= i <= j+(A->ml).             *
 *                                                                *
 ******************************************************************/

void BandScale(double c, BandMat A);

/******************************************************************
 *                                                                *
 * Function : BandAddI                                            *
 * Usage    : BandAddI(A);                                        *
 *----------------------------------------------------------------*
 * A(j,j) <- A(j,j)+1.0,   0 <= j <= (A->size)-1.                 *
 *                                                                *
 ******************************************************************/

void BandAddI(BandMat A);

/******************************************************************
 *                                                                *
 * Function : BandFreeMat                                         *
 * Usage    : BandFreeMat(A);                                     *
 *----------------------------------------------------------------*
 * BandFreeMat frees the memory allocated by BandAllocMat for     *
 * the band matrix A.                                             *
 *                                                                *
 ******************************************************************/

void BandFreeMat(BandMat A);

/******************************************************************
 *                                                                *
 * Function : BandFreePiv                                         *
 * Usage    : BandFreePiv(p);                                     *
 *----------------------------------------------------------------*
 * BandFreePiv frees the memory allocated by BandAllocPiv for     *
 * the pivot information array p.                                 *
 *                                                                *
 ******************************************************************/

void BandFreePiv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : BandPrint                                           *
 * Usage    : BandPrint(A);                                       *
 *----------------------------------------------------------------*
 * This routine prints the N by N band matrix A (upper and lower  *
 * bandwidths A->mu and A->ml, respectively) to standard output   *
 * as it would normally appear on paper. It is intended as a      *
 * debugging tool with small values of N. The elements are        *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/

void BandPrint(BandMat A);

/* Functions that use the double ** representation for a band matrix */

/******************************************************************
 *                                                                *
 * Function : bandalloc                                           *
 * Usage    : double **a;                                           *
 *            a = bandalloc(n, smu, ml);                          *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * bandalloc(n, smu, ml) allocates storage for an n by n band     *
 * matrix A with storage upper bandwidth smu and lower bandwidth  *
 * ml. It returns a pointer to the newly allocated storage if     *
 * successful. If the memory request cannot be satisfied, then    *
 * bandalloc returns NULL. If, mathematically, A has upper and    *
 * lower bandwidths mu and ml, respectively, then the value       *
 * passed to bandalloc for smu may need to be greater than mu.    *
 * The gbfa routine writes the LU factors into the storage (named *
 * "a" in the above usage documentation) for A (thus destroying   *
 * the original elements of A). The upper triangular factor U,    *
 * however, may have a larger upper bandwidth than the upper      *
 * bandwidth mu of A. Thus some "extra" storage for A must be     *
 * allocated if A is to be factored by gbfa. Pass smu as follows: *
 *                                                                *
 * (1) Pass smu = mu if A will not be factored.                   *
 *                                                                *
 * (2) Pass smu = MIN(n-1,mu+ml) if A will be factored.           *
 *                                                                *
 * The underlying type of the band matrix returned is double **. If *
 * we allocate a band matrix A in double **a by                     *
 * a = bandalloc(n,smu,ml), then a[0] is a pointer to             *
 * n * (smu + ml + 1) contiguous storage locations and a[j] is a  *
 * pointer to the uppermost element in the storage for the jth    *
 * column. The expression a[j][i-j+smu] references the (i,j)th    *
 * element of A, where 0 <= i,j <= n-1 and j-mu <= i <= j+ml.     *
 * (The elements a[j][0], a[j][1], ..., a[j][smu-mu-1] are used   *
 * by gbfa and gbsl.)                                             *
 *                                                                *
 ******************************************************************/

double **bandalloc(int64 n, int64 smu, int64 ml);

/******************************************************************
 *                                                                *
 * Function : bandallocpiv                                        *
 * Usage    : int64 *pivot;                                     *
 *            pivot = bandallocpiv(n);                            *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * bandallocpiv(n) allocates an array of n integers. It returns a *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 *                                                                *
 ******************************************************************/

int64 *bandallocpiv(int64 n);

/******************************************************************
 *                                                                *
 * Function : gbfa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * gbfa(a,n,mu,ml,smu,p) factors the n by n band matrix A (upper  *
 * and lower bandwidths mu and ml, storage upper bandwidth smu)   *
 * stored in "a". It overwrites the elements of A with the LU     *
 * factors and it keeps track of the pivot rows chosen in the     *
 * pivot array p.                                                 *
 *                                                                *
 * A successful LU factorization leaves a and pivot array p with  *
 * the following information:                                     *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., n-1.  *
 *                                                                *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * gbfa returns 0 if successful. Otherwise it encountered a zero  *
 * diagonal element during the factorization. In this case it     *
 * returns the column index (numbered from one) at which it       *
 * encountered the zero.                                          *
 *                                                                *
 * IMPORTANT NOTE: Suppose A is a band matrix with upper          *
 * bandwidth mu and lower bandwidth ml, then the upper triangular *
 * factor U can have upper bandwidth as big as MIN(n-1,mu+ml)     *
 * because of partial pivoting. The lower triangular factor L has *
 * lower bandwidth ml. Thus, if A is to be factored and           *
 * backsolved using gbfa and gbsl, then it should be allocated    *
 * as a = bandalloc(n,smu,ml), where smu = MIN(n-1,mu+ml). The    *
 * call to gbfa is ier = gbfa(a,n,mu,ml,smu,p). The corresponding *
 * call to gbsl is gbsl(a,n,smu,ml,p,b). The user does not need   *
 * to zero the "extra" storage allocated for the purpose of       *
 * factorization. This is handled by the gbfa routine. If A is    *
 * not going to be factored and backsolved, then it can be        *
 * allocated as a = bandalloc(n,smu,ml). In either case, all      *
 * routines in this section use the parameter name smu for a      *
 * parameter which must be the "storage upper bandwidth" which    *
 * was passed to bandalloc.                                       *
 *                                                                *
 ******************************************************************/

int64 gbfa(double **a, int64 n, int64 mu, int64 ml, int64 smu, int64 *p);

/******************************************************************
 *                                                                *
 * Function : gbsl                                                *
 * Usage    : double *b;                                            *
 *            ier = gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier == 0) gbsl(a,n,smu,ml,p,b);                 *
 *----------------------------------------------------------------*
 * gbsl(a,n,smu,ml,p,b) solves the n by n linear system           *
 * Ax = b, where A is band matrix stored in "a" with storage      *
 * upper bandwidth smu and lower bandwidth ml. It assumes that A  *
 * has been LU factored and the pivot array p has been set by a   *
 * successful call gbfa(a,n,mu,ml,smu,p). The solution x is       *
 * written into the b array.                                      *
 *                                                                *
 ******************************************************************/

void gbsl(double **a, int64 n, int64 smu, int64 ml, int64 *p, double *b);

/******************************************************************
 *                                                                *
 * Function : bandzero                                            *
 * Usage    : bandzero(a,n,mu,ml,smu);                            *
 *----------------------------------------------------------------*
 * a(i,j) <- 0.0,   0 <= i,j <= n-1, j-mu <= i <= j+ml.           *
 *                                                                *
 ******************************************************************/

void bandzero(double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 *                                                                *
 * Function : bandcopy                                            *
 * Usage    : bandcopy(a,b,n,a_smu,b_smu,copymu,copyml);          *
 *----------------------------------------------------------------*
 * b(i,j) <- a(i,j), 0 <= i,j <= n-1, j-copymu <= i <= j+copyml.  *
 *                                                                *
 ******************************************************************/

void bandcopy(double **a, double **b, int64 n, int64 a_smu, int64 b_smu,
              int64 copymu, int64 copyml);

/******************************************************************
 *                                                                *
 * Function : bandscale                                           *
 * Usage    : bandscale(c,a,n,mu,ml);                             *
 *----------------------------------------------------------------*
 * a(i,j) <- c*a(i,j),   0 <= i,j <= n-1, j-mu <= i <= j+ml.      *
 *                                                                *
 ******************************************************************/

void bandscale(double c, double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 *                                                                *
 * Function : bandaddI                                            *
 * Usage    : bandaddI(a,n,smu);                                  *
 *----------------------------------------------------------------*
 * a(j,j) <- a(j,j)+1.0,   0 <= j <= n-1.                         *
 *                                                                *
 ******************************************************************/

void bandaddI(double **a, int64 n, int64 smu);

/******************************************************************
 *                                                                *
 * Function : bandfreepiv                                         *
 * Usage    : bandfreepiv(p);                                     *
 *----------------------------------------------------------------*
 * bandfreepiv(p) frees the pivot array p allocated by            *
 * bandallocpiv.                                                  *
 *                                                                *
 ******************************************************************/

void bandfreepiv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : bandfree                                            *
 * Usage    : bandfree(a);                                        *
 *----------------------------------------------------------------*
 * bandfree(a) frees the band matrix a allocated by bandalloc.    *
 *                                                                *
 ******************************************************************/

void bandfree(double **a);

/******************************************************************
 *                                                                *
 * Function : bandprint                                           *
 * Usage    : bandprint(a,n,mu,ml,smu);                           *
 *----------------------------------------------------------------*
 * bandprint(a,n,mu,ml,smu) prints the n by n band matrix stored  *
 * in a (with upper bandwidth mu and lower bandwidth ml) to       *
 * standard output as it would normally appear on paper. It is    *
 * intended as a debugging tool with small values of n. The       *
 * elements are printed using the %g option. A blank line is      *
 * printed before and after the matrix.                           *
 *                                                                *
 ******************************************************************/

void bandprint(double **a, int64 n, int64 mu, int64 ml, int64 smu);

#endif

#ifndef BROWSE_H
#define BROWSE_H

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
    double **data;
    int32 istart, iend;
} BROWSER;

extern BROWSER my_browser;

double **get_browser_data(void);
void set_browser_data(double **data, int32 col0);
double *get_data_col(int32 c);
int32 gettimenow(void);
void waitasec(int32 msec);
int32 get_maxrow_browser(void);
void write_mybrowser_data(FILE *fp);
void write_browser_data(FILE *fp, BROWSER *b);
int32 check_for_stor(double **data);
void del_stor_col(char *var, BROWSER *b);
void data_del_col(BROWSER *b);
void data_add_col(BROWSER *b);
int32 add_stor_col(char *name, char *formula, BROWSER *b);
void chk_seq(char *f, int32 *seq, double *a1, double *a2);
void replace_column(char *var, char *form, double **dat, int32 n);
void wipe_rep(void);
void unreplace_column(void);
void make_d_table(double xlo, double xhi, int32 col, char *filename, BROWSER b);
void find_value(int32 col, double val, int32 *row, BROWSER b);
void find_variable(char *s, int32 *col);
void browse_but_on(BROWSER *b, int32 i, Window w, int32 yn);
void enter_browser(XEvent ev, BROWSER *b, int32 yn);
void display_browser(Window w, BROWSER b);
void redraw_browser(BROWSER b);
void new_browse_dat(double **new_dat, int32 dat_len);
void refresh_browser(int32 length);
void reset_browser(void);
void draw_data(BROWSER b);
void init_browser(void);
void kill_browser(BROWSER *b);
void make_new_browser(void);
Window br_button(Window root, int32 row, int32 col, int32 iflag);
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
void get_data_xyz(double *x, double *y, double *z, int32 i1, int32 i2, int32 i3,
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
void write_mybrowser_data(FILE *fp);
void data_get_mybrowser(int32 row);

#endif

#ifndef CALC_H
#define CALC_H

#include <X11/Xlib.h>

void draw_calc(Window w);
void make_calc(double z);
void quit_calc(void);
void ini_calc_string(char *name, char *value, int32 *pos, int32 *col);
void q_calc(void);
int32 do_calc(char *temp, double *z);
int32 has_eq(char *z, char *w, int32 *where);
double calculate(char *expr, int32 *ok);

#endif

#ifndef CHOICE_BOX_H
#define CHOICE_BOX_H

#include <X11/Xlib.h>
#include "struct.h"

void destroy_choice(CHOICE_BOX p);
void display_choice(Window w, CHOICE_BOX p);
void do_checks(CHOICE_BOX p);
void base_choice(char *wname, int32 n, int32 mcc, char **names, int32 *check,
                 int32 type);
int32 do_choice_box(Window root, char *wname, int32 n, int32 mcc, char **names,
                    int32 *check, int32 type);
int32 choice_box_event_loop(CHOICE_BOX p);

#endif

#ifndef COLOR_H
#define COLOR_H

#include <X11/Xlib.h>

extern int32 custom_color;

extern int32 color_min, color_total, color_max;
extern int32 COLOR;

void tst_color(Window w);
void set_scolor(int32 col);
void set_color(int32 col);
void make_cmaps(int32 *r, int32 *g, int32 *b, int32 n, int32 type);
int32 rfun(double y, int32 per);
int32 gfun(double y, int32 per);
int32 bfun(double y, int32 per);
void NewColormap(int32 type);
void get_ps_color(int32 i, double *r, double *g, double *b);
void get_svg_color(int32 i, int32 *r, int32 *g, int32 *b);
void MakeColormap(void);
uint32 ColorMap(int32 i);

#endif

#ifndef COMLINE_H
#define COMLINE_H

typedef struct {
    char *name;
    char *does;
    uint32 use;
} INTERN_SET;

typedef struct {
    char *name;
    struct SET_NAME *next;
} SET_NAME;

extern int32 NincludedFiles;
extern int32 Nintern_2_use;

int32 is_set_name(SET_NAME *set, char *nam);
SET_NAME *add_set(SET_NAME *set, char *nam);
SET_NAME *rm_set(SET_NAME *set, char *nam);
void do_comline(int32 argc, char **argv);
int32 if_needed_select_sets(void);
int32 if_needed_load_set(void);
int32 if_needed_load_par(void);
int32 if_needed_load_ic(void);
int32 if_needed_load_ext_options(void);
int32 parse_it(char *com);

#endif

#ifndef DAE_FUN_H
#define DAE_FUN_H

int32 add_svar(char *name, char *rhs);
int32 add_svar_names(void);
int32 add_aeqn(char *rhs);
int32 compile_svars(void);
void reset_dae(void);
void set_init_guess(void);
void err_dae(void);
void init_dae_work(void);
void get_dae_fun(double *y, double *f);
void do_daes(void);
int32 solve_dae(void);
void get_new_guesses(void);

#endif

#ifndef DELAY_HANDLE_H
#define DELAY_HANDLE_H

double delay_stab_eval(double delay, int32 var);
int32 alloc_delay(double big);
void free_delay(void);
void stor_delay(double *y);
double get_delay_old(int32 in, double tau);
void polint(double *xa, double *ya, int32 n, double x, double *y, double *dy);
double get_delay(int32 in, double tau);
int32 do_init_delay(double big);

#endif

#ifndef DEL_STAB_H
#define DEL_STAB_H

typedef struct {
    double r, i;
} COMPLEX;

void do_delay_sing(double *x, double eps, double err, double big, int32 maxit,
                   int32 n, int32 *ierr, double *stabinfo);
COMPLEX csum(COMPLEX z, COMPLEX w);
COMPLEX cdif(COMPLEX z, COMPLEX w);
COMPLEX cmlt(COMPLEX z, COMPLEX w);
COMPLEX cdivv(COMPLEX z, COMPLEX w);
COMPLEX cexp2(COMPLEX z);
void switch_rows(COMPLEX *z, int32 i1, int32 i2, int32 n);
COMPLEX rtoc(double x, double y);
void cprintn(COMPLEX z);
void cprint(COMPLEX z);
void cprintarr(COMPLEX *z, int32 n, int32 m);
double c_abs(COMPLEX z);
COMPLEX cdeterm(COMPLEX *z, int32 n);
COMPLEX cxdeterm(COMPLEX *z, int32 n);
void make_z(COMPLEX *z, double *delay, int32 n, int32 m, double *coef,
            COMPLEX lambda);
int32 find_positive_root(double *coef, double *delay, int32 n, int32 m,
                         double err, double eps, double big, int32 maxit,
                         double *rr);
void process_root(double real, double im);
double get_arg(double *delay, double *coef, int32 m, int32 n, COMPLEX lambda);
int32 test_sign(double old, double new);
int32 plot_args(double *coef, double *delay, int32 n, int32 m, int32 npts,
                double almax, double wmax);

#endif

/******************************************************************
 *                                                                *
 * File          : dense.h                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a generic DENSE linear solver      *
 * package. There are two sets of dense solver routines listed in *
 * this file: one set uses type DenseMat defined below and the    *
 * other set uses the type double ** for dense matrix arguments.    *
 * The two sets of dense solver routines make it easy to work     *
 * with two types of dense matrices:                              *
 *                                                                *
 * (1) The DenseMat type is intended for use with large dense     *
 *     matrices whose elements/columns may be stored in           *
 *     non-contiguous memory locations or even distributed into   *
 *     different processor memories. This type may be modified to *
 *     include such distribution information. If this is done,    *
 *     then all the routines that use DenseMat must be modified   *
 *     to reflect the new data structure.                         *
 *                                                                *
 * (2) The set of routines that use double ** (and NOT the DenseMat *
 *     type) is intended for use with small matrices which can    *
 *     easily be allocated within a contiguous block of memory    *
 *     on a single processor.                                     *
 *                                                                *
 * Routines that work with the type DenseMat begin with "Dense".  *
 * The DenseAllocMat function allocates a dense matrix for use in *
 * the other DenseMat routines listed in this file. Matrix        *
 * storage details are given in the documentation for the type    *
 * DenseMat. The DenseAllocPiv function allocates memory for      *
 * pivot information. The storage allocated by DenseAllocMat and  *
 * DenseAllocPiv is deallocated by the routines DenseFreeMat and  *
 * DenseFreePiv, respectively. The DenseFactor and DenseBacksolve *
 * routines perform the actual solution of a dense linear system. *
 * Note that the DenseBacksolve routine has a parameter b of type *
 * N_Vector. The current implementation makes use of a machine    *
 * environment-specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type N_Vector. Thus, the          *
 * implementation of DenseBacksolve may need to change if the     *
 * type N_Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "den" (except for   *
 * the factor and solve routines which are called gefa and gesl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for denalloc.                                *
 *                                                                *
 ******************************************************************/

#ifndef DENSE_H
#define DENSE_H

#include "llnltyps.h"
#include "vector.h"

/******************************************************************
 *                                                                *
 * Type: DenseMat                                                 *
 *----------------------------------------------------------------*
 * The type DenseMat is defined to be a pointer to a structure    *
 * with a size and a data field. The size field indicates the     *
 * number of columns (== number of rows) of a dense matrix, while *
 * the data field is a two dimensional array used for component   *
 * storage. The elements of a dense matrix are stored columnwise  *
 * (i.e columns are stored one on top of the other in memory). If *
 * A is of type DenseMat, then the (i,j)th element of A (with     *
 * 0 <= i,j <= size-1) is given by the expression (A->data)[j][i] *
 * or by the expression (A->data)[0][j*n+i]. The macros below     *
 * allow a user to access efficiently individual matrix           *
 * elements without writing out explicit data structure           *
 * references and without knowing too much about the underlying   *
 * element storage. The only storage assumption needed is that    *
 * elements are stored columnwise and that a pointer to the jth   *
 * column of elements can be obtained via the DENSE_COL macro.    *
 * Users should use these macros whenever possible.               *
 *                                                                *
 ******************************************************************/

typedef struct {
    int64 size;
    double **data;
} *DenseMat;

/* DenseMat accessor macros */

/******************************************************************
 *                                                                *
 * Macro : DENSE_ELEM                                             *
 * Usage : DENSE_ELEM(A,i,j) = a_ij;  OR                          *
 *         a_ij = DENSE_ELEM(A,i,j);                              *
 *----------------------------------------------------------------*
 * DENSE_ELEM(A,i,j) references the (i,j)th element of the N by N *
 * DenseMat A, 0 <= i,j <= N-1.                                   *
 *                                                                *
 ******************************************************************/

#define DENSE_ELEM(A, i, j) ((A->data)[j][i])

/******************************************************************
 *                                                                *
 * Macro : DENSE_COL                                              *
 * Usage : col_j = DENSE_COL(A,j);                                *
 *----------------------------------------------------------------*
 * DENSE_COL(A,j) references the jth column of the N by N         *
 * DenseMat A, 0 <= j <= N-1. The type of the expression          *
 * DENSE_COL(A,j) is double *. After the assignment in the usage    *
 * above, col_j may be treated as an array indexed from 0 to N-1. *
 * The (i,j)th element of A is referenced by col_j[i].            *
 *                                                                *
 ******************************************************************/

#define DENSE_COL(A, j) ((A->data)[j])

/* Functions that use the DenseMat representation for a dense matrix */

/******************************************************************
 *                                                                *
 * Function : DenseAllocMat                                       *
 * Usage    : A = DenseAllocMat(N);                               *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocMat allocates memory for an N by N dense matrix and  *
 * returns the storage allocated (type DenseMat). DenseAllocMat   *
 * returns NULL if the request for matrix storage cannot be       *
 * satisfied. See the above documentation for the type DenseMat   *
 * for matrix storage details.                                    *
 *                                                                *
 ******************************************************************/

DenseMat DenseAllocMat(int64 N);

/******************************************************************
 *                                                                *
 * Function : DenseAllocPiv                                       *
 * Usage    : p = DenseAllocPiv(N);                               *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocPiv allocates memory for pivot information to be     *
 * filled in by the DenseFactor routine during the factorization  *
 * of an N by N dense matrix. The underlying type for pivot       *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, DenseAllocPiv returns NULL. *
 *                                                                *
 ******************************************************************/

int64 *DenseAllocPiv(int64 N);

/******************************************************************
 *                                                                *
 * Function : DenseFactor                                         *
 * Usage    : ier = DenseFactor(A, p);                            *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * DenseFactor performs the LU factorization of the N by N dense  *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the matrix A and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 *                                                                *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * DenseFactor returns 0 if successful. Otherwise it encountered  *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 ******************************************************************/

int64 DenseFactor(DenseMat A, int64 *p);

/******************************************************************
 *                                                                *
 * Function : DenseBacksolve                                      *
 * Usage    : DenseBacksolve(A, p, b);                            *
 *----------------------------------------------------------------*
 * DenseBacksolve solves the N-dimensional system A x = b using   *
 * the LU factorization in A and the pivot information in p       *
 * computed in DenseFactor. The solution x is returned in b. This *
 * routine cannot fail if the corresponding call to DenseFactor   *
 * did not fail.                                                  *
 *                                                                *
 ******************************************************************/

void DenseBacksolve(DenseMat A, int64 *p, N_Vector b);

/******************************************************************
 *                                                                *
 * Function : DenseZero                                           *
 * Usage    : DenseZero(A);                                       *
 *----------------------------------------------------------------*
 * DenseZero sets all the elements of the N by N matrix A to 0.0. *
 *                                                                *
 ******************************************************************/

void DenseZero(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : DenseCopy                                           *
 * Usage    : DenseCopy(A, B);                                    *
 *----------------------------------------------------------------*
 * DenseCopy copies the contents of the N by N matrix A into the  *
 * N by N matrix B.                                               *
 *                                                                *
 ******************************************************************/

void DenseCopy(DenseMat A, DenseMat B);

/******************************************************************
 *                                                                *
 * Function: DenseScale                                           *
 * Usage   : DenseScale(c, A);                                    *
 *----------------------------------------------------------------*
 * DenseScale scales the elements of the N by N matrix A by the   *
 * constant c and stores the result back in A.                    *
 *                                                                *
 ******************************************************************/

void DenseScale(double c, DenseMat A);

/******************************************************************
 *                                                                *
 * Function : DenseAddI                                           *
 * Usage    : DenseAddI(A);                                       *
 *----------------------------------------------------------------*
 * DenseAddI adds the identity matrix to A and stores the result  *
 * back in A.                                                     *
 *                                                                *
 ******************************************************************/

void DenseAddI(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : DenseFreeMat                                        *
 * Usage    : DenseFreeMat(A);                                    *
 *----------------------------------------------------------------*
 * DenseFreeMat frees the memory allocated by DenseAllocMat for   *
 * the N by N matrix A.                                           *
 *                                                                *
 ******************************************************************/

void DenseFreeMat(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : DenseFreePiv                                        *
 * Usage    : DenseFreePiv(p);                                    *
 *----------------------------------------------------------------*
 * DenseFreePiv frees the memory allocated by DenseAllocPiv for   *
 * the pivot information array p.                                 *
 *                                                                *
 ******************************************************************/

void DenseFreePiv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : DensePrint                                          *
 * Usage    : DensePrint(A);                                      *
 *----------------------------------------------------------------*
 * This routine prints the N by N dense matrix A to standard      *
 * output as it would normally appear on paper. It is intended    *
 * as a debugging tool with small values of N. The elements are   *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/

void DensePrint(DenseMat A);

/* Functions that use the double ** representation for a dense matrix */

/******************************************************************
 *                                                                *
 * Function : denalloc                                            *
 * Usage    : double **a;                                           *
 *            a = denalloc(n);                                    *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * denalloc(n) allocates storage for an n by n dense matrix. It   *
 * returns a pointer to the newly allocated storage if            *
 * successful. If the memory request cannot be satisfied, then    *
 * denalloc returns NULL. The underlying type of the dense matrix *
 * returned is double **. If we allocate a dense matrix double **a by *
 * a = denalloc(n), then a[j][i] references the (i,j)th element   *
 * of the matrix a, 0 <= i,j <= n-1, and a[j] is a pointer to the *
 * first element in the jth column of a. The location a[0]        *
 * contains a pointer to n^2 contiguous locations which contain   *
 * the elements of a.                                             *
 *                                                                *
 ******************************************************************/

double **denalloc(int64 n);

/******************************************************************
 *                                                                *
 * Function : denallocpiv                                         *
 * Usage    : int64 *pivot;                                     *
 *            pivot = denallocpiv(n);                             *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * denallocpiv(n) allocates an array of n integers. It returns a  *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 *                                                                *
 ******************************************************************/

int64 *denallocpiv(int64 n);

/******************************************************************
 *                                                                *
 * Function : gefa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = gefa(a,n,p);                                  *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * gefa(a,n,p) factors the n by n dense matrix a. It overwrites   *
 * the elements of a with its LU factors and keeps track of the   *
 * pivot rows chosen in the pivot array p.                        *
 *                                                                *
 * A successful LU factorization leaves the matrix a and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., n-1.  *
 *                                                                *
 * (2) If the unique LU factorization of a is given by Pa = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of a     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of a contains the multipliers, I-L.        *
 *                                                                *
 * gefa returns 0 if successful. Otherwise it encountered a zero  *
 * diagonal element during the factorization. In this case it     *
 * returns the column index (numbered from one) at which it       *
 * encountered the zero.                                          *
 *                                                                *
 ******************************************************************/

int64 gefa(double **a, int64 n, int64 *p);

/******************************************************************
 *                                                                *
 * Function : gesl                                                *
 * Usage    : double *b;                                            *
 *            ier = gefa(a,n,p);                                  *
 *            if (ier == 0) gesl(a,n,p,b);                        *
 *----------------------------------------------------------------*
 * gesl(a,n,p,b) solves the n by n linear system ax = b. It       *
 * assumes that a has been LU factored and the pivot array p has  *
 * been set by a successful call to gefa(a,n,p). The solution x   *
 * is written into the b array.                                   *
 *                                                                *
 ******************************************************************/

void gesl(double **a, int64 n, int64 *p, double *b);

/******************************************************************
 *                                                                *
 * Function : denzero                                             *
 * Usage    : denzero(a,n);                                       *
 *----------------------------------------------------------------*
 * denzero(a,n) sets all the elements of the n by n dense matrix  *
 * a to be 0.0.                                                   *
 *                                                                *
 ******************************************************************/

void denzero(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : dencopy                                             *
 * Usage    : dencopy(a,b,n);                                     *
 *----------------------------------------------------------------*
 * dencopy(a,b,n) copies the n by n dense matrix a into the       *
 * n by n dense matrix b.                                         *
 *                                                                *
 ******************************************************************/

void dencopy(double **a, double **b, int64 n);

/******************************************************************
 *                                                                *
 * Function : denscale                                            *
 * Usage    : denscale(c,a,n);                                    *
 *----------------------------------------------------------------*
 * denscale(c,a,n) scales every element in the n by n dense       *
 * matrix a by c.                                                 *
 *                                                                *
 ******************************************************************/

void denscale(double c, double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : denaddI                                             *
 * Usage    : denaddI(a,n);                                       *
 *----------------------------------------------------------------*
 * denaddI(a,n) increments the n by n dense matrix a by the       *
 * identity matrix.                                               *
 *                                                                *
 ******************************************************************/

void denaddI(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : denfreepiv                                          *
 * Usage    : denfreepiv(p);                                      *
 *----------------------------------------------------------------*
 * denfreepiv(p) frees the pivot array p allocated by             *
 * denallocpiv.                                                   *
 *                                                                *
 ******************************************************************/

void denfreepiv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : denfree                                             *
 * Usage    : denfree(a);                                         *
 *----------------------------------------------------------------*
 * denfree(a) frees the dense matrix a allocated by denalloc.     *
 *                                                                *
 ******************************************************************/

void denfree(double **a);

/******************************************************************
 *                                                                *
 * Function : denprint                                            *
 * Usage    : denprint(a,n);                                      *
 *----------------------------------------------------------------*
 * denprint(a,n) prints the n by n dense matrix a to standard     *
 * output as it would normally appear on paper. It is intended as *
 * a debugging tool with small values of n. The elements are      *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/

void denprint(double **a, int64 n);

#endif

#ifndef DERIVED_H
#define DERIVED_H

void free_derived(void);
int32 compile_derived(void);
void evaluate_derived(void);
int32 add_derived(char *name, char *rhs);

#endif

#ifndef DIAGRAM_H
#define DIAGRAM_H
#include "auto_nox.h"

extern Diagram *bifd;

#include <stdio.h>

extern int32 NBifs;

void start_diagram(int32 n);
void edit_start(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar,
                double a, double *uhi, double *ulo, double *u0, double *ubar,
                double *par, double per, int32 n, int32 icp1, int32 icp2,
                int32 icp3, int32 icp4, double *evr, double *evi);
void edit_diagram(Diagram *d, int32 ibr, int32 ntot, int32 itp, int32 lab,
                  int32 nfpar, double a, double *uhi, double *ulo, double *u0,
                  double *ubar, double *par, double per, int32 n, int32 icp1,
                  int32 icp2, int32 icp3, int32 icp4, int32 flag2, double *evr,
                  double *evi, double tp);
void add_diagram(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar,
                 double a, double *uhi, double *ulo, double *u0, double *ubar,
                 double *par, double per, int32 n, int32 icp1, int32 icp2,
                 int32 icp3, int32 icp4, int32 flag2, double *evr, double *evi);
void kill_diagrams(void);
void redraw_diagram(void);
void write_info_out(void);
void write_init_data_file(void);
void write_pts(void);
void post_auto(void);
void svg_auto(void);
void bound_diagram(double *xlo, double *xhi, double *ylo, double *yhi);
int32 save_diagram(FILE *fp, int32 n);
int32 load_diagram(FILE *fp, int32 node);
void load_browser_with_branch(int32 ibr, int32 pts, int32 pte);

#endif

#ifndef DIALOG_BOX_H
#define DIALOG_BOX_H

#include <X11/Xlib.h>
#include "struct.h"

int32 get_dialog(char *wname, char *name, char *value, char *ok, char *cancel,
                 int32 max);
int32 dialog_event_loop(DIALOG *d, int32 *pos, int32 *col);
void display_dialog(Window w, DIALOG d, int32 col);

#endif

#ifndef DO_FIT_H
#define DO_FIT_H

typedef struct {
    char file[25];
    char varlist[25], collist[25];
    char parlist1[25], parlist2[25];
    int32 dim, npars, nvars, npts, maxiter;
    int32 icols[50], ipar[50], ivar[50];
    double tol, eps;
} FITINFO;

void init_fit_info(void);
void get_fit_info(double *y, double *a, double *t0, int32 *flag, double eps,
                  double *yfit, double **yderv, int32 npts, int32 npars,
                  int32 nvars, int32 *ivar, int32 *ipar);
void printem(double **yderv, double *yfit, double *t0, int32 npars, int32 nvars,
             int32 npts);
int32 one_step_int(double *y, double t0, double t1, int32 *istart);
void print_fit_info(void);
void test_fit(void);
int32 run_fit(char *filename, int32 npts, int32 npars, int32 nvars,
              int32 maxiter, int32 ndim, double eps, double tol, int32 *ipar,
              int32 *ivar, int32 *icols, double *y0, double *a, double *yfit);
int32 marlevstep(double *t0, double *y0, double *y, double *sig, double *a,
                 int32 npts, int32 nvars, int32 npars, int32 *ivar, int32 *ipar,
                 double *covar, double *alpha, double *chisq, double *alambda,
                 double *work, double **yderv, double *yfit, double *ochisq,
                 int32 ictrl, double eps);
int32 mrqcof(double *t0, double *y0, double *y, double *sig, double *a,
             int32 npts, int32 nvars, int32 npars, int32 *ivar, int32 *ipar,
             double *alpha, double *chisq, double *beta, double **yderv,
             double *yfit, double eps);
int32 get_fit_params(void);
void parse_collist(char *collist, int32 *icols, int32 *n);
void parse_varlist(char *varlist, int32 *ivars, int32 *n);
void parse_parlist(char *parlist, int32 *ipars, int32 *n);

#endif

/*      DOP853
        ------

This code computes the numerical solution of a system of first order ordinary
differential equations y'=f(x,y). It uses an explicit Runge-Kutta method of
order 8(5,3) due to Dormand & Prince with step size control and dense output.

Authors : E. Hairer & G. Wanner
          Universite de Geneve, dept. de Mathematiques
          CH-1211 GENEVE 4, SWITZERLAND
          E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).

Version of Mai 2, 1994.

Remarks about the C version : this version allocates memory by itself, the
iwork array (among the initial FORTRAN parameters) has been splitted into
independant initial parameters, the statistical variables and last step size
and x have been encapsulated in the module and are now accessible through
dedicated functions; the variable names have been kept to maintain a kind
of reading compatibility between the C and FORTRAN codes; adaptation made by
J.Colinge (COLINGE@DIVSUN.UNIGE.CH).

INPUT PARAMETERS
----------------

n        Dimension of the system (n < UINT_MAX).

fcn      A pointer the the function definig the differential equation, this
         function must have the following prototype

           void fcn (uint32 n, double x, double *y, double *f)

         where the array f will be filled with the function result.

x        Initial x value.

*y       Initial y values (double y[n]).

xend     Final x value (xend-x may be positive or negative).

*rtoler  Relative and absolute error tolerances. They can be both scalars or
*atoler  vectors of length n (in the scalar case pass the addresses of
         variables where you have placed the tolerance values).

itoler   Switch for atoler and rtoler :
           itoler=0 : both atoler and rtoler are scalars, the code keeps
                      roughly the local error of y[i] below
                      rtoler*abs(y[i])+atoler.
           itoler=1 : both rtoler and atoler are vectors, the code keeps
                      the local error of y[i] below
                      rtoler[i]*abs(y[i])+atoler[i].

solout   A pointer to the output function called during integration.
         If iout >= 1, it is called after every successful step. If iout = 0,
         pass a pointer equal to NULL. solout must must have the following
         prototype

           solout (long nr, double xold, double x, double* y, uint32 n, int32*
irtrn)

         where y is the solution the at nr-th grid point x, xold is the
         previous grid point and irtrn serves to interrupt the integration
         (if set to a negative value).

         Continuous output : during the calls to solout, a continuous solution
         for the interval (xold,x) is available through the function

           contd8(i,s)

         which provides an approximation to the i-th component of the solution
         at the point s (s must lie in the interval (xold,x)).

iout     Switch for calling solout :
           iout=0 : no call,
           iout=1 : solout only used for output,
           iout=2 : dense output is performed in solout (in this case nrdens
                    must be greater than 0).

fileout  A pointer to the stream used for messages, if you do not want any
         message, just pass NULL.

icont    An array containing the indexes of components for which dense
         output is required. If no dense output is required, pass NULL.

licont   The number of cells in icont.

Sophisticated setting of parameters
-----------------------------------

         Several parameters have a default value (if set to 0) but, to better
         adapt the code to your problem, you can specify particular initial
         values.

uround   The rounding unit, default 2.3E-16 (this default value can be
         replaced in the code by DBL_EPSILON providing double.h defines it
         in your system).

safe     Safety factor in the step size prediction, default 0.9.

fac1     Parameters for step size selection; the new step size is chosen
fac2     subject to the restriction  fac1 <= hnew/hold <= fac2.
         Default values are fac1=0.333 and fac2=6.0.

beta     The "beta" for stabilized step size control (see section IV.2 of our
         book). Larger values for beta ( <= 0.1 ) make the step size control
         more stable. Negative initial value provoke beta=0; default beta=0.

hmax     Maximal step size, default xend-x.

h        Initial step size, default is a guess computed by the function hinit.

nmax     Maximal number of allowed steps, default 100000.

meth     Switch for the choice of the method coefficients; at the moment the
         only possibility and default value are 1.

nstiff   Test for stiffness is activated when the current step number is a
         multiple of nstiff. A negative value means no test and the default
         is 1000.

nrdens   Number of components for which dense outpout is required, default 0.
         For 0 < nrdens < n, the components have to be specified in icont[0],
         icont[1], ... icont[nrdens-1]. Note that if nrdens=0 or nrdens=n, no
         icont is needed, pass NULL.

Memory requirements
-------------------

         The function dop853 allocates dynamically 11*n doubles for the method
         stages, 8*nrdens doubles for the interpolation if dense output is
         performed and n uint32 if 0 < nrdens < n.

OUTPUT PARAMETERS
-----------------

y       numerical solution at x=xRead() (see below).

dopri5 returns the following values

         1 : computation successful,
         2 : computation successful interrupted by solout,
        -1 : input is not consistent,
        -2 : larger nmax is needed,
        -3 : step size becomes too small,
        -4 : the problem is probably stff (interrupted).

Several functions provide access to different values :

xRead   x value for which the solution has been computed (x=xend after
        successful return).

hRead   Predicted step size of the last accepted step (useful for a subsequent
        call to dop853).

nstepRead   Number of used steps.
naccptRead  Number of accepted steps.
nrejctRead  Number of rejected steps.
nfcnRead    Number of function calls.

*/

/*      DOPRI5
        ------

This code computes the numerical solution of a system of first order ordinary
differential equations y'=f(x,y). It uses an explicit Runge-Kutta method of
order (4)5 due to Dormand & Prince with step size control and dense output.

Authors : E. Hairer & G. Wanner
          Universite de Geneve, dept. de Mathematiques
          CH-1211 GENEVE 4, SWITZERLAND
          E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).

Version of April 28, 1994.

Remarks about the C version : this version allocates memory by itself, the
iwork array (among the initial FORTRAN parameters) has been splitted into
independant initial parameters, the statistical variables and last step size
and x have been encapsulated in the module and are now accessible through
dedicated functions, the variable names have been kept to maintain a kind
of reading compatibility between the C and FORTRAN codes; adaptation made by
J.Colinge (COLINGE@DIVSUN.UNIGE.CH).

INPUT PARAMETERS
----------------

n        Dimension of the system (n < UINT_MAX).

fcn      A pointer the the function definig the differential equation, this
         function must have the following prototype

           void fcn (uint32 n, double x, double *y, double *f)

         where the array f will be filled with the function result.

x        Initial x value.

*y       Initial y values (double y[n]).

xend     Final x value (xend-x may be positive or negative).

*rtoler  Relative and absolute error tolerances. They can be both scalars or
*atoler  vectors of length n (in the scalar case pass the addresses of
         variables where you have placed the tolerance values).

itoler   Switch for atoler and rtoler :
           itoler=0 : both atoler and rtoler are scalars, the code keeps
                      roughly the local error of y[i] below
                      rtoler*abs(y[i])+atoler.
           itoler=1 : both rtoler and atoler are vectors, the code keeps
                      the local error of y[i] below
                      rtoler[i]*abs(y[i])+atoler[i].

solout   A pointer to the output function called during integration.
         If iout >= 1, it is called after every successful step. If iout = 0,
         pass a pointer equal to NULL. solout must must have the following
         prototype

           solout (long nr, double xold, double x, double* y, uint32 n, int32*
irtrn)

         where y is the solution the at nr-th grid point x, xold is the
         previous grid point and irtrn serves to interrupt the integration
         (if set to a negative value).

         Continuous output : during the calls to solout, a continuous solution
         for the interval (xold,x) is available through the function

           contd5(i,s)

         which provides an approximation to the i-th component of the solution
         at the point s (s must lie in the interval (xold,x)).

iout     Switch for calling solout :
           iout=0 : no call,
           iout=1 : solout only used for output,
           iout=2 : dense output is performed in solout (in this case nrdens
                    must be greater than 0).

fileout  A pointer to the stream used for messages, if you do not want any
         message, just pass NULL.

icont    An array containing the indexes of components for which dense
         output is required. If no dense output is required, pass NULL.

licont   The number of cells in icont.

Sophisticated setting of parameters
-----------------------------------

         Several parameters have a default value (if set to 0) but, to better
         adapt the code to your problem, you can specify particular initial
         values.

uround   The rounding unit, default 2.3E-16 (this default value can be
         replaced in the code by DBL_EPSILON providing double.h defines it
         in your system).

safe     Safety factor in the step size prediction, default 0.9.

fac1     Parameters for step size selection; the new step size is chosen
fac2     subject to the restriction  fac1 <= hnew/hold <= fac2.
         Default values are fac1=0.2 and fac2=10.0.

beta     The "beta" for stabilized step size control (see section IV.2 of our
         book). Larger values for beta ( <= 0.1 ) make the step size control
         more stable. dopri5 needs a larger beta than Higham & Hall. Negative
         initial value provoke beta=0; default beta=0.04.

hmax     Maximal step size, default xend-x.

h        Initial step size, default is a guess computed by the function hinit.

nmax     Maximal number of allowed steps, default 100000.

meth     Switch for the choice of the method coefficients; at the moment the
         only possibility and default value are 1.

nstiff   Test for stiffness is activated when the current step number is a
         multiple of nstiff. A negative value means no test and the default
         is 1000.

nrdens   Number of components for which dense outpout is required, default 0.
         For 0 < nrdens < n, the components have to be specified in icont[0],
         icont[1], ... icont[nrdens-1]. Note that if nrdens=0 or nrdens=n, no
         icont is needed, pass NULL.

Memory requirements
-------------------

         The function dopri5 allocates dynamically 8*n doubles for the method
         stages, 5*nrdens doubles for the interpolation if dense output is
         performed and n uint32 if 0 < nrdens < n.

OUTPUT PARAMETERS
-----------------

y       numerical solution at x=xRead() (see below).

dopri5 returns the following values

         1 : computation successful,
         2 : computation successful interrupted by solout,
        -1 : input is not consistent,
        -2 : larger nmax is needed,
        -3 : step size becomes too small,
        -4 : the problem is probably stff (interrupted).

Several functions provide access to different values :

xRead   x value for which the solution has been computed (x=xend after
        successful return).

hRead   Predicted step size of the last accepted step (useful for a
        subsequent call to dopri5).

nstepRead   Number of used steps.
naccptRead  Number of accepted steps.
nrejctRead  Number of rejected steps.
nfcnRead    Number of function calls.

*/

#include <stdio.h>
#include <limits.h>

typedef void (*FcnEqDiff)(uint32 n, double x, double *y, double *f);
typedef void (*SolTrait)(long nr, double xold, double x, double *y, uint32 n,
                         int32 *irtrn);

extern int32
dop853(uint32 n,        /* dimension of the system <= UINT_MAX-1*/
       FcnEqDiff fcn,   /* function computing the value of f(x,y) */
       double x,        /* initial x-value */
       double *y,       /* initial values for y */
       double xend,     /* final x-value (xend-x may be positive or negative) */
       double *rtoler,  /* relative error tolerance */
       double *atoler,  /* absolute error tolerance */
       int32 itoler,    /* switch for rtoler and atoler */
       SolTrait solout, /* function providing the numerical solution during
                           integration */
       int32 iout,      /* switch for calling solout */
       FILE *fileout,   /* messages stream */
       double uround,   /* rounding unit */
       double safe,     /* safety factor */
       double fac1,     /* parameters for step size selection */
       double fac2, double beta, /* for stabilized step size control */
       double hmax,              /* maximal step size */
       double h,                 /* initial step size */
       long nmax,                /* maximal number of allowed steps */
       int32 meth,               /* switch for the choice of the coefficients */
       long nstiff,              /* test for stiffness */
       uint32 nrdens, /* number of components for which dense outpout is
                           required */
       uint32 *icont, /* indexes of components for which dense output is
                           required, >= nrdens */
       uint32 licont, /* declared length of icon */
       double *work);

extern double contd8(uint32 ii, /* index of desired component */
                     double x   /* approximation at x */
);

extern int32
dopri5(uint32 n,        /* dimension of the system <= UINT_MAX-1*/
       FcnEqDiff fcn,   /* function computing the value of f(x,y) */
       double x,        /* initial x-value */
       double *y,       /* initial values for y */
       double xend,     /* final x-value (xend-x may be positive or negative) */
       double *rtoler,  /* relative error tolerance */
       double *atoler,  /* absolute error tolerance */
       int32 itoler,    /* switch for rtoler and atoler */
       SolTrait solout, /* function providing the numerical solution during
                           integration */
       int32 iout,      /* switch for calling solout */
       FILE *fileout,   /* messages stream */
       double uround,   /* rounding unit */
       double safe,     /* safety factor */
       double fac1,     /* parameters for step size selection */
       double fac2, double beta, /* for stabilized step size control */
       double hmax,              /* maximal step size */
       double h,                 /* initial step size */
       long nmax,                /* maximal number of allowed steps */
       int32 meth,               /* switch for the choice of the coefficients */
       long nstiff,              /* test for stiffness */
       uint32 nrdens, /* number of components for which dense outpout is
                           required */
       uint32 *icont, /* indexes of components for which dense output is
                           required, >= nrdens */
       uint32 licont, /* declared length of icon */
       double *work);

extern double contd5(uint32 ii, /* index of desired component */
                     double x   /* approximation at x */
);

void dprhs(uint32 n, double t, double *y, double *f);
void dp_err(int32 k);
int32 dp(int32 *istart, double *y, double *t, int32 n, double tout, double *tol,
         double *atol, int32 flag, int32 *kflag);
int32 dormprin(int32 *istart, double *y, double *t, int32 n, double tout,
               double *tol, double *atol, int32 flag, int32 *kflag);
long nfcnRead(void);
long nstepRead(void);
long naccptRead(void);
long nrejctRead(void);
double hRead(void);
double xRead(void);
int32 dop853(uint32 n, FcnEqDiff fcn, double x, double *y, double xend,
             double *rtoler, double *atoler, int32 itoler, SolTrait solout,
             int32 iout, FILE *fileout, double uround, double safe, double fac1,
             double fac2, double beta, double hmax, double h, long nmax,
             int32 meth, long nstiff, uint32 nrdens, uint32 *icont,
             uint32 licont, double *work);
double contd8(uint32 ii, double x);
int32 dopri5(uint32 n, FcnEqDiff fcn, double x, double *y, double xend,
             double *rtoler, double *atoler, int32 itoler, SolTrait solout,
             int32 iout, FILE *fileout, double uround, double safe, double fac1,
             double fac2, double beta, double hmax, double h, long nmax,
             int32 meth, long nstiff, uint32 nrdens, uint32 *icont,
             uint32 licont, double *work);
double contd5(uint32 ii, double x);

extern long nfcnRead(void); /* encapsulation of statistical data */
extern long nstepRead(void);
extern long naccptRead(void);
extern long nrejctRead(void);
extern double hRead(void);
extern double xRead(void);

#ifndef EDIT_RHS_H
#define EDIT_RHS_H

#include <X11/Xlib.h>
#include "xpplim.h"
#include <stdio.h>

#define NEQMAXFOREDIT 20
#define MAXARG 20
#define MAX_N_EBOX MAX_ODE
#define MAX_LEN_EBOX 86
#define FORGET_ALL 0
#define DONE_ALL 2
#define FORGET_THIS 3
#define DONE_THIS 1
#define RESET_ALL 4

#define MAX_UFUN 50

/*typedef struct {
  int32 narg;
  char args[MAXARG][11];
} UFUN_ARG;
*/

/*  This is a edit box widget which handles a list of
        editable strings
 */

typedef struct {
    Window base, ok, cancel, reset;
    Window win[MAX_N_EBOX];
    char name[MAX_N_EBOX][MAX_LEN_EBOX], value[MAX_N_EBOX][MAX_LEN_EBOX],
        rval[MAX_N_EBOX][MAX_LEN_EBOX];
    int32 n, hot;
} EDIT_BOX;

void reset_ebox(EDIT_BOX *sb, int32 *pos, int32 *col);
int32 do_edit_box(int32 n, char *title, char **names, char **values);
void expose_ebox(EDIT_BOX *sb, Window w, int32 pos);
void ereset_hot(int32 inew, EDIT_BOX *sb);
void enew_editable(EDIT_BOX *sb, int32 inew, int32 *pos, int32 *col,
                   int32 *done, Window *w);
int32 e_box_event_loop(EDIT_BOX *sb, int32 *pos, int32 *col);
void make_ebox_windows(EDIT_BOX *sb, char *title);
void edit_menu(void);
void edit_rhs(void);
void user_fun_info(FILE *fp);
void edit_functions(void);
int32 save_as(void);

#endif

#ifndef EIG_LIST_H
#define EIG_LIST_H

#include <X11/Xlib.h>

void draw_eq_list(Window w);
void create_eq_list(void);
void eq_list_keypress(XEvent ev, int32 *used);
void enter_eq_stuff(Window w, int32 b);
void eq_list_button(XEvent ev);
void eq_list_up(void);
void eq_list_down(void);
void eq_box_import(void);
void get_new_size(Window win, uint32 *wid, uint32 *hgt);
void resize_eq_list(Window win);
void eq_box_button(Window w);
void create_eq_box(int32 cp, int32 cm, int32 rp, int32 rm, int32 im, double *y,
                   int32 n);
void draw_eq_box(Window w);

#endif

#ifndef EXTRA_H
#define EXTRA_H

#define MAXW 50

void load_new_dll(void);
int32 my_fun(double *in, double *out, int32 nin, int32 nout, double *v,
             double *c);
void auto_load_dll(void);
void do_in_out(void);
void add_export_list(char *in, char *out);
void check_inout(void);
int32 get_export_count(char *s);
void do_export_list(void);
void parse_inout(char *l, int32 flag);
void get_import_values(int32 n, double *ydot, char *soname, char *sofun,
                       int32 ivar, double *wgt[MAXW], double *var, double *con);

#endif

/*--------------------------------*-C-*---------------------------------*
 * File:
 *	fftn.h
 *
 * Singleton's multivariate floatcomplex Fourier transform, computed in
 * place using mixed-radix Fast Fourier Transform algorithm.
 *
 * Called here `fftn' since it does a radix-n FFT on n-dimensional data
 *
 * Copyright(c)1995,97 Mark Olesen <olesen@me.QueensU.CA>
 *		Queen's Univ at Kingston (Canada)
 *
 * Permission to use, copy, modify, and distribute this software for
 * any purpose without fee is hereby granted, provided that this
 * entire notice is included in all copies of any software which is
 * or includes a copy or modification of this software and in all
 * copies of the supporting documentation for such software.
 *
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR QUEEN'S
 * UNIVERSITY AT KINGSTON MAKES ANY REPRESENTATION OR WARRANTY OF ANY
 * KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS
 * FITNESS FOR ANY PARTICULAR PURPOSE.
 *
 * All of which is to say that you can do what you like with this
 * source code provided you don't try to sell it as your own and you
 * include an unaltered copy of this message (including the
 * copyright).
 *
 * It is also implicitly understood that bug fixes and improvements
 * should make their way back to the general Internet community so
 * that everyone benefits.
 *
 * Brief overview of parameters:
 * ---------------------------------------------------------------------*
 * Re[]:	real value array
 * Im[]:	imaginary value array
 * nTotal:	total number of floatcomplex values
 * nPass:	number of elements involved in this pass of transform
 * nSpan:	nspan/nPass = number of bytes to increment pointer
 *		in Re[] and Im[]
 * isign:	exponent: +1 = forward  -1 = reverse
 * scaling:	normalizing constant by which the final result is DIVIDED
 *	scaling == -1, normalize by total dimension of the transform
 *	scaling <  -1, normalize by the square-root of the total dimension
 *
 *
 * Slightly more detailed information:
 * ----------------------------------------------------------------------*
 * void fft_free (void);
 *
 * free-up allocated temporary storage after finished all the Fourier
 * transforms.
 *
 * ----------------------------------------------------------------------*
 *
 * int32 fftn (int32 ndim,  int32 dims[], REAL Re[], REAL Im[],
 *	    int32 iSign, double scaling);
 *
 * NDIM = the total number dimensions
 * DIMS = a vector of array sizes
 *	if NDIM is zero then DIMS must be zero-terminated
 *
 * RE and IM hold the real and imaginary components of the data, and
 * return the resulting real and imaginary Fourier coefficients.
 * Multidimensional data *must* be allocated contiguously.  There is
 * no limit on the number of dimensions.
 *
 * ISIGN = the sign of the floatcomplex exponential
 *	(ie, forward or inverse FFT)
 *	the magnitude of ISIGN (normally 1) is used to determine
 *	the correct indexing increment (see below).
 *
 * SCALING = normalizing constant by which the final result is DIVIDED
 *	if SCALING == -1, normalize by total dimension of the transform
 *	if SCALING <  -1, normalize by the square-root of the total dimension
 *
 * example:
 * tri-variate transform with Re[n3][n2][n1], Im[n3][n2][n1]
 *
 *	int32 dims[3] = {n1,n2,n3}
 *	fftn (3, dims, Re, Im, 1, scaling);
 *
 * or, using a null terminated dimension list
 *	int32 dims[4] = {n1,n2,n3,0}
 *	fftn (0, dims, Re, Im, 1, scaling);
 * ----------------------------------------------------------------------*/
#ifndef FFTN_H
#define FFTN_H

extern void fft_free(void);

/* double precision routine */
extern int32 fftn(int32 /* ndim */, int32 /* dims */[], double /* Re */[],
                  double /* Im */[], int32 /* isign */, double /* scaling */);

/* double precision routine */
extern int32 fftnf(int32 /* ndim */, int32 /* dims */[], double /* Re */[],
                   double /* Im */[], int32 /* isign */, double /* scaling */);

#endif /* _FFTN_H */

/*----------------------- end-of-file (C header) -----------------------*/
#ifndef FLAGS_H
#define FLAGS_H

int32 add_global(char *cond, int32 sign, char *rest);
void show_flags(void);
int32 compile_flags(void);
int32 one_flag_step(double *yold, double *ynew, int32 *istart, double told,
                    double *tnew, int32 neq, double *s);
int32 one_flag_step_symp(double *y, double dt, double *work, int32 neq,
                         double *tim, int32 *istart);
int32 one_flag_step_euler(double *y, double dt, double *work, int32 neq,
                          double *tim, int32 *istart);
int32 one_flag_step_discrete(double *y, double dt, double *work, int32 neq,
                             double *tim, int32 *istart);
int32 one_flag_step_heun(double *y, double dt, double *yval[2], int32 neq,
                         double *tim, int32 *istart);
int32 one_flag_step_rk4(double *y, double dt, double *yval[3], int32 neq,
                        double *tim, int32 *istart);
void printflaginfo(void);
int32 one_flag_step_gear(int32 neq, double *t, double tout, double *y,
                         double hmin, double hmax, double eps, int32 mf,
                         double *error, int32 *kflag, int32 *jstart,
                         double *work, int32 *iwork);
int32 one_flag_step_rosen(double *y, double *tstart, double tfinal,
                          int32 *istart, int32 n, double *work, int32 *ierr);
int32 one_flag_step_dp(int32 *istart, double *y, double *t, int32 n,
                       double tout, double *tol, double *atol, int32 flag,
                       int32 *kflag);
int32 one_flag_step_adap(double *y, int32 neq, double *t, double tout,
                         double eps, double *hguess, double hmin, double *work,
                         int32 *ier, double epjac, int32 iflag, int32 *jstart);
int32 one_flag_step_backeul(double *y, double *t, double dt, int32 neq,
                            double *yg, double *yp, double *yp2, double *ytemp,
                            double *errvec, double *jac, int32 *istart);
int32 one_flag_step_cvode(int32 *command, double *y, double *t, int32 n,
                          double tout, int32 *kflag, double *atol,
                          double *rtol);

#endif

#ifndef FORM_ODE_H
#define FORM_ODE_H

#include "xpplim.h"
#include "newpars.h"
#include <stdio.h>

#define MAXVNAM 33
#define MAXLINES 5000

/*void break_up_list(char *rhs);
void compile_em();
void free_varinfo();
void remove_blanks(char *s1);
void read_a_line(FILE *fp,char *s);

void subsk(char *big,char *new,int32 k,int32 flag);
void free_comments();

void add_comment(char *s);
void init_varinfo();
void add_varinfo(int32 type,char *lhs,char *rhs,int32 nargs,char
args[MAXARG][NAMLEN+1]); void stor_internopts(char *s1);
*/

typedef struct {
    char *name, *value;
} FIXINFO;

int32 make_eqn(void);
void strip_saveqn(void);
int32 disc(char *string);
void dump_src(void);
void dump_comments(void);
void format_list(char **s, int32 n);
int32 get_a_filename(char *filename, char *wild);
void list_em(char *wild);
int32 read_eqn(void);
int32 get_eqn(FILE *fptr);
int32 compiler(char *bob, FILE *fptr);
void list_upar(void);
void welcome(void);
void show_syms(void);
void take_apart(char *bob, double *value, char *name);
char *get_first(char *string, char *src);
char *get_next(char *src);
void find_ker(char *string, int32 *alt);
void pos_prn(char *s);
void clrscr(void);
int32 getuch(void);
int32 getchi(void);
int32 if_include_file(char *old, char *nf);
int32 if_end_include(char *old);
int32 do_new_parser(FILE *fp, char *first, int32 nnn);
void create_plot_list(void);
void add_only(char *s);
void break_up_list(char *rhs);
int32 find_the_name(char list[MAX_ODE1][MAXVNAM], int32 n, char *name);
void compile_em(void);
int32 formula_or_number(char *expr, double *z);
void strpiece(char *dest, char *src, int32 i0, int32 ie);
int32 parse_a_string(char *s1, VAR_INFO *v);
void init_varinfo(void);
void add_varinfo(int32 type, char *lhs, char *rhs, int32 nargs,
                 char args[20][13 + 1]);
void free_varinfo(void);
int32 extract_ode(char *s1, int32 *ie, int32 i1);
int32 strparse(char *s1, char *s2, int32 i0, int32 *i1);
int32 extract_args(char *s1, int32 i0, int32 *ie, int32 *narg,
                   char args[20][13 + 1]);
int32 find_char(char *s1, char *s2, int32 i0, int32 *i1);
int32 next_nonspace(char *s1, int32 i0, int32 *i1);
void remove_blanks(char *s1);
void read_a_line(FILE *fp, char *s);
int32 search_array(char *old, char *new, int32 *i1, int32 *i2, int32 *flag);
int32 check_if_ic(char *big);
int32 not_ker(char *s, int32 i);
int32 is_comment(char *s);
void subsk(char *big, char *new, int32 k, int32 flag);
void keep_orig_comments(void);
void default_comments(void);
void free_comments(void);
void new_comment(FILE *f);
void add_comment(char *s);

/* for parsing par, init with whitespace correctly */
char *new_string2(char *old, int32 length);
void advance_past_first_word(char **sptr);
char *get_next2(char **tokens_ptr);
void strcpy_trim(char *dest, char *source);
void strncpy_trim(char *dest, char *source, int32 n);
#endif

#ifndef GEAR_H
#define GEAR_H

void silent_fixpt(double *x, double eps, double err, double big, int32 maxit,
                  int32 n, double *er, double *em, int32 *ierr);
void do_sing(double *x, double eps, double err, double big, int32 maxit,
             int32 n, int32 *ierr, double *stabinfo);
void do_sing_info(double *x, double eps, double err, double big, int32 maxit,
                  int32 n, double *er, double *em, int32 *ierr);

void shoot_this_now(void);
void pr_evec(double *x, double *ev, int32 n, int32 type);
void get_complex_evec(double *m, double evr, double evm, double *br, double *bm,
                      int32 n, int32 maxit, double err, int32 *ierr);
void get_evec(double *a, double *anew, double *b, double *bp, int32 n,
              int32 maxit, double err, int32 *ipivot, double eval, int32 *ierr);
void eigen(int32 n, double *a, double *ev, double *work, int32 *ierr);
void hqrx(int32 n, int32 low, int32 igh, double *h, double *ev, int32 *ierr);
void orthesx(int32 n, int32 low, int32 igh, double *a, double *ort);
double sign(double x, double y);
int32 imin(int32 x, int32 y);
double amax(double u, double v);
void getjactrans(double *x, double *y, double *yp, double *xp, double eps,
                 double *d, int32 n);
void getjac(double *x, double *y, double *yp, double *xp, double eps,
            double *dermat, int32 n);
void rooter(double *x, double err, double eps, double big, double *work,
            int32 *ierr, int32 maxit, int32 n);
double sqr2(double z);
int32 gear(int32 n, double *t, double tout, double *y, double hmin, double hmax,
           double eps, int32 mf, double *error, int32 *kflag, int32 *jstart,
           double *work, int32 *iwork);
int32 ggear(int32 n, double *t, double tout, double *y, double hmin,
            double hmax, double eps, int32 mf, double *error, int32 *kflag,
            int32 *jstart, double *work, int32 *iwork);
double sgnum(double x, double y);
double Max(double x, double y);
double Min(double x, double y);
void sgefa(double *a, int32 lda, int32 n, int32 *ipvt, int32 *info);
void sgesl(double *a, int32 lda, int32 n, int32 *ipvt, double *b, int32 job);
void saxpy(int32 n, double sa, double *sx, int32 incx, double *sy, int32 incy);
int32 isamax(int32 n, double *sx, int32 incx);
double sdot(int32 n, double *sx, int32 incx, int32 incy);
void sscal(int32 n, double sa, double *sx, int32 incx);
void save_batch_shoot(void);

#endif

#ifndef GGETS_H
#define GGETS_H

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
void bottom_msg(char *msg);
void gputs(char *string, Window win);
void err_msg(char *string);
int32 plintf(char *fmt, ...);
int32 show_position(XEvent ev);
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
void display_command(char *name, char *value, int32 pos);
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

int32 go_go_auto(void);

#ifndef GRAF_PAR_H
#define GRAF_PAR_H

#define RUBBOX 0
#define RUBLINE 1

#define SCRNFMT 0
#define PSFMT 1
#define SVGFMT 2

#define REAL_SMALL 1.e-6
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
    double *x[MAXBIFCRV], *y[MAXBIFCRV];
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
void add_bd_crv(double *x, double *y, int32 len, int32 type, int32 ncrv);
void frz_bd(void);
void read_bd(FILE *fp);
int32 get_frz_index(Window w);
void export_graf_data(void);
void add_a_curve_com(int32 c);
void default_window(void);
void dump_ps(int32 i);

#endif

#ifndef GRAPHICS_H
#define GRAPHICS_H

void get_scale(double *x1, double *y1, double *x2, double *y2);
void set_scale(double x1, double y1, double x2, double y2);
void get_draw_area_flag(int32 flag);
void get_draw_area(void);
void change_current_linestyle(int32 new, int32 *old);
void set_normal_scale(void);
void point(int32 x, int32 y);
void line(int32 x1, int32 y1, int32 x2, int32 y2);
void bead(int32 x1, int32 y1);
void frect(int32 x1, int32 y1, int32 w, int32 h);
void put_text(int32 x, int32 y, char *str);
void init_x11(void);
void init_ps(void);
void init_svg(void);
void point_x11(int32 xp, int32 yp);
void set_linestyle(int32 ls);
void set_line_style_x11(int32 ls);
void bead_x11(int32 x, int32 y);
void rect_x11(int32 x, int32 y, int32 w, int32 h);
void line_x11(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void put_text_x11(int32 x, int32 y, char *str);
void special_put_text_x11(int32 x, int32 y, char *str, int32 size);
void fancy_put_text_x11(int32 x, int32 y, char *str, int32 size, int32 font);
void scale_dxdy(double x, double y, double *i, double *j);
void scale_to_screen(double x, double y, int32 *i, int32 *j);
void scale_to_real(int32 i, int32 j, double *x, double *y);
void init_all_graph(void);
void set_extra_graphs(void);
void reset_graph(void);
void get_graph(void);
void init_graph(int32 i);
void copy_graph(int32 i, int32 l);
void make_rot(double theta, double phi);
void scale3d(double x, double y, double z, double *xp, double *yp, double *zp);
double proj3d(double theta, double phi, double x, double y, double z, int32 in);
int32 threedproj(double x2p, double y2p, double z2p, double *xp, double *yp);
void text3d(double x, double y, double z, char *s);
void text_3d(double x, double y, double z, char *s);
int32 threed_proj(double x, double y, double z, double *xp, double *yp);
void point_3d(double x, double y, double z);
void line3dn(double xs1, double ys1, double zs1, double xsp1, double ysp1,
             double zsp1);
void line3d(double x01, double y01, double z01, double x02, double y02,
            double z02);
void line_3d(double x, double y, double z, double xp, double yp, double zp);
void pers_line(double x, double y, double z, double xp, double yp, double zp);
void rot_3dvec(double x, double y, double z, double *xp, double *yp,
               double *zp);
void point_abs(double x1, double y1);
void line_nabs(double x1_out, double y1_out, double x2_out, double y2_out);
void bead_abs(double x1, double y1);
void frect_abs(double x1, double y1, double w, double h);
void line_abs(double x1, double y1, double x2, double y2);
void text_abs(double x, double y, char *text);
void fillintext(char *old, char *new);
void fancy_text_abs(double x, double y, char *old, int32 size);
int32 clip3d(double x1, double y1, double z1, double x2, double y2, double z2,
             double *x1p, double *y1p, double *z1p, double *x2p, double *y2p,
             double *z2p);
int32 clip(double x1, double x2, double y1, double y2, double *x1_out,
           double *y1_out, double *x2_out, double *y2_out);
void eq_symb(double *x, int32 type);
void draw_symbol(double x, double y, double size, int32 my_symb);
void reset_all_line_type(void);
void draw_many_lines(void);

#endif

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

int32 two_d_hist(int32 col1, int32 col2, int32 ndat, int32 n1, int32 n2,
                 double xlo, double xhi, double ylo, double yhi);
void four_back(void);
void hist_back(void);
void new_four(int32 nmodes, int32 col);
int32 new_2d_hist(void);
void new_hist(int32 nbins, double zlo, double zhi, int32 col, int32 col2,
              char *condition, int32 which);
void column_mean(void);
int32 get_col_info(int32 *col, char *prompt);
void compute_power(void);
int32 spectrum(double *data, int32 nr, int32 win, int32 w_type, double *pow);
int32 cross_spectrum(double *data, double *data2, int32 nr, int32 win,
                     int32 w_type, double *pow, int32 type);
void compute_sd(void);
void compute_fourier(void);
void compute_correl(void);
void compute_stacor(void);
void mycor(double *x, double *y, int32 n, double zlo, double zhi, int32 nbins,
           double *z, int32 flag);
void mycor2(double *x, double *y, int32 n, int32 nbins, double *z, int32 flag);
void compute_hist(void);
void sft(double *data, double *ct, double *st, int32 nmodes, int32 grid);
void fftxcorr(double *data1, double *data2, int32 length, int32 nlag,
              double *cr, int32 flag);
void fft(double *data, double *ct, double *st, int32 nmodes, int32 length);
void just_fourier(int32 flag);
void just_sd(int32 flag);
void post_process_stuff(void);

#endif

#ifndef INIT_CONDS_H
#define INIT_CONDS_H

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
void do_slide_button(uint32 w, PAR_SLIDER *p);
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
void do_slide_release(uint32 w, PAR_SLIDER *p);
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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdio.h>

void init_ar_ic(void);
void dump_range(FILE *fp, int32 f);
void init_range(void);
int32 set_up_eq_range(void);
void cont_integ(void);
int32 range_item(void);
int32 range_item2(void);
int32 set_up_range(void);
int32 set_up_range2(void);
void init_monte_carlo(void);
void monte_carlo(void);
void do_monte_carlo_search(int32 append, int32 stuffbrowse, int32 ishoot);
void do_eq_range(double *x);
void swap_color(int32 *col, int32 rorw);
void set_cycle(int32 flag, int32 *icol);
int32 do_range(double *x, int32 flag);
void find_equilib_com(int32 com);
void batch_integrate(void);
void do_batch_dry_run(void);
void batch_integrate_once(void);
int32 write_this_run(char *file, int32 i);
void do_init_data(int32 com);
void run_from_x(double *x);
void run_now(void);
void do_start_flags(double *x, double *t);
void usual_integrate_stuff(double *x);
void do_new_array_ic(char *new, int32 j1, int32 j2);
void store_new_array_ic(char *new, int32 j1, int32 j2, char *formula);
void evaluate_ar_ic(char *v, char *f, int32 j1, int32 j2);
int32 extract_ic_data(char *big);
void arr_ic_start(void);
int32 set_array_ic(void);
int32 form_ic(void);
void get_ic(int32 it, double *x);
int32 ode_int(double *y, double *t, int32 *istart, int32 ishow);
int32 integrate(double *t, double *x, double tend, double dt, int32 count,
                int32 nout, int32 *start);
void send_halt(void);
void send_output(double *y, double t);
void do_plot(double *oldxpl, double *oldypl, double *oldzpl, double *xpl,
             double *ypl, double *zpl);
void export_data(FILE *fp);
void plot_the_graphs(double *xv, double *xvold, double ddt, int32 *tc,
                     int32 flag);
void plot_one_graph(double *xv, double *xvold, double ddt, int32 *tc);
void restore(int32 i1, int32 i2);
void comp_color(double *v1, double *v2, int32 n, double dt);
void shoot(double *x, double *xg, double *evec, int32 sgn);
void shoot_easy(double *x);
void stop_integration(void);
int32 stor_full(void);
int32 do_auto_range_go(void);
void silent_equilibria(void);

#endif

/******************************************************************
 *                                                                *
 * File          : iterativ.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This header file contains declarations intended for use by     *
 * generic iterative solvers of Ax = b. The enumeration gives     *
 * symbolic names for the type  of preconditioning to be used.    *
 * The function type declarations give the prototypes for the     *
 * functions to be called within an iterative linear solver, that *
 * are responsible for                                            *
 *    multiplying A by a given vector v (ATimesFn), and           *
 *    solving the preconditioner equation Pz = r (PSolveFn).      *
 *                                                                *
 ******************************************************************/

#ifndef ITERATIV_H
#define ITERATIV_H

#include "llnltyps.h"
#include "vector.h"

/******************************************************************
 *                                                                *
 * enum : types of preconditioning                                *
 *----------------------------------------------------------------*
 * NONE  : The iterative linear solver should not use             *
 *         preconditioning.                                       *
 *                                                                *
 * LEFT  : The iterative linear solver uses preconditioning on    *
 *         the left only.                                         *
 *                                                                *
 * RIGHT : The iterative linear solver uses preconditioning on    *
 *         the right only.                                        *
 *                                                                *
 * BOTH  : The iterative linear solver uses preconditioning on    *
 *         both the left and the right.                           *
 *                                                                *
 ******************************************************************/

enum {
    NONE,
    LEFT,
    RIGHT,
    BOTH
};

/******************************************************************
 *                                                                *
 * enum : types of Gram-Schmidt routines                          *
 *----------------------------------------------------------------*
 * MODIFIED_GS  : The iterative solver uses the modified          *
 *                Gram-Schmidt routine ModifiedGS listed in this  *
 *                file.                                           *
 *                                                                *
 * CLASSICAL_GS : The iterative solver uses the classical         *
 *                Gram-Schmidt routine ClassicalGS listed in this *
 *                file.                                           *
 *                                                                *
 ******************************************************************/

enum {
    MODIFIED_GS,
    CLASSICAL_GS
};

/******************************************************************
 *                                                                *
 * Type: ATimesFn                                                 *
 *----------------------------------------------------------------*
 * An ATimesFn multiplies Av and stores the result in z. The      *
 * caller is responsible for allocating memory for the z vector.  *
 * The parameter A_data is a pointer to any information about A   *
 * which the function needs in order to do its job. The vector v  *
 * is unchanged. An ATimesFn returns 0 if successful and a        *
 * non-zero value if unsuccessful.                                *
 *                                                                *
 ******************************************************************/

typedef int32 (*ATimesFn)(void *A_data, N_Vector v, N_Vector z);

/******************************************************************
 *                                                                *
 * Type: PSolveFn                                                 *
 *----------------------------------------------------------------*
 * A PSolveFn solves the preconditioner equation Pz = r for the   *
 * vector z. The caller is responsible for allocating memory for  *
 * the z vector. The parameter P_data is a pointer to any         *
 * information about P which the function needs in order to do    *
 * its job. The parameter lr indicates whether P is a left        *
 * preconditioner or a right preconditioner. The vector r is      *
 * unchanged.  A PSolveFn returns 0 if successful and a non-zero  *
 * value if unsuccessful.  On a failure, a negative return value  *
 * indicates an unrecoverable condition, while a positive value   *
 * indicates a recoverable one, in which the calling routine may  *
 * reattempt the solution after updating preconditioner data.     *
 *                                                                *
 ******************************************************************/

typedef int32 (*PSolveFn)(void *P_data, N_Vector r, N_Vector z, int32 lr);

/******************************************************************
 *                                                                *
 * Function: ModifiedGS                                           *
 *----------------------------------------------------------------*
 * ModifiedGS performs a modified Gram-Schmidt orthogonalization  *
 * of the N_Vector v[k] against the p unit N_Vectors at           *
 * v[k-1], v[k-2], ..., v[k-p].                                   *
 *                                                                *
 * v is an array of (k+1) N_Vectors v[i], i=0, 1, ..., k.         *
 * v[k-1], v[k-2], ..., v[k-p] are assumed to have L2-norm        *
 * equal to 1.                                                    *
 *                                                                *
 * h is the output k by k Hessenberg matrix of inner products.    *
 * This matrix must be allocated row-wise so that the (i,j)th     *
 * entry is h[i][j]. The inner products (v[i],v[k]),              *
 * i=i0, i0+1, ..., k-1, are stored at h[i][k-1]. Here            *
 * i0=MAX(0,k-p).                                                 *
 *                                                                *
 * k is the index of the vector in the v array that needs to be   *
 * orthogonalized against previous vectors in the v array.        *
 *                                                                *
 * p is the number of previous vectors in the v array against     *
 * which v[k] is to be orthogonalized.                            *
 *                                                                *
 * new_vk_norm is a pointer to memory allocated by the caller to  *
 * hold the Euclidean norm of the orthogonalized vector v[k].     *
 *                                                                *
 * If (k-p) < 0, then ModifiedGS uses p=k. The orthogonalized     *
 * v[k] is NOT normalized and is stored over the old v[k]. Once   *
 * the orthogonalization has been performed, the Euclidean norm   *
 * of v[k] is stored in (*new_vk_norm).                           *
 *                                                                *
 * ModifiedGS returns 0 to indicate success. It cannot fail.      *
 *                                                                *
 ******************************************************************/

int32 ModifiedGS(N_Vector *v, double **h, int32 k, int32 p,
                 double *new_vk_norm);

/******************************************************************
 *                                                                *
 * Function: ClassicalGS                                          *
 *----------------------------------------------------------------*
 * ClassicalGS performs a classical Gram-Schmidt                  *
 * orthogonalization of the N_Vector v[k] against the p unit      *
 * N_Vectors at v[k-1], v[k-2], ..., v[k-p]. The parameters v, h, *
 * k, p, and new_vk_norm are as described in the documentation    *
 * for ModifiedGS.                                                *
 *                                                                *
 * temp is an N_Vector which can be used as workspace by the      *
 * ClassicalGS routine.                                           *
 *                                                                *
 * s is a length k array of reals which can be used as workspace  *
 * by the ClassicalGS routine.                                    *
 *                                                                *
 * ClassicalGS returns 0 to indicate success. It cannot fail.     *
 *                                                                *
 ******************************************************************/

int32 ClassicalGS(N_Vector *v, double **h, int32 k, int32 p,
                  double *new_vk_norm, N_Vector temp, double *s);

/******************************************************************
 *                                                                *
 * Function: QRfact                                               *
 *----------------------------------------------------------------*
 * QRfact performs a QR factorization of the Hessenberg matrix H. *
 *                                                                *
 * n is the problem size; the matrix H is (n+1) by n.             *
 *                                                                *
 * h is the (n+1) by n Hessenberg matrix H to be factored. It is  *
 * stored row-wise.                                               *
 *                                                                *
 * q is an array of length 2*n containing the Givens rotations    *
 * computed by this function. A Givens rotation has the form:     *
 * | c  -s |                                                      *
 * | s   c |.                                                     *
 * The components of the Givens rotations are stored in q as      *
 * (c, s, c, s, ..., c, s).                                       *
 *                                                                *
 * job is a control flag. If job==0, then a new QR factorization  *
 * is performed. If job!=0, then it is assumed that the first     *
 * n-1 columns of h have already been factored and only the last  *
 * column needs to be updated.                                    *
 *                                                                *
 * QRfact returns 0 if successful. If a zero is encountered on    *
 * the diagonal of the triangular factor R, then QRfact returns   *
 * the equation number of the zero entry, where the equations are *
 * numbered from 1, not 0. If QRsol is subsequently called in     *
 * this situation, it will return an error because it could not   *
 * divide by the zero diagonal entry.                             *
 *                                                                *
 ******************************************************************/

int32 QRfact(int32 n, double **h, double *q, int32 job);

/******************************************************************
 *                                                                *
 * Function: QRsol                                                *
 *----------------------------------------------------------------*
 * QRsol solves the linear least squares problem                  *
 *                                                                *
 * min (b - H*x, b - H*x), x in R^n,                              *
 *                                                                *
 * where H is a Hessenberg matrix, and b is in R^(n+1).           *
 * It uses the QR factors of H computed by QRfact.                *
 *                                                                *
 * n is the problem size; the matrix H is (n+1) by n.             *
 *                                                                *
 * h is a matrix (computed by QRfact) containing the upper        *
 * triangular factor R of the original Hessenberg matrix H.       *
 *                                                                *
 * q is an array of length 2*n (computed by QRfact) containing    *
 * the Givens rotations used to factor H.                         *
 *                                                                *
 * b is the (n+1)-vector appearing in the least squares problem   *
 * above.                                                         *
 *                                                                *
 * On return, b contains the solution x of the least squares      *
 * problem, if QRsol was successful.                              *
 *                                                                *
 * QRsol returns a 0 if successful.  Otherwise, a zero was        *
 * encountered on the diagonal of the triangular factor R.        *
 * In this case, QRsol returns the equation number (numbered      *
 * from 1, not 0) of the zero entry.                              *
 *                                                                *
 ******************************************************************/

int32 QRsol(int32 n, double **h, double *q, double *b);

#endif

#ifndef KINESCOPE_H
#define KINESCOPE_H

void do_movie_com(int32 c);
void reset_film(void);
int32 film_clip(void);
int32 show_frame(int32 i, int32 h, int32 w);
void play_back(void);
void save_kine(void);
void make_anigif(void);
void save_movie(char *basename, int32 fmat);
void auto_play(void);
void too_small(void);

#endif

/******************************************************************
 *                                                                *
 * File          : llnlmath.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a C math library. The routines     *
 * listed here work with the type double as defined in llnltyps.h.  *
 * To do single precision floating point arithmetic, set the type *
 * double to be double. To do double precision arithmetic, set the   *
 * type double to be double. The default implementations for        *
 * RPowerR and RSqrt call standard math library functions which   *
 * do double precision arithmetic. If this is unacceptable when   *
 * double is double, then the user should re-implement these two     *
 * routines by calling single precision routines available on     *
 * his/her machine.                                               *
 *                                                                *
 ******************************************************************/

#ifndef LLNLMATH_H
#define LLNLMATH_H

#include "llnltyps.h"

/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 *                                                                *
 * MAX(A, B) returns the maximum of A and B.                      *
 *                                                                *
 * ABS(A) returns the absolute value of A.                        *
 *                                                                *
 * SQR(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/

#define MIN(A, B) ((A) < (B) ? (A) : (B))

#define MAX(A, B) ((A) > (B) ? (A) : (B))

#define ABS(A) ((A > 0) ? (A) : -(A))

#define SQR(A) ((A) * (A))

/******************************************************************
 *                                                                *
 * Function : UnitRoundoff                                        *
 * Usage    : double uround;                                        *
 *            uround = UnitRoundoff();                            *
 *----------------------------------------------------------------*
 * UnitRoundoff returns the unit roundoff u for double floating     *
 * point arithmetic, where u is defined to be the largest         *
 * positive double such that 1.0 + u != 1.0.                        *
 *                                                                *
 ******************************************************************/

double UnitRoundoff(void);

/******************************************************************
 *                                                                *
 * Function : RPowerI                                             *
 * Usage    : int32 exponent;                                       *
 *            double base, ans;                                     *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerI returns the value base^exponent, where base is a double  *
 * and exponent is an int32.                                        *
 *                                                                *
 ******************************************************************/

double RPowerI(double base, int32 exponent);

/******************************************************************
 *                                                                *
 * Function : RPowerR                                             *
 * Usage    : double base, exponent, ans;                           *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerR returns the value base^exponent, where both base and   *
 * exponent are reals. If base < 0.0, then RPowerR returns 0.0.   *
 *                                                                *
 ******************************************************************/

double RPowerR(double base, double exponent);

/******************************************************************
 *                                                                *
 * Function : RSqrt                                               *
 * Usage    : double sqrt_x;                                        *
 *            sqrt_x = RSqrt(x);                                  *
 *----------------------------------------------------------------*
 * RSqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
 * returns 0.0.                                                   *
 *                                                                *
 ******************************************************************/

double RSqrt(double x);

#endif

#ifndef LOAD_EQN_H
#define LOAD_EQN_H

#include <stdio.h>

/*
The acutual max filename length is determined by the
FILENAME_MAX (see <stdio.h>), and usually 4096 -- but
this is huge and usually overkill.  On the otherhand
the old Xpp default string buffer size of 100 is a bit
restricitive for lengths of filenames. You could also
set this define in the Makefile or at compile time to
override the below definition.
*/

#ifndef XPP_MAX_NAME
#define XPP_MAX_NAME 300
#if (XPP_MAX_NAME > FILENAME_MAX)
#define XPP_MAX_NAME FILENAME_MAX
#endif

#endif

/*
Options are set accroding to an order of precedence

command line < mfile < .xpprc < default.opt

Add any options here that you might want to track.
*/
typedef struct {
    int32 BIG_FONT_NAME;
    int32 SMALL_FONT_NAME;
    int32 BACKGROUND;
    int32 IXPLT;
    int32 IYPLT;
    int32 IZPLT;
    int32 AXES;
    int32 NMESH;
    int32 METHOD;
    int32 TIMEPLOT;
    int32 MAXSTOR;
    int32 TEND;
    int32 DT;
    int32 T0;
    int32 TRANS;
    int32 BOUND;
    int32 TOLER;
    int32 DELAY;
    int32 XLO;
    int32 XHI;
    int32 YLO;
    int32 YHI;
    int32 UserBlack;
    int32 UserWhite;
    int32 UserMainWinColor;
    int32 UserDrawWinColor;
    int32 UserGradients;
    int32 UserBGBitmap;
    int32 UserMinWidth;
    int32 UserMinHeight;
    int32 YNullColor;
    int32 XNullColor;
    int32 StableManifoldColor;
    int32 UnstableManifoldColor;
    int32 START_LINE_TYPE;
    int32 RandSeed;
    int32 PaperWhite;
    int32 COLORMAP;
    int32 NPLOT;
    int32 DLL_LIB;
    int32 DLL_FUN;
    int32 XP;
    int32 YP;
    int32 ZP;
    int32 NOUT;
    int32 VMAXPTS;
    int32 TOR_PER;
    int32 JAC_EPS;
    int32 NEWT_TOL;
    int32 NEWT_ITER;
    int32 FOLD;
    int32 DTMIN;
    int32 DTMAX;
    int32 ATOL;
    int32 TOL;
    int32 BANDUP;
    int32 BANDLO;
    int32 PHI;
    int32 THETA;
    int32 XMIN;
    int32 XMAX;
    int32 YMIN;
    int32 YMAX;
    int32 ZMIN;
    int32 ZMAX;
    int32 POIVAR;
    int32 OUTPUT;
    int32 POISGN;
    int32 POIEXT;
    int32 POISTOP;
    int32 STOCH;
    int32 POIPLN;
    int32 POIMAP;
    int32 RANGEOVER;
    int32 RANGESTEP;
    int32 RANGELOW;
    int32 RANGEHIGH;
    int32 RANGERESET;
    int32 RANGEOLDIC;
    int32 RANGE;
    int32 NTST;
    int32 NMAX;
    int32 NPR;
    int32 NCOL;
    int32 DSMIN;
    int32 DSMAX;
    int32 DS;
    int32 PARMAX;
    int32 NORMMIN;
    int32 NORMMAX;
    int32 EPSL;
    int32 EPSU;
    int32 EPSS;
    int32 RUNNOW;
    int32 SEC;
    int32 UEC;
    int32 SPC;
    int32 UPC;
    int32 AUTOEVAL;
    int32 AUTOXMAX;
    int32 AUTOYMAX;
    int32 AUTOXMIN;
    int32 AUTOYMIN;
    int32 AUTOVAR;
    int32 PS_FONT;
    int32 PS_LW;
    int32 PS_FSIZE;
    int32 PS_COLOR;
    int32 FOREVER;
    int32 BVP_TOL;
    int32 BVP_EPS;
    int32 BVP_MAXIT;
    int32 BVP_FLAG;
    int32 SOS;
    int32 FFT;
    int32 HIST;
    int32 PltFmtFlag;
    int32 ATOLER;
    int32 MaxEulIter;
    int32 EulTol;
    int32 EVEC_ITER;
    int32 EVEC_ERR;
    int32 NULL_ERR;
    int32 NEWT_ERR;
    int32 NULL_HERE;
    int32 TUTORIAL;
    int32 SLIDER1;
    int32 SLIDER2;
    int32 SLIDER3;
    int32 SLIDER1LO;
    int32 SLIDER2LO;
    int32 SLIDER3LO;
    int32 SLIDER1HI;
    int32 SLIDER2HI;
    int32 SLIDER3HI;
    int32 POSTPROCESS;
    int32 HISTCOL;
    int32 HISTLO;
    int32 HISTHI;
    int32 HISTBINS;
    int32 SPECCOL;
    int32 SPECCOL2;
    int32 SPECWIDTH;
    int32 SPECWIN;
    int32 PLOTFORMAT;
    int32 DFGRID;
    int32 DFBATCH;
    int32 NCBATCH;
    int32 COLORVIA;
    int32 COLORIZE;
    int32 COLORLO;
    int32 COLORHI;
    int32 HISTCOL2;
    int32 HISTLO2;
    int32 HISTHI2;
    int32 HISTBINS2;
} OptionsSet;

void dump_torus(FILE *fp, int32 f);
void load_eqn(void);
void set_X_vals(void);
void set_all_vals(void);
void read_defaults(FILE *fp);
void fil_flt(FILE *fpt, double *val);
void fil_int(FILE *fpt, int32 *val);
void add_intern_set(char *name, char *does);
void extract_action(char *ptr);
void extract_internset(int32 j);
void do_intern_set(char *name1, char *value);
int32 msc(char *s1, char *s2);
void set_internopts(OptionsSet *mask);
void set_internopts_xpprc_and_comline(void);
void split_apart(char *bob, char *name, char *value);
void check_for_xpprc(void);
void stor_internopts(char *s1);
void set_option(char *s1, char *s2, int32 force, OptionsSet *mask);

#endif

#ifndef LUNCH_NEW_H
#define LUNCH_NEW_H

#include <stdio.h>

void file_inf(void);
void ps_write_pars(FILE *fp);
void do_info(FILE *fp);
int32 read_lunch(FILE *fp);
void do_lunch(int32 f);
void dump_eqn(FILE *fp);
void io_numerics(int32 f, FILE *fp);
void io_parameter_file(char *fn, int32 flag);
void io_ic_file(char *fn, int32 flag);
void io_parameters(int32 f, FILE *fp);
void io_exprs(int32 f, FILE *fp);
void io_graph(int32 f, FILE *fp);
void io_int(int32 *i, FILE *fp, int32 f, char *ss);
void io_double(double *z, FILE *fp, int32 f, char *ss);
void io_float(double *z, FILE *fp, int32 f, char *ss);
void io_string(char *s, int32 len, FILE *fp, int32 f);

#endif

#ifndef main_h_

#define main_h_

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
                uint32 min_wid, uint32 min_hgt, int32 argc, char **argv);
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
int32 get_command_width(void);

#endif

#ifndef MANY_POPS_H
#define MANY_POPS_H

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
int32 man_xy(double *xe, double *ye);
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

#ifndef MARKOV_H
#define MARKOV_H

#include <stdio.h>

void add_wiener(int32 index);
void set_wieners(double dt, double *x, double t);
void add_markov(int32 nstate, char *name);
int32 build_markov(char **ma, char *name);
int32 old_build_markov(FILE *fptr, char *name);
void extract_expr(char *source, char *dest, int32 *i0);
void create_markov(int32 nstates, double *st, int32 type, char *name);
void add_markov_entry(int32 index, int32 j, int32 k, char *expr);
void compile_all_markov(void);
int32 compile_markov(int32 index, int32 j, int32 k);
void update_markov(double *x, double t, double dt);
double new_state(double old, int32 index, double dt);
void make_gill_nu(double *nu, int32 n, int32 m, double *v);
void one_gill_step(int32 meth, int32 nrxn, int32 *rxn, double *v);
void do_stochast_com(int32 i);
void mean_back(void);
void variance_back(void);
void compute_em(void);
void free_stoch(void);
void init_stoch(int32 len);
void append_stoch(int32 first, int32 length);
void do_stats(int32 ierr);
double gammln(double xx);
double poidev(double xm);
double ndrand48(void);
void nsrand48(int32 seed);
double ran1(long *idum);

#endif

#ifndef MENUDRIVE_H
#define MENUDRIVE_H

#include <X11/Xlib.h>

#define M_IR 0
#define M_I2 1
#define M_IL 2
#define M_IO 3
#define M_IG 4
#define M_IM 5
#define M_IS 6
#define M_IN 7
#define M_IH 8
#define M_IF 9
#define M_IU 10
#define M_II 11
#define M_ID 12
#define M_IB 13

#define MAX_M_I 13

#define M_C 20

#define M_NN 31
#define M_NR 32
#define M_NA 33
#define M_NM 34
#define M_NF 35
#define M_NS 36

#define M_NFF 37
#define M_NFD 38
#define M_NFR 39
#define M_NFA 40

#define M_SG 50
#define M_SM 51
#define M_SR 52
#define M_SC 53

#define M_DD 60
#define M_DF 61
#define M_DN 62
#define M_DC 63
#define M_DS 64

#define M_WW 70
#define M_WZ 71
#define M_WO 72
#define M_WF 73
#define M_WD 74
#define M_WS 75

#define M_AA 80
#define M_AN 81
#define M_AC 82

#define M_KC 90
#define M_KR 91
#define M_KP 92
#define M_KA 93
#define M_KS 94
#define M_KM 95
#define M_KX 96

#define M_GA 100
#define M_GD 101
#define M_GR 102
#define M_GE 103
#define M_GP 104
#define M_GV 105
#define M_GF 106
#define M_GX 107
#define M_GO 108
#define M_GC 109

#define M_GFF 110
#define M_GFD 111
#define M_GFE 112
#define M_GFR 113
#define M_GFK 114
#define M_GFB 115
#define M_GFC 116
#define M_GFO 117

#define M_GFKN 120
#define M_GFKK 121

#define M_GCN 130
#define M_GCP 131
#define M_GCH 132
#define M_GCC 133
#define M_GCB 134
#define M_GCG 135
#define M_GCU 136

#define M_P 140

#define M_EE 141

#define M_R 142

#define M_X 143

#define M_3 144

#define M_MC 150
#define M_MK 151
#define M_MD 152
#define M_MB 153
#define M_MA 154
#define M_MM 155
#define M_MS 156

#define M_TT 160
#define M_TA 161
#define M_TP 162
#define M_TM 163
#define M_TE 164
#define M_TD 165
#define M_TS 166
#define M_TEM 170
#define M_TEC 171
#define M_TED 172

#define M_V2 180
#define M_V3 181
#define M_VA 182
#define M_VT 183

#define M_BR 190
#define M_BN 191
#define M_BS 192
#define M_BP 193
#define M_BH 194

#define M_FP 200
#define M_FW 201
#define M_FR 202
#define M_FA 203
#define M_FC 204
#define M_FS 205
#define M_FB 206
#define M_FH 207
#define M_FQ 208
#define M_FT 209
#define M_FI 210
#define M_FG 211

#define M_FER 212
#define M_FEF 213
#define M_FES 214
#define M_FEL 215

#define M_FX 216
#define M_FU 217

/* CLONE change ! */
#define M_FL 218

/*  some numerics commands */

#define M_UAN 300
#define M_UAM 301
#define M_UAA 302
#define M_UAO 303
#define M_UAH 304
#define M_UAP 305
#define M_UAR 306

#define M_UCN 310
#define M_UCV 311
#define M_UCA 312

#define M_UPN 320
#define M_UPS 321
#define M_UPM 322
#define M_UPP 323

#define M_UHN 330
#define M_UHC 331
#define M_UHD 332
#define M_UHM 333
#define M_UHV 334
#define M_UHH 335
#define M_UHO 336
#define M_UHF 337
#define M_UHP 338
#define M_UHI 339
#define M_UHS 340
#define M_UHL 341
#define M_UHA 342
#define M_UHX 343
#define M_UHE 344
#define M_UH2 345

#define M_UKE 350
#define M_UKV 351

/*  one shot numerics commands */

#define M_UT 400
#define M_US 401
#define M_UR 402
#define M_UD 403
#define M_UN 404
#define M_UV 405
#define M_UI 406
#define M_UO 407
#define M_UB 408
#define M_UE 409
#define M_UC 410

void xpp_hlp(void);
void MessageBox(char *m);
void RedrawMessageBox(Window w);
void KillMessageBox(void);
int32 TwoChoice(char *c1, char *c2, char *q, char *key);
int32 GetMouseXY(int32 *x, int32 *y);
void FlushDisplay(void);
void clear_draw_window(void);
void drw_all_scrns(void);
void clr_all_scrns(void);
void run_the_commands(int32 com);
void do_stochast(void);
void get_pmap_pars(void);
void set_col_par(void);
void make_adj(void);
void do_file_com(int32 com);
void do_gr_objs(void);
void new_lookup(void);
void find_bvp(void);
void change_view(void);
void do_windows(void);
void add_a_curve(void);
void do_movie(void);
void do_torus(void);
void window_zoom(void);
void direct_field(void);
void new_clines(void);
void froz_cline_stuff(void);
void find_equilibrium(void);
void ini_data_menu(void);
void new_param(void);
void clear_screens(void);
void x_vs_t(void);
void redraw_them_all(void);
void get_3d_par(void);
void edit_xpprc(void);
void do_tutorial(void);

#endif

#ifndef XPPMENU_H
#define XPPMENU_H

#include <X11/Xlib.h>

void flash(int32 num);
void add_menu(Window base, int32 j, int32 n, char **names, char *key,
              char **hint);
void create_the_menus(Window base);
void show_menu(int32 j);
void unshow_menu(int32 j);
void help(void);
void help_num(void);
void help_file(void);
void menu_crossing(Window win, int32 yn);
void menu_expose(Window win);
void menu_button(Window win);
void draw_help(void);

#endif

#ifndef MY_PS_H
#define MY_PS_H

int32 ps_init(char *filename, int32 color);
void ps_stroke(void);
void ps_do_color(int32 color);
void ps_setcolor(int32 color);
void ps_end(void);
void ps_frect(int32 x, int32 y, int32 w, int32 h);
void ps_last_pt_off(void);
void ps_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void chk_ps_lines(void);
void ps_linetype(int32 linetype);
void ps_point(int32 x, int32 y);
void ps_write(char *str);
void ps_fnt(int32 cf, int32 scale);
void ps_show(char *str, int32 type);
void ps_abs(int32 x, int32 y);
void ps_rel(int32 x, int32 y);
void special_put_text_ps(int32 x, int32 y, char *str, int32 size);
void fancy_ps_text(int32 x, int32 y, char *str, int32 size, int32 font);
void ps_text(int32 x, int32 y, char *str);

#endif

#ifndef MY_RHS_H
#define MY_RHS_H

int32 MAIN__(void);
int32 main(int32 argc, char **argv);
void extra(double *y__y, double t, int32 nod, int32 neq);
void set_fix_rhs(double t, double *y);
int32 my_rhs(double t, double *y, double *ydot, int32 neq);
void update_based_on_current(void);
void fix_only(void);
void rhs_only(double *ydot);
void vec_rhs(double t, double *y, double *ydot, int32 neq);

#endif

#ifndef MY_SVG_H
#define MY_SVG_H

int32 svg_init(char *filename);
void svg_stroke(void);
void svg_do_color(int32 color);
void svg_end(void);
void svg_bead(void);
void svg_frect(int32 x, int32 y, int32 w, int32 h);
void svg_last_pt_off(void);
void svg_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void chk_svg_lines(void);
void svg_linetype(int32 linetype);
void svg_point(int32 x, int32 y);
void svg_write(char *str);
void special_put_text_svg(int32 x, int32 y, char *str, int32 size);
void svg_text(int32 x, int32 y, char *str);

#endif

#ifndef NULLCLINE_H
#define NULLCLINE_H

#include <stdio.h>

typedef struct {
    double x, y, z;
} Pt;

typedef struct nclines {
    double *xn, *yn;
    int32 nmx, nmy;
    int32 n_ix, n_iy;
    struct nclines *n, *p;
} NCLINES;

void create_new_cline(void);
void froz_cline_stuff_com(int32 i);
void do_range_clines(void);
void start_ncline(void);
void clear_froz_cline(void);
int32 get_nullcline_floats(double **v, int32 *n, int32 who, int32 type);
void save_frozen_clines(char *fn);
void redraw_froz_cline(int32 flag);
void add_froz_cline(double *xn, int32 nmx, int32 n_ix, double *yn, int32 nmy,
                    int32 n_iy);
void get_max_dfield(double *y, double *ydot, double u0, double v0, double du,
                    double dv, int32 n, int32 inx, int32 iny, double *mdf);
void redraw_dfield(void);
void direct_field_com(int32 c);
void save_the_nullclines(void);
void restore_nullclines(void);
void dump_clines(FILE *fp, double *x, int32 nx, double *y, int32 ny);
void dump_clines_old(FILE *fp, double *x, int32 nx, double *y, int32 ny);
void restor_null(double *v, int32 n, int32 d);
void new_clines_com(int32 c);
void new_nullcline(int32 course, double xlo, double ylo, double xhi, double yhi,
                   double *stor, int32 *npts);
void stor_null(double x1, double y1, double x2, double y2);
double fnull(double x, double y);
int32 interpolate(Pt p1, Pt p2, double z, double *x, double *y);
void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4);
void triangle_contour(Pt p1, Pt p2, Pt p3);
void do_cline(int32 ngrid, double x1, double y1, double x2, double y2);
void do_batch_nclines(void);
void do_batch_dfield(void);
void set_colorization_stuff(void);
void silent_nullclines(void);
void silent_dfields(void);

#endif

#ifndef NUMERICS_H
#define NUMERICS_H

/*       Numerics.h   */

extern double DELTA_T, TEND, T0, TRANS, NULL_ERR, EVEC_ERR, NEWT_ERR;
extern double BOUND, DELAY, TOLER, HMIN, HMAX;
extern double *fft_data, *hist_data, color_scale, min_scale;
extern double POIPLN;

extern int32 NMESH, NJMP, METHOD, color_flag, NC_ITER;
extern int32 EVEC_ITER, FOREVER;

extern int32 POIMAP, POIVAR, POISGN, SOS;

extern int32 HIST, HVAR, hist_ind;

extern int32 XSHFT, YSHFT, ZSHFT;

void chk_volterra(void);
void check_pos(int32 *j);
void quick_num(int32 com);
void get_num_par(int32 ch);
void chk_delay(void);
void set_delay(void);
void ruelle(void);
void init_numerics(void);
void meth_dialog(void);
void get_pmap_pars_com(int32 l);
void get_method(void);
void set_col_par_com(int32 i);
void do_meth(void);
void set_total(double total);
void user_set_color_par(int32 flag, char *via, double lo, double hi);
void compute_one_period(double period, double *x, char *name);
#endif

#ifndef odesol2_h
#define odesol2_h

int32 symplect3(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 discrete(double *y, double *tim, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 bak_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 one_bak_step(double *y, double *t, double dt, int32 neq, double *yg,
                   double *yp, double *yp2, double *ytemp, double *errvec,
                   double *jac);
void one_step_discrete(double *y, double dt, double *yp, int32 neq, double *t);
void one_step_symp(double *y, double h, double *f, int32 n, double *t);
void one_step_euler(double *y, double dt, double *yp, int32 neq, double *t);
void one_step_rk4(double *y, double dt, double *yval[3], int32 neq,
                  double *tim);
void one_step_heun(double *y, double dt, double *yval[2], int32 neq,
                   double *tim);
int32 euler(double *y, double *tim, double dt, int32 nt, int32 neq,
            int32 *istart, double *work);
int32 mod_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 rung_kut(double *y, double *tim, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 adams(double *y, double *tim, double dt, int32 nstep, int32 neq,
            int32 *ist, double *work);
int32 abmpc(double *y, double *t, double dt, int32 neq);
int32 rb23(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
           double *work, int32 *ierr);
int32 rosen(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
            double *work, int32 *ierr);
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy,
                 int32 neq, double eps, double scal);
void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold,
                  int32 n, double eps, double scal);
int32 bandfac(double *a, int32 ml, int32 mr, int32 n);
void bandsol(double *a, double *b, int32 ml, int32 mr, int32 n);

#endif

#ifndef POP_LIST_H
#define POP_LIST_H

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include "xpplim.h"

#define MAX_N_SBOX 22

#define FORGET_ALL 0
#define DONE_ALL 2
#define FORGET_THIS 3
#define DONE_THIS 1

#define EV_MASK                                                                \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

#define BUT_MASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     EnterWindowMask | LeaveWindowMask)

extern Display *display;
extern int32 DisplayWidth, DisplayHeight;
extern int32 screen;
extern Atom deleteWindowAtom;
extern Window main_win, info_pop, draw_win;
extern int32 DCURY, DCURX, CURY_OFF, DCURXs, DCURYs, CURY_OFFs, xor_flag;
extern GC gc, small_gc;
extern uint32 MyBackColor, MyForeColor;

extern int32 TipsFlag;
extern char UserBlack[8];
extern char UserWhite[8];
extern int32 UserGradients;

/*  This is a string box widget which handles a list of
        editable strings
 */

typedef struct {
    Window base, ok, cancel;
    Window win[MAX_N_SBOX];
    char name[MAX_N_SBOX][MAX_LEN_SBOX], value[MAX_N_SBOX][MAX_LEN_SBOX];
    int32 n, hot;
    int32 hgt, wid;
    int32 hh[MAX_N_SBOX];
} STRING_BOX;

typedef struct {
    char **list;
    int32 n;
} SCRBOX_LIST;

extern int32 NUPAR, NEQ, NODE, NMarkov;
extern char upar_names[MAX_PAR][14], uvar_names[MAX_ODE][12];
extern char *color_names[];
extern SCRBOX_LIST scrbox_list[10];

/*  This is a new improved pop_up widget */
typedef struct {
    Window base, tit;
    Window *w;
    char *title;
    char **entries;
    char **hints;
    int32 n, max;
    char *key;
    int32 hot;
} POP_UP;

typedef struct {
    Window base, slide, close, text;
    int32 i0;
    int32 exist, len, nlines;
    char **list;
} TEXTWIN;

typedef struct {
    Window base, slide;
    Window *w;
    int32 nw, nent, i0;
    int32 len, exist;
    char **list;
} SCROLLBOX;

extern TEXTWIN mytext;
#define SB_PLOTTABLE 0
#define SB_VARIABLE 1
#define SB_PARAMETER 2
#define SB_PARVAR 3
#define SB_COLOR 4
#define SB_MARKER 5
#define SB_METHOD 6

void set_window_title(Window win, char *string);
void make_scrbox_lists(void);
int32 get_x_coord_win(Window win);
void destroy_scroll_box(SCROLLBOX *sb);
void create_scroll_box(Window root, int32 x0, int32 y0, int32 nent, int32 nw,
                       char **list, SCROLLBOX *sb);
void expose_scroll_box(Window w, SCROLLBOX sb);
void redraw_scroll_box(SCROLLBOX sb);
void crossing_scroll_box(Window w, int32 c, SCROLLBOX sb);
int32 scroll_box_motion(XEvent ev, SCROLLBOX *sb);
int32 select_scroll_item(Window w, SCROLLBOX sb);
void scroll_popup(STRING_BOX *sb, SCROLLBOX *scrb);
int32 do_string_box(int32 n, int32 row, int32 col, char *title, char **names,
                    char values[][MAX_LEN_SBOX], int32 maxchar);
void expose_sbox(STRING_BOX sb, Window w, int32 pos);
void do_hilite_text(char *name, char *value, int32 flag, Window w, int32 pos);
void reset_hot(int32 inew, STRING_BOX *sb);
void new_editable(STRING_BOX *sb, int32 inew, int32 *pos, int32 *col,
                  int32 *done, Window *w);
void set_sbox_item(STRING_BOX *sb, int32 item);
int32 s_box_event_loop(STRING_BOX *sb, int32 *pos, int32 *col, SCROLLBOX *scrb);
void make_sbox_windows(STRING_BOX *sb, int32 row, int32 col, char *title,
                       int32 maxchar);
Window make_fancy_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw);
Window make_unmapped_window(Window root, int32 x, int32 y, int32 width,
                            int32 height, int32 bw);
Window make_plain_unmapped_window(Window root, int32 x, int32 y, int32 width,
                                  int32 height, int32 bw);
Window make_window(Window root, int32 x, int32 y, int32 width, int32 height,
                   int32 bw);
Window make_plain_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw);
void expose_resp_box(char *button, char *message, Window wb, Window wm,
                     Window w);
void respond_box(char *button, char *message);
void message_box(Window *w, int32 x, int32 y, char *message);
void expose_choice(char *choice1, char *choice2, char *msg, Window c1,
                   Window c2, Window wm, Window w);
int32 two_choice(char *choice1, char *choice2, char *string, char *key, int32 x,
                 int32 y, Window w, char *title);
int32 yes_no_box(void);
int32 pop_up_list(Window *root, char *title, char **list, char *key, int32 n,
                  int32 max, int32 def, int32 x, int32 y, char **hints,
                  Window hwin, char *httxt);
void draw_pop_up(POP_UP p, Window w);
Window make_unmapped_icon_window(Window root, int32 x, int32 y, int32 width,
                                 int32 height, int32 bw, uchar *icdata);
Window make_icon_window(Window root, int32 x, int32 y, int32 width,
                        int32 height, int32 bw, uchar *icdata);

#endif

#ifndef PP_SHOOT_H
#define PP_SHOOT_H

#include <stdio.h>

void do_bc(double *y__0, double t0, double *y__1, double t1, double *f,
           int32 n);
void compile_bvp(void);
void reset_bvp(void);
void init_shoot_range(char *s);
void dump_shoot_range(FILE *fp, int32 f);
void bad_shoot(int32 iret);
void do_sh_range(double *ystart, double *yend);
int32 set_up_homoclinic(void);
int32 set_up_periodic(int32 *ipar, int32 *ivar, double *sect, int32 *ishow);
void find_bvp_com(int32 com);
void last_shot(int32 flag);
int32 set_up_sh_range(void);
void bvshoot(double *y, double *yend, double err, double eps, int32 maxit,
             int32 *iret, int32 n, int32 ishow, int32 iper, int32 ipar,
             int32 ivar, double sect);

#endif

#ifndef RUBBER_H
#define RUBBER_H

#include <X11/Xlib.h>

int32 rubber(int32 *x1, int32 *y1, int32 *x2, int32 *y2, Window w, int32 f);
void rbox(int32 i1, int32 j1, int32 i2, int32 j2, Window w, int32 f);

#endif

#ifndef SCRNGIF_H
#define SCRNGIF_H

#include <stdio.h>
#include <X11/Xlib.h>

typedef struct GifTree {
    char typ;   /* terminating, lookup, or search */
    int32 code; /* the code to be output */
    uchar ix;   /* the color map index */
    struct GifTree **node, *nxt, *alt;
} GifTree;

typedef struct {
    uchar r, g, b;
} GIFCOL;

void set_global_map(int32 flag);
int32 ppmtopix(uchar r, uchar g, uchar b, int32 *n);
void end_ani_gif(FILE *fp);
void add_ani_gif(Window win, FILE *fp, int32 count);
void screen_to_gif(Window win, FILE *fp);
void get_global_colormap(Window win);
void local_to_global(void);
int32 use_global_map(uchar *pixels, uchar *ppm, int32 h, int32 w);
int32 make_local_map(uchar *pixels, uchar *ppm, int32 h, int32 w);
void gif_stuff(Window win, FILE *fp, int32 task);
void write_global_header(int32 cols, int32 rows, FILE *dst);
void GifLoop(FILE *fout, uint32 repeats);
void write_local_header(int32 cols, int32 rows, FILE *fout, int32 colflag,
                        int32 delay);
void make_gif(uchar *pixels, int32 cols, int32 rows, FILE *dst);
int32 GifEncode(FILE *fout, uchar *pixels, int32 depth, int32 siz);
void ClearTree(int32 cc, GifTree *root);
uchar *AddCodeToBuffer(int32 code, int16 n, uchar *buf);

#endif

#ifndef SIMPLENET_H
#define SIMPLENET_H

double net_interp(double x, int32 i);
double network_value(double x, int32 i);
void init_net(double *v, int32 n);
int32 add_spec_fun(char *name, char *rhs);
void add_special_name(char *name, char *rhs);
int32 add_vectorizer(char *name, char *rhs);
void add_vectorizer_name(char *name, char *rhs);
int32 is_network(char *s);
void eval_all_nets(void);
void evaluate_network(int32 ind);
void update_all_ffts(void);
void update_fft(int32 ind);
void fft_conv(int32 it, int32 n, double *values, double *yy, double *fftr,
              double *ffti, double *dr, double *di);
int32 gilparse(char *s, int32 *ind, int32 *nn);
int32 g_namelist(char *s, char *root, int32 *flag, int32 *i1, int32 *i2);
double vector_value(double x, int32 i);

#endif

/******************************************************************
 *                                                                *
 * File          : spgmr.h                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for the implementation of the scaled   *
 * preconditioned GMRES (SPGMR) iterative linear solver.          *
 *                                                                *
 ******************************************************************/

#ifndef SPGMR_H
#define SPGMR_H

#include "llnltyps.h"
#include "vector.h"

/******************************************************************
 *                                                                *
 * Types: SpgmrMemRec, SpgmrMem                                   *
 *----------------------------------------------------------------*
 * SpgmrMem is a pointer to an SpgmrMemRec which contains         *
 * the memory needed by SpgmrSolve. The SpgmrMalloc routine       *
 * returns a pointer of type SpgmrMem which should then be passed *
 * in subsequent calls to SpgmrSolve. The SpgmrFree routine frees *
 * the memory allocated by SpgmrMalloc.                           *
 *                                                                *
 * N is the linear system size.                                   *
 *                                                                *
 * l_max is the maximum Krylov dimension that SpgmrSolve will be  *
 * permitted to use.                                              *
 *                                                                *
 * V is the array of Krylov basis vectors v_1, ..., v_(l_max+1),  *
 * stored in V[0], ..., V[l_max], where l_max is the second       *
 * parameter to SpgmrMalloc. Each v_i is a length N vector of     *
 * type N_Vector. (N is the first parameter to SpgmrMalloc and    *
 * represents the size of the linear system.)                     *
 *                                                                *
 * Hes is the (l_max+1) x l_max Hessenberg matrix. It is stored   *
 * row-wise so that the (i,j)th element is given by Hes[i][j].    *
 *                                                                *
 * givens is a length 2*l_max array which represents the          *
 * Givens rotation matrices that arise in the algorithm. The      *
 * Givens rotation matrices F_0, F_1, ..., F_j, where F_i is      *
 *                                                                *
 *             1                                                  *
 *               1                                                *
 *                 c_i  -s_i      <--- row i                      *
 *                 s_i   c_i                                      *
 *                           1                                    *
 *                             1                                  *
 *                                                                *
 * are represented in the givens vector as                        *
 * givens[0]=c_0, givens[1]=s_0, givens[2]=c_1, givens[3]=s_1,    *
 * ..., givens[2j]=c_j, givens[2j+1]=s_j.                         *
 *                                                                *
 * xcor is a length N vector (type N_Vector) which holds the      *
 * scaled, preconditioned correction to the initial guess.        *
 *                                                                *
 * yg is a length (l_max+1) array of reals used to hold "int16"   *
 * vectors (e.g. y and g).                                        *
 *                                                                *
 * vtemp is a length N vector (type N_Vector) used as temporary   *
 * vector storage during calculations.                            *
 *                                                                *
 ******************************************************************/

typedef struct {
    int64 N;
    int32 l_max;

    N_Vector *V;
    double **Hes;
    double *givens;
    N_Vector xcor;
    double *yg;
    N_Vector vtemp;
} SpgmrMemRec, *SpgmrMem;

/******************************************************************
 *                                                                *
 * Function : SpgmrMalloc                                         *
 *----------------------------------------------------------------*
 * SpgmrMalloc allocates the memory used by SpgmrSolve. It        *
 * returns a pointer of type SpgmrMem which the user of the       *
 * SPGMR package should pass to SpgmrSolve. The parameter N       *
 * is the size of the system to be solved by SpgmrSolve and l_max *
 * is the maximum Krylov dimension that SpgmrSolve will be        *
 * permitted to use. The parameter machEnv is a pointer to        *
 * machine environment-specific information. Pass NULL in the     *
 * ordinary sequential case (see vector.h). This routine returns  *
 * NULL if there is a memory request failure.                     *
 *                                                                *
 ******************************************************************/

SpgmrMem SpgmrMalloc(int64 N, int32 l_max);

/******************************************************************
 *                                                                *
 * Function : SpgmrSolve                                          *
 *----------------------------------------------------------------*
 * SpgmrSolve solves the linear system Ax = b using the SPGMR     *
 * method. The return values are given by the symbolic constants  *
 * below. The first SpgmrSolve parameter is a pointer to memory   *
 * allocated by a prior call to SpgmrMalloc. The system size N    *
 * passed in the call to SpgmrMalloc should be the same as the    *
 * length of all N_Vector arguments passed to SpgmrSolve.         *
 *                                                                *
 * mem is the pointer returned by SpgmrMalloc to the structure    *
 * containing the memory needed by SpgmrSolve.                    *
 *                                                                *
 * A_data is a pointer to information about the coefficient       *
 * matrix A. This pointer is passed to the user-supplied function *
 * atimes.                                                        *
 *                                                                *
 * x is the initial guess x_0 upon entry and the solution         *
 * N_Vector upon exit with return value SPGMR_SUCCESS or          *
 * SPGMR_RES_REDUCED. For all other return values, the output x   *
 * is undefined.                                                  *
 *                                                                *
 * b is the right hand side N_Vector. It is undisturbed by this   *
 * function.                                                      *
 *                                                                *
 * pretype is the type of preconditioning to be used. Its         *
 * legal possible values are enumerated in iterativ.h. These      *
 * values are NONE=0, LEFT=1, RIGHT=2, and BOTH=3.                *
 *                                                                *
 * gstype is the type of Gram-Schmidt orthogonalization to be     *
 * used. Its legal values are enumerated in iterativ.h. These     *
 * values are MODIFIED_GS=0 and CLASSICAL_GS=1.                   *
 *                                                                *
 * delta is the tolerance on the L2 norm of the scaled,           *
 * preconditioned residual. On return with value SPGMR_SUCCESS,   *
 * this residual satisfies || sb P1_inv (b - Ax) ||_2 <= delta.   *
 *                                                                *
 * max_restarts is the maximum number of times the algorithm is   *
 * allowed to restart.                                            *
 *                                                                *
 * P_data is a pointer to preconditioner information. This        *
 * pointer is passed to the user-supplied function psolve.        *
 *                                                                *
 * sx is the N_Vector of positive scale factors for x (not        *
 * tested). Pass NULL if x scaling not required.                  *
 *                                                                *
 * sb is the N_Vector of positive scale factors for b (not        *
 * tested). Pass NULL if b scaling not required.                  *
 *                                                                *
 * atimes is the user-supplied function which performs the        *
 * operation of multiplying A by a given vector. Its description  *
 * is given in iterativ.h.                                        *
 *                                                                *
 * psolve is the user-supplied function which solves a            *
 * preconditioned equation Pz = r. Its description is also given  *
 * in iterativ.h. The psolve function will not be called if       *
 * pretype is NONE. In this case, the user should pass NULL for   *
 * psolve.                                                        *
 *                                                                *
 * res_norm is a pointer to the L2 norm of the scaled,            *
 * preconditioned residual. On return with value SPGMR_SUCCESS or *
 * SPGMR_RES_REDUCED, (*res_norm) contains the value              *
 * || sb P1_inv (b - Ax) ||_2. Here x is the computed solution,   *
 * sb is the diagonal scaling matrix for the right hand side b,   *
 * and P1_inv is the inverse of the left preconditioner matrix.   *
 * For all other return values, (*res_norm) is undefined. The     *
 * caller is responsible for allocating the memory (*res_norm)    *
 * to be filled in by SpgmrSolve.                                 *
 *                                                                *
 * nli is a pointer to the number of linear iterations done in    *
 * the execution of SpgmrSolve. The caller is responsible for     *
 * allocating the memory (*nli) to be filled in by SpgmrSolve.    *
 *                                                                *
 * nps is a pointer to the number of calls made to psolve during  *
 * the execution of SpgmrSolve. The caller is responsible for     *
 * allocating the memory (*nps) to be filled in by SpgmrSolve.    *
 *                                                                *
 * Note.. Repeated calls can be made to SpgmrSolve with varying   *
 * input arguments. If, however, the problem size N or the        *
 * maximum Krylov dimension l_max changes, then a call to         *
 * SpgmrMalloc must be made to obtain new memory for SpgmrSolve   *
 * to use.                                                        *
 *                                                                *
 ******************************************************************/

int32 SpgmrSolve(SpgmrMem mem, void *A_data, N_Vector x, N_Vector b,
                 int32 pretype, int32 gstype, double delta, int32 max_restarts,
                 void *P_data, N_Vector sx, N_Vector sb, ATimesFn atimes,
                 PSolveFn psolve, double *res_norm, int32 *nli, int32 *nps);

/* Return values for SpgmrSolve */

#define SPGMR_SUCCESS 0 /* Converged                    */
#define SPGMR_RES_REDUCED                                                      \
    1                              /* Did not converge, but reduced            \
                                      norm of residual             */
#define SPGMR_CONV_FAIL 2          /* Failed to converge           */
#define SPGMR_QRFACT_FAIL 3        /* QRfact found singular matrix */
#define SPGMR_PSOLVE_FAIL_REC 4    /* psolve failed recoverably    */
#define SPGMR_MEM_NULL -1          /* mem argument is NULL         */
#define SPGMR_ATIMES_FAIL -2       /* atimes returned failure flag */
#define SPGMR_PSOLVE_FAIL_UNREC -3 /* psolve failed unrecoverably  */
#define SPGMR_GS_FAIL                                                          \
    -4                      /* Gram-Schmidt routine                            \
                               returned failure flag        */
#define SPGMR_QRSOL_FAIL -5 /* QRsol found singular R       */

/******************************************************************
 *                                                                *
 * Function : SpgmrFree                                           *
 *----------------------------------------------------------------*
 * SpgmrMalloc frees the memory allocated by SpgmrMalloc. It is   *
 * illegal to use the pointer mem after a call to SpgmrFree.      *
 *                                                                *
 ******************************************************************/

void SpgmrFree(SpgmrMem mem);

#endif

#ifndef STIFF_H
#define STIFF_H

void jacobn(double x, double *y, double *dfdx, double *dermat, double eps,
            double *work, int32 n);
int32 adaptive(double *ystart, int32 nvar, double *xs, double x2, double eps,
               double *hguess, double hmin, double *work, int32 *ier,
               double epjac, int32 iflag, int32 *jstart);
int32 gadaptive(double *ystart, int32 nvar, double *xs, double x2, double eps,
                double *hguess, double hmin, double *work, int32 *ier,
                double epjac, int32 iflag);
int32 stiff(double y[], double dydx[], int32 n, double *x, double htry,
            double eps, double yscal[], double *hdid, double *hnext,
            double *work, double epjac, int32 *ier);
int32 rkqs(double *y, double *dydx, int32 n, double *x, double htry, double eps,
           double *yscal, double *hdid, double *hnext, double *work,
           int32 *ier);
void rkck(double *y, double *dydx, int32 n, double x, double h, double *yout,
          double *yerr, double *work);

#endif

#ifndef STORAGE_H
#define STORAGE_H

void init_alloc_info(void);
void alloc_meth(void);
void init_stor(int32 nrow, int32 ncol);
void free_storage(int32 ncol);
int32 reallocstor(int32 ncol, int32 nrow);

#endif

#ifndef TABULAR_H
#define TABULAR_H

void set_auto_eval_flags(int32 f);
void set_table_name(char *name, int32 index);
void view_table(int32 index);
void new_lookup_com(int32 i);
void new_lookup_ok(void);
double lookupxy(double x, int32 n, double *xv, double *yv);
double tab_interp(double xlo, double h, double x, double *y, int32 i);
double lookup(double x, int32 index);
void init_table(void);
void redo_all_fun_tables(void);
int32 eval_fun_table(int32 n, double xlo, double xhi, char *formula, double *y);
int32 create_fun_table(int32 npts, double xlo, double xhi, char *formula,
                       int32 index);
int32 load_table(char *filename, int32 index);
int32 get_lookup_len(int32 i);

#endif

#ifndef TORUS_H
#define TORUS_H

#include <X11/Xlib.h>

void do_torus_com(int32 c);
void draw_tor_var(int32 i);
void draw_torus_box(Window win);
void choose_torus(void);
void make_tor_box(char *title);
void do_torus_events(void);

#endif

#ifndef TXTREAD_H
#define TXTREAD_H

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

#ifndef USERBUT_H
#define USERBUT_H

#include <X11/Xlib.h>

typedef struct {
    Window w;
    char bname[10];
    int32 com;
} USERBUT;

void user_button_events(XEvent report);
void user_button_press(Window w);
void user_button_draw(Window w);
void user_button_cross(Window w, int32 b);
int32 get_button_info(char *s, char *bname, char *sc);
int32 find_kbs(char *sc);
void add_user_button(char *s);
void create_user_buttons(int32 x0, int32 y0, Window base);

#endif

/****************************************************************
 *                                                              *
 * File          : vector.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL  *
 * Last Modified : 1 September 1994                             *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the header file for a generic VECTOR package. It     *
 * exports the type N_Vector.                                   *
 *                                                              *
 * Part I of this file contains declarations which are specific *
 * to the particular machine environment in which this version  *
 * of the vector package is to be used. This includes the       *
 * typedef for the type N_Vector, as well as accessor macros    *
 * that allow the user to use efficiently the type N_Vector     *
 * without making explicit references to its underlying         *
 * representation. The underlying type of N_Vector will always  *
 * be some pointer type.                                        *
 *                                                              *
 * Part II of this file contains the prototypes for the vector  *
 * kernels which operate on the type N_Vector. These prototypes *
 * are fixed for all implementations of the vector package. The *
 * definitions of the types double and int64 are in the header  *
 * file llnltyps.h and these may be changed according to the    *
 * user's needs. The llnltyps.h file also contains the          *
 * definition for the type bool (int16 for boolean) that is the *
 * return type for the routine N_VInvTest.                      *
 *                                                              *
 * Important Note: N_Vector arguments to arithmetic kernels     *
 * need not be distinct. Thus, for example, the call            *
 *                                                              *
 *         N_VLinearSum(a,x,b,y,y);    y <- ax+by               *
 *                                                              *
 * is legal.                                                    *
 *                                                              *
 * This version of vector.h is for the ordinary sequential      *
 * machine environment. In the documentation given below, N is  *
 * the length of all N_Vector parameters and x[i] denotes the   *
 * ith component of the N_Vector x, where 0 <= i <= N-1.        *
 *                                                              *
 ****************************************************************/

#ifndef volterra2_h
#define volterra2_h

double ker_val(int32 in);
void alloc_v_memory(void);
void allocate_volterra(int32 npts, int32 flag);
void re_evaluate_kernels(void);
void alloc_kernels(int32 flag);
void init_sums(double t0, int32 n, double dt, int32 i0, int32 iend,
               int32 ishift);
double alpha1n(double mu, double dt, double t, double t0);
double alpbetjn(double mu, double dt, int32 l);
double betnn(double mu, double dt);
void get_kn(double *y, double t);
int32 volterra(double *y, double *t, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 volt_step(double *y, double t, double dt, int32 neq, double *yg,
                double *yp, double *yp2, double *errvec, double *jac);

#endif

#endif
