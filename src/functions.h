#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "stdbool.h"
#include "integers.h"
#include "X11/Xlib.h"

extern Display *display;

extern char *timesfonts[];
extern char *symbolfonts[];

extern char *null_hint[];
extern char *null_freeze[];
extern char *ic_hint[];
extern char *wind_hint[];
extern char *flow_hint[];
extern char *phas_hint[];
extern char *kin_hint[];
extern char *graf_hint[];
extern char *cmap_hint[];
extern char *frz_hint[];
extern char *stoch_hint[];
extern char *bvp_hint[];
extern char *adj_hint[];
extern char *map_hint[];
extern char *view_hint[];
extern char *half_hint[];
extern char *text_hint[];
extern char *edit_hint[];
extern char *sing_hint[];
extern char *meth_hint[];
extern char *color_hint[];
extern char *tab_hint[];
extern char *edrh_hint[];
extern char *auto_hint[];
extern char *no_hint[];
extern char *aaxes_hint[];
extern char *afile_hint[];
extern char *aspecial_hint[];
extern char *arun_hint[];
extern char *browse_hint[];

#define MAIN_MENU 0
#define FILE_MENU 1
#define NUM_MENU 2

#define LENGTH(X) (sizeof(X) / sizeof(*X))

#define SETVAR(i, x)                                                           \
    do {                                                                       \
        if ((i) < NVAR)                                                        \
            variables[(i)] = (x);                                              \
    } while (0)
#define GETVAR(i) (i) < NVAR ? variables[(i)] : 0.0

#define MAX_LEN_SBOX 30

#ifndef ADJ2_H
#define ADJ2_H

#include <stdio.h>

extern bool adj_range;

int32 adj2_do_transpose(void);
void adj2_alloc_h_stuff(void);
void adj2_alloc_liap(int32 n);
void adj2_do_liapunov(void);
void adj2_do_this_liaprun(int32 i, double p);
void adj2_dump_h_stuff(FILE *fp, int32 f);
void adj2_dump_transpose_info(FILE *fp, int32 f);
void adj2_init_trans(void);
void adj2_make_adj_com(int32 com);
void adj2_new_adjoint(void);
void adj2_new_h_fun(int32 silent);
void adj2_data_back(void);

#endif

#ifndef ANI_H
#define ANI_H

#include <X11/Xlib.h>
#include <stdio.h>

/**************  New stuff for the Grabber ***************************/
#define MAX_GEVENTS 20  /* maximum variables you can change per grabbable */
#define MAX_ANI_GRAB 50 /* max grabbable objects  */

/* tasks have the form {name1=formula1;name2=formula2;...} */
typedef struct GrabTask {
    double vrhs[MAX_GEVENTS];
    char lhsname[MAX_GEVENTS][11];
    int32 lhsivar[MAX_GEVENTS];
    int32 *comrhs[MAX_GEVENTS];
    int32 runnow;
    int32 n; /* number of tasks <= MAX_GEVENTS */
} GrabTask;

typedef struct AniGrab {
    int32 ok;
    double zx;
    double zy;
    double tol;
    int32 *x;
    int32 *y;
    GrabTask start;
    GrabTask end;
} AniGrab;

/***************  End of grabber stuff  in header **************/

typedef struct AniCom {
    struct {
        int32 n;
        int32 *x;
        int32 *y;
        int32 *col;
        int32 i;
    } c;
    int32 type;
    int32 flag;
    int32 *col;
    int32 *x1;
    int32 *y1;
    int32 *x2;
    int32 *y2;
    int32 *who;
    double zcol;
    double zx1;
    double zy1;
    double zx2;
    double zy2;
    double zrad;
    double zval;
    int32 zthick;
    int32 tfont;
    int32 tsize;
    int32 tcolor;
} AniCom;

extern int32 animation_on_the_fly;

void ani_new_vcr(void);
void ani_do_events(XEvent event);
void ani_expose(Window window);
void ani_on_the_fly(int32 task);
int32 ani_get_ppm_bits(Window window, int32 *wid, int32 *hgt, uchar *out);
int32 ani_write_frame(char *filename, Window window, int32 wid, int32 hgt);
void ani_zero(void);
void ani_get_file(char *fname);
void ani_de_space(char *s);
int32 ani_check_pause(XEvent event);

#endif

#ifndef ARRAY_PLOT_H
#define ARRAY_PLOT_H

#include <X11/Xlib.h>
#include <stdio.h>

extern int32 array_plot_range;

void array_plot_close_files(void);
void array_plot_draw_one(char *);
void array_plot_optimize(int32 *plist);
void array_plot_make_my(char *name);
void array_plot_expose(Window window);
void array_plot_do_events(XEvent event);
void array_plot_init_my(void);
void array_plot_edit(void);
void array_plot_dump(FILE *fp, int32 f);

#endif

#ifndef ARRAY_PRINT_H
#define ARRAY_PRINT_H

int32 array_print(char *filename, char *xtitle, char *ytitle, char *bottom,
                  int32 nacross, int32 ndown, int32 col0, int32 row0,
                  int32 nskip, int32 ncskip, int32 maxrow, int32 maxcol,
                  double **data, double zmin, double zmax, double tlo,
                  double thi, int32 type);

#endif

#ifndef AUTO_X11_H
#define AUTO_X11_H

#include <X11/Xlib.h>

extern int32 auto_redraw_flag;

extern int32 mark_flag;
extern int32 mark_ibrs;
extern int32 mark_ibre;
extern int32 mark_ipts;
extern int32 mark_ipte;
extern int32 mark_ixs;
extern int32 mark_ixe;
extern int32 mark_iys;
extern int32 mark_iye;

void auto_x11_line(int32 a, int32 b, int32 c, int32 d);
void auto_x11_line_trans(double a, double b, double c, double d);
void auto_x11_text(int32 a, int32 b, char *c);
void auto_x11_clr_stab(void);
void auto_x11_stab_line(int32 x, int32 y, int32 xp, int32 yp);
void auto_x11_clear_plot(void);
void auto_x11_redraw_menus(void);
void auto_x11_traverse_diagram(void);
void auto_x11_clear_info(void);
void auto_x11_draw_info(char *bob, int32 x, int32 y);
void auto_x11_refresh_display(void);
int32 auto_x11_bye(int32 *iflag);
void auto_x11_circle(int32 x, int32 y, int32 r);
void auto_x11_col(int32 col);
void auto_x11_bw(void);
void auto_x11_scroll(void);
int32 auto_x11_rubber(int32 *i1, int32 *j1, int32 *i2, int32 *j2, int32 flag);
int32 auto_x11_pop_up_list(char *title, char **list, char *key, int32 n,
                           int32 max, int32 def, int32 x, int32 y, char **hints,
                           char *httxt);
void auto_x11_xor_cross(int32 x, int32 y);
void auto_x11_fill_circle(int32 x, int32 y, int32 r);
void auto_x11_line_width(int32 wid);
void auto_x11_motion(XEvent event);
void auto_x11_display(Window window);
void auto_x11_make(char *wname, char *iname);
void auto_x11_resize_window(XEvent event);
void auto_x11_enter(Window window, int32 v);
void auto_x11_button(XEvent event);
void auto_x11_keypress(XEvent event, int32 *used);
void auto_x11_get_info(int32 *n, char *pname);
void auto_x11_set_mark(int32 i);
void auto_x11_do_range(void);

#endif

#ifndef AXES2_H
#define AXES2_H

#include <X11/Xlib.h>

extern int32 axes2_doing;
extern int32 axes2_doing_box;

void axes2_redraw_cube_pt(double theta, double phi);
void axes2_do(void);
void axes2_redraw_cube(double theta, double phi);
void axes2_box(double x_min, double x_max, double y_min, double y_max, char *sx,
               char *sy, int32 flag);

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
 * The band_alloc_mat function allocates a band matrix for use in   *
 * the other matrix routines listed in this file. Matrix storage  *
 * details are given in the documentation for the type BandMat.   *
 * The BandAllocPiv function allocates memory for pivot           *
 * information. The storage allocated by band_alloc_mat and         *
 * BandAllocPiv is deallocated by the routines band_free_mat and    *
 * band_free_piv, respectively. The band_factor and band_back_solve    *
 * routines perform the actual solution of a band linear system.  *
 * Note that the band_back_solve routine has a parameter b of type  *
 * Vector. The current implementation makes use of a machine    *
 * environment specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type Vector. Thus, the          *
 * implementation of band_back_solve may need to change if the      *
 * type Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "band" (except for  *
 * the factor and solve routines which are called gbfa and gbsl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for bandalloc.                               *
 *                                                                *
 ******************************************************************/

#ifndef BAND_H
#define BAND_H

#include "vector.h"

/******************************************************************
 * Type: BandMat                                                  *
 *----------------------------------------------------------------*
 * The type BandMat is the type of a large (possibly distributed) *
 * band matrix. It is defined to be a pointer to a structure      *
 * with the following fields:                                     *
 *                                                                *
 * size is the number of columns (== number of rows)              *
 * mu   is the upper bandwidth, 0 <= mu <= size-1                 *
 * ml   is the lower bandwidth, 0 <= ml <= size-1                 *
 * smu  is the storage upper bandwidth, mu <= smu <= size-1.      *
 *         The band_factor routine writes the LU factors           *
 *         into the storage for A. The upper triangular factor U, *
 *         however, may have an upper bandwidth as big as         *
 *         MIN(size-1,mu+ml) because of partial pivoting. The smu *
 *         field holds the upper bandwidth allocated for A.       *
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
 * data[j] is a pointer to the uppermost element within the band  *
 *            in the jth column. This pointer may be treated as   *
 *            an array indexed from smu-mu (to access the         *
 *            uppermost element within the band in the jth        *
 *            column) to smu+ml (to access the lowest element     *
 *            within the band in the jth column). (Indices from 0 *
 *            to smu-mu-1 give access to extra storage elements   *
 *            required by band_factor.)                            *
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
    int64 mu;
    int64 ml;
    int64 smu;
    double **data;
} *BandMat;

/* BandMat accessor macros */

/******************************************************************
 * Macro : BAND_ELEM                                              *
 * Usage : BAND_ELEM(A,i,j) = a_ij;  OR                           *
 *         a_ij = BAND_ELEM(A,i,j);                               *
 *----------------------------------------------------------------*
 * BAND_ELEM(A,i,j) references the (i,j)th element of the         *
 * N by N band matrix A, where 0 <= i,j <= N-1. The location      *
 * (i,j) should further satisfy j-(A->mu) <= i <= j+(A->ml).      *
 ******************************************************************/
#define BAND_ELEM(A, i, j) ((A->data)[j][i - j + (A->smu)])

/******************************************************************
 * Macro : BAND_COL                                               *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *----------------------------------------------------------------*
 * BAND_COL(A,j) references the diagonal element of the jth       *
 * column of the N by N band matrix A, 0 <= j <= N-1. The type of *
 * the expression BAND_COL(A,j) is double *. The pointer returned   *
 * by the call BAND_COL(A,j) can be treated as an array which is  *
 * indexed from -(A->mu) to (A->ml).                              *
 ******************************************************************/
#define BAND_COL(A, j) (((A->data)[j]) + (A->smu))

/******************************************************************
 * Macro : BAND_COL_ELEM                                          *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *         BAND_COL_ELEM(col_j,i,j) = a_ij;  OR                   *
 *         a_ij = BAND_COL_ELEM(col_j,i,j);                       *
 *----------------------------------------------------------------*
 * This macro references the (i,j)th entry of the band matrix A   *
 * when used in conjunction with BAND_COL as shown above. The     *
 * index (i,j) should satisfy j-(A->mu) <= i <= j+(A->ml).        *
 ******************************************************************/
#define BAND_COL_ELEM(col_j, i, j) (col_j[i - j])

/* Functions that use the BandMat representation for a band matrix */

/******************************************************************
 * Function : band_alloc_mat                                        *
 * Usage    : A = band_alloc_mat(N, mu, ml, smu);                   *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * band_alloc_mat allocates memory for an N by N band matrix with   *
 * upper bandwidth mu, lower bandwidth ml, and storage upper      *
 * bandwidth smu. Pass smu as follows depending on whether A will *
 * be factored by band_factor:                                     *
 *                                                                *
 * (1) Pass smu = mu if A will not be factored.                   *
 * (2) Pass smu = MIN(N-1,mu+ml) if A will be factored.           *
 *                                                                *
 * band_alloc_mat returns the storage allocated (type BandMat) or   *
 * NULL if the request for matrix storage cannot be satisfied.    *
 * See the documentation for the type BandMat for matrix storage  *
 * details.                                                       *
 ******************************************************************/
BandMat band_alloc_mat(int64 N, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : BandAllocPiv                                        *
 * Usage    : p = band_alloc_piv(N);                                *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * BandAllocPiv allocates memory for pivot information to be      *
 * filled in by the band_factor routine during the factorization   *
 * of an N by N band matrix. The underlying type for pivot        *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, BandAllocPiv returns NULL.  *
 ******************************************************************/
int64 *band_alloc_piv(int64 N);

/******************************************************************
 * Function : band_factor                                          *
 * Usage    : ier = band_factor(A, p);                             *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * band_factor performs the LU factorization of the N by N band    *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the "matrix" A and the    *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * band_factor returns 0 if successful. Otherwise it encountered   *
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
 * call A = band_alloc_mat(N,mu,ml,smu), where mu, ml, and smu are  *
 * as defined above. The user does not have to zero the "extra"   *
 * storage allocated for the purpose of factorization. This will  *
 * handled by the band_factor routine.                             *
 ******************************************************************/
int64 band_factor(BandMat A, int64 *p);

/******************************************************************
 * Function : band_back_solve                                       *
 * Usage    : band_back_solve(A, p, b);                             *
 *----------------------------------------------------------------*
 * band_back_solve solves the N-dimensional system A x = b using    *
 * the LU factorization in A and the pivot information in p       *
 * computed in band_factor. The solution x is returned in b. This  *
 * routine cannot fail if the corresponding call to band_factor    *
 * did not fail.                                                  *
 ******************************************************************/
void band_back_solve(BandMat A, int64 *p, Vector b);

/******************************************************************
 * Function : band_zero                                            *
 * Usage    : band_zero(A);                                        *
 *----------------------------------------------------------------*
 * A(i,j) <- 0.0,    j-(A->mu) <= i <= j+(A->ml).                 *
 ******************************************************************/
void band_zero(BandMat A);

/******************************************************************
 * Function : band_copy                                            *
 * Usage    : band_copy(A, B, copymu, copyml);                     *
 *----------------------------------------------------------------*
 * band_copy copies the submatrix with upper and lower bandwidths  *
 * copymu, copyml of the N by N band matrix A into the N by N     *
 * band matrix B.                                                 *
 ******************************************************************/
void band_copy(BandMat A, BandMat B, int64 copymu, int64 copyml);

/******************************************************************
 * Function: band_scale                                            *
 * Usage   : band_scale(c, A);                                     *
 *----------------------------------------------------------------*
 * A(i,j) <- c*A(i,j),   j-(A->mu) <= i <= j+(A->ml).             *
 ******************************************************************/
void band_scale(double c, BandMat A);

/******************************************************************
 * Function : band_add_i                                            *
 * Usage    : band_add_i(A);                                        *
 *----------------------------------------------------------------*
 * A(j,j) <- A(j,j)+1.0,   0 <= j <= (A->size)-1.                 *
 ******************************************************************/
void band_add_i(BandMat A);

/******************************************************************
 * Function : band_free_mat                                         *
 * Usage    : band_free_mat(A);                                     *
 *----------------------------------------------------------------*
 * band_free_mat frees the memory allocated by band_alloc_mat for     *
 * the band matrix A.                                             *
 ******************************************************************/
void band_free_mat(BandMat A);

/******************************************************************
 * Function : band_free_piv                                         *
 * Usage    : band_free_piv(p);                                     *
 *----------------------------------------------------------------*
 * band_free_piv frees the memory allocated by BandAllocPiv for     *
 * the pivot information array p.                                 *
 ******************************************************************/
void band_free_piv(int64 *p);

/******************************************************************
 * Function : band_print                                           *
 * Usage    : band_print(A);                                       *
 *----------------------------------------------------------------*
 * This routine prints the N by N band matrix A (upper and lower  *
 * bandwidths A->mu and A->ml, respectively) to standard output   *
 * as it would normally appear on paper. It is intended as a      *
 * debugging tool with small values of N. The elements are        *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 ******************************************************************/
void band_print(BandMat A);

/* Functions that use the double ** representation for a band matrix */

/******************************************************************
 * Function : bandalloc                                           *
 * Usage    : double **a;                                           *
 *            a = band_alloc2(n, smu, ml);                          *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * band_alloc2(n, smu, ml) allocates storage for an n by n band     *
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
 * (2) Pass smu = MIN(n-1,mu+ml) if A will be factored.           *
 *                                                                *
 * The underlying type of the band matrix returned is double **. If *
 * we allocate a band matrix A in double **a by                     *
 * a = band_alloc2(n,smu,ml), then a[0] is a pointer to             *
 * n*(smu + ml + 1) contiguous storage locations and a[j] is a  *
 * pointer to the uppermost element in the storage for the jth    *
 * column. The expression a[j][i-j+smu] references the (i,j)th    *
 * element of A, where 0 <= i,j <= n-1 and j-mu <= i <= j+ml.     *
 * (The elements a[j][0], a[j][1], ..., a[j][smu-mu-1] are used   *
 * by gbfa and gbsl.)                                             *
 ******************************************************************/
double **band_alloc2(int64 n, int64 smu, int64 ml);

/******************************************************************
 * Function : bandallocpiv                                        *
 * Usage    : int64 *pivot;                                     *
 *            pivot = band_alloc_piv2(n);                            *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * band_alloc_piv2(n) allocates an array of n integers. It returns a *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 ******************************************************************/
int64 *band_alloc_piv2(int64 n);

/******************************************************************
 * Function : gbfa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = band_gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * band_gbfa(a,n,mu,ml,smu,p) factors the n by n band matrix A (upper  *
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
 * as a = band_alloc2(n,smu,ml), where smu = MIN(n-1,mu+ml). The    *
 * call to gbfa is ier = band_gbfa(a,n,mu,ml,smu,p). The corresponding *
 * call to gbsl is band_gbsl(a,n,smu,ml,p,b). The user does not need   *
 * to zero the "extra" storage allocated for the purpose of       *
 * factorization. This is handled by the gbfa routine. If A is    *
 * not going to be factored and backsolved, then it can be        *
 * allocated as a = band_alloc2(n,smu,ml). In either case, all      *
 * routines in this section use the parameter name smu for a      *
 * parameter which must be the "storage upper bandwidth" which    *
 * was passed to bandalloc.                                       *
 ******************************************************************/
int64 band_gbfa(double **a, int64 n, int64 mu, int64 ml, int64 smu, int64 *p);

/******************************************************************
 * Function : gbsl                                                *
 * Usage    : double *b;                                            *
 *            ier = band_gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier == 0) band_gbsl(a,n,smu,ml,p,b);                 *
 *----------------------------------------------------------------*
 * band_gbsl(a,n,smu,ml,p,b) solves the n by n linear system           *
 * Ax = b, where A is band matrix stored in "a" with storage      *
 * upper bandwidth smu and lower bandwidth ml. It assumes that A  *
 * has been LU factored and the pivot array p has been set by a   *
 * successful call band_gbfa(a,n,mu,ml,smu,p). The solution x is       *
 * written into the b array.                                      *
 ******************************************************************/
void band_gbsl(double **a, int64 n, int64 smu, int64 ml, int64 *p, double *b);

/******************************************************************
 * Function : band_zero2                                            *
 * Usage    : band_zero2(a,n,mu,ml,smu);                            *
 *----------------------------------------------------------------*
 * a(i,j) <- 0.0,   0 <= i,j <= n-1, j-mu <= i <= j+ml.           *
 ******************************************************************/
void band_zero2(double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : band_copy2                                            *
 * Usage    : band_copy2(a,b,n,a_smu,b_smu,copymu,copyml);          *
 *----------------------------------------------------------------*
 * b(i,j) <- a(i,j), 0 <= i,j <= n-1, j-copymu <= i <= j+copyml.  *
 ******************************************************************/
void band_copy2(double **a, double **b, int64 n, int64 a_smu, int64 b_smu,
                int64 copymu, int64 copyml);

/******************************************************************
 * Function : band_scale2                                           *
 * Usage    : band_scale2(c,a,n,mu,ml);                             *
 *----------------------------------------------------------------*
 * a(i,j) <- c*a(i,j),   0 <= i,j <= n-1, j-mu <= i <= j+ml.      *
 ******************************************************************/
void band_scale2(double c, double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : band_add_i2                                            *
 * Usage    : band_add_i2(a,n,smu);                                  *
 *----------------------------------------------------------------*
 * a(j,j) <- a(j,j)+1.0,   0 <= j <= n-1.                         *
 ******************************************************************/
void band_add_i2(double **a, int64 n, int64 smu);

/******************************************************************
 * Function : band_free_pvi2                                         *
 * Usage    : band_free_pvi2(p);                                     *
 *----------------------------------------------------------------*
 * band_free_pvi2(p) frees the pivot array p allocated by            *
 * bandallocpiv.                                                  *
 ******************************************************************/
void band_free_pvi2(int64 *p);

/******************************************************************
 * Function : band_free2                                            *
 * Usage    : band_free2(a);                                        *
 *----------------------------------------------------------------*
 * band_free2(a) frees the band matrix a allocated by bandalloc.    *
 ******************************************************************/
void band_free2(double **a);

/******************************************************************
 * Function : band_print2                                           *
 * Usage    : band_print2(a,n,mu,ml,smu);                           *
 *----------------------------------------------------------------*
 * band_print2(a,n,mu,ml,smu) prints the n by n band matrix stored  *
 * in a (with upper bandwidth mu and lower bandwidth ml) to       *
 * standard output as it would normally appear on paper. It is    *
 * intended as a debugging tool with small values of n. The       *
 * elements are printed using the %g option. A blank line is      *
 * printed before and after the matrix.                           *
 ******************************************************************/
void band_print2(double **a, int64 n, int64 mu, int64 ml, int64 smu);

#endif

#ifndef BROWSE_H
#define BROWSE_H

#define BMAXCOL 20

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct Browser {
    Window base;
    Window upper;
    Window find, up, down, pgup, pgdn, home, end, left, right;
    Window first, last, restore, write, get, close;
    Window load, repl, unrepl, table, addcol, delcol;
    Window main;
    Window label[BMAXCOL];
    Window time;
    Window hint;
    char hinttxt[256];
    int32 dataflag;
    int32 xflag;
    int32 col0, row0, ncol, nrow;
    int32 maxrow;
    int32 maxcol;
    double **data;
    int32 istart;
    int32 iend;
} Browser;

extern Browser my_browser;

double **get_browser_data(void);
void set_browser_data(double **data, int32 col0);
double *browse_get_data_col(int32 c);
int32 browse_get_time_now(void);
void browse_wait_a_sec(int32 msec);
int32 get_max_row_browser(void);
void write_my_browser_data(FILE *fp);
void browse_wipe_rep(void);
void browse_find_variable(char *s, int32 *col);
void refresh_browser(int32 length);
void reset_browser(void);
void init_browser(void);
void make_new_browser(void);
Window browse_button2(Window root, int32 row, int32 col, int32 iflag);
Window browse_button_data(Window root, int32 row, int32 col, char *name,
                          int32 iflag);
void expose_my_browser(XEvent event);
void enter_my_browser(XEvent event, int32 yn);
void my_browse_button(XEvent event);
void my_browse_keypress(XEvent event, int32 *used);
void resize_my_browser(Window win);
void browse_get_data_xyz(double *x, double *y, double *z, int32 i1, int32 i2,
                         int32 i3, int32 off);
void browse_open_write_file(FILE **fp, char *fil, int32 *ok);
void write_my_browser_data(FILE *fp);
void data_get_my_browser(int32 row);

#endif

#ifndef CALC_H
#define CALC_H

#include <X11/Xlib.h>

void calc_q_calc(void);
int32 calc_do_calc(char *temp, double *z);
double calc(char *expr, int32 *ok);

#endif

#ifndef CHOICE_BOX_H
#define CHOICE_BOX_H

#include <X11/Xlib.h>
#include "struct.h"

void choice_box_base(char *wname, int32 n, int32 mcc, char **names,
                     int32 *check, int32 type);

#endif

#ifndef COLOR_H
#define COLOR_H

#include <X11/Xlib.h>

extern int32 custom_color;

extern int32 color_min;
extern int32 color_total;
extern int32 color_max;
extern int32 COLOR;

void color_set_s(int32 col);
void color_set(int32 col);
void color_new_map(int32 type);
void color_get_ps(int32 i, double *r, double *g, double *b);
void color_get_svg(int32 i, int32 *r, int32 *g, int32 *b);
void color_map_make(void);
uint32 color_map(int32 i);

#endif

#ifndef GGETS_H
#define GGETS_H

#include <X11/Xlib.h>

#define MAX_INCLUDE_FILES 10
#define ClickTime 200

extern int32 MSStyle;
extern char *info_message;
extern int32 CURS_X;
extern int32 CURS_Y;
extern int32 PS_Color;

void ggets_ping(void);
void ggets_reset_graphics(void);
void ggets_blank_screen(Window window);
void ggets_set_fore(void);
void ggets_set_back(void);
void ggets_show_char(int32 ch, int32 col, int32 row, Window or);
void ggets_chk_xor(void);
void ggets_clr_command(void);
void ggets_draw_info_pop(Window win);
void ggets_bottom_msg(char *msg);
void ggets_err_msg(char *string);
int32 ggets_plintf(char *fmt, ...);
int32 ggets_show_position(XEvent event);
void ggets_put_command(char *string);
int32 ggets_get_key_press(XEvent *event);
void ggets_cput_text(void);
int32 ggets_mouse_xy(int32 *x, int32 *y, Window window);
void ggets_f_text(int32 x, int32 y, char *string, Window o);
void ggets_bar(int32 x, int32 y, int32 x2, int32 y2, Window window);
void ggets_rectangle(int32 x, int32 y, int32 x2, int32 y2, Window window);
void ggets_circle(int32 x, int32 y, int32 radius, Window window);
void ggets_xline(int32 x0, int32 y0, int32 x1, int32 y1, Window window);
int32 ggets_new_float(char *name, double *value);
int32 ggets_new_int(char *name, int32 *value);
void ggets_display_command(char *name, char *value, int32 pos);
void ggets_put_cursor_at(Window window, int32 col0, int32 pos);
void ggets_mov_mem(char *s1, char *s2, int32 len);
void ggets_mem_mov(char *s1, char *s2, int32 len);
void ggets_edit_window(Window window, int32 *pos, char *value, int32 *col,
                       int32 *done, int32 ch);
void ggets_edit_command_string(XEvent event, char *name, char *value,
                               int32 *done, int32 *pos, int32 *col);
int32 ggets_new_string(char *name, char *value);

#endif

#ifndef LOAD_EQN_H
#define LOAD_EQN_H

#include <stdio.h>

typedef struct InternSet {
    char *name;
    char *does;
    uint32 use;
} InternSet;

extern int32 RunImmediately;
extern int32 MultiWin;
extern int32 START_LINE_TYPE;
extern int32 Nintern_set;
extern double TOR_PERIOD;
extern int32 TORUS;
extern int32 IX_PLT[10];
extern int32 IY_PLT[10];
extern int32 IZ_PLT[10];
extern int32 NPltV;
extern double X_LO[10];
extern double Y_LO[10];
extern double X_HI[10];
extern double Y_HI[10];
extern InternSet intern_set[MAX_INTERN_SET];
extern double last_ic[MAX_ODE];
extern char delay_string[MAX_ODE][80];
extern int32 itor[MAX_ODE];

extern int32 POIEXT;
extern int32 HIST;

extern int32 ENDSING;
extern int32 SHOOT;
extern int32 PAR_FOL;
extern int32 xorfix;
extern int32 silent;
extern int32 got_file;

/* The acutual max filename length is determined by the
 * FILENAME_MAX (see <stdio.h>), and usually 4096 -- but
 * this is huge and usually overkill.  On the otherhand
 * the old Xpp default string buffer size of 100 is a bit
 * restricitive for lengths of filenames. You could also
 * set this define in the Makefile or at compile time to
 * override the below definition. */

#ifndef XPP_MAX_NAME
#define XPP_MAX_NAME 512

extern char this_file[XPP_MAX_NAME];
extern char this_internset[XPP_MAX_NAME];

extern int32 mov_ind;
extern int32 storind;
extern int32 STORFLAG;
extern int32 INFLAG;
extern int32 MAXSTOR;

extern double x_3d[2];
extern double y_3d[2];
extern double z_3d[2];

extern int32 IXPLT;
extern int32 IYPLT;
extern int32 IZPLT;

extern int32 AXES;
extern int32 TIMPLOT;

extern int32 TIMPLOT;
extern int32 PLOT_3D;
extern double MY_XLO;
extern double MY_YLO;
extern double MY_XHI;
extern double MY_YHI;

extern char options[100];
extern double ATOLER;
extern double HMIN;
extern double HMAX;
extern double BVP_EPS;
extern double BVP_TOL;
extern int32 euler_max_iter;
extern double euler_tol;
extern int32 BVP_MAXIT;
extern int32 BVP_FLAG;
extern int32 FFT;
extern int32 NULL_HERE;

#endif

#ifndef COMLINE_H
#define COMLINE_H

extern int32 NincludedFiles;
extern int32 Nintern_2_use;
extern int32 loadincludefile;
extern int32 querysets;
extern int32 querypars;
extern int32 queryics;
extern int32 dryrun;
extern int32 noicon;
extern int32 newseed;

extern char includefilename[MAX_INCLUDE_FILES][XPP_MAX_NAME];

void comline_do(int32 argc, char **argv);
int32 comline_if_needed_select_sets(void);
int32 comline_if_needed_load_set(void);
int32 comline_if_needed_load_par(void);
int32 comline_if_needed_load_ic(void);
int32 comline_if_needed_load_ext_options(void);

#endif

#ifndef DAE_FUN_H
#define DAE_FUN_H

int32 dae_fun_add_svar(char *name, char *rhs);
int32 dae_fun_add_svar_names(void);
int32 dae_fun_add_aeqn(char *rhs);
int32 dae_fun_compile_svars(void);
void dae_fun_reset_dae(void);
void dae_fun_set_init_guess(void);
void dae_fun_err_dae(void);
void dae_fun_do_daes(void);
void dae_fun_get_new_guesses(void);

#endif

#ifndef DELAY_HANDLE_H
#define DELAY_HANDLE_H

extern double AlphaMax;
extern double OmegaMax;
extern int32 DelayFlag;
extern int32 DelayGrid;
extern int32 NDelay;
extern int32 del_stab_flag;
extern int32 WhichDelay;

extern double variable_shift[2][MAX_ODE];
extern double delay_list[MAX_DELAY];

double delay_handle_stab_eval(double delay, int32 var);
int32 delay_handle_alloc_delay(double big);
void delay_handle_free_delay(void);
void delay_handle_stor_delay(double *y);
double delay_handle_get_delay(int32 in, double tau);
int32 delay_handle_do_init_delay(double big);

#endif

#ifndef DEL_STAB_H
#define DEL_STAB_H

typedef struct COMPLEX {
    double r;
    double i;
} COMPLEX;

void del_stab_do_delay_sing(double *x, double eps, double err, double big,
                            int32 maxit, int32 n, int32 *ierr,
                            double *stabinfo);
int32 del_stab_find_positive_root(double *coef, double *delay, int32 n, int32 m,
                                  double err, double eps, double big,
                                  int32 maxit, double *rr);

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
 * The dense_alloc_mat function allocates a dense matrix for use in *
 * the other DenseMat routines listed in this file. Matrix        *
 * storage details are given in the documentation for the type    *
 * DenseMat. The DenseAllocPiv function allocates memory for      *
 * pivot information. The storage allocated by dense_alloc_mat and  *
 * DenseAllocPiv is deallocated by the routines dense_free_mat and  *
 * dense_free_piv, respectively. The dense_factor and dense_back_solve *
 * routines perform the actual solution of a dense linear system. *
 * Note that the dense_back_solve routine has a parameter b of type *
 * Vector. The current implementation makes use of a machine    *
 * environment-specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type Vector. Thus, the          *
 * implementation of dense_back_solve may need to change if the     *
 * type Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "den" (except for   *
 * the factor and solve routines which are called gefa and gesl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for denalloc.                                *
 *                                                                *
 ******************************************************************/

#ifndef DENSE_H
#define DENSE_H

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
 * Function : dense_alloc_mat                                       *
 * Usage    : A = dense_alloc_mat(N);                               *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * dense_alloc_mat allocates memory for an N by N dense matrix and  *
 * returns the storage allocated (type DenseMat). dense_alloc_mat   *
 * returns NULL if the request for matrix storage cannot be       *
 * satisfied. See the above documentation for the type DenseMat   *
 * for matrix storage details.                                    *
 *                                                                *
 ******************************************************************/
DenseMat dense_alloc_mat(int64 N);

/******************************************************************
 *                                                                *
 * Function : DenseAllocPiv                                       *
 * Usage    : p = dense_alloc_piv(N);                               *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocPiv allocates memory for pivot information to be     *
 * filled in by the dense_factor routine during the factorization  *
 * of an N by N dense matrix. The underlying type for pivot       *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, DenseAllocPiv returns NULL. *
 *                                                                *
 ******************************************************************/
int64 *dense_alloc_piv(int64 N);

/******************************************************************
 *                                                                *
 * Function : dense_factor                                         *
 * Usage    : ier = dense_factor(A, p);                            *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * dense_factor performs the LU factorization of the N by N dense  *
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
 * dense_factor returns 0 if successful. Otherwise it encountered  *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 ******************************************************************/
int64 dense_factor(DenseMat A, int64 *p);

/******************************************************************
 *                                                                *
 * Function : dense_back_solve                                      *
 * Usage    : dense_back_solve(A, p, b);                            *
 *----------------------------------------------------------------*
 * dense_back_solve solves the N-dimensional system A x = b using   *
 * the LU factorization in A and the pivot information in p       *
 * computed in dense_factor. The solution x is returned in b. This *
 * routine cannot fail if the corresponding call to dense_factor   *
 * did not fail.                                                  *
 *                                                                *
 ******************************************************************/
void dense_back_solve(DenseMat A, int64 *p, Vector b);

/******************************************************************
 *                                                                *
 * Function : dense_zero                                           *
 * Usage    : dense_zero(A);                                       *
 *----------------------------------------------------------------*
 * dense_zero sets all the elements of the N by N matrix A to 0.0. *
 *                                                                *
 ******************************************************************/
void dense_zero(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_copy                                           *
 * Usage    : dense_copy(A, B);                                    *
 *----------------------------------------------------------------*
 * dense_copy copies the contents of the N by N matrix A into the  *
 * N by N matrix B.                                               *
 *                                                                *
 ******************************************************************/
void dense_copy(DenseMat A, DenseMat B);

/******************************************************************
 *                                                                *
 * Function: dense_scal                                           *
 * Usage   : dense_scal(c, A);                                    *
 *----------------------------------------------------------------*
 * dense_scal scales the elements of the N by N matrix A by the   *
 * constant c and stores the result back in A.                    *
 *                                                                *
 ******************************************************************/
void dense_scal(double c, DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_add_i                                           *
 * Usage    : dense_add_i(A);                                       *
 *----------------------------------------------------------------*
 * dense_add_i adds the identity matrix to A and stores the result  *
 * back in A.                                                     *
 *                                                                *
 ******************************************************************/
void dense_add_i(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_free_mat                                        *
 * Usage    : dense_free_mat(A);                                    *
 *----------------------------------------------------------------*
 * dense_free_mat frees the memory allocated by dense_alloc_mat for   *
 * the N by N matrix A.                                           *
 *                                                                *
 ******************************************************************/
void dense_free_mat(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_free_piv                                        *
 * Usage    : dense_free_piv(p);                                    *
 *----------------------------------------------------------------*
 * dense_free_piv frees the memory allocated by DenseAllocPiv for   *
 * the pivot information array p.                                 *
 *                                                                *
 ******************************************************************/
void dense_free_piv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : dense_print                                          *
 * Usage    : dense_print(A);                                      *
 *----------------------------------------------------------------*
 * This routine prints the N by N dense matrix A to standard      *
 * output as it would normally appear on paper. It is intended    *
 * as a debugging tool with small values of N. The elements are   *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/
void dense_print(DenseMat A);

/* Functions that use the double ** representation for a dense matrix */

/******************************************************************
 *                                                                *
 * Function : denalloc                                            *
 * Usage    : double **a;                                           *
 *            a = dense_alloc2(n);                                    *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * dense_alloc2(n) allocates storage for an n by n dense matrix. It   *
 * returns a pointer to the newly allocated storage if            *
 * successful. If the memory request cannot be satisfied, then    *
 * denalloc returns NULL. The underlying type of the dense matrix *
 * returned is double **. If we allocate a dense matrix double **a by *
 * a = dense_alloc2(n), then a[j][i] references the (i,j)th element   *
 * of the matrix a, 0 <= i,j <= n-1, and a[j] is a pointer to the *
 * first element in the jth column of a. The location a[0]        *
 * contains a pointer to n^2 contiguous locations which contain   *
 * the elements of a.                                             *
 *                                                                *
 ******************************************************************/
double **dense_alloc2(int64 n);

/******************************************************************
 *                                                                *
 * Function : denallocpiv                                         *
 * Usage    : int64 *pivot;                                     *
 *            pivot = dense_alloc_piv2(n);                             *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * dense_alloc_piv2(n) allocates an array of n integers. It returns a  *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 *                                                                *
 ******************************************************************/
int64 *dense_alloc_piv2(int64 n);

/******************************************************************
 *                                                                *
 * Function : gefa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = dense_gefa(a,n,p);                                  *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * dense_gefa(a,n,p) factors the n by n dense matrix a. It overwrites   *
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
int64 dense_gefa(double **a, int64 n, int64 *p);

/******************************************************************
 *                                                                *
 * Function : gesl                                                *
 * Usage    : double *b;                                            *
 *            ier = dense_gefa(a,n,p);                                  *
 *            if (ier == 0) dense_gesl(a,n,p,b);                        *
 *----------------------------------------------------------------*
 * dense_gesl(a,n,p,b) solves the n by n linear system ax = b. It       *
 * assumes that a has been LU factored and the pivot array p has  *
 * been set by a successful call to dense_gefa(a,n,p). The solution x   *
 * is written into the b array.                                   *
 *                                                                *
 ******************************************************************/
void dense_gesl(double **a, int64 n, int64 *p, double *b);

/******************************************************************
 *                                                                *
 * Function : den_zero                                             *
 * Usage    : dense_zero2(a,n);                                       *
 *----------------------------------------------------------------*
 * dense_zero2(a,n) sets all the elements of the n by n dense matrix  *
 * a to be 0.0.                                                   *
 *                                                                *
 ******************************************************************/
void dense_zero2(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_copy                                             *
 * Usage    : dense_copy2(a,b,n);                                     *
 *----------------------------------------------------------------*
 * dense_copy2(a,b,n) copies the n by n dense matrix a into the       *
 * n by n dense matrix b.                                         *
 *                                                                *
 ******************************************************************/
void dense_copy2(double **a, double **b, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_scale                                            *
 * Usage    : dense_scale2(c,a,n);                                    *
 *----------------------------------------------------------------*
 * dense_scale2(c,a,n) scales every element in the n by n dense       *
 * matrix a by c.                                                 *
 *                                                                *
 ******************************************************************/
void dense_scale2(double c, double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_add_i                                             *
 * Usage    : dense_add_i2(a,n);                                       *
 *----------------------------------------------------------------*
 * dense_add_i2(a,n) increments the n by n dense matrix a by the       *
 * identity matrix.                                               *
 *                                                                *
 ******************************************************************/
void dense_add_i2(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_free_piv                                          *
 * Usage    : dense_free_piv2(p);                                      *
 *----------------------------------------------------------------*
 * dense_free_piv2(p) frees the pivot array p allocated by             *
 * denallocpiv.                                                   *
 *                                                                *
 ******************************************************************/
void dense_free_piv2(int64 *p);

/******************************************************************
 *                                                                *
 * Function : den_free                                             *
 * Usage    : dense_free2(a);                                         *
 *----------------------------------------------------------------*
 * dense_free2(a) frees the dense matrix a allocated by denalloc.     *
 *                                                                *
 ******************************************************************/
void dense_free2(double **a);

/******************************************************************
 *                                                                *
 * Function : den_print                                            *
 * Usage    : dense_print2(a,n);                                      *
 *----------------------------------------------------------------*
 * dense_print2(a,n) prints the n by n dense matrix a to standard     *
 * output as it would normally appear on paper. It is intended as *
 * a debugging tool with small values of n. The elements are      *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/
void dense_print2(double **a, int64 n);

#endif

#ifndef DERIVED_H
#define DERIVED_H

int32 derived_compile(void);
void derived_evaluate(void);
int32 derived_add(char *name, char *rhs);

#endif

#ifndef DIAGRAM_H
#define DIAGRAM_H
#include "auto_nox.h"

extern Diagram *bifd;

#include <stdio.h>

extern int32 NBifs;

void diagram_start(int32 n);
void diagram_edit_start(int32 ibr, int32 ntot, int32 itp, int32 lab,
                        int32 nfpar, double a, double *uhi, double *ulo,
                        double *u0, double *ubar, double *par, double per,
                        int32 n, int32 icp1, int32 icp2, int32 icp3, int32 icp4,
                        double *evr, double *evi);
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
void diagram_redraw(void);
void diagram_write_info_out(void);
void diagram_write_init_data_file(void);
void diagram_write_pts(void);
void diagram_post_auto(void);
void diagram_svg_auto(void);
void diagram_bound(double *xlo, double *xhi, double *ylo, double *yhi);
int32 diagram_save(FILE *fp, int32 n);
int32 diagram_load(FILE *fp, int32 node);
void diagram_load_browser_with_branch(int32 ibr, int32 pts, int32 pte);

#endif

#ifndef DIALOG_BOX_H
#define DIALOG_BOX_H

#include <X11/Xlib.h>
#include "struct.h"

int32 dialog_box_get(char *wname, char *name, char *value, char *ok,
                     char *cancel, int32 max);

#endif

#ifndef DO_FIT_H
#define DO_FIT_H

void do_fit_init_info(void);
void do_fit_get_info(double *y, double *a, double *t0, int32 *flag, double eps,
                     double *yfit, double **yderv, int32 npts, int32 npars,
                     int32 nvars, int32 *ivar, int32 *ipar);
void do_fit_printem(double **yderv, double *yfit, double *t0, int32 npars,
                    int32 nvars, int32 npts);
int32 do_fit_one_step_int(double *y, double t0, double t1, int32 *istart);
void do_fit_test(void);
int32 do_fit_run(char *filename, int32 npts, int32 npars, int32 nvars,
                 int32 maxiter, int32 ndim, double eps, double tol, int32 *ipar,
                 int32 *ivar, int32 *icols, double *y0, double *a,
                 double *yfit);
int32 do_fit_marlev_step(double *t0, double *y0, double *y, double *sig,
                         double *a, int32 npts, int32 nvars, int32 npars,
                         int32 *ivar, int32 *ipar, double *covar, double *alpha,
                         double *chisq, double *alambda, double *work,
                         double **yderv, double *yfit, double *ochisq,
                         int32 ictrl, double eps);
int32 do_fit_mrqcof(double *t0, double *y0, double *y, double *sig, double *a,
                    int32 npts, int32 nvars, int32 npars, int32 *ivar,
                    int32 *ipar, double *alpha, double *chisq, double *beta,
                    double **yderv, double *yfit, double eps);

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

int32
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

void dormpri_dp_err(int32 k);
int32 dp(int32 *istart, double *y, double *t, int32 n, double tout, double *tol,
         double *atol, int32 flag, int32 *kflag);
int32 dormprin(int32 *istart, double *y, double *t, int32 n, double tout,
               double *tol, double *atol, int32 flag, int32 *kflag);
int32 dopri5(uint32 n, FcnEqDiff fcn, double x, double *y, double xend,
             double *rtoler, double *atoler, int32 itoler, SolTrait solout,
             int32 iout, FILE *fileout, double uround, double safe, double fac1,
             double fac2, double beta, double hmax, double h, long nmax,
             int32 meth, long nstiff, uint32 nrdens, uint32 *icont,
             uint32 licont, double *work);

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

void edit_rhs_menu(void);
void edit_rhs(void);
void edit_rhs_user_fun_info(FILE *fp);
void edit_rhs_functions(void);
int32 edit_rhs_save_as(void);

#endif

#ifndef EIG_LIST_H
#define EIG_LIST_H

#include <X11/Xlib.h>

void eig_list_draw_eq_list(Window window);
void eig_list_create_eq_list(void);
void eig_list_eq_list_keypress(XEvent event, int32 *used);
void eig_list_enter_eq_stuff(Window window, int32 b);
void eig_list_eq_list_button(XEvent event);
void eig_list_get_new_size(Window win, uint32 *wid, uint32 *hgt);
void eig_list_resize_eq_list(Window win);
void eig_list_create_eq_box(int32 cp, int32 cm, int32 rp, int32 rm, int32 im,
                            double *y, int32 n);
void eig_list_draw_eq_box(Window window);

#endif

#ifndef EXTRA_H
#define EXTRA_H

#define MAXW 50

extern int32 dll_flag;
extern int32 dll_loaded;
extern char dll_lib[256];
extern char dll_fun[256];

void extra_load_new_dll(void);
int32 extra_my_fun(double *in, double *out, int32 nin, int32 nout, double *v,
                   double *c);
void extra_auto_load_dll(void);
void extra_do_in_out(void);
void extra_add_export_list(char *in, char *out);
void extra_do_export_list(void);
void extra_get_import_values(int32 n, double *ydot, char *soname, char *sofun,
                             int32 ivar, double *wgt[MAXW], double *var,
                             double *con);

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

extern int32 NFlags;
extern double STOL;

int32 flags_add_global(char *cond, int32 sign, char *rest);
void flags_show(void);
int32 flags_compile(void);
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
#define MAXCOMMENTS 500

typedef struct FixInfo {
    char *name;
    char *value;
} FixInfo;

typedef struct Action {
    char *text;
    char *action;
    int32 aflag;
} Action;

typedef struct {
    int32 *com;
    char *string;
    char *name;
    int32 side;
} BcStruct;

extern char *ode_names[MAX_ODE];
extern char *save_eqn[MAXLINES];
extern double default_val[MAX_PAR];
extern double default_ic[MAX_ODE];
extern int32 PrimeStart;
extern int32 NCON_START;
extern int32 NSYM_START;
extern int32 BVP_N;
extern FILE *convertf;
extern int32 FIX_VAR;
extern int32 EqType[MAX_ODE];

extern FixInfo fixinfo[MAX_ODE];

extern int32 *my_ode[MAX_ODE];
extern int32 *plotlist;
extern int32 N_plist;
extern Action comments[MAXCOMMENTS];
extern BcStruct my_bc[MAX_ODE];

extern int32 n_comments;
extern int32 ConvertStyle;

extern int32 NODE;
extern int32 NUPAR;
extern int32 NLINES;
extern int32 IN_VARS;
extern int32 leng[MAX_ODE];

int32 form_ode_make_eqn(void);
void form_ode_strip_saveqn(void);
int32 form_ode_idsc(char *string);
int32 form_ode_get_eqn(FILE *fptr);
int32 form_ode_compiler(char *bob, FILE *fptr);
char *form_ode_get_first(char *string, char *src);
char *form_ode_do_fit_get_next(char *src);
void form_ode_create_plot_list(void);
void form_ode_add_varinfo(int32 type, char *lhs, char *rhs, int32 nargs,
                          char args[20][13 + 1]);
int32 form_ode_extract_args(char *s1, int32 i0, int32 *ie, int32 *narg,
                            char args[20][13 + 1]);
int32 form_ode_find_char(char *s1, char *s2, int32 i0, int32 *i1);
int32 form_ode_search_array(char *old, char *new, int32 *i1, int32 *i2,
                            int32 *flag);
void form_ode_subsk(char *big, char *new, int32 k, int32 flag);

/* for parsing par, init with whitespace correctly */
char *form_ode_new_string2(char *old, int32 length);
char *form_ode_get_next2(char **tokens_ptr);
#endif

#ifndef GEAR_H
#define GEAR_H

extern int32 UnstableManifoldColor;
extern int32 StableManifoldColor;
extern double ShootIC[8][MAX_ODE];
extern int32 ShootICFlag;
extern int32 ShootICFlag;
extern int32 ShootIndex;

void gear_do_sing(double *x, double eps, double err, double big, int32 maxit,
                  int32 n, int32 *ierr, double *stabinfo);
void gear_do_sing_info(double *x, double eps, double err, double big,
                       int32 maxit, int32 n, double *er, double *em,
                       int32 *ierr);

void gear_shoot_this_now(void);
void gear_get_complex_evec(double *m, double evr, double evm, double *br,
                           double *bm, int32 n, int32 maxit, double err,
                           int32 *ierr);
void gear_get_evec(double *a, double *anew, double *b, double *bp, int32 n,
                   int32 maxit, double err, int32 *ipivot, double eval,
                   int32 *ierr);
void gear_eigen(int32 n, double *a, double *ev, double *work, int32 *ierr);
double gear_sign(double x, double y);
int32 gear_imin(int32 x, int32 y);
double gear_amax(double u, double v);
void gear_jac_trans(double *x, double *y, double *yp, double *xp, double eps,
                    double *d, int32 n);
void gear_get_jac(double *x, double *y, double *yp, double *xp, double eps,
                  double *dermat, int32 n);
void gear_rooter(double *x, double err, double eps, double big, double *work,
                 int32 *ierr, int32 maxit, int32 n);
int32 gear(int32 n, double *t, double tout, double *y, double hmin, double hmax,
           double eps, int32 mf, double *error, int32 *kflag, int32 *jstart,
           double *work, int32 *iwork);
int32 ggear(int32 n, double *t, double tout, double *y, double hmin,
            double hmax, double eps, int32 mf, double *error, int32 *kflag,
            int32 *jstart, double *work, int32 *iwork);
double gear_min(double x, double y);
void gear_sgefa(double *a, int32 lda, int32 n, int32 *ipvt, int32 *info);
void gear_sgesl(double *a, int32 lda, int32 n, int32 *ipvt, double *b,
                int32 job);
void gear_save_batch_shoot(void);

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
#define LMAX(a, b) ((a > b) ? a : b)

#include <X11/Xlib.h>
#include <stdio.h>

extern int32 AutoFreezeFlag;
extern int32 colorline[];

void graf_par_change_view_com(int32 com);
void graf_par_ind_to_sym(int32 ind, char *str);
void graf_par_get_max(int32 index, double *vmin, double *vmax);
void graf_par_check_windows(void);
void graf_par_xi_vs_t(void);
void graf_par_redraw_the_graph(void);
void graf_par_get_3d_com(void);
void graf_par_window_zoom_com(int32 c);
void graf_par_graph_all(int32 *list, int32 n, int32 type);
void graf_par_change_cmap_com(int32 i);
void graf_par_freeze_com(int32 c);
void graf_par_key_frz_com(int32 c);
void graf_par_auto_freeze_it(void);
void graf_par_draw_freeze(Window window);
void graf_par_init_bd(void);
void graf_par_add_a_curve_com(int32 c);
void graf_par_default_window(void);
void graf_par_dump_ps(int32 i);

#endif

#ifndef GRAPHICS_H
#define GRAPHICS_H

extern double THETA0;
extern double PHI0;
extern int32 PS_Port;
extern int32 PointRadius;
extern int32 PointType;

extern XFontStruct *symfonts[5];
extern XFontStruct *romfonts[5];
extern int32 avsymfonts[5];
extern int32 avromfonts[5];

extern int32 DLeft;
extern int32 DRight;
extern int32 DTop;
extern int32 DBottom;
extern int32 VTic;
extern int32 HTic;
extern int32 VChar;
extern int32 HChar;
extern int32 XDMax;
extern int32 YDMax;
extern double XMin;
extern double YMin;
extern double XMax;
extern double YMax;
extern int32 TextJustify;
extern int32 TextAngle;

void graphics_get_scale(double *x1, double *y1, double *x2, double *y2);
void graphics_set_scale(double x1, double y1, double x2, double y2);
void graphics_get_draw_area_flag(int32 flag);
void graphics_get_draw_area(void);
void graphics_change_current_linestyle(int32 new, int32 *old);
void graphics_set_normal_scale(void);
void graphics_point(int32 x, int32 y);
void graphics_line(int32 x1, int32 y1, int32 x2, int32 y2);
void graphics_frect(int32 x1, int32 y1, int32 w, int32 h);
void graphics_put_text(int32 x, int32 y, char *str);
void graphics_init_x11(void);
void graphics_init_ps(void);
void graphics_init_svg(void);
void graphics_set_linestyle(int32 ls);
void graphics_put_text_x11(int32 x, int32 y, char *str);
void graphics_special_put_text_x11(int32 x, int32 y, char *str, int32 size);
void graphics_fancy_put_text_x11(int32 x, int32 y, char *str, int32 size,
                                 int32 font);
void graphics_scale_dxdy(double x, double y, double *i, double *j);
void graphics_scale_to_screen(double x, double y, int32 *i, int32 *j);
void graphics_scale_to_real(int32 i, int32 j, double *x, double *y);
void graphics_init_all(void);
void graphics_set_extra(void);
void graphics_reset_graph(void);
void graphics_get_graph(void);
void graphics_copy_graph(int32 i, int32 l);
void graphics_make_rot(double theta, double phi);
void graphics_scale3d(double x, double y, double z, double *xp, double *yp,
                      double *zp);
int32 graphics_threedproj(double x2p, double y2p, double z2p, double *xp,
                          double *yp);
void graphics_text3d(double x, double y, double z, char *s);
int32 graphics_threed_proj(double x, double y, double z, double *xp,
                           double *yp);
void graphics_point_3d(double x, double y, double z);
void graphics_line3dn(double xs1, double ys1, double zs1, double xsp1,
                      double ysp1, double zsp1);
void graphics_line3d(double x01, double y01, double z01, double x02, double y02,
                     double z02);
void graphics_line_3d(double x, double y, double z, double xp, double yp,
                      double zp);
void graphics_point_abs(double x1, double y1);
void graphics_bead_abs(double x1, double y1);
void graphics_frect_abs(double x1, double y1, double w, double h);
void graphics_line_abs(double x1, double y1, double x2, double y2);
void graphics_text_abs(double x, double y, char *text);
void graphics_fillin_text(char *old, char *new);
void graphics_fancy_text_abs(double x, double y, char *old, int32 size);
int32 graphics_clip3d(double x1, double y1, double z1, double x2, double y2,
                      double z2, double *x1p, double *y1p, double *z1p,
                      double *x2p, double *y2p, double *z2p);
int32 graphics_clip(double x1, double x2, double y1, double y2, double *x1_out,
                    double *y1_out, double *x2_out, double *y2_out);
void graphics_eq_symb(double *x, int32 type);
void graphics_reset_all_line_type(void);
void graphics_draw_many_lines(void);

#endif

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

typedef struct BoxList {
    int32 use, type, xuse;
    int32 n;
    int32 n0;
    int32 nwin, minwid, minhgt;
    Window up;
    Window dn;
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

typedef struct HistInfo {
    int32 nbins, nbins2, type, col, col2, fftc;
    double xlo;
    double xhi;
    double ylo, yhi;
    char cond[80];
} HistInfo;

extern HistInfo hist_inf;
extern int32 spec_col;
extern int32 spec_wid;
extern int32 spec_win;
extern int32 spec_col2;
extern int32 post_process;

extern int32 FOUR_HERE;

extern BoxList ParamBox;

int32 histogram_two_d(int32 col1, int32 col2, int32 ndat, int32 n1, int32 n2,
                      double xlo, double xhi, double ylo, double yhi);
void histogram_back(void);
int32 histogram_new_2d(void);
void histogram_new(int32 nbins, double zlo, double zhi, int32 col, int32 col2,
                   char *condition, int32 which);
void histogram_column_mean(void);
void histogram_compute_power(void);
int32 histogram_cross_spectrum(double *data, double *data2, int32 nr, int32 win,
                               int32 w_type, double *pow, int32 type);
void histogram_compute_sd(void);
void histogram_compute_fourier(void);
void histogram_compute_correl(void);
void histogram_compute_stacor(void);
void histogram_mycor(double *x, double *y, int32 n, double zlo, double zhi,
                     int32 nbins, double *z, int32 flag);
void histogram_compute(void);
void histogram_fft_xcorr(double *data1, double *data2, int32 length, int32 nlag,
                         double *cr, int32 flag);
void histogram_post_process_stuff(void);

#endif

#ifndef INIT_CONDS_H
#define INIT_CONDS_H

#include <X11/Xlib.h>
#include "read_dir.h"

typedef struct BoxListold {
    int32 use;
    int32 type;
    int32 n;
    Window base;
    Window cancel, ok, def, go;
    Window *w;
    Window *we;
    char **value;
    int32 mc, *off, *pos;
} BoxListold;

void init_conds_clone_ode(void);
int32 init_conds_find_user_name(int32 type, char *oname);
void init_conds_create_par_sliders(Window base, int32 x0, int32 h0);
void init_conds_resize_par_slides(int32 h);
void init_conds_slide_button_press(Window window);
int32 init_conds_file_selector(char *title, char *file, char *wild);
void init_conds_reset_sliders(void);
void init_conds_slide_release(Window window);

void init_conds_expose_slides(Window window);
void init_conds_enter_slides(Window window, int32 val);
void init_conds_slider_motion(XEvent event);
void init_conds_make_new_ic_box(void);
void init_conds_make_new_bc_box(void);
void init_conds_make_new_delay_box(void);
void init_conds_make_new_param_box(void);
void init_conds_initialize_box(void);
void init_conds_resize_par_box(Window win);
void init_conds_make_box_list(BoxList *b, char *wname, char *iname, int32 n, int32 type,
                   int32 use);
void init_conds_do_box_expose(Window window);
void init_conds_draw_one_box(BoxList b, int32 index);
void init_conds_redraw_params(void);
void init_conds_redraw_delays(void);
void init_conds_redraw_ics(void);
void init_conds_redraw_bcs(void);
void init_conds_box_enter_events(Window window, int32 yn);
void init_conds_box_buttons(Window window);
void init_conds_box_keypress(XEvent event, int32 *used);
void init_conds_man_ic(void);
void init_conds_new_parameter(void);
void init_conds_redo_stuff(void);
void init_conds_set_edit_params(BoxList *b, int32 i, char *string);

#endif

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdio.h>

extern int32 MakePlotFlag;
extern int32 SuppressOut;
extern int32 SuppressBounds;

typedef struct Range {
    char item[30];
    char item2[30];
    int32 steps, steps2, reset, oldic, index, index2, cycle, type, type2, movie;
    double plow, phigh, plow2, phigh2;
    int32 rtype;
} Range;

extern Range range;

extern int32 (*solver)(double *y, double *tim, double dt, int32 nt, int32 neq,
                       int32 *istart, double *work);

typedef struct XppVec {
    int32 nvec;
    int32 node;
    double *x;
} XppVec;

extern XppVec xpv;

extern int32 DelayErr;
extern double MyData[MAX_ODE];
extern double MyTime;
extern double MyData[MAX_ODE];
extern double MyTime;
extern int32 MyStart;
extern int32 RANGE_FLAG;
extern double LastTime;

void integrate_init_ar_ic(void);
void integrate_dump_range(FILE *fp, int32 f);
void integrate_init_range(void);
void integrate_cont_integ(void);
void integrate_swap_color(int32 *col, int32 rorw);
void integrate_set_cycle(int32 flag, int32 *icol);
int32 integrate_do_range(double *x, int32 flag);
void integrate_find_equilib_com(int32 com);
void integrate_batch(void);
void integrate_do_init_data(int32 com);
void integrate_run_now(void);
void usual_integrate_stuff(double *x);
int32 integrate_extract_ic_data(char *big);
void integrate_arr_ic_start(void);
void integrate_get_ic(int32 it, double *x);
int32 integrate_ode_int(double *y, double *t, int32 *istart, int32 ishow);
int32 integrate(double *t, double *x, double tend, double dt, int32 count,
                int32 nout, int32 *start);
void integrate_send_halt(void);
void integrate_send_output(double *y, double t);
void integrate_plot(double *oldxpl, double *oldypl, double *oldzpl, double *xpl,
                    double *ypl, double *zpl);
void integrate_export_data(FILE *fp);
void integrate_plot_the_graphs(double *xv, double *xvold, double ddt, int32 *tc,
                               int32 flag);
void integrate_restore(int32 i1, int32 i2);
void integrate_comp_color(double *v1, double *v2, int32 n, double dt);
void integrate_shoot(double *x, double *xg, double *evec, int32 sgn);
void integrate_shoot_easy(double *x);
void integrate_stop_integration(void);
int32 integrate_do_auto_range_go(void);
void integrate_silent_equilibria(void);

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
    PRE_NONE,
    PRE_LEFT,
    PRE_RIGHT,
    PRE_BOTH
};

/******************************************************************
 *                                                                *
 * enum : types of Gram-Schmidt routines                          *
 *----------------------------------------------------------------*
 * MODIFIED_GS  : The iterative solver uses the modified          *
 *                Gram-Schmidt routine modified_gs listed in this  *
 *                file.                                           *
 *                                                                *
 * CLASSICAL_GS : The iterative solver uses the classical         *
 *                Gram-Schmidt routine classical_gs listed in this *
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
typedef int32 (*ATimesFn)(void *A_data, Vector v, Vector z);

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
typedef int32 (*PSolveFn)(void *P_data, Vector r, Vector z, int32 lr);

/******************************************************************
 *                                                                *
 * Function: modified_gs                                           *
 *----------------------------------------------------------------*
 * modified_gs performs a modified Gram-Schmidt orthogonalization  *
 * of the Vector v[k] against the p unit N_Vectors at           *
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
 * If (k-p) < 0, then modified_gs uses p=k. The orthogonalized     *
 * v[k] is NOT normalized and is stored over the old v[k]. Once   *
 * the orthogonalization has been performed, the Euclidean norm   *
 * of v[k] is stored in (*new_vk_norm).                           *
 *                                                                *
 * modified_gs returns 0 to indicate success. It cannot fail.      *
 *                                                                *
 ******************************************************************/
int32 iterativ_modified_gs(Vector *v, double **h, int32 k, int32 p,
                           double *new_vk_norm);

/******************************************************************
 *                                                                *
 * Function: classical_gs                                          *
 *----------------------------------------------------------------*
 * classical_gs performs a classical Gram-Schmidt                  *
 * orthogonalization of the Vector v[k] against the p unit      *
 * N_Vectors at v[k-1], v[k-2], ..., v[k-p]. The parameters v, h, *
 * k, p, and new_vk_norm are as described in the documentation    *
 * for modified_gs.                                                *
 *                                                                *
 * temp is an Vector which can be used as workspace by the      *
 * classical_gs routine.                                           *
 *                                                                *
 * s is a length k array of reals which can be used as workspace  *
 * by the classical_gs routine.                                    *
 *                                                                *
 * classical_gs returns 0 to indicate success. It cannot fail.     *
 *                                                                *
 ******************************************************************/
int32 iterativ_classical_gs(Vector *v, double **h, int32 k, int32 p,
                            double *new_vk_norm, Vector temp, double *s);

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
int32 iterativ_qr_fact(int32 n, double **h, double *q, int32 job);

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
int32 iterativ_qr_sol(int32 n, double **h, double *q, double *b);

#endif

#ifndef KINESCOPE_H
#define KINESCOPE_H

void kinescope_do_movie_com(int32 c);
void kinescope_reset_film(void);
int32 kinescope_film_clip(void);

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

/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 * MAX(A, B) returns the maximum of A and B.                      *
 * ABS(A) returns the absolute value of A.                        *
 * SQR(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define ABS(A) ((A > 0) ? (A) : -(A))
#define SQR(A) ((A)*(A))

/******************************************************************
 *                                                                *
 * Function : UnitRoundoff                                        *
 * Usage    : double uround;                                        *
 *            uround = llnlmath_unit_roundoff();                            *
 *----------------------------------------------------------------*
 * UnitRoundoff returns the unit roundoff u for double floating     *
 * point arithmetic, where u is defined to be the largest         *
 * positive double such that 1.0 + u != 1.0.                        *
 *                                                                *
 ******************************************************************/
double llnlmath_unit_roundoff(void);

/******************************************************************
 *                                                                *
 * Function : RPowerI                                             *
 * Usage    : int32 exponent;                                       *
 *            double base, ans;                                     *
 *            ans = llnlmath_rpower_i(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerI returns the value base^exponent, where base is a double  *
 * and exponent is an int32.                                        *
 *                                                                *
 ******************************************************************/
double llnlmath_rpower_i(double base, int32 exponent);

/******************************************************************
 *                                                                *
 * Function : RPowerR                                             *
 * Usage    : double base, exponent, ans;                           *
 *            ans = llnlmath_rpower_i(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerR returns the value base^exponent, where both base and   *
 * exponent are reals. If base < 0.0, then RPowerR returns 0.0.   *
 *                                                                *
 ******************************************************************/
double llnlmath_rpower_r(double base, double exponent);
/******************************************************************
 *                                                                *
 * Function : RSqrt                                               *
 * Usage    : double sqrt_x;                                        *
 *            sqrt_x = llnlmath_rsqrt(x);                                  *
 *----------------------------------------------------------------*
 * llnlmath_rsqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
 * returns 0.0.                                                   *
 *                                                                *
 ******************************************************************/
double llnlmath_rsqrt(double x);

#endif

/*
Options are set accroding to an order of precedence

command line < mfile < .xpprc < default.opt

Add any options here that you might want to track.
*/
typedef struct OptionsSet {
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
    int32 euler_max_iter;
    int32 euler_tol;
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

void load_eqn_dump_torus(FILE *fp, int32 f);
void load_eqn(void);
void load_eqn_set_x_vals(void);
void load_eqn_set_all_vals(void);
void load_eqn_add_intern_set(char *name, char *does);
void load_eqn_extract_action(char *ptr);
void load_eqn_extract_internset(int32 j);
int32 load_eqn_msc(char *s1, char *s2);
void load_eqn_set_internopts(OptionsSet *mask);
void load_eqn_set_internopts_xpprc_and_comline(void);
void load_eqn_check_for_xpprc(void);
void load_eqn_stor_internopts(char *s1);
void load_eqn_set_option(char *s1, char *s2, int32 force, OptionsSet *mask);

#endif

#ifndef LUNCH_NEW_H
#define LUNCH_NEW_H

#include <stdio.h>

void lunch_file_inf(void);
void lunch_ps_write_pars(FILE *fp);
int32 lunch_read(FILE *fp);
void do_lunch(int32 f);
void lunch_io_parameter_file(char *fn, int32 flag);
void lunch_io_ic_file(char *fn, int32 flag);
void lunch_io_int(int32 *i, FILE *fp, int32 f, char *ss);
void lunch_io_double(double *z, FILE *fp, int32 f, char *ss);
void lunch_io_string(char *s, int32 len, FILE *fp, int32 f);

#endif

#ifndef MAIN_H
#define MAIN_H

#include <X11/Xlib.h>

extern int32 all_win_vis;
extern int32 use_ani_file;
extern int32 XPPBatch;
extern int32 batch_range;
extern int32 BatchEquil;
extern int32 PaperWhite;

extern int32 UserGradients;
extern int32 UserMinWidth;
extern int32 UserMinHeight;
extern int32 periodic;
extern int32 XPPVERBOSE;
extern int32 OVERRIDE_QUIET;
extern int32 OVERRIDE_LOGFILE;

extern int32 SLIDER1;
extern int32 SLIDER2;
extern int32 SLIDER3;

extern double SLIDER1LO;
extern double SLIDER2LO;
extern double SLIDER3LO;
extern double SLIDER1HI;
extern double SLIDER2HI;
extern double SLIDER3HI;

extern int32 do_tutorial;
extern char anifile[XPP_MAX_NAME];

extern double xpp_version_maj;
extern double xpp_version_min;
extern int32 Xup;

extern char batchout[XPP_MAX_NAME];
extern char UserOUTFILE[XPP_MAX_NAME];
extern int32 TrueColorFlag;
extern char font_name_big[100];
extern char font_name_small[100];
extern char PlotFormat[10];

extern Window command_pop;
extern GC gc_graph;
extern GC font_gc;

void *XMALLOC(usize size, const char *function, int32 line);
#ifdef MALLOC_DEBUG
#define xmalloc(X) XMALLOC(X, __func__, __LINE__)
#else
void *xmalloc(usize size);
#endif
void main_plot_command(int32 nit, int32 icount, int32 cwidth);
int32 main_my_abort(void);
void do_main(int32 argc, char **argv) __attribute__((noreturn));
void main_do_vis_env(void);
void main_do_events(uint32 min_wid, uint32 min_hgt) __attribute__((noreturn));
void main_bye_bye(void) __attribute__((noreturn));
void main_clr_scrn(void);
void main_redraw_all(void);
void main_commander(int32 ch);
Window main_init_win(uint32 bw, char *icon_name, char *win_name, int32 x,
                     int32 y, uint32 min_wid, uint32 min_hgt, int32 argc,
                     char **argv);
void main_top_button_draw(Window window);
void main_fix_window_size(Window window, int32 width, int32 height, int32 flag);
int32 main_get_command_width(void);

#endif

#ifndef MANY_POPS_H
#define MANY_POPS_H

#include <X11/Xlib.h>

extern int32 manual_expose;
extern int32 SimulPlotFlag;
extern Graph graph[MAXPOP];
extern Curve frz[MAXFRZ];
extern Graph *MyGraph;
extern int32 current_pop;
extern int32 num_pops;
extern int32 ActiveWinList[MAXPOP];

int32 many_pops_select_table(void);
void many_pops_get_intern_set(void);
void many_pops_make_icon(char *icon, int32 wid, int32 hgt, Window window);
void many_pops_title_text(char *string);
void many_pops_gtitle_text(char *string, Window win);
void many_pops_restore_off(void);
void many_pops_restore_on(void);
void many_pops_add_label(char *s, int32 x, int32 y, int32 size, int32 font);
void many_pops_draw_label(Window window);
void many_pops_add_grob(double xs, double ys, double xe, double ye, double size,
                        int32 type, int32 color);
void many_pops_edit_object_com(int32 com);
void many_pops_do_gr_objs_com(int32 com);
void many_pops_do_windows_com(int32 c);
void many_pops_init_grafs(int32 x, int32 y, int32 w, int32 h);
void many_pops_ps_restore(void);
void many_pops_svg_restore(void);
int32 many_pops_rotate_3dcheck(XEvent event);
void many_pops_do_motion_events(XEvent event);
void many_pops_do_expose(XEvent event);
void many_pops_resize_all(int32 wid, int32 hgt);
void many_pops_create_a_pop(void);
void many_pops_gr_col(void);
void many_pops_base_col(void);
void many_pops_small_gr(void);
void many_pops_small_base(void);
void many_pops_change_plot_vars(int32 k);
int32 many_pops_check_active_plot(int32 k);
void many_pops_make_active(int32 i, int32 flag);
void many_pops_hi_lite(Window wi);
void many_pops_canvas_xy(char *buf);
void many_pops_check_draw_button(XEvent event);
void many_pops_set_active_windows(void);

#endif

#ifndef MARKOV_H
#define MARKOV_H

#include <stdio.h>

extern int32 STOCH_FLAG;
extern int32 NWiener;

void markov_add_wiener(int32 index);
void markov_set_wieners(double dt, double *x, double t);
void markov_add(int32 nstate, char *name);
int32 build_markov(char **ma, char *name);
int32 markov_old_build(FILE *fptr, char *name);
void markov_compile_all(void);
void markov_make_gill_nu(double *nu, int32 n, int32 m, double *v);
void markov_one_gill_step(int32 meth, int32 nrxn, int32 *rxn, double *v);
void markov_do_stochast_com(int32 i);
void markov_mean_back(void);
void markov_variance_back(void);
void markov_append_stoch(int32 first, int32 length);
void markov_do_stats(int32 ierr);
double markov_poidev(double xm);
double markov_ndrand48(void);
void markov_nsrand48(int32 seed);

#endif

#ifndef MENU_H
#define MENU_H

extern int32 help_menu;

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

void menudrive_xpp_hlp(void);
void menudrive_message_box(char *m);
void menudrive_message_box_redraw(Window window);
void menudrive_message_box_kill(void);
int32 menudrive_two_choice(char *c1, char *c2, char *q, char *key);
int32 menudrive_get_mouse_xy(int32 *x, int32 *y);
void menudrive_flush_display(void);
void menudrive_clear_draw_window(void);
void menudrive_drw_all_scrns(void);
void menudrive_clr_all_scrns(void);
void menudrive_run_the_commands(int32 com);
void menudrive_do_stochast(void);
void menudrive_get_pmap_pars(void);
void menudrive_set_col_par(void);
void menudrive_make_adj(void);
void menudrive_do_gr_objs(void);
void menudrive_new_lookup(void);
void menudrive_find_bvp(void);
void menudrive_change_view(void);
void menudrive_do_windows(void);
void menudrive_add_a_curve(void);
void menudrive_do_movie(void);
void menudrive_do_torus(void);
void menudrive_window_zoom(void);
void menudrive_direct_field(void);
void menudrive_new_clines(void);
void menudrive_froz_cline_stuff(void);
void menudrive_find_equilibrium(void);
void menudrive_ini_data_menu(void);
void menudrive_new_param(void);
void menudrive_clear_screens(void);
void menudrive_x_vs_t(void);
void menudrive_redraw_them_all(void);
void menudrive_get_3d_par(void);
void menudrive_edit_xpprc(void);
void menudrive_do_tutorial(void);

#endif

#ifndef XPPMENU_H
#define XPPMENU_H

#include <X11/Xlib.h>

void add_menu(Window base, int32 j, int32 n, char **names, char *key,
              char **hint);
void create_the_menus(Window base);
void menu_help(void);
void menu_help_num(void);
void menu_help_file(void);
void menu_crossing(Window win, int32 yn);
void menu_expose(Window win);
void menu_button(Window win);
void menu_draw_help(void);

#endif

#ifndef MY_PS_H
#define MY_PS_H

extern int32 NoBreakLine;
extern int32 PS_FONTSIZE;

extern double PS_LW;
extern char PS_FONT[200];
extern FILE *psfile;
extern int32 PltFmtFlag;
extern int32 PSColorFlag;
extern int32 PSLines;
extern int32 LastPSX;
extern int32 LastPSY;
extern int32 LastPtLine;

int32 ps_init(char *filename, int32 color);
void ps_stroke(void);
void ps_do_color(int32 color);
void ps_end(void);
void ps_frect(int32 x, int32 y, int32 w, int32 h);
void ps_last_pt_off(void);
void ps_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void ps_linetype(int32 linetype);
void ps_point(int32 x, int32 y);
void ps_fnt(int32 cf, int32 scale);
void ps_show(char *str, int32 type);
void ps_abs(int32 x, int32 y);
void ps_special_put_text(int32 x, int32 y, char *str, int32 size);
void ps_text(int32 x, int32 y, char *str);

#endif

#ifndef MY_RHS_H
#define MY_RHS_H

int32 main(int32 argc, char **argv);
void main_rhs_extra(double *y__y, double t, int32 nod, int32 neq);
void main_rhs_set_fix(double t, double *y);
int32 main_rhs(double t, double *y, double *ydot, int32 neq);
void main_rhs_update_based_on_current(void);
void main_rhs_fix_only(void);
void main_rhs_only(double *ydot);

#endif

#ifndef MY_SVG_H
#define MY_SVG_H

extern FILE *svgfile;

int32 svg_init(char *filename);
void svg_do_color(int32 color);
void svg_end(void);
void svg_bead(void);
void svg_frect(int32 x, int32 y, int32 w, int32 h);
void svg_last_pt_off(void);
void svg_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void svg_linetype(int32 linetype);
void svg_point(int32 x, int32 y);
void special_put_text_svg(int32 x, int32 y, char *str, int32 size);
void svg_text(int32 x, int32 y, char *str);

#endif

#ifndef NULLCLINE_H
#define NULLCLINE_H

extern int32 DFBatch;
extern int32 NCBatch;
extern int32 XNullColor;
extern int32 YNullColor;
extern int32 DF_GRID;
extern int32 DF_FLAG;
extern int32 DOING_DFIELD;

extern char ColorVia[15];
extern double ColorViaLo;
extern double ColorViaHi;
extern int32 ColorizeFlag;

#include <stdio.h>

typedef struct NullClines {
    double *xn;
    double *yn;
    int32 nmx, nmy;
    int32 n_ix;
    int32 n_iy;
    struct NullClines *n, *p;
} NullClines;

void nullcline_create_new_cline(void);
void nullcline_froz_cline_stuff_com(int32 i);
int32 get_nullcline_floats(double **v, int32 *n, int32 who, int32 type);
void nullcline_add_froz(double *xn, int32 nmx, int32 n_ix, double *yn,
                        int32 nmy, int32 n_iy);
void nullcline_get_max_dfield(double *y, double *ydot, double u0, double v0,
                              double du, double dv, int32 n, int32 inx,
                              int32 iny, double *mdf);
void nullcline_redraw_dfield(void);
void nullcline_direct_field_com(int32 c);
void restore_nullclines(void);
void nullcline_new_clines_com(int32 c);
void new_nullcline(int32 course, double xlo, double ylo, double xhi, double yhi,
                   double *stor, int32 *npts);
void nullcline_do_batch_nclines(void);
void nullcline_do_batch_dfield(void);
void nullcline_set_colorization_stuff(void);
void silent_nullclines(void);
void nullcline_silent_dfields(void);

#endif

#ifndef NUMERICS_H
#define NUMERICS_H

extern double DELTA_T;
extern double TEND;
extern double T0;
extern double TRANS;
extern double NULL_ERR;
extern double EVEC_ERR;
extern double NEWT_ERR;
extern double BOUND;
extern double DELAY;
extern double TOLER;
extern double HMIN;
extern double HMAX;
extern double *fft_data;
extern double *hist_data;
extern double POIPLN;

extern int32 cv_bandflag;
extern int32 cv_bandupper;
extern int32 cv_bandlower;

extern int32 NMESH;
extern int32 NJMP;
extern int32 METHOD;
extern int32 color_flag;
extern int32 NC_ITER;
extern int32 EVEC_ITER;
extern int32 FOREVER;

extern int32 POIMAP;
extern int32 POIVAR;
extern int32 POISGN;
extern int32 SOS;

extern int32 HIST;
extern int32 HVAR;
extern int32 hist_ind;

void numerics_chk_volterra(void);
void numerics_quick_num(int32 com);
void numerics_get_num_par(int32 ch);
void numerics_chk_delay(void);
void numerics_set_delay(void);
void numerics_get_pmap_pars_com(int32 l);
void numerics_set_col_par_com(int32 i);
void numerics_do_meth(void);
void numerics_set_total(double total);
void numerics_user_set_color_par(int32 flag, char *via, double lo, double hi);
void numerics_compute_one_period(double period, double *x, char *name);
#endif

#ifndef ODESOL2_H
#define ODESOL2_H

extern int32 (*rhs_function)(double t, double *y, double *ydot, int32 neq);

int32 symplect3(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 discrete(double *y, double *tim, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 bak_euler(double *y, double *tim, double dt, int32 nt, int32 neq,
                int32 *istart, double *work);
int32 one_bak_step(double *y, double *t, double dt, int32 neq, double *yg,
                   double *yp, double *yp2, double *ytemp, double *errvec,
                   double *jac);
void odesol2_one_step_discrete(double *y, double dt, double *yp, int32 neq,
                               double *t);
void odesol2_one_step_symp(double *y, double h, double *f, int32 n, double *t);
void odesol2_one_step_euler(double *y, double dt, double *yp, int32 neq,
                            double *t);
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
int32 rb23(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
           double *work, int32 *ierr);
int32 rosen(double *y, double *tstart, double tfinal, int32 *istart, int32 n,
            double *work, int32 *ierr);
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy,
                 int32 neq, double eps, double scal);
void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold,
                  int32 n, double eps, double scal);

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

#define EV_MASK                                                                \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

#define BUT_MASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     EnterWindowMask | LeaveWindowMask)

extern int32 DisplayWidth;
extern int32 DisplayHeight;
extern int32 screen;
extern Atom deleteWindowAtom;
extern Window main_win;
extern Window info_pop;
extern Window draw_win;
extern int32 DCURY;
extern int32 DCURX;
extern int32 CURY_OFF;
extern int32 DCURXs;
extern int32 DCURYs;
extern int32 CURY_OFFs;
extern int32 xor_flag;
extern GC gc;
extern GC small_gc;
extern uint32 MyBackColor;
extern uint32 MyForeColor;

extern int32 TipsFlag;
extern char UserBlack[8];
extern char UserWhite[8];
extern int32 UserGradients;
extern char UserMainWinColor[8];
extern char UserDrawWinColor[8];
extern char UserBGBitmap[XPP_MAX_NAME];
extern uint32 MyMainWinColor;
extern uint32 MyDrawWinColor;
extern uint32 GrFore;
extern uint32 GrBack;

extern int32 SCALEY;
extern int32 DCURYb;
extern int32 CURY_OFFb;
extern int32 DCURYs;
extern int32 DCURXs;
extern FILE *logfile;
extern int32 tfBell;
extern char SLIDER1VAR[20];
extern char SLIDER2VAR[20];
extern char SLIDER3VAR[20];

extern OptionsSet notAlreadySet;
extern XFontStruct *font_small;

/*  This is a string box widget which handles a list of
        editable strings
 */

typedef struct StringBox {
    Window base, ok, cancel;
    Window win[MAX_N_SBOX];
    char name[MAX_N_SBOX][MAX_LEN_SBOX], value[MAX_N_SBOX][MAX_LEN_SBOX];
    int32 n;
    int32 hot;
    int32 hgt, wid;
    int32 hh[MAX_N_SBOX];
} StringBox;

extern int32 NUPAR;
extern int32 NEQ;
extern int32 NODE;
extern int32 NMarkov;
extern char upar_names[MAX_PAR][MAX_ODE_NAME_LENGTH];
extern char uvar_names[MAX_ODE][MAX_ODE_NAME_LENGTH];
extern char *color_names[];

#define SB_PLOTTABLE 0
#define SB_VARIABLE 1
#define SB_PARAMETER 2
#define SB_PARVAR 3
#define SB_COLOR 4
#define SB_MARKER 5
#define SB_METHOD 6

void pop_list_set_window_title(Window win, char *string);
void pop_list_make_scrbox_lists(void);
int32 pop_list_do_string_box(int32 n, int32 row, int32 col, char *title, char **names,
                    char values[][MAX_LEN_SBOX], int32 maxchar);
void pop_list_do_hilite_text(char *name, char *value, int32 flag, Window window,
                    int32 pos);
void pop_list_new_editable(StringBox *sb, int32 inew, int32 *pos, int32 *col,
                  int32 *done, Window *w);
Window pop_list_make_fancy_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw);
Window pop_list_make_unmapped_window(Window root, int32 x, int32 y, int32 width,
                            int32 height, int32 bw);
Window pop_list_make_plain_unmapped_window(Window root, int32 x, int32 y, int32 width,
                                  int32 height, int32 bw);
Window pop_list_make_window(Window root, int32 x, int32 y, int32 width, int32 height,
                   int32 bw);
Window pop_list_make_plain_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw);
void pop_list_respond_box(char *button, char *message);
void pop_list_message_box(Window *w, int32 x, int32 y, char *message);
void pop_list_expose_choice(char *choice1, char *choice2, char *msg, Window c1,
                   Window c2, Window wm, Window window);
int32 pop_list_two_choice(char *choice1, char *choice2, char *string, char *key, int32 x,
                 int32 y, Window window, char *title);
int32 pop_list_yes_no_box(void);
int32 pop_up_list(Window *root, char *title, char **list, char *key, int32 n,
                  int32 max, int32 def, int32 x, int32 y, char **hints,
                  Window hwin, char *httxt);
Window pop_list_make_unmapped_icon_window(Window root, int32 x, int32 y, int32 width,
                                 int32 height, int32 bw, uchar *icdata);
Window pop_list_make_icon_window(Window root, int32 x, int32 y, int32 width,
                        int32 height, int32 bw, uchar *icdata);

#endif

#ifndef PP_SHOOT_H
#define PP_SHOOT_H

#include <stdio.h>

void pp_shoot_do_bc(double *y__0, double t0, double *y__1, double t1, double *f,
                    int32 n);
void pp_shoot_compile_bvp(void);
void pp_shoot_reset_bvp(void);
void pp_shoot_init_shoot_range(char *s);
void pp_shoot_dump_shoot_range(FILE *fp, int32 f);
void pp_shoot_find_bvp_com(int32 com);
void pp_shoot_bv(double *y, double *yend, double err, double eps, int32 maxit,
             int32 *iret, int32 n, int32 ishow, int32 iper, int32 ipar,
             int32 ivar, double sect);

#endif

#ifndef RUBBER_H
#define RUBBER_H

#include <X11/Xlib.h>

int32 rubber(int32 *x1, int32 *y1, int32 *x2, int32 *y2, Window window,
             int32 f);

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

void scrngif_set_global_map(int32 flag);
void scrngif_end_ani_gif(FILE *fp);
void scrngif_add_ani_gif(Window win, FILE *fp, int32 count);
void scrngif_screen_to_gif(Window win, FILE *fp);
void scrngif_get_global_colormap(Window win);

#endif

#ifndef SIMPLENET_H
#define SIMPLENET_H

double simplenet_network_value(double x, int32 i);
int32 simplenet_add_spec_fun(char *name, char *rhs);
void simplenet_add_special_name(char *name, char *rhs);
int32 simplenet_add_vectorizer(char *name, char *rhs);
void simplenet_add_vectorizer_name(char *name, char *rhs);
void simplenet_eval_all_nets(void);
void simplenet_update_all_ffts(void);
void simplenet_fft_conv(int32 it, int32 n, double *values, double *yy, double *fftr,
              double *ffti, double *dr, double *di);
double simplenet_vector_value(double x, int32 i);

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

#include "vector.h"

/******************************************************************
 *                                                                *
 * Types: SpgmrMemRec, SpgmrMem                                   *
 *----------------------------------------------------------------*
 * SpgmrMem is a pointer to an SpgmrMemRec which contains         *
 * the memory needed by spgmr_solve. The SpgmrMalloc routine       *
 * returns a pointer of type SpgmrMem which should then be passed *
 * in subsequent calls to spgmr_solve. The SpgmrFree routine frees *
 * the memory allocated by SpgmrMalloc.                           *
 *                                                                *
 * N is the linear system size.                                   *
 *                                                                *
 * l_max is the maximum Krylov dimension that spgmr_solve will be  *
 * permitted to use.                                              *
 *                                                                *
 * V is the array of Krylov basis vectors v_1, ..., v_(l_max+1),  *
 * stored in V[0], ..., V[l_max], where l_max is the second       *
 * parameter to SpgmrMalloc. Each v_i is a length N vector of     *
 * type Vector. (N is the first parameter to SpgmrMalloc and    *
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
 * xcor is a length N vector (type Vector) which holds the      *
 * scaled, preconditioned correction to the initial guess.        *
 *                                                                *
 * yg is a length (l_max+1) array of reals used to hold "int16"   *
 * vectors (e.g. y and g).                                        *
 *                                                                *
 * vtemp is a length N vector (type Vector) used as temporary   *
 * vector storage during calculations.                            *
 *                                                                *
 ******************************************************************/
typedef struct {
    int64 N;
    int32 l_max;

    Vector *V;
    double **Hes;
    double *givens;
    Vector xcor;
    double *yg;
    Vector vtemp;
} SpgmrMemRec, *SpgmrMem;

/******************************************************************
 *                                                                *
 * Function : SpgmrMalloc                                         *
 *----------------------------------------------------------------*
 * SpgmrMalloc allocates the memory used by spgmr_solve. It        *
 * returns a pointer of type SpgmrMem which the user of the       *
 * SPGMR package should pass to spgmr_solve. The parameter N       *
 * is the size of the system to be solved by spgmr_solve and l_max *
 * is the maximum Krylov dimension that spgmr_solve will be        *
 * permitted to use. The parameter machEnv is a pointer to        *
 * machine environment-specific information. Pass NULL in the     *
 * ordinary sequential case (see vector.h). This routine returns  *
 * NULL if there is a memory request failure.                     *
 *                                                                *
 ******************************************************************/
SpgmrMem spgmr_malloc(int64 N, int32 l_max);

/******************************************************************
 *                                                                *
 * Function : spgmr_solve                                          *
 *----------------------------------------------------------------*
 * spgmr_solve solves the linear system Ax = b using the SPGMR     *
 * method. The return values are given by the symbolic constants  *
 * below. The first spgmr_solve parameter is a pointer to memory   *
 * allocated by a prior call to SpgmrMalloc. The system size N    *
 * passed in the call to SpgmrMalloc should be the same as the    *
 * length of all Vector arguments passed to spgmr_solve.         *
 *                                                                *
 * mem is the pointer returned by SpgmrMalloc to the structure    *
 * containing the memory needed by spgmr_solve.                    *
 *                                                                *
 * A_data is a pointer to information about the coefficient       *
 * matrix A. This pointer is passed to the user-supplied function *
 * atimes.                                                        *
 *                                                                *
 * x is the initial guess x_0 upon entry and the solution         *
 * Vector upon exit with return value SPGMR_SUCCESS or          *
 * SPGMR_RES_REDUCED. For all other return values, the output x   *
 * is undefined.                                                  *
 *                                                                *
 * b is the right hand side Vector. It is undisturbed by this   *
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
 * sx is the Vector of positive scale factors for x (not        *
 * tested). Pass NULL if x scaling not required.                  *
 *                                                                *
 * sb is the Vector of positive scale factors for b (not        *
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
 * to be filled in by spgmr_solve.                                 *
 *                                                                *
 * nli is a pointer to the number of linear iterations done in    *
 * the execution of spgmr_solve. The caller is responsible for     *
 * allocating the memory (*nli) to be filled in by spgmr_solve.    *
 *                                                                *
 * nps is a pointer to the number of calls made to psolve during  *
 * the execution of spgmr_solve. The caller is responsible for     *
 * allocating the memory (*nps) to be filled in by spgmr_solve.    *
 *                                                                *
 * Note.. Repeated calls can be made to spgmr_solve with varying   *
 * input arguments. If, however, the problem size N or the        *
 * maximum Krylov dimension l_max changes, then a call to         *
 * SpgmrMalloc must be made to obtain new memory for spgmr_solve   *
 * to use.                                                        *
 *                                                                *
 ******************************************************************/
int32 spgmr_solve(SpgmrMem mem, void *A_data, Vector x, Vector b, int32 pretype,
                  int32 gstype, double delta, int32 max_restarts, void *P_data,
                  Vector sx, Vector sb, ATimesFn atimes, PSolveFn psolve,
                  double *res_norm, int32 *nli, int32 *nps);

/* Return values for spgmr_solve */

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
void spgmr_free(SpgmrMem mem);

#endif

#ifndef STIFF_H
#define STIFF_H

void stiff_jacobn(double x, double *y, double *dfdx, double *dermat, double eps,
                  double *work, int32 n);
int32 stiff_adaptive(double *ystart, int32 nvar, double *xs, double x2,
                     double eps, double *hguess, double hmin, double *work,
                     int32 *ier, double epjac, int32 iflag, int32 *jstart);
int32 stiff_gadaptive(double *ystart, int32 nvar, double *xs, double x2,
                      double eps, double *hguess, double hmin, double *work,
                      int32 *ier, double epjac, int32 iflag);
int32 stiff(double y[], double dydx[], int32 n, double *x, double htry,
            double eps, double yscal[], double *hdid, double *hnext,
            double *work, double epjac, int32 *ier);
int32 stiff_rkqs(double *y, double *dydx, int32 n, double *x, double htry,
                 double eps, double *yscal, double *hdid, double *hnext,
                 double *work, int32 *ier);
void stiff_rkck(double *y, double *dydx, int32 n, double x, double h,
                double *yout, double *yerr, double *work);

#endif

#ifndef STORAGE_H
#define STORAGE_H

extern double **storage;
extern double *WORK;
extern int32 IWORK[10000];

void storage_init_alloc_info(void);
void storage_alloc_meth(void);
void storage_init_stor(int32 nrow, int32 ncol);
int32 storage_realloc(int32 ncol, int32 nrow);

#endif

#ifndef TABULAR_H
#define TABULAR_H

typedef struct Tabular {
    double xlo, xhi, dx;
    double *y;
    double *x;
    int32 n, flag, interp, autoeval;
    int32 xyvals;
    /* flag=0 if virgin array, flag=1 if already allocated; flag=2 for function
                             interp=0 for normal interpolation, interp=1 for
       'step' interp=2 for cubic spline table   and finally, xyvals=1 if both x
       and y vals are needed (xyvals=0 is faster lookup )*/
    char filename[128], name[12];
} Tabular;

extern Tabular my_table[MAX_TAB];

void tabular_set_auto_eval_flags(int32 f);
void tabular_set_table_name(char *name, int32 index);
void tabular_new_lookup_com(int32 i);
double tabular_lookup(double x, int32 index);
void tabular_init_table(void);
void tabular_redo_all_fun_tables(void);
int32 tabular_create_fun(int32 npts, double xlo, double xhi, char *formula,
                         int32 index);
int32 tabular_load_table(char *filename, int32 index);
int32 tabular_get_lookup_len(int32 i);

#endif

#ifndef TORUS_H
#define TORUS_H

#include <X11/Xlib.h>

void do_torus_com(int32 c);

#endif

#ifndef TXTREAD_H
#define TXTREAD_H

void txt_view_events(XEvent event);
void txt_redraw_view(Window window);
void txt_init_view(void);
void txt_make_view(void);

#endif

#ifndef VOLTERRA2_H
#define VOLTERRA2_H

extern int32 AutoEvaluate;

double volterra_ker_val(int32 in);
void volterra_alloc_memory(void);
void volterra_allocate(int32 npts, int32 flag);
void volterra_re_evaluate_kernels(void);
void volterra_alloc_kernels(int32 flag);
void volterra_init_sums(double t0, int32 n, double dt, int32 i0, int32 iend,
                        int32 ishift);
int32 volterra(double *y, double *t, double dt, int32 nt, int32 neq,
               int32 *istart, double *work);
int32 volterra_step(double *y, double t, double dt, int32 neq, double *yg,
                    double *yp, double *yp2, double *errvec, double *jac);

#endif

#endif
