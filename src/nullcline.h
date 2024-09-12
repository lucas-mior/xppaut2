#ifndef _nullcline_h_
#define _nullcline_h_
#include "integers.h"

#include <stdio.h>

typedef struct {
    float x, y, z;
} Pt;

typedef struct nclines {
    float *xn, *yn;
    int32 nmx, nmy;
    int32 n_ix, n_iy;
    struct nclines *n, *p;
} NCLINES;

void create_new_cline();
void froz_cline_stuff_com(int32 i);
void do_range_clines(void);
void start_ncline(void);
void clear_froz_cline(void);
int32 get_nullcline_floats(float **v, int32 *n, int32 who, int32 type);
void save_frozen_clines(char *fn);
void redraw_froz_cline(int32 flag);
void add_froz_cline(float *xn, int32 nmx, int32 n_ix, float *yn, int32 nmy,
                    int32 n_iy);
void get_max_dfield(double *y, double *ydot, double u0, double v0, double du,
                    double dv, int32 n, int32 inx, int32 iny, double *mdf);
void redraw_dfield(void);
void direct_field_com(int32 c);
void save_the_nullclines(void);
void restore_nullclines(void);
void dump_clines(FILE *fp, float *x, int32 nx, float *y, int32 ny);
void dump_clines_old(FILE *fp, float *x, int32 nx, float *y, int32 ny);
void restor_null(float *v, int32 n, int32 d);
void new_clines_com(int32 c);
void new_nullcline(int32 course, double xlo, double ylo, double xhi, double yhi,
                   float *stor, int32 *npts);
void stor_null(double x1, double y1, double x2, double y2);
float fnull(double x, double y);
int32 interpolate(Pt p1, Pt p2, double z, float *x, float *y);
void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4);
void triangle_contour(Pt p1, Pt p2, Pt p3);
void do_cline(int32 ngrid, double x1, double y1, double x2, double y2);
void do_batch_nclines();
void do_batch_dfield();

#endif
