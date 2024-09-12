#ifndef _graphics_h
#define _graphics_h
#include "integers.h"

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
void scale_to_real(int32 i, int32 j, float *x, float *y);
void init_all_graph(void);
void set_extra_graphs(void);
void reset_graph(void);
void get_graph(void);
void init_graph(int32 i);
void copy_graph(int32 i, int32 l);
void make_rot(double theta, double phi);
void scale3d(double x, double y, double z, float *xp, float *yp, float *zp);
double proj3d(double theta, double phi, double x, double y, double z, int32 in);
int32 threedproj(double x2p, double y2p, double z2p, float *xp, float *yp);
void text3d(double x, double y, double z, char *s);
void text_3d(double x, double y, double z, char *s);
int32 threed_proj(double x, double y, double z, float *xp, float *yp);
void point_3d(double x, double y, double z);
void line3dn(double xs1, double ys1, double zs1, double xsp1, double ysp1,
             double zsp1);
void line3d(double x01, double y01, double z01, double x02, double y02,
            double z02);
void line_3d(double x, double y, double z, double xp, double yp, double zp);
void pers_line(double x, double y, double z, double xp, double yp, double zp);
void rot_3dvec(double x, double y, double z, float *xp, float *yp, float *zp);
void point_abs(double x1, double y1);
void line_nabs(double x1_out, double y1_out, double x2_out, double y2_out);
void bead_abs(double x1, double y1);
void frect_abs(double x1, double y1, double w, double h);
void line_abs(double x1, double y1, double x2, double y2);
void text_abs(double x, double y, char *text);
void fillintext(char *old, char *new);
void fancy_text_abs(double x, double y, char *old, int32 size, int32 font);
int32 clip3d(double x1, double y1, double z1, double x2, double y2, double z2,
             float *x1p, float *y1p, float *z1p, float *x2p, float *y2p,
             float *z2p);
int32 clip(double x1, double x2, double y1, double y2, float *x1_out,
           float *y1_out, float *x2_out, float *y2_out);
void eq_symb(double *x, int32 type);
void draw_symbol(double x, double y, double size, int32 my_symb);
void reset_all_line_type(void);

#endif
