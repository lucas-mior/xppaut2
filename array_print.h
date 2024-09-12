#ifndef _array_print_h_
#define _array_print_h_
#include "integers.h"

/* array_print.c */
int32 array_print(char *filename, char *xtitle, char *ytitle, char *bottom,
                int32 nacross, int32 ndown, int32 col0, int32 row0, int32 nskip,
                int32 ncskip, int32 maxrow, int32 maxcol, float **data, double zmin,
                double zmax, double tlo, double thi, int32 type);
void ps_replot(float **z, int32 col0, int32 row0, int32 nskip, int32 ncskip, int32 maxrow,
               int32 maxcol, int32 nacross, int32 ndown, double zmin, double zmax,
               int32 type);
void ps_begin(double xlo, double ylo, double xhi, double yhi, double sx,
              double sy);
void ps_convert(double x, double y, float *xs, float *ys);
void ps_col_scale(double y0, double x0, double dy, double dx, int32 n, double zlo,
                  double zhi, int32 type, double mx);
void ps_boxit(double tlo, double thi, double jlo, double jhi, double zlo,
              double zhi, char *sx, char *sy, char *sb, int32 type);
void ps_close(void);
void ps_setline(double fill, int32 thick);
void ps_put_char(int32 ch, float *x, float *y);
void ps_text2(char *str, double xr, double yr, int32 icent);
void ps_line2(double x1r, double y1r, double x2r, double y2r);
void ps_set_text(double angle, double slant, double x_size, double y_size);
void ps_rect(double x, double y, double wid, double len);
void ps_bar(double x, double y, double wid, double len, double fill, int32 flag);
void ps_rgb_bar(double x, double y, double wid, double len, double fill,
                int32 flag, int32 rgb);
void ps_hsb_bar(double x, double y, double wid, double len, double fill,
                int32 flag);

#endif
