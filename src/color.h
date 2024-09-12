#ifndef _color_h_
#define _color_h_
#include "integers.h"

#include <X11/Xlib.h>

void tst_color(Window w);
void set_scolor(int32 col);
void set_color(int32 col);
void make_cmaps(int32 *r, int32 *g, int32 *b, int32 n, int32 type);
int32 rfun(double y, int32 per);
int32 gfun(double y, int32 per);
int32 bfun(double y, int32 per);
void NewColormap(int32 type);
void get_ps_color(int32 i, float *r, float *g, float *b);
void get_svg_color(int32 i, int32 *r, int32 *g, int32 *b);
void MakeColormap(void);
int32 ColorMap(int32 i);

#endif
