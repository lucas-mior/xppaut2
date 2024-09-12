#ifndef _calc_h_
#define _calc_h_
#include "integers.h"

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
