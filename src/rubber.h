
#ifndef _rubber_h
#define _rubber_h
#include "integers.h"

#include <X11/Xlib.h>

int32 rubber(int32 *x1, int32 *y1, int32 *x2, int32 *y2, Window w, int32 f);
void rbox(int32 i1, int32 j1, int32 i2, int32 j2, Window w, int32 f);

#endif
