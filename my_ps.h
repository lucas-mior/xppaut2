#ifndef _my_ps_h_
#define _my_ps_h_
#include "integers.h"

int32 ps_init(char *filename, int32 color);
void ps_stroke(void);
void ps_do_color(int32 color);
void ps_setcolor(int32 color);
void ps_end(void);
void ps_bead(int32 x, int32 y);
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
