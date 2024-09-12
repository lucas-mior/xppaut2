
#ifndef _my_svg_h_
#define _my_svg_h_
#include "integers.h"

int32 svg_init(char *filename, int32 color);
void svg_stroke(void);
void svg_do_color(int32 color);
void svg_setcolor(int32 color);
void svg_end(void);
void svg_bead(int32 x, int32 y);
void svg_frect(int32 x, int32 y, int32 w, int32 h);
void svg_last_pt_off(void);
void svg_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2);
void chk_svg_lines(void);
void svg_linetype(int32 linetype);
void svg_point(int32 x, int32 y);
void svg_write(char *str);
void svg_fnt(int32 cf, int32 scale);
void svg_show(char *str, int32 type);
void svg_abs(int32 x, int32 y);
void svg_rel(int32 x, int32 y);
void special_put_text_svg(int32 x, int32 y, char *str, int32 size);
void fancy_svg_text(int32 x, int32 y, char *str, int32 size, int32 font);
void svg_text(int32 x, int32 y, char *str);

#endif
