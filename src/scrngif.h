#ifndef _scrngif_h_
#define _scrngif_h_
#include "integers.h"

#include <stdio.h>
#include <X11/Xlib.h>

typedef struct GifTree {
    char typ;         /* terminating, lookup, or search */
    int32 code;       /* the code to be output */
    unsigned char ix; /* the color map index */
    struct GifTree **node, *nxt, *alt;
} GifTree;

typedef struct {
    unsigned char r, g, b;
} GIFCOL;

void set_global_map(int32 flag);
int32 ppmtopix(unsigned char r, unsigned char g, unsigned char b, int32 *n);
void end_ani_gif(FILE *fp);
void add_ani_gif(Window win, FILE *fp, int32 count);
void screen_to_gif(Window win, FILE *fp);
void get_global_colormap(Window win);
void local_to_global(void);
int32 use_global_map(unsigned char *pixels, unsigned char *ppm, int32 h,
                     int32 w);
int32 make_local_map(unsigned char *pixels, unsigned char *ppm, int32 h,
                     int32 w);
void gif_stuff(Window win, FILE *fp, int32 task);
void write_global_header(int32 cols, int32 rows, FILE *dst);
void GifLoop(FILE *fout, uint32 repeats);
void write_local_header(int32 cols, int32 rows, FILE *fout, int32 colflag,
                        int32 delay);
void make_gif(unsigned char *pixels, int32 cols, int32 rows, FILE *dst);
int32 GifEncode(FILE *fout, unsigned char *pixels, int32 depth, int32 siz);
void ClearTree(int32 cc, GifTree *root);
unsigned char *AddCodeToBuffer(int32 code, short n, unsigned char *buf);

#endif
