#ifndef struct_h
#define struct_h
#include "integers.h"

#include "xpplim.h"
#define MAXCHAR 60
#define MAXENTRY 20
#define RADIO 0
#define CHOICE 1
#define ICMAX 25

#define MAXPERPLOT 10
#define MAXFRZ 26
#define MAXPOP 21

#define MAXNCLINE 26

#define ICLENGTH 30
#define NAMELENGTH 10

typedef struct {
    double xlo;
    double xhi;
    char rv[10];
    int32 nstep, ic, stor;
} RANGE_INFO;

typedef struct {
    Window base, ok, cancel, old, last, more, range;
    Window wrlo, wrhi, wstep, wreset, woldic;
    RANGE_INFO *rinf;
    double *yold, *y, *ylast;
    int32 n;
    int32 node;
    char **name;
    char ascval[MAX_ODE][ICLENGTH];
    Window wname[ICMAX], wval[ICMAX];
} IC_BOX;

typedef struct {
    Window window;
    Window w_info;

    int32 Use;
    int32 state;
    int32 Restore;
    int32 Nullrestore;
    int32 x0;
    int32 y0;
    int32 Width;
    int32 Height;
    int32 x11Wid;
    int32 x11Hgt;
    int32 nvars;
    double rm[3][3];
    double min_scale;
    double color_scale;
    double xmin, ymin, zmin, xmax, ymax, zmax, xorg, yorg, zorg;
    double xbar, ybar, zbar, dx, dy, dz;
    int32 xv[MAXPERPLOT], yv[MAXPERPLOT], zv[MAXPERPLOT];
    int32 line[MAXPERPLOT], color[MAXPERPLOT];
    double Theta;
    double Phi;
    double ZPlane, ZView;
    double xlo, ylo, xhi, yhi, oldxlo, oldxhi, oldylo, oldyhi;
    int32 grtype, ThreeDFlag, TimeFlag, PerspFlag;
    int32 xshft, yshft, zshft;
    int32 xorgflag, yorgflag, zorgflag;
    int32 ColorFlag;
    int32 ColorValue;
    char xlabel[30], ylabel[30], zlabel[30];
    char gr_info[256];
} GRAPH;

typedef struct {
    GC gc;
    int32 dx, dy, yoff;
    uint32 fcol;
    uint32 bcol;
} TEXTGC;

typedef struct Curve {
    Window window;
    char key[20], name[10];
    int16 use;
    int16 type;
    double *xv, *yv, *zv;
    int32 len;
    int32 color;
} Curve;

typedef struct {
    Window window;
    char name[10];
    int16 use;
    double *x_n;
    double *y_n;
    int32 ix, iy, num_x, num_y;
} NCLINE;

typedef struct {
    Window mes;
    Window ok;
    Window cancel;
    Window input;
    Window base;
    char mes_s[MAXCHAR];
    char input_s[MAXCHAR];
    char ok_s[MAXCHAR];
    char cancel_s[MAXCHAR];
} DIALOG;

typedef struct {
    char title[MAXCHAR];
    int32 n;
    Window base;
    Window ok;
    Window cancel;
    int16 type;
    int32 mc;
    Window cw[MAXENTRY];
    char **name;
    int32 *flag;
} CHOICE_BOX;

typedef struct {
    Window window;
    char name[MAXCHAR];
    char value[MAXCHAR];
} PARAM;

typedef struct {
    Window base;
    char title[MAXCHAR];
    PARAM *p;
    int32 n;
    Window ok;
    Window cancel;
} PARAM_BOX;

typedef struct {
    char name[10];
    char value[80];
    Window window;
} TCHOICE;

typedef struct {
    char title[100];
    Window who, what, cancel, ok;
    TCHOICE tc[100];
} TXTCHOICE;

#endif
