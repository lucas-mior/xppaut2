#ifndef STRUCT_H
#define STRUCT_H
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

typedef struct RangeInfo {
    double xlo;
    double xhi;
    char rv[10];
    int32 nstep;
    int32 ic;
    int32 stor;
} RangeInfo;

typedef struct IcBox {
    Window base;
    Window ok;
    Window cancel;
    Window old;
    Window last;
    Window more;
    Window range;
    Window wrlo;
    Window wrhi;
    Window wstep;
    Window wreset;
    Window woldic;
    RangeInfo *rinf;
    double *yold, *y, *ylast;
    int32 n;
    int32 node;
    char **name;
    char ascval[MAX_ODE][ICLENGTH];
    Window wname[ICMAX];
    Window wval[ICMAX];
} IcBox;

typedef struct Graph {
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
    double xmin;
    double ymin;
    double zmin;
    double xmax;
    double ymax;
    double zmax;
    double xorg;
    double yorg;
    double zorg;
    double xbar;
    double ybar;
    double zbar;
    double dx;
    double dy;
    double dz;
    int32 xv[MAXPERPLOT];
    int32 yv[MAXPERPLOT];
    int32 zv[MAXPERPLOT];
    int32 line[MAXPERPLOT];
    int32 color[MAXPERPLOT];
    double Theta;
    double Phi;
    double ZPlane;
    double ZView;
    double xlo;
    double ylo;
    double xhi;
    double yhi;
    double oldxlo;
    double oldxhi;
    double oldylo;
    double oldyhi;
    int32 grtype;
    int32 ThreeDFlag;
    int32 TimeFlag;
    int32 PerspFlag;
    int32 xshft;
    int32 yshft;
    int32 zshft;
    int32 xorgflag;
    int32 yorgflag;
    int32 zorgflag;
    int32 ColorFlag;
    int32 ColorValue;
    char xlabel[30];
    char ylabel[30];
    char zlabel[30];
    char gr_info[256];
} Graph;

typedef struct TextGC {
    GC gc;
    int32 dx;
    int32 dy;
    int32 yoff;
    uint32 fcol;
    uint32 bcol;
} TextGC;

typedef struct Curve {
    Window window;
    char key[20];
    char name[10];
    int16 use;
    int16 type;
    double *xv, *yv, *zv;
    int32 len;
    int32 color;
} Curve;

typedef struct NullCline {
    Window window;
    char name[10];
    int16 use;
    double *x_n;
    double *y_n;
    int32 ix;
    int32 iy;
    int32 num_x;
    int32 num_y;
} NullCline;

typedef struct Dialog {
    Window mes;
    Window ok;
    Window cancel;
    Window input;
    Window base;
    char mes_s[MAXCHAR];
    char input_s[MAXCHAR];
    char ok_s[MAXCHAR];
    char cancel_s[MAXCHAR];
} Dialog;

typedef struct ChoiceBox {
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
} ChoiceBox;

typedef struct Param {
    Window window;
    char name[MAXCHAR];
    char value[MAXCHAR];
} Param;

typedef struct ParamBox2 {
    Window base;
    char title[MAXCHAR];
    Param *p;
    int32 n;
    Window ok;
    Window cancel;
} ParamBox2;

typedef struct TChoice {
    char name[10];
    char value[80];
    Window window;
} TChoice;

typedef struct TxtChoice {
    char title[100];
    Window who;
    Window what;
    Window cancel;
    Window ok;
    TChoice tc[100];
} TxtChoice;

#endif
