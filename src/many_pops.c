#include "xpplim.h"
#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include "graph.bitmap"
#include "struct.h"

#define MAXLAB 50
#define MAXGROB 400
#define POINTER 0
#define ARROW 1
#define MARKER 2 /* markers start at 2  there are several of them */

#define WDMARK .001
#define HTMARK .0016

typedef struct {
    double xs, ys, xe, ye;
    double size;
    int16 use;
    Window window;
    int32 type;
    int32 color;
} GROB;

typedef struct {
    int32 type;
    int32 color;
    int32 number, start, skip;
    double size;
} MARKINFO;

MARKINFO markinfo = {2, 0, 1, 0, 1, 1.0};

int32 manual_expose = 0;
extern char *info_message;
extern Atom deleteWindowAtom;
LABEL lb[MAXLAB];
GROB grob[MAXGROB];
GRAPH graph[MAXPOP];
CURVE frz[MAXFRZ];
GRAPH *MyGraph;
extern int32 screen;
extern int32 SCALEY;
extern int32 CURY_OFF;
extern int32 CURY_OFFs;
extern int32 DCURYs;
extern int32 DCURXs;
extern int32 DCURYb;
int32 SimulPlotFlag = 0;
extern int32 storind;
extern int32 PltFmtFlag;
extern char *text_hint[];
extern char *edit_hint[];
extern char *no_hint[];
extern Window main_win;
extern Window draw_win;
extern Window command_pop;
extern Window info_pop;
int32 current_pop;
extern uint32 MyBackColor;
extern uint32 MyForeColor;
extern uint32 MyMainWinColor;
extern uint32 MyDrawWinColor;
extern uint32 GrFore;
extern uint32 GrBack;
extern GC gc;
extern GC gc_graph;
extern GC small_gc;
extern int32 xor_flag;
extern int32 DCURX;
extern int32 DCURY;
int32 num_pops;
int32 MINI_H = 300;
int32 MINI_W = 450;

extern int32 Xup;
int32 ActiveWinList[MAXPOP];

typedef struct {
    double xlo, xhi, dx;
    double *y;
    double *x;
    int32 n, flag, interp, autoeval;
    int32 xyvals;
    /* flag=0 if virgin array, flag=1 if already allocated; flag=2 for function
                             interp=0 for normal interpolation, interp=1 for
       'step' interp=2 for cubic spline table   and finally, xyvals=1 if both x
       and y vals are needed (xyvals=0 is faster lookup )*/
    char filename[128], name[12];
} TABULAR;

extern TABULAR my_table[MAX_TAB];

extern int32 NTable;

extern InternSet intern_set[MAX_INTERN_SET];
extern int32 Nintern_set;

static void select_sym(Window window);
static void lo_lite(Window wi);
static void set_gr_back(void);
static void set_gr_fore(void);
static void select_window(Window window);
static void kill_all_pops(void);
static void destroy_a_pop(void);
static void set_restore(int32 flag);
static void add_pntarr(int32 type);
static void add_markers(void);
static void add_marker(void);
static int32 get_markers_info(void);
static int32 get_marker_info(void);
static int32 select_marker_type(int32 *type);
static void destroy_label(Window window);
static void destroy_grob(Window window);
static void arrow_head(double xs, double ys, double xe, double ye, double size);
static void draw_grob(int32 i);
static void draw_marker(double x, double y, double size, int32 type);

int32
select_table(void) {
    int32 i;
    int32 j;
    Window temp = main_win;
    char *n[MAX_TAB], key[MAX_TAB], ch;
    for (i = 0; i < NTable; i++) {
        n[i] = xmalloc(25);
        key[i] = 'a' + (char)i;
        sprintf(n[i], "%c: %s", key[i], my_table[i].name);
    }
    key[NTable] = 0;
    ch = (char)pop_up_list(&temp, "Table", n, key, NTable, 12, 0, 10, 0,
                           no_hint, info_pop, info_message);
    for (i = 0; i < NTable; i++)
        free(n[i]);
    j = (int32)(ch - 'a');
    if (j < 0 || j >= NTable) {
        err_msg("Not a valid table");
        return -1;
    }
    return j;
}

void
get_intern_set(void) {
    char *n[MAX_INTERN_SET], key[MAX_INTERN_SET], ch;
    int32 i;
    int32 j;
    int32 count = Nintern_set;
    Window temp = main_win;
    if (count == 0)
        return;
    for (i = 0; i < Nintern_set; i++) {
        n[i] = xmalloc(256);
        key[i] = 'a' + (char)i;
        sprintf(n[i], "%c: %s", key[i], intern_set[i].name);
    }
    key[count] = 0;
    ch = (char)pop_up_list(&temp, "Param set", n, key, count, 12, 0, 10, 0,
                           no_hint, info_pop, info_message);
    for (i = 0; i < count; i++)
        free(n[i]);
    j = (int32)(ch - 'a');
    if (j < 0 || j >= Nintern_set) {
        err_msg("Not a valid set");
        return;
    }
    get_graph();
    extract_internset(j);
    chk_delay();
    redraw_params();
    redraw_ics();
    reset_graph();
}

void
make_icon(char *icon, int32 wid, int32 hgt, Window window) {
    Pixmap icon_map;
    XWMHints wm_hints;
    icon_map =
        XCreateBitmapFromData(display, window, icon, (uint)wid, (uint)hgt);
    wm_hints.initial_state = NormalState;
    wm_hints.input = True;
    wm_hints.icon_pixmap = icon_map;
    wm_hints.flags = StateHint | IconPixmapHint | InputHint;

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";
        XSetWMProperties(display, window, NULL, NULL, NULL, 0, NULL, &wm_hints,
                         &class_hints);
    }
    return;
}

void
title_text(char *string) {
    gtitle_text(string, draw_win);
    return;
}

void
gtitle_text(char *string, Window win) {
    XTextProperty wname;
    XTextProperty iname;
    GrCol();
    if (win != graph[0].window) {
        XStringListToTextProperty(&string, 1, &wname);
        XStringListToTextProperty(&string, 1, &iname);
        XSetWMProperties(display, win, &wname, &iname, NULL, 0, NULL, NULL,
                         NULL);
    } else {
        int32 len = (int32)strlen(string);
        int32 x, y;
        uint32 w, h, bw, de;
        int32 xs, ys = 2;
        Window root;
        XGetGeometry(display, win, &root, &x, &y, &w, &h, &bw, &de);
        xs = ((int32)w - len*DCURX) / 2;
        if (xs < 0)
            xs = 0;
        Ftext(xs, ys, string, win);
        set_color(0);
        xline(0, 18, (int32)w, 18, win);
    }
    BaseCol();
    return;
}

void
restore_off(void) {
    MyGraph->Restore = 0;
    /* MyGraph->Nullrestore=0; */
    return;
}

void
restore_on(void) {
    MyGraph->Restore = 1;
    /*  MyGraph->Nullrestore=1; */
    return;
}

void
add_label(char *s, int32 x, int32 y, int32 size, int32 font) {
    int32 i;
    double xp;
    double yp;
    scale_to_real(x, y, &xp, &yp);
    for (i = 0; i < MAXLAB; i++) {
        if (lb[i].use == 0) {
            lb[i].use = 1;
            lb[i].x = xp;
            lb[i].y = yp;
            lb[i].window = draw_win;
            lb[i].font = font;
            lb[i].size = size;
            strcpy(lb[i].s, s);
            return;
        }
    }
    return;
}

void
draw_marker(double x, double y, double size, int32 type) {
    int32 pen = 0;
    double x1 = x, y1 = y, x2, y2;
    int32 ind = 0;
    int32 offset;

    static int32 sym_dir[] = {
        /*          box              */
        0,
        -6,
        -6,
        1,
        12,
        0,
        1,
        0,
        12,
        1,
        -12,
        0,
        1,
        0,
        -12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          diamond             */
        0,
        8,
        0,
        1,
        -8,
        -8,
        1,
        8,
        -8,
        1,
        8,
        8,
        1,
        -8,
        8,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        /*          triangle         */
        0,
        -6,
        -6,
        1,
        12,
        0,
        1,
        -6,
        12,
        1,
        -6,
        -12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          plus            */
        0,
        -6,
        0,
        1,
        12,
        0,
        0,
        -6,
        -6,
        1,
        0,
        12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          cross            */
        0,
        -6,
        6,
        1,
        12,
        -12,
        0,
        -12,
        0,
        1,
        12,
        12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          circle           */
        0,
        6,
        0,
        1,
        -1,
        3,
        1,
        -2,
        2,
        1,
        -3,
        1,
        1,
        -3,
        -1,
        1,
        -2,
        -2,
        1,
        -1,
        -3,
        1,
        1,
        -3,
        1,
        2,
        -2,
        1,
        3,
        -1,
        1,
        3,
        1,
        1,
        2,
        2,
        1,
        1,
        3,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

    };
    double dx = (MyGraph->xhi - MyGraph->xlo)*WDMARK*size;
    double dy = (MyGraph->yhi - MyGraph->ylo)*HTMARK*size;
    while (true) {
        offset = 48*type + 3*ind;
        pen = sym_dir[offset];
        if (pen == 3)
            break;
        x2 = dx*sym_dir[offset + 1] + x1;
        y2 = dy*sym_dir[offset + 2] + y1;
        if (pen == 1)
            line_abs(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
        ind++;
    }
    return;
}

void
draw_grob(int32 i) {
    double xs = grob[i].xs, ys = grob[i].ys, xe = grob[i].xe, ye = grob[i].ye;
    set_linestyle(grob[i].color);
    if (grob[i].type == POINTER)
        line_abs(xs, ys, xe, ye);
    if (grob[i].type == ARROW || grob[i].type == POINTER)
        arrow_head(xs, ys, xe, ye, grob[i].size);
    if (grob[i].type >= MARKER)
        draw_marker(xs, ys, grob[i].size, grob[i].type - 2);
    return;
}

void
arrow_head(double xs, double ys, double xe, double ye, double size) {
    double l = xe - xs, h = ye - ys;
    double ar = (MyGraph->xhi - MyGraph->xlo) / (MyGraph->yhi - MyGraph->ylo);
    double x0 = xs + size*l, y0 = ys + size*h;
    /* double tot=(double)sqrt((double)(l*l+h*h)); */

    double xp = x0 + .5*size*h*ar, yp = y0 - .5*size*l / ar;
    double xm = x0 - .5*size*h*ar, ym = y0 + .5*size*l / ar;
    line_abs(xs, ys, xp, yp);
    line_abs(xs, ys, xm, ym);
    return;
}

void
destroy_grob(Window window) {
    int32 i;
    for (i = 0; i < MAXGROB; i++) {
        if ((grob[i].use == 1) && (grob[i].window == window)) {
            grob[i].use = 0;
            grob[i].window = (Window)0;
        }
    }
    return;
}

void
destroy_label(Window window) {
    int32 i;
    for (i = 0; i < MAXLAB; i++) {
        if ((lb[i].use == 1) && (lb[i].window == window)) {
            lb[i].use = 0;
            lb[i].window = (Window)0;
        }
    }
    return;
}

void
draw_label(Window window) {
    int32 i;
    GrCol();
    for (i = 0; i < MAXLAB; i++) {
        if ((lb[i].use == 1) && (lb[i].window == window))
            fancy_text_abs(lb[i].x, lb[i].y, lb[i].s, lb[i].size);
    }
    for (i = 0; i < MAXGROB; i++) {
        if ((grob[i].use == 1) && (grob[i].window == window))
            draw_grob(i);
    }
    BaseCol();
    return;
}

void
add_grob(double xs, double ys, double xe, double ye, double size, int32 type,
         int32 color) {
    int32 i;
    for (i = 0; i < MAXGROB; i++) {
        if (grob[i].use == 0) {
            grob[i].use = 1;
            grob[i].xs = xs;
            grob[i].xe = xe;
            grob[i].ys = ys;
            grob[i].ye = ye;
            grob[i].window = draw_win;
            grob[i].size = size;
            grob[i].color = color;
            grob[i].type = type;
            /*     redraw_all(); */
            return;
        }
    }
    return;
}

int32
select_marker_type(int32 *type) {
    int32 ival = *type - MARKER;
    int32 i;
    char *list[] = {"Box", "Diamond", "Triangle", "Plus", "X", "Circle"};
    static char key[] = "bdtpxc";
    Window temp = main_win;
    char ch;
    ch = (char)pop_up_list(&temp, "Markers", list, key, 6, 9, ival, 10,
                           4*DCURY + 8, no_hint, info_pop, info_message);
    if (ch == 27)
        return 0;
    for (i = 0; i < 6; i++) {
        if (ch == key[i])
            ival = i;
    }
    if (ival < 6)
        *type = MARKER + ival;

    return 1;
}

int32
get_marker_info(void) {
    static char *n[] = {"*5Type", "*4Color", "Size"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", markinfo.type);
    snprintf(values[1], sizeof(values[1]), "%d", markinfo.color);
    snprintf(values[2], sizeof(values[2]), "%g", markinfo.size);
    status = do_string_box(3, 3, 1, "Add Marker", n, values, 25);
    if (status != 0) {
        markinfo.type = atoi(values[0]);
        markinfo.size = atof(values[2]);
        markinfo.color = atoi(values[1]);
        return 1;
    }
    return 0;
}

int32
get_markers_info(void) {
    static char *n[] = {"*5Type", "*4Color", "Size", "Number", "Row1", "Skip"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", markinfo.type);
    snprintf(values[1], sizeof(values[1]), "%d", markinfo.color);
    snprintf(values[2], sizeof(values[2]), "%g", markinfo.size);
    snprintf(values[3], sizeof(values[3]), "%d", markinfo.number);
    snprintf(values[4], sizeof(values[4]), "%d", markinfo.start);
    snprintf(values[5], sizeof(values[5]), "%d", markinfo.skip);
    status = do_string_box(6, 6, 1, "Add Markers", n, values, 25);
    if (status != 0) {
        markinfo.type = atoi(values[0]);
        markinfo.size = atof(values[2]);
        markinfo.color = atoi(values[1]);
        markinfo.number = atoi(values[3]);
        markinfo.start = atoi(values[4]);
        markinfo.skip = atoi(values[5]);

        return 1;
    }
    return 0;
}

void
add_marker(void) {
    int32 flag, i1, j1, status;
    double xe = 0.0, ye = 0.0, xs, ys;
    status = get_marker_info();
    if (status == 0)
        return;
    MessageBox("Position");
    flag = GetMouseXY(&i1, &j1);
    KillMessageBox();
    XFlush(display);
    if (flag == 0)
        return;
    scale_to_real(i1, j1, &xs, &ys);
    add_grob(xs, ys, xe, ye, markinfo.size, markinfo.type, markinfo.color);
    redraw_all();
}

void
add_markers(void) {
    int32 i;
    double xe = 0.0, ye = 0.0, xs, ys, x, y, z;

    if (get_markers_info() == 0)
        return;
    for (i = 0; i < markinfo.number; i++) {
        get_data_xyz(&x, &y, &z, MyGraph->xv[0], MyGraph->yv[0], MyGraph->zv[0],
                     markinfo.start + i*markinfo.skip);
        if (MyGraph->ThreeDFlag == 0) {
            xs = x;
            ys = y;
        } else {
            threed_proj(x, y, z, &xs, &ys);
        }
        add_grob(xs, ys, xe, ye, markinfo.size, markinfo.type, markinfo.color);
    }
    redraw_all();
}

void
add_pntarr(int32 type) {
    double size = .1;
    int32 i1, j1, i2, j2, color = 0;
    double xe, ye, xs, ys;
    /*Window temp;*/
    int32 flag;
    /*temp=main_win;*/
    if (new_float("Size: ", &size))
        return;
    if (new_int("Color: ", &color))
        return;
    /* message_box(&temp,0,SCALEY-5*DCURY,"Choose start/end"); */
    MessageBox("Choose start/end");
    flag = rubber(&i1, &j1, &i2, &j2, draw_win, 1);
    /*  XDestroyWindow(display,temp); */
    KillMessageBox();
    XFlush(display);
    if (flag) {
        scale_to_real(i1, j1, &xs, &ys);
        scale_to_real(i2, j2, &xe, &ye);
        if (i1 == i2 && j1 == j2)
            return;
        add_grob(xs, ys, xe, ye, size, type, color);
        redraw_all();
    }
    return;
}

void
edit_object_com(int32 com) {
    char ans, str[80];
    int32 i, j, ilab = -1, flag, type;
    double x;
    double y;
    double dist = 1e20, dd;

    MessageBox("Choose Object");
    flag = GetMouseXY(&i, &j);

    KillMessageBox();
    XFlush(display);
    if (flag) {
        scale_to_real(i, j, &x, &y);
        /* now search all labels to find the best */
        type = 0; /* label =  0, arrows, etc =1 */
        for (i = 0; i < MAXLAB; i++) {
            if (lb[i].use == 1 && lb[i].window == draw_win) {
                dd = (x - lb[i].x)*(x - lb[i].x) +
                     (y - lb[i].y)*(y - lb[i].y);
                if (dd < dist) {
                    ilab = i;
                    dist = dd;
                }
            }
        }
        for (i = 0; i < MAXGROB; i++) {
            if (grob[i].use == 1 && grob[i].window == draw_win) {
                dd = (x - grob[i].xs)*(x - grob[i].xs) +
                     (y - grob[i].ys)*(y - grob[i].ys);
                if (dd < dist) {
                    ilab = i;
                    dist = dd;
                    type = 1;
                }
            }
        }
        if (ilab >= 0 && type == 0) {
            switch (com) {
            case 0:
                snprintf(str, sizeof(str), "Move %s ?", lb[ilab].s);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    MessageBox("Click on new position");
                    flag = GetMouseXY(&i, &j);

                    KillMessageBox();
                    XFlush(display);
                    if (flag) {
                        scale_to_real(i, j, &x, &y);
                        lb[ilab].x = x;
                        lb[ilab].y = y;
                        clr_scrn();
                        redraw_all();
                    }
                }
                break;
            case 1:
                snprintf(str, sizeof(str), "Change %s ?", lb[ilab].s);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    new_string("Text: ", lb[ilab].s);
                    new_int("Size 0-4 :", &lb[ilab].size);
                    /* new_int("Font  0-times/1-symbol :",&lb[ilab].font); */
                    if (lb[ilab].size > 4)
                        lb[ilab].size = 4;
                    if (lb[ilab].size < 0)
                        lb[ilab].size = 0;
                    clr_scrn();
                    redraw_all();
                }
                break;
            case 2:
                snprintf(str, sizeof(str), "Delete %s ?", lb[ilab].s);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    lb[ilab].window = 0;
                    lb[ilab].use = 0;
                    clr_scrn();
                    redraw_all();
                }
                break;
            default:
                break;
            }
        }
        if (ilab >= 0 && type == 1) {
            switch (com) {
            case 0:
                snprintf(str, sizeof(str), "Move graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    MessageBox("Reposition");
                    flag = GetMouseXY(&i, &j);

                    KillMessageBox();
                    XFlush(display);
                    if (flag) {
                        scale_to_real(i, j, &x, &y);
                        grob[ilab].xe = grob[ilab].xe - grob[ilab].xs + x;
                        grob[ilab].ye = grob[ilab].ye - grob[ilab].ys + y;
                        grob[ilab].xs = x;
                        grob[ilab].ys = y;
                        clr_scrn();
                        redraw_all();
                    }
                }
                break;
            case 1:
                snprintf(str, sizeof(str), "Change graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    if (grob[ilab].type >= MARKER)
                        select_marker_type(&grob[ilab].type);
                    new_float("Size ", &grob[ilab].size);
                    new_int("Color :", &grob[ilab].color);
                    clr_scrn();
                    redraw_all();
                }
                break;
            case 2:
                snprintf(str, sizeof(str), "Delete graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)TwoChoice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    grob[ilab].window = 0;
                    grob[ilab].use = 0;
                    clr_scrn();
                    redraw_all();
                }
                break;
            default:
                break;
            }
        }
    }
    return;
}

void
do_gr_objs_com(int32 com) {
    switch (com) {
    case 0:
        cput_text();
        break;
    case 1:
        add_pntarr(ARROW);
        break;
    case 2:
        add_pntarr(POINTER);
        break;
    case 3:
        add_marker();
        break;
    case 6:
        add_markers();
        break;
    case 5:
        destroy_label(draw_win);
        destroy_grob(draw_win);
        clr_scrn();
        redraw_all();
        break;
    default:
        break;
    }
    return;
}

void
set_active_windows(void) {
    int32 i, np = 0;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use == 1) {
            ActiveWinList[np] = i;
            np++;
        }
    }
    num_pops = np;
    return;
}

void
do_windows_com(int32 c) {
    switch (c) {
    case 0:
        create_a_pop();
        break;
    case 1:
        if (yes_no_box())
            kill_all_pops();
        break;
    case 3:
        XLowerWindow(display, draw_win);
        break;
    case 2:
        destroy_a_pop();
        break;
    case 5:
        set_restore(0);
        break;
    case 4:
        set_restore(1);
        break;
    case 6:
        SimulPlotFlag = 1 - SimulPlotFlag;
        break;
    default:
        break;
    }

    set_active_windows();
    return;
}

void
set_restore(int32 flag) {
    int32 i;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].window == draw_win) {
            graph[i].Restore = flag;
            graph[i].Nullrestore = flag;
            return;
        }
    }
    return;
}

void
destroy_a_pop(void) {
    int32 i;
    if (draw_win == graph[0].window) {
        respond_box("Okay", "Can't destroy big window!");
        /*respond_box(main_win,0,0,"Okay","Can't destroy big window!");*/
        return;
    }
    for (i = 1; i < MAXPOP; i++)
        if (graph[i].window == draw_win)
            break;
    if (i >= MAXPOP)
        return;
    select_window(graph[0].window);
    graph[i].Use = 0;
    destroy_label(graph[i].window);
    destroy_grob(graph[i].window);
    waitasec(ClickTime);
    XDestroySubwindows(display, graph[i].window);
    XDestroyWindow(display, graph[i].window);
    num_pops--;
    return;
}

void
init_grafs(int32 x, int32 y, int32 w, int32 h) {
    int32 i;
    GrCol();
    for (i = 0; i < MAXLAB; i++) {
        lb[i].use = 0;
        lb[i].window = (Window)0;
    }
    for (i = 0; i < MAXGROB; i++) {
        grob[i].window = (Window)0;
        grob[i].use = 0;
    }
    init_bd();
    for (i = 0; i < MAXFRZ; i++)
        frz[i].use = 0;
    /* for(i=0;i<MAXNCLINE ... */
    for (i = 0; i < MAXPOP; i++)
        graph[i].Use = 0;
    ActiveWinList[0] = 0;
    init_all_graph();

    graph[0].window = XCreateSimpleWindow(display, main_win, x, y + 4, (uint)w,
                                          (uint)h, 2, GrFore, MyDrawWinColor);
    graph[0].w_info = info_pop;

    info_message = graph[0].gr_info;
    graph[0].Use = 1;
    graph[0].Restore = 1;
    graph[0].Nullrestore = 1;
    graph[0].x0 = x;
    graph[0].y0 = y;
    graph[0].Height = h;
    graph[0].Width = w;
    XSelectInput(display, graph[0].window,
                 KeyPressMask | ButtonPressMask | ExposureMask |
                     ButtonReleaseMask | ButtonMotionMask |
                     StructureNotifyMask);
    num_pops = 1;
    XMapWindow(display, graph[0].window);
    draw_win = graph[0].window;
    current_pop = 0;
    get_draw_area();
    select_sym(graph[0].window);
    BaseCol();
}

void
ps_restore(void) {
    if (Xup) {
        redraw_dfield();
        ps_do_color(0);
        if (MyGraph->Nullrestore) {
            restore_nullclines();
            ps_stroke();
        }
    }
    ps_last_pt_off();

    restore(0, my_browser.maxrow);

    do_batch_nclines();
    do_batch_dfield();
    axes2_do();

    ps_do_color(0);
    if (Xup) {
        draw_label(draw_win);
        draw_freeze(draw_win);
    }
    ps_end();
    return;
}

void
svg_restore(void) {
    /* restore(0,my_browser.maxrow);
     */
    /*ps_do_color(0);
    if(MyGraph->Nullrestore){restore_nullclines();ps_stroke();}
     */

    redraw_dfield();
    if (MyGraph->Nullrestore)
        restore_nullclines();
    svg_last_pt_off();
    /*ps_do_color(0);*/
    restore(0, my_browser.maxrow);
    axes2_do();
    if (Xup) {
        draw_label(draw_win);
        draw_freeze(draw_win);
    }
    do_batch_nclines();
    do_batch_dfield();
    svg_end();
    return;
}

int32
rotate3dcheck(XEvent event) {
    Window window = event.xbutton.window;
    XEvent z;
    int32 xini, yini, dx, dy;
    double theta;
    double phi;
    if (window == draw_win && MyGraph->ThreeDFlag) {
        xini = event.xbutton.x;
        yini = event.xbutton.y;
        phi = MyGraph->Phi;
        theta = MyGraph->Theta;
        while (true) {
            XNextEvent(display, &z);
            if (z.type == ButtonRelease) {
                axes2_do();
                redraw_all();
                hi_lite(draw_win);
                return 1;
            }
            if (z.type == MotionNotify) {
                dx = z.xmotion.x - xini;
                dy = z.xmotion.y - yini;
                MyGraph->Phi = phi - (double)dy;
                MyGraph->Theta = theta - (double)dx;
                axes2_redraw_cube_pt(MyGraph->Theta, MyGraph->Phi);
            }
        }
    }
    return 0;
}

void
do_motion_events(XEvent event) {
    int32 i = event.xmotion.x;
    int32 j = event.xmotion.y;
    double x;
    double y;
    char buf[256];
    slider_motion(event);
#ifdef AUTO
    auto_motion(event);
#endif
    if (event.xmotion.window == draw_win) {
        scale_to_real(i, j, &x, &y);
        snprintf(buf, sizeof(buf), "x=%f y=%f ", x, y);
        canvas_xy(buf);
    }
    return;
}

void
do_expose(XEvent event) {
    int32 i;
    int32 cp = current_pop;
    Window temp;

    temp = draw_win;
    top_button_draw(event.xany.window);
    array_plot_expose(event.xany.window);
    /* redraw_txtview(ev.xany.window);  */
    ani_expose(event.xany.window);
    expose_my_browser(event);
    /* draw_info_pop(ev.xany.window); */
    RedrawMessageBox(event.xany.window);
    draw_eq_list(event.xany.window);
    draw_eq_box(event.xany.window);
    do_box_expose(event.xany.window);
    expose_slides(event.xany.window);
    menu_expose(event.xany.window);
#ifdef AUTO
    display_auto(event.xany.window);
#endif
    /* if(ev.xexpose.window==menu_pop){
          draw_help();

     return;
      }
      */
    if (manual_expose == 0) {
        GrCol();

        for (i = 0; i < MAXPOP; i++) {
            if ((graph[i].Use) && (event.xexpose.window == graph[i].w_info)) {
                XClearWindow(display, graph[i].w_info);
                if (i == 0) {
                    BaseCol();
                    XDrawString(display, graph[i].w_info, gc, 5, CURY_OFF,
                                graph[i].gr_info,
                                (int)strlen(graph[i].gr_info));
                } else {
                    SmallBase();
                    XDrawString(display, graph[i].w_info, small_gc, 0,
                                CURY_OFFs, graph[i].gr_info,
                                (int)strlen(graph[i].gr_info));
                    SmallGr();
                }
            }
            if ((event.type == Expose) && (graph[i].Use) &&
                (event.xexpose.window == graph[i].window)) {
                /* redraw_dfield(); */

                current_pop = i;
                MyGraph = &graph[i];
                draw_win = graph[i].window;
                get_draw_area();
                axes2_do();
                if (graph[i].Restore)
                    restore(0, my_browser.maxrow);
                draw_label(graph[i].window);
                draw_freeze(graph[i].window);
                if (graph[i].Nullrestore)
                    restore_nullclines();
            }
        }
    } /* namual expose */
    draw_win = temp;
    MyGraph = &graph[cp];
    current_pop = cp;
    hi_lite(draw_win);
    get_draw_area();
    BaseCol();
    SmallBase();
    return;
}

void
resize_all_pops(int32 wid, int32 hgt) {
    int32 nw = wid - 16 - 16*DCURX + 7,
          nh = hgt - 3*DCURYb - 4*DCURYs - 24;
    nw = 4*((nw / 4));
    nh = 4*((nh / 4));
    XResizeWindow(display, graph[0].window, (uint)nw, (uint)nh);
    graph[0].Width = nw;
    graph[0].Height = nh;
    get_draw_area();
    return;
}

void
kill_all_pops(void) {
    int32 i;
    select_window(graph[0].window);

    for (i = 1; i < MAXPOP; i++)
        if (graph[i].Use) {
            graph[i].Use = 0;
            destroy_label(graph[i].window);
            destroy_grob(graph[i].window);

            XDestroySubwindows(display, graph[i].window);
            XDestroyWindow(display, graph[i].window);
        }
    num_pops = 1;
    return;
}

void
create_a_pop(void) {
    int32 i;
    int32 index;

    for (i = 1; i < MAXPOP; i++)
        if (graph[i].Use == 0)
            break;
    if (i >= MAXPOP) {
        /*respond_box(main_win,0,0,"Okay","Too many windows!");*/
        respond_box("Okay", "Too many windows!");
        return;
    }
    index = i;

    graph[index].window =
        XCreateSimpleWindow(display, RootWindow(display, screen), 0, 0,
                            (uint)MINI_W, (uint)MINI_H, 2, GrFore, GrBack);
    graph[index].w_info =
        make_window(graph[index].window, 10, 0, 40*DCURXs, DCURYs, 0);
    XSetWindowBackground(display, graph[i].window, MyDrawWinColor);

    copy_graph(index, current_pop);
    graph[index].Width = MINI_W;
    graph[index].Height = MINI_H;
    graph[index].x0 = 0;
    graph[index].y0 = 0;
    num_pops++;
    make_icon((char *)graph_bits, graph_width, graph_height,
              graph[index].window);
    XSelectInput(display, graph[index].window,
                 KeyPressMask | ButtonPressMask | ExposureMask |
                     ButtonReleaseMask | ButtonMotionMask);
    XMapWindow(display, graph[index].window);
    XRaiseWindow(display, graph[index].window);
    XSetWMProtocols(display, graph[index].window, &deleteWindowAtom, 1);
    select_window(graph[index].window);
    /*  select_window(graph[0].window);
        select_window(graph[index].window); */
    XRaiseWindow(display, graph[0].window);
    /*  XDestroyWindow(display,temp); */
    return;
}

void
GrCol(void) {
    XSetForeground(display, gc, GrFore);
    XSetBackground(display, gc, GrBack);
    return;
}

void
BaseCol(void) {
    XSetForeground(display, gc, MyForeColor);
    XSetBackground(display, gc, MyBackColor);
    return;
}

void
SmallGr(void) {
    XSetForeground(display, small_gc, GrFore);
    XSetBackground(display, small_gc, GrBack);
    return;
}

void
SmallBase(void) {
    XSetForeground(display, small_gc, MyForeColor);
    XSetBackground(display, small_gc, MyBackColor);
    return;
}

void
change_plot_vars(int32 k) {
    int32 i;
    int32 ip;
    int32 np;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use) {
            np = graph[i].nvars;
            for (ip = 0; ip < np; ip++) {
                if (graph[i].xv[ip] > k)
                    graph[i].xv[ip] = graph[i].xv[ip] - 1;
                if (graph[i].yv[ip] > k)
                    graph[i].yv[ip] = graph[i].yv[ip] - 1;
                if (graph[i].zv[ip] > k)
                    graph[i].zv[ip] = graph[i].zv[ip] - 1;
            }
        }
    }
    return;
}

int32
check_active_plot(int32 k) {
    int32 i;
    int32 ip;
    int32 np;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use) {
            np = graph[i].nvars;
            for (ip = 0; ip < np; ip++) {
                if (graph[i].xv[ip] == k || graph[i].yv[ip] == k ||
                    graph[i].zv[ip] == k)
                    return 1;
            }
        }
    }
    return 0;
}

void
make_active(int32 i, int32 flag) {
    current_pop = i;
    MyGraph = &graph[current_pop];
    draw_win = MyGraph->window;
    get_draw_area_flag(flag);
    return;
}

void
select_window(Window window) {
    int32 i;

    if (window == draw_win)
        return;
    GrCol();
    if (window == graph[0].window)
        current_pop = 0;
    else {
        for (i = 1; i < MAXPOP; i++)
            if ((graph[i].Use) && (window == graph[i].window))
                current_pop = i;
    }
    MyGraph = &graph[current_pop];
    lo_lite(draw_win);
    draw_win = window;
    hi_lite(window);
    XRaiseWindow(display, window);
    get_draw_area();
    BaseCol();
    return;
}

void
set_gr_fore(void) {
    XSetForeground(display, gc, GrFore);
    return;
}

void
set_gr_back(void) {
    XSetForeground(display, gc, GrBack);
    return;
}

void
hi_lite(Window wi) {
    set_gr_fore();
    select_sym(wi);
    return;
}

void
lo_lite(Window wi) {
    set_gr_back();
    bar(0, 0, 5, 5, wi);
    return;
}

void
select_sym(Window window) {
    bar(0, 0, 5, 5, window);
    return;
}

void
canvas_xy(char *buf) {
    XClearWindow(display, MyGraph->w_info);
    strcpy(MyGraph->gr_info, buf);
    if (MyGraph->w_info == info_pop) {
        BaseCol();
        XDrawString(display, MyGraph->w_info, gc, 5, CURY_OFF, buf,
                    (int)strlen(buf));
    } else {
        SmallBase();
        XDrawString(display, MyGraph->w_info, small_gc, 0, CURY_OFFs, buf,
                    (int)strlen(buf));
        /* SmallGr(); */
    }
    return;
}

void
check_draw_button(XEvent event) {
    int32 k;
    char buf[256];

    int32 button;
    int32 i;
    int32 j;
    double x, y;
    int32 flag = 0;
    Window window;
    button = (int32)event.xbutton.button;
    window = event.xbutton.window;
    i = event.xbutton.x;
    j = event.xbutton.y;
    if (button == 1) { /* select window   */

        for (k = 1; k < MAXPOP; k++)
            if ((graph[k].Use) && (window == graph[k].window))
                flag = 1;
        if ((window == graph[1].window) || (flag == 1))
            select_window(window);
    } else /* any other button   */
    {
        if (window != draw_win)
            return;
        scale_to_real(i, j, &x, &y);
        snprintf(buf, sizeof(buf), "x=%f y=%f ", x, y);
        canvas_xy(buf);
    }
}
